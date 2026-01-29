#include "ThermalModel.h"
#include "../../Physics/physics.h"

/// @brief Registers the model into the model creator
void Register_ThermalModel(){
    ModelNames.push_back("Thermal/ThermalModel");
    ModelCreators.push_back(New_ThermalModel);
}

/// @brief Creates a new phasefield model
/// @param My_Physics input: physics object
/// @param MyNameIn input: name of this model, from where to use input properties
/// @return pointer to the newly created model
BaseModel* New_ThermalModel(Physics& My_Physics, std::string MyNameIn){
    return new ThermalModel(My_Physics, MyNameIn);
}

/// @brief Initializer, forwards to baseModel
/// @param My_Physics //input: physics object, reference copied into model
/// @param MyName   //name of this model, from where to use input properties
ThermalModel::ThermalModel(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "Thermal/ThermalModel";
}

ThermalModel::~ThermalModel(){

}

/// @brief Initialize this model from scratch
/// @param inputs input: JSON object representing the input file
void ThermalModel::init(inputData& inputs){
    Setup(inputs);

    std::vector<size_t> UniqueNodes_T = mesh->ElementGroups[ElemGroupIndex_T].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes_T, dofType_T);
}

/// @brief Loads model from a previous restart savefile
/// @param inputs input: JSON object representing the input file
/// @param data  input: restart data from where to load variables
void ThermalModel::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);
}

/// @brief Save model to a restart file
/// @param data input: pointer to object in which to save restart data
void ThermalModel::save(SaveDataFile& data){

}

/// @brief Performs set-up of this model based on input file
/// @param inputs input: Input json data file
void ThermalModel::Setup(inputData& inputs){

    //get indices for relevant element groups
    std::string ElemGroupName;

    inputs.GetRequired(ElemGroupName, {"Models", MyName, "ElementGroup_T"});
    ElemGroupIndex_T = mesh->GetElementGroupIdx(ElemGroupName);
    dofs->getDofTypesSteps(DofNames_T, dofType_T, Step_T);

    //material parameters
    inputs.GetRequired(SolidMatName, {"Models", MyName, "Material"});
    M = new ThermalMaterial(inputs, SolidMatName);
}

/// @brief Assemble the tangent stiffness matrix and residual force vector
/// @param K Tangent matrix
/// @param f residual force vector
/// @param cons constraints applied by this model [none]
/// @param step input: Step for which to assemble the tangential matrix
void ThermalModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    if (step == Step_T){
        size_t ipcount = mesh->ElementGroups[ElemGroupIndex_T].BaseElem->ipcount;   //number of integration points
        size_t nNodes_T = mesh->ElementGroups[ElemGroupIndex_T].NNodes_per_elem;        //number of nodes per Temperature element
    
        std::vector<double> w(ipcount); //integration weights

        std::vector<Eigen::RowVectorXd> NT(ipcount); for (size_t ip=0; ip<ipcount;  ip++) NT[ip].resize(nNodes_T);  //displacement shape functions
        std::vector<Eigen::MatrixXd> GT(ipcount); for (size_t ip=0; ip<ipcount;  ip++) GT[ip].resize(mesh->dim, nNodes_T);  //displacement shape function gradients

        std::vector<size_t> El_T(nNodes_T); //nodes contained within this element
        std::vector<PetscInt> dofsT(nNodes_T);        //displacement degree of freedom indices

        Eigen::VectorXd T(nNodes_T), dT(nNodes_T);  

        Eigen::VectorXd F_T(nNodes_T);    //element force vector
        Eigen::MatrixXd K_TT(nNodes_T, nNodes_T);   //tangential matrix, displacement/displacement component

        Eigen::VectorXd WLumped(nNodes_T);

        double dXdt = physics->timeScheme->du_dt;

        for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex_T].NElems; el++){   //loop over all elements

            //get shape functions and gradients
            mesh->getShapeGrads(ElemGroupIndex_T, el, El_T, NT, GT, w);

            //obtain degrees of freedom and state vectors
            dofs->getDofForNodes(El_T, dofType_T, dofsT);

            physics->StateVectors[Step_T].GetValues(dofsT, T);
            physics->dStateVectors[Step_T].GetValues(dofsT, dT);

            //zero matrices before assembly
            F_T.setZero();
            K_TT.setZero();
            WLumped.setZero();

            for (size_t ip = 0; ip < ipcount; ip++){

                //capacity 
                // double cap = poros*Sw*cp_w*density_w + (1-poros)*cp_s*density_s;
                // double dcap_dp = poros*dSw*cp_w*density_w;

                // F_T  += w[ip]*cap*NT[ip].transpose()*NT[ip]*dT;
                // K_TT += w[ip]*cap*dXdt*NT[ip].transpose()*NT[ip];
                // K_TP += -w[ip]*dcap_dp*NT[ip].transpose()*(NT[ip]*dT)*Np[ip];

                //diffusion
                F_T += w[ip]*M->k*GT[ip].transpose()*GT[ip]*T;
                K_TT+= w[ip]*M->k*GT[ip].transpose()*GT[ip];

                WLumped += w[ip]*NT[ip];
            }

            for (size_t i = 0; i < nNodes_T; i++){
                F_T(i)   += WLumped(i)*M->Density*M->cp*dT(i);
                K_TT(i,i)+= WLumped(i)*M->Density*M->cp*dXdt;
            }

            VecAdd(f, dofsT, F_T);
            MatAdd(K, dofsT, dofsT, K_TT);
        }
    }
}

/// @brief Save integration-point specific values to output files (for later post-processing)
/// @param SaveLoc input: location of the variables to save (ip/nodes)
/// @param DataName input: name of parameters to save (not necesarily coming from this model)
/// @param ElemGroup input: element group to save outputs for
/// @param Data in/outputs: object to save data into
/// @return output: has this saved any data?
bool ThermalModel::SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>>& Data){

    return false;
}