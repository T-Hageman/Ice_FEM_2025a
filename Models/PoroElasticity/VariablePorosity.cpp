#include "VariablePorosity.h"
#include "../../Physics/physics.h"

/// @brief Registers the model into the model creator
void Register_VariablePorosity(){
    ModelNames.push_back("PoroElasticity/VariablePorosity");
    ModelCreators.push_back(New_VariablePorosity);
}

/// @brief Creates a new phasefield model
/// @param My_Physics input: physics object
/// @param MyNameIn input: name of this model, from where to use input properties
/// @return pointer to the newly created model
BaseModel* New_VariablePorosity(Physics& My_Physics, std::string MyNameIn){
    return new VariablePorosity(My_Physics, MyNameIn);
}

/// @brief Initializer, forwards to baseModel
/// @param My_Physics //input: physics object, reference copied into model
/// @param MyName   //name of this model, from where to use input properties
VariablePorosity::VariablePorosity(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "PoroElasticity/VariablePorosity";
}

VariablePorosity::~VariablePorosity(){

}

/// @brief Initialize this model from scratch
/// @param inputs input: JSON object representing the input file
void VariablePorosity::init(inputData& inputs){
    Setup(inputs);

    std::vector<size_t> UniqueNodes_T = mesh->ElementGroups[ElemGroupIndex_T].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes_T, dofType_T);

    if (Fluid){
        std::vector<size_t> UniqueNodes_p = mesh->ElementGroups[ElemGroupIndex_p].GetUniqueNodes();
        dofs->AddDofs(UniqueNodes_p, dofType_p);
    }

    std::vector<size_t> UniqueNodes_Por = mesh->ElementGroups[ElemGroupIndex_Por].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes_Por, dofType_Por);
}

/// @brief Loads model from a previous restart savefile
/// @param inputs input: JSON object representing the input file
/// @param data  input: restart data from where to load variables
void VariablePorosity::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);
}

/// @brief Save model to a restart file
/// @param data input: pointer to object in which to save restart data
void VariablePorosity::save(SaveDataFile& data){

}

/// @brief Performs set-up of this model based on input file
/// @param inputs input: Input json data file
void VariablePorosity::Setup(inputData& inputs){

    //material parameters
    inputs.GetRequired(SolidMatName, {"Models", MyName, "Material"});
    inputs.GetRequired(FluidMatName, {"Models", MyName, "Fluid"});

    S = new ThermalMaterial(inputs, SolidMatName);
    F = new ThermalMaterial(inputs, FluidMatName);

    
    if (dofs->hasDofType(DofNames_p)){
        Fluid = true;
        PoroRelations = new PoroElasticRelations(inputs);
    }

    //get indices for relevant element groups
    std::string ElemGroupName;
    inputs.GetRequired(ElemGroupName, {"Models", MyName, "ElementGroup_T"});
    ElemGroupIndex_T = mesh->GetElementGroupIdx(ElemGroupName);

    inputs.GetRequired(ElemGroupName, {"Models", MyName, "ElementGroup_poros"});
    ElemGroupIndex_Por = mesh->GetElementGroupIdx(ElemGroupName);    

    if (Fluid){
        inputs.GetRequired(ElemGroupName, {"Models", MyName, "ElementGroup_p"});
        ElemGroupIndex_p = mesh->GetElementGroupIdx(ElemGroupName);
        dofs->getDofTypesSteps(DofNames_p, dofType_p, Step_p);
    } else {
        Step_p = 99;
    }
    
    //get relevant degree of freedmo steps and indices
    dofs->getDofTypesSteps(DofNames_T, dofType_T, Step_T);
    dofs->getDofTypesSteps(DofNames_Por, dofType_Por, Step_Por);

    PhaseChange = false;
    inputs.GetOptional(PhaseChange,  {"Models", MyName, "PhaseChange"});
}

/// @brief Assemble the tangent stiffness matrix and residual force vector
/// @param K Tangent matrix
/// @param f residual force vector
/// @param cons constraints applied by this model [none]
/// @param step input: Step for which to assemble the tangential matrix
void VariablePorosity::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    if (step == Step_T || step==Step_Por || step == Step_p){
        size_t ipcount = mesh->ElementGroups[ElemGroupIndex_T].BaseElem->ipcount;   //number of integration points
        size_t nNodes_T = mesh->ElementGroups[ElemGroupIndex_T].NNodes_per_elem;        //number of nodes per Temperature element
        size_t nNodes_p = 0;//number of nodes per pressure element
        size_t nNodes_por = mesh->ElementGroups[ElemGroupIndex_Por].NNodes_per_elem;//number of nodes per pressure element

        if (Fluid){
            nNodes_p = mesh->ElementGroups[ElemGroupIndex_p].NNodes_per_elem;
        }
    
        std::vector<double> w(ipcount); //integration weights

        std::vector<Eigen::RowVectorXd> NT(ipcount); for (size_t ip=0; ip<ipcount;  ip++) NT[ip].resize(nNodes_T);  //temperature shape functions
        std::vector<Eigen::MatrixXd> GT(ipcount); for (size_t ip=0; ip<ipcount;  ip++) GT[ip].resize(2, nNodes_T);  //temperature shape function gradients

        std::vector<Eigen::RowVectorXd> Np(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Np[ip].resize(nNodes_p); //pressure shape functions
        std::vector<Eigen::MatrixXd> Gp(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gp[ip].resize(2, nNodes_p); //pressure shape function gradients

        std::vector<Eigen::RowVectorXd> Npor(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Npor[ip].resize(nNodes_por); //porosity shape functions
        std::vector<Eigen::MatrixXd> Gpor(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gpor[ip].resize(2, nNodes_por); //porosity shape function gradients


        std::vector<size_t> El_T(nNodes_T), El_p(nNodes_p), El_por(nNodes_por); //nodes contained within this element

        std::vector<PetscInt> dofsT(nNodes_T);  //temperature degree of freedom indices
        std::vector<PetscInt> dofsP(nNodes_p);  //phase field degree of freedom indices
        std::vector<PetscInt> dofsPor(nNodes_por);  //porosity degree of freedom indices

        Eigen::VectorXd P(nNodes_p), T(nNodes_T), dT(nNodes_T), Por(nNodes_por), PorOld(nNodes_por);  

        Eigen::VectorXd F_T(nNodes_T);    //element force vector
        Eigen::MatrixXd K_TT(nNodes_T, nNodes_T);   //tangential matrix, displacement/displacement component
        Eigen::MatrixXd K_TP(nNodes_T, nNodes_p);   //tangential matrix, displacement/displacement component
        Eigen::MatrixXd K_TPor(nNodes_T, nNodes_por);   //tangential matrix, displacement/displacement component

        Eigen::VectorXd F_Por(nNodes_por);    //element force vector
        Eigen::MatrixXd K_PorT(nNodes_por, nNodes_T);   //tangential matrix, displacement/displacement component
        Eigen::MatrixXd K_PorP(nNodes_por, nNodes_p);   //tangential matrix, displacement/displacement component
        Eigen::MatrixXd K_PorPor(nNodes_por, nNodes_por);   //tangential matrix, displacement/displacement component

        Eigen::VectorXd F_P(nNodes_por);    //element force vector
        Eigen::MatrixXd K_PT(nNodes_por, nNodes_T);   //tangential matrix, displacement/displacement component
        Eigen::MatrixXd K_PP(nNodes_por, nNodes_p);   //tangential matrix, displacement/displacement component
        Eigen::MatrixXd K_PPor(nNodes_por, nNodes_por);   //tangential matrix, displacement/displacement component

        Eigen::VectorXd WLumped(nNodes_T);

        double dXdt = physics->timeScheme->du_dt;

        for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex_T].NElems; el++){   //loop over all elements

            //get shape functions and gradients
            mesh->getShapeGrads(ElemGroupIndex_T, el, El_T, NT, GT, w);
            mesh->getShapeGrads(ElemGroupIndex_Por, el, El_por, Npor, Gpor, w);

            //obtain degrees of freedom and state vectors
            dofs->getDofForNodes(El_T, dofType_T, dofsT);
            dofs->getDofForNodes(El_por, dofType_Por, dofsPor);

            physics->StateVectors[Step_T].GetValues(dofsT, T);
            physics->dStateVectors[Step_T].GetValues(dofsT, dT);
            physics->StateVectors[Step_Por].GetValues(dofsPor, Por);
            physics->StateVectorsOld[Step_Por].GetValues(dofsPor, PorOld);

            if (Fluid){
                mesh->getShapeGrads(ElemGroupIndex_p, el, El_p, Np, Gp, w);
                dofs->getDofForNodes(El_p, dofType_p, dofsP);
                physics->StateVectors[Step_p].GetValues(dofsP, P);
            }

            //zero matrices before assembly
            F_T.setZero();
            K_TT.setZero();
            K_TP.setZero();
            K_TPor.setZero();

            F_Por.setZero();
            K_PorT.setZero();
            K_PorP.setZero();
            K_PorPor.setZero();

            F_P.setZero();
            K_PT.setZero();
            K_PP.setZero();
            K_PPor.setZero();

            WLumped.setZero();

            for (size_t ip = 0; ip < ipcount; ip++){
                //double p=Np[ip]*P;
                //double dSw, ddSw;
                //double Sw = PoroRelations->GetSaturation(p ,dSw, ddSw);

                WLumped += w[ip]*Npor[ip];
            }

            for (size_t i = 0; i < nNodes_T; i++){
                //capacity
                F_Por(i)      += WLumped(i)*(Por(i)-PorOld(i))/physics->timeScheme->dt;
                K_PorPor(i,i) += WLumped(i)/physics->timeScheme->dt;

                if (Por(i)<0.0){
                    double KDummy = 1e2;
                    F_Por(i) += WLumped(i)*KDummy*Por(i);
                    K_PorPor(i,i) += WLumped(i)*KDummy;
                }

                if (PhaseChange && Fluid){
                    double p=P(i);
                    double dSw, ddSw;
                    double Sw = PoroRelations->GetSaturation(p ,dSw, ddSw);
                    double tau = 1.0;

                    double q=0.0, dq_dpor=0.0, dq_dT=0.0, dq_dp=0.0;
                    if (T(i)>=0){ //melting
                        q       = -tau*F->Density*S->LatentHeat*T(i)*(1-Por(i))/physics->timeScheme->dt;
                        dq_dpor = tau*F->Density*S->LatentHeat*T(i)/physics->timeScheme->dt;
                        dq_dT   = -tau*F->Density*S->LatentHeat*(1-Por(i))/physics->timeScheme->dt;
                        dq_dp   = 0.0;
                    } else { //freezing   (positive q = freezing)
                        q       = -tau*F->Density*S->LatentHeat*T(i)*(Sw-PoroRelations->S_0)*Por(i)/physics->timeScheme->dt;
                        dq_dpor = -tau*F->Density*S->LatentHeat*T(i)*(Sw-PoroRelations->S_0)/physics->timeScheme->dt;
                        dq_dT   = -tau*F->Density*S->LatentHeat*(Sw-PoroRelations->S_0)*Por(i)/physics->timeScheme->dt;
                        dq_dp   = -tau*F->Density*S->LatentHeat*T(i)*dSw*Por(i)/physics->timeScheme->dt;
                    }
                    F_Por(i)     += WLumped(i)/S->LatentHeat/S->Density*q;
                    K_PorPor(i,i)+= WLumped(i)/S->LatentHeat/S->Density*dq_dpor;
                    K_PorT(i,i)  += WLumped(i)/S->LatentHeat/S->Density*dq_dT;
                    K_PorP(i,i)  += WLumped(i)/S->LatentHeat/S->Density*dq_dp;

                    F_T(i)       += -WLumped(i)*q;
                    K_TT(i,i)    += -WLumped(i)*dq_dT;
                    K_TP(i,i)    += -WLumped(i)*dq_dp;
                    K_TPor(i,i)  += -WLumped(i)*dq_dpor;

                    F_P(i)     += WLumped(i)/S->LatentHeat/F->Density*q;
                    K_PPor(i,i)+= WLumped(i)/S->LatentHeat/F->Density*dq_dpor;
                    K_PT(i,i)  += WLumped(i)/S->LatentHeat/F->Density*dq_dT;
                    K_PP(i,i)  += WLumped(i)/S->LatentHeat/F->Density*dq_dp + WLumped(i)*0.0e-11; //added offset for stabilisation
                    
                }
            }

            if (step==Step_T){
                VecAdd(f, dofsT, F_T);
                MatAdd(K, dofsT, dofsT, K_TT);
                if (step==Step_Por){
                    MatAdd(K, dofsT, dofsPor, K_TPor);
                }
                if (step==Step_p){
                    MatAdd(K, dofsT, dofsP, K_TP);
                }
            }
            if (step==Step_Por){
                VecAdd(f, dofsPor, F_Por);
                MatAdd(K, dofsPor, dofsPor, K_PorPor);
                if (step==Step_T){
                    MatAdd(K, dofsPor, dofsT, K_PorT);
                }
                if (step==Step_p){
                    MatAdd(K, dofsPor, dofsP, K_PorP);
                }
            }
            if (step==Step_p){
                VecAdd(f, dofsP, F_P);
                MatAdd(K, dofsP, dofsP, K_PP);
                if (step==Step_Por){
                    MatAdd(K, dofsP, dofsPor, K_PPor);
                }
                if (step==Step_T){
                    MatAdd(K, dofsP, dofsT, K_PT);
                }
            }
        }
    }
}

/// @brief Save integration-point specific values to output files (for later post-processing)
/// @param SaveLoc input: location of the variables to save (ip/nodes)
/// @param DataName input: name of parameters to save (not necesarily coming from this model)
/// @param ElemGroup input: element group to save outputs for
/// @param Data in/outputs: object to save data into
/// @return output: has this saved any data?
bool VariablePorosity::SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>>& Data){

    return false;
}