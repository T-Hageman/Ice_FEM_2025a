#include "GeneralSolidModel.h"
#include "../../Physics/physics.h"
#include <Eigen/Dense>

/// @brief Registers the model into the model creator
void Register_GeneralSolidModel(){
    ModelNames.push_back("GeneralSolidModel");
    ModelCreators.push_back(New_GeneralSolidModel);
}

/// @brief Creates a new general solid model
/// @param My_Physics input: physics object
/// @param MyNameIn input: name of this model, from where to use input properties
/// @return pointer to the newly created model
BaseModel* New_GeneralSolidModel(Physics& My_Physics, std::string MyNameIn){
    return new GeneralSolidModel(My_Physics, MyNameIn);
}

/// @brief Initializer, forwards to baseModel
/// @param My_Physics //input: physics object, reference copied into model
/// @param MyName   //name of this model, from where to use input properties
GeneralSolidModel::GeneralSolidModel(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "GeneralSolidModel";
}

GeneralSolidModel::~GeneralSolidModel(){

}

/// @brief Initialize this model from scratch
/// @param inputs input: JSON object representing the input file
void GeneralSolidModel::init(inputData& inputs){
    Setup(inputs);

    std::vector<size_t> UniqueNodes_u = mesh->ElementGroups[ElemGroupIndex_u].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes_u, dofTypes_u);

    size_t NElems = mesh->ElementGroups[ElemGroupIndex_u].NElems;
    size_t NIP = mesh->ElementGroups[ElemGroupIndex_u].BaseElem->ipcount;
    Plastic_Strain.resize(NElems);
    Plastic_StrainOld.resize(NElems);
    for (size_t i = 0; i < NElems; i++){
        Plastic_Strain[i].resize(NIP);
        Plastic_StrainOld[i].resize(NIP);
        for (size_t j = 0; j < NIP; j++){
            Plastic_Strain[i][j].setZero();
            Plastic_StrainOld[i][j].setZero();
        }
    }
    if (M->HistSize>0){
        MatHist.resize(NElems);
        MatHistOld.resize(NElems);
        for (size_t i = 0; i < NElems; i++){
            MatHist[i].resize(NIP);
            MatHistOld[i].resize(NIP);
            for (size_t j = 0; j < NIP; j++){
                MatHist[i][j].resize(M->HistSize);
                MatHistOld[i][j].resize(M->HistSize);
                for (size_t k = 0; k < M->HistSize; k++){
                    MatHist[i][j][k] = 0.0;
                    MatHistOld[i][j][k] = 0.0;
                }
                
            }
        }
    }
}

/// @brief Loads model from a previous restart savefile
/// @param inputs input: JSON object representing the input file
/// @param data  input: restart data from where to load variables
void GeneralSolidModel::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);
    //data.Load("Energy", Energy_Hist);
    //data.Load("PlasticStrains", Plastic_Strain);
    Commit(CommitTypes::TIMEDEP_COMMIT_TYPE);
}

/// @brief Save model to a restart file
/// @param data input: pointer to object in which to save restart data
void GeneralSolidModel::save(SaveDataFile& data){
    //data.Save("Energy", Energy_Hist);
    //data.Save("PlasticStrains", Plastic_Strain);
}

/// @brief Commits time dependent history variables (Energy history and plastic strains)
/// @param CommitType
void GeneralSolidModel::Commit(int CommitType){
    if (CommitType==CommitTypes::TIMEDEP_COMMIT_TYPE){
        size_t NElems = mesh->ElementGroups[ElemGroupIndex_u].NElems;
        size_t NIP = mesh->ElementGroups[ElemGroupIndex_u].BaseElem->ipcount;
        for (size_t i = 0; i < NElems; i++){
            for (size_t j = 0; j < NIP; j++){
                Plastic_StrainOld[i][j] = 1.0*Plastic_Strain[i][j];
                if (M->HistSize>0){
                    for (size_t k = 0; k < M->HistSize; k++){
                        MatHistOld[i][j][k] = MatHist[i][j][k];
                    }
                };
            }
        }
    }
}

void GeneralSolidModel::ResetStep(){
    size_t NElems = mesh->ElementGroups[ElemGroupIndex_u].NElems;
    size_t NIP = mesh->ElementGroups[ElemGroupIndex_u].BaseElem->ipcount;
    for (size_t i = 0; i < NElems; i++){
        for (size_t j = 0; j < NIP; j++){
            Plastic_Strain[i][j] = 1.0*Plastic_StrainOld[i][j];
            if (M->HistSize>0){
                for (size_t k = 0; k < M->HistSize; k++){
                    MatHist[i][j][k] = MatHistOld[i][j][k];
                }
            }
        }
    }
}

/// @brief Performs set-up of this model based on input file
/// @param inputs input: Input json data file
void GeneralSolidModel::Setup(inputData& inputs){
    if (dim==3){
        DofNames_u.push_back("uz");
        g = {0, 0, -9.81};
    }

    //get indices for relevant element groups
    std::string ElemGroupName;
    inputs.GetRequired(ElemGroupName, {"Models", MyName, "ElementGroup"});
    ElemGroupIndex_u = mesh->GetElementGroupIdx(ElemGroupName);

    //get relevant degree of freedmo steps and indices
    std::vector<size_t> dofSteps;
    dofs->getDofTypesSteps(DofNames_u, dofTypes_u, dofSteps);
    if (dofSteps[0] != dofSteps[1]) throw std::invalid_argument(ModelName+" requires "+ DofNames_u[0] + " and " + DofNames_u[1] + " to be in the same solver step\n");
    Step_u = dofSteps[0];


    //material parameters
    inputs.GetRequired(SolidMatName, {"Models", MyName, "Material"});
    M = new SolidMaterial(inputs, SolidMatName);

    IncludeInertia = false;
    inputs.GetOptional(IncludeInertia, {"Models", MyName, "Intertia"});

    bool HasGravity = false;
    inputs.GetOptional(HasGravity, {"Models", MyName, "Gravity"});
    if (HasGravity == false){
        g.setZero();
    }
}

void GeneralSolidModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    double t = physics->time;
    double dt = physics->timeScheme->dt;

    if (step == Step_u){    //assemble stiffness matrix for momentum balance
        Assemble_U(f, K);
    }

}

void GeneralSolidModel::Assemble_U(Vec &f, Mat &K){
    size_t ipcount = mesh->ElementGroups[ElemGroupIndex_u].BaseElem->ipcount;    // number of integration points
    size_t nNodes_u = mesh->ElementGroups[ElemGroupIndex_u].NNodes_per_elem;         // number of nodes per displacement element

    std::vector<double> w(ipcount); // integration weights
    std::vector<Eigen::RowVectorXd> Nu(ipcount); for (size_t ip = 0; ip < ipcount; ip++) Nu[ip].resize(nNodes_u); // displacement shape functions
    std::vector<Eigen::MatrixXd> Gu(ipcount); for (size_t ip = 0; ip < ipcount; ip++) Gu[ip].resize(dim, nNodes_u);        // displacement shape function gradients
    Eigen::MatrixXd B(6, dim * nNodes_u);    // stress to strain mapping matrix
    Eigen::MatrixXd N_uu(3, dim * nNodes_u); // state to displacement mapping matrix

    std::vector<size_t> El_u(nNodes_u); // nodes contained within this element
    std::vector<PetscInt> dofsU(dim * nNodes_u);     // displacement degree of freedom indices

    Eigen::VectorXd U(dim * nNodes_u), dU(dim * nNodes_u), ddU(dim * nNodes_u), UOld(dim * nNodes_u); // current displacement, velocity, acceleration, and old displacement
    Vector6d Stress;                                             
    Vector6d Strain, StrainOld, Strain_Pl, Strain_Pl_Old;                 // total and plastic strain vectors
    Matrix6d DPlast;                                                   // tangential matrices (damageable, undamageable, plastic)

    double ddXdt2 = physics->timeScheme->ddu_dt;                       // state to acceleration derivative
    double dXdt = physics->timeScheme->du_dt;
    double dt= physics->timeScheme->dt;

    Eigen::VectorXd F_U(dim * nNodes_u);                // element force vector
    Eigen::MatrixXd K_UU(dim * nNodes_u, dim * nNodes_u); // tangential matrix, displacement/displacement component
    Eigen::VectorXd WLumped(nNodes_u);

    for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex_u].NElems; el++){ // loop over all elements
        // get shape functions and gradients
        mesh->getShapeGrads(ElemGroupIndex_u, el, El_u, Nu, Gu, w);
        dofs->getDofForNodes(El_u, dofTypes_u, dofsU);

        physics->StateVectors[Step_u].GetValues(dofsU, U);
        physics->StateVectorsOld[Step_u].GetValues(dofsU, UOld);
        physics->dStateVectors[Step_u].GetValues(dofsU, dU);
        physics->ddStateVectors[Step_u].GetValues(dofsU, ddU);

        // zero matrices before assembly
        F_U.setZero();
        K_UU.setZero();

        WLumped.setZero();

        // integration over element
        for (size_t ip = 0; ip < ipcount; ip++){
            // special matrices
            ConstructBMat(B, Gu[ip], dim);

            N_uu.setZero();
            N_uu(0, Eigen::seq(0, nNodes_u - 1)) = Nu[ip];
            N_uu(1, Eigen::seq(nNodes_u, 2 * nNodes_u - 1)) = Nu[ip];
            if (dim==3){
                N_uu(2, Eigen::seq(2*nNodes_u, 3 * nNodes_u - 1)) = Nu[ip];
            }

            // stiffness
            Strain = B * U;
            StrainOld = B * UOld;
            Strain_Pl_Old = Plastic_StrainOld[el][ip];
            Strain_Pl = Plastic_Strain[el][ip];

            double PlasticDiss, DamPlasticDiss;
            if (M->HistSize>0){
                M->UpdatePlasticStrains(DPlast, Strain_Pl, Strain, Strain_Pl_Old, StrainOld, 0.0, 0.0, dt, MatHist[el][ip], MatHistOld[el][ip], PlasticDiss, DamPlasticDiss);
            } else {
                std::vector<double> Unused, Unused2;
                M->UpdatePlasticStrains(DPlast, Strain_Pl, Strain, Strain_Pl_Old, StrainOld, 0.0, 0.0, dt, Unused, Unused2, PlasticDiss, DamPlasticDiss);
            }
            Plastic_Strain[el][ip] = Strain_Pl;


            Strain = Strain - Strain_Pl;
            Stress = M->D * Strain;

            F_U += w[ip] * B.transpose() * Stress;
            K_UU += w[ip] * B.transpose() * M->D * DPlast * B;

            // Lumped weights
            WLumped += w[ip] * Nu[ip];// * DamI;
        }
        for (size_t n = 0; n < nNodes_u; n++){
            int ny = n + nNodes_u;
            int nz = n + 2*nNodes_u;

            // inertia
            if (IncludeInertia){
                F_U(n) += WLumped(n) * M->Density * ddU(n);
                K_UU(n, n) += WLumped(n)  * M->Density * ddXdt2;

                F_U(ny) += WLumped(n) * M->Density * ddU(ny);
                K_UU(ny, ny) += WLumped(n)  * M->Density * ddXdt2;

                if (dim==3){
                    F_U(nz) += WLumped(n) * M->Density * ddU(nz);
                    K_UU(nz, nz) += WLumped(n)  * M->Density * ddXdt2;
                }
            }

            // Damping
            if (M->Damping>0.0){
                F_U(n) += WLumped(n) * M->Damping * dU(n);
                K_UU(n, n) += WLumped(n) * M->Damping * dXdt;

                F_U(ny) += WLumped(n) * M->Damping * dU(ny);
                K_UU(ny, ny) += WLumped(n) * M->Damping * dXdt;

                if (dim==3){
                    F_U(nz) += WLumped(n) * M->Damping * dU(nz);
                    K_UU(nz, nz) += WLumped(n) * M->Damping * dXdt;
                }
            }

            //gravity
            F_U(n) += -WLumped(n) * M->Density * g(0);
            F_U(ny)+= -WLumped(n) * M->Density * g(1);
            if (dim==3){
                F_U(nz)+= -WLumped(n) * M->Density * g(2);
            }
        }

        // add to total vector and matrix
        VecAdd(f, dofsU, F_U);
        MatAdd(K, dofsU, dofsU, K_UU);
    }
}


/// @brief Save integration-point specific values to output files (for later post-processing)
/// @param SaveLoc input: location of the variables to save (ip/nodes)
/// @param DataName input: name of parameters to save (not necesarily coming from this model)
/// @param ElemGroup input: element group to save outputs for
/// @param Data in/outputs: object to save data into
/// @return output: has this saved any data?
bool GeneralSolidModel::SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>>& Data){

    //save displacement norm (mainly for plotting)
    if (SaveLoc=="Nodes" && DataName == "unorm" && ElemGroup == ElemGroupIndex_u){
        size_t nNodes_u = mesh->ElementGroups[ElemGroupIndex_u].NNodes_per_elem;//number of nodes per pressure element

        std::vector<Eigen::RowVectorXd> NExport(mesh->ElementGroups[ElemGroupIndex_u].BaseElem->NExport.size());
        for (size_t i = 0; i < NExport.size(); i++) NExport[i].resize(nNodes_u);

        std::vector<size_t> El_u(nNodes_u); //nodes contained within this element
        std::vector<PetscInt> dofsU(dim*nNodes_u);  //phase field degree of freedom indices

        Eigen::MatrixXd N_uu(3,dim*nNodes_u); //state to displacement mapping matrix
        Eigen::VectorXd U(dim*nNodes_u);  //current phasefield parameter and phasefield time derivative

        for (size_t i = 0; i < mesh->ElementGroups[ElemGroupIndex_u].NElems; i++){
            mesh->getExportShape(ElemGroupIndex_u, i, El_u, NExport);

            mesh->GetNodesForElem(El_u, ElemGroupIndex_u, i);
            dofs->getDofForNodes(El_u, dofTypes_u, dofsU);
            physics->StateVectors[Step_u].GetValues(dofsU, U);

            for (size_t j = 0; j < NExport.size(); j++){
                N_uu.setZero();
                N_uu(0,Eigen::seq(0,nNodes_u-1)) = NExport[j];
                N_uu(1,Eigen::seq(nNodes_u,2*nNodes_u-1)) = NExport[j];
                if (dim==3){
                    N_uu(2,Eigen::seq(2*nNodes_u,3*nNodes_u-1)) = NExport[j];
                }

                Eigen::Vector3d disp = N_uu*U;
                double dispNorm = disp.norm();
                Data[i][j] = dispNorm;
            }
        }
        return true;
    }

    //save stresses
    if (SaveLoc == "ip" && ElemGroup == ElemGroupIndex_u && (DataName == "sxx" || DataName == "syy" || DataName == "szz" || DataName == "sxy" || DataName == "syz" || DataName == "sxz" || DataName == "s1" || DataName == "s2" || DataName == "s3")){
        SaveStressComponent(DataName, Data, true);
        return true;
    }

    //save plastic strains
    if (SaveLoc == "ip" && ElemGroup == ElemGroupIndex_u && (DataName == "exxp" || DataName == "eyyp" || DataName == "ezzp" || DataName == "exyp" || DataName == "eyzp" || DataName == "exzp")){
        size_t comp;
        if (DataName == "exxp"){comp=0;}
        if (DataName == "eyyp"){comp=1;}
        if (DataName == "ezzp"){comp=2;}
        if (DataName == "exyp"){comp=3;}
        if (DataName == "eyzp"){comp=4;}
        if (DataName == "exzp"){comp=5;}
        for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex_u].NElems; el++){
            for (size_t ip = 0; ip < mesh->ElementGroups[ElemGroupIndex_u].BaseElem->ipcount; ip++){
                Data[el][ip] = Plastic_Strain[el][ip](comp);
            }
        }

        return true;
    }

    //save total strains
    if (SaveLoc == "ip" && ElemGroup == ElemGroupIndex_u && (DataName == "exx" || DataName == "eyy" || DataName == "ezz" || DataName == "exy" || DataName == "eyz" || DataName == "exz")){
        SaveStrainComponent(DataName, Data);
        return true;
    }

    return false;
}


void GeneralSolidModel::SaveStrainComponent(std::string DataName, std::vector<std::vector<double>> &Data){
    size_t Comp;
    if (DataName == "exx"){
        Comp = 0;
    }
    if (DataName == "eyy"){
        Comp = 1;
    }
    if (DataName == "ezz"){
        Comp = 2;
    }
    if (DataName == "exy"){
        Comp = 3;
    }
    if (DataName == "eyz"){
        Comp = 4;
    }
    if (DataName == "exz"){
        Comp = 5;
    }

    size_t ipcount = mesh->ElementGroups[ElemGroupIndex_u].BaseElem->ipcount;
    size_t nNodes_u = mesh->ElementGroups[ElemGroupIndex_u].NNodes_per_elem;

    std::vector<double> w(ipcount);

    std::vector<Eigen::RowVectorXd> Nu(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Nu[ip].resize(nNodes_u);
    std::vector<Eigen::MatrixXd> Gu(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gu[ip].resize(dim, nNodes_u);
    Eigen::MatrixXd B(6, dim*nNodes_u);

    std::vector<size_t> El_u(nNodes_u);

    std::vector<PetscInt> dofsU(dim*nNodes_u);

    Eigen::VectorXd U(dim*nNodes_u);
    Vector6d Strain;

    for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex_u].NElems; el++){
        mesh->getShapeGrads(ElemGroupIndex_u, el, El_u, Nu, Gu, w);
        dofs->getDofForNodes(El_u, dofTypes_u, dofsU);
        physics->StateVectors[Step_u].GetValues(dofsU, U);

        for (size_t ip = 0; ip < ipcount; ip++){
            //special matrices
            ConstructBMat(B, Gu[ip], dim);

            Strain = B*U;
            Data[el][ip] = Strain(Comp);
        }
    }
}


/// @brief Saves the stress components to a file
/// @param DataName
/// @param Data
/// @param Damaged
void GeneralSolidModel::SaveStressComponent(std::string DataName, std::vector<std::vector<double>>& Data, bool Damaged){
    size_t Comp;
    bool Pri;
    if (DataName == "sxx"){
        Comp = 0;
        Pri=false;
    }
    if (DataName == "syy"){
        Comp = 1;
        Pri=false;
    }
    if (DataName == "szz"){
        Comp = 2;
        Pri=false;
    }
    if (DataName == "sxy"){
        Comp = 3;
        Pri=false;
    }
    if (DataName == "syz"){
        Comp = 4;
        Pri=false;
    }
    if (DataName == "sxz"){
        Comp = 5;
        Pri=false;
    }
    if (DataName == "s1"){
        Comp = 0;
        Pri=true;
    }
    if (DataName == "s2"){
        Comp = 1;
        Pri=true;
    }
    if (DataName == "s3"){
        Comp = 2;
        Pri=true;
    }

    size_t ipcount = mesh->ElementGroups[ElemGroupIndex_u].BaseElem->ipcount;
    size_t nNodes_u = mesh->ElementGroups[ElemGroupIndex_u].NNodes_per_elem;
    size_t nNodes_Phase = mesh->ElementGroups[ElemGroupIndex_u].NNodes_per_elem;

    std::vector<double> w(ipcount);

    std::vector<Eigen::RowVectorXd> Nu(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Nu[ip].resize(nNodes_u);
    std::vector<Eigen::MatrixXd> Gu(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gu[ip].resize(dim, nNodes_u);
    Eigen::MatrixXd B(6, dim*nNodes_u);

    std::vector<size_t> El_u(nNodes_u);

    std::vector<PetscInt> dofsU(dim*nNodes_u);

    Eigen::VectorXd U(dim*nNodes_u);
    Vector6d Stress, Strain;
    Matrix6d D;
    Eigen::Matrix3d Stress_Total;
    Eigen::Vector3d PriStresses;

    double Dam, dDamdPhi, DamI, dDamIdPhi, NotNeeded;

    for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex_u].NElems; el++){
        mesh->getShapeGrads(ElemGroupIndex_u, el, El_u, Nu, Gu, w);
        dofs->getDofForNodes(El_u, dofTypes_u, dofsU);
        physics->StateVectors[Step_u].GetValues(dofsU, U);

        for (size_t ip = 0; ip < ipcount; ip++){

            //special matrices
            ConstructBMat(B, Gu[ip], dim);

            //stiffness
            Strain = B*U-Plastic_Strain[el][ip];
            Stress = M->D * Strain;

            if (Pri==false){
                if (Damaged){
                    Data[el][ip] = Stress(Comp);
                } else {
                    Data[el][ip] = Stress(Comp);
                }
            } else {
                Stress_Total.setZero();
                Stress_Total(0,0) = Stress(0) ;
                Stress_Total(1,1) = Stress(1) ;
                Stress_Total(2,2) = Stress(2) ;
                Stress_Total(0,1) = Stress(3) ;
                Stress_Total(1,0) = Stress(3) ;
                Stress_Total(1,2) = Stress(4) ;
                Stress_Total(2,1) = Stress(4) ;
                Stress_Total(0,2) = Stress(5) ;
                Stress_Total(2,0) = Stress(5) ;

                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(Stress_Total);
                PriStresses = eigensolver.eigenvalues();
                Data[el][ip] = PriStresses(Comp);
            }
        }
    }
}
