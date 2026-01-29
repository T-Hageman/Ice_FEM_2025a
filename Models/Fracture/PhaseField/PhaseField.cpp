
#include "PhaseField.h"
#include "../../../Physics/physics.h"

#include <Eigen/Dense>

/// @brief Registers the model into the model creator
void Register_PhaseFieldDamage(){
    ModelNames.push_back("Fracture/PhaseField");
    ModelCreators.push_back(New_PhaseFieldDamage);
}

/// @brief Creates a new phasefield model
/// @param My_Physics input: physics object
/// @param MyNameIn input: name of this model, from where to use input properties
/// @return pointer to the newly created model
BaseModel* New_PhaseFieldDamage(Physics& My_Physics, std::string MyNameIn){
    return new PhaseFieldDamage(My_Physics, MyNameIn);
}

/// @brief Initializer, forwards to baseModel
/// @param My_Physics //input: physics object, reference copied into model
/// @param MyName   //name of this model, from where to use input properties
PhaseFieldDamage::PhaseFieldDamage(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "Fracture/PhaseField";
}

PhaseFieldDamage::~PhaseFieldDamage(){

}

/// @brief Initialize this model from scratch
/// @param inputs input: JSON object representing the input file
void PhaseFieldDamage::init(inputData& inputs){
    Setup(inputs);

    std::vector<size_t> UniqueNodes_u = mesh->ElementGroups[ElemGroupIndex_u].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes_u, dofTypes_u);

    std::vector<size_t> UniqueNodes_phase = mesh->ElementGroups[ElemGroupIndex_phase].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes_phase, dofType_phase);

    if (HasT){
        std::vector<size_t> UniqueNodes_T = mesh->ElementGroups[ElemGroupIndex_T].GetUniqueNodes();
        dofs->AddDofs(UniqueNodes_T, dofType_T);   
    }

    size_t NElems = mesh->ElementGroups[ElemGroupIndex_phase].NElems;
    size_t NIP = mesh->ElementGroups[ElemGroupIndex_phase].BaseElem->ipcount;
    Energy_Hist.resize(NElems);
    Energy_HistOld.resize(NElems);
    Plastic_Strain.resize(NElems);
    Plastic_StrainOld.resize(NElems);
    ViscDiss.resize(NElems);
    ViscDissOld.resize(NElems);
    DamViscDiss.resize(NElems);
    DamViscDissOld.resize(NElems);
    for (size_t i = 0; i < NElems; i++){
        Energy_Hist[i].resize(NIP);
        Energy_HistOld[i].resize(NIP);
        Plastic_Strain[i].resize(NIP);
        Plastic_StrainOld[i].resize(NIP);
        ViscDiss[i].resize(NIP);
        ViscDissOld[i].resize(NIP);
        DamViscDiss[i].resize(NIP);
        DamViscDissOld[i].resize(NIP);
        for (size_t j = 0; j < NIP; j++){
            Energy_Hist[i][j] = 0.0;
            Energy_HistOld[i][j] = 0.0;
            ViscDiss[i][j] = 0.0;
            ViscDissOld[i][j] = 0.0;
            DamViscDiss[i][j] = 0.0;
            DamViscDissOld[i][j] = 0.0;
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

    if (PF_Util->IrreversibleMethod==PF_Util->LagrangeMultiplier || PF_Util->IrreversibleMethod==PF_Util->AugLagrangeMultiplier){
        dofs->AddDofs(UniqueNodes_phase, dofType_lMult);
    }
    if (PF_Util->IrreversibleMethod==PF_Util->AugLagrangeMultiplier2 || PF_Util->IrreversibleMethod==PF_Util->AugmentedLagrangeMult3){
        dofs->AddDofs(UniqueNodes_phase, dofType_lMult2);
    }
}

/// @brief Loads model from a previous restart savefile
/// @param inputs input: JSON object representing the input file
/// @param data  input: restart data from where to load variables
void PhaseFieldDamage::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);
    data.Load("Energy", Energy_Hist);
    //data.Load("PlasticStrains", Plastic_Strain);
    Commit(CommitTypes::TIMEDEP_COMMIT_TYPE);
}

/// @brief Save model to a restart file
/// @param data input: pointer to object in which to save restart data
void PhaseFieldDamage::save(SaveDataFile& data){
    data.Save("Energy", Energy_Hist);
    //data.Save("PlasticStrains", Plastic_Strain);
}

/// @brief Commits time dependent history variables (Energy history and plastic strains)
/// @param CommitType
void PhaseFieldDamage::Commit(int CommitType){
    if (CommitType==CommitTypes::TIMEDEP_COMMIT_TYPE){
        double LCrackRed = LCrack;
        MPI_Allreduce(MPI_IN_PLACE, &LCrackRed, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        LCrackOld = LCrackRed;
        
        size_t NElems = mesh->ElementGroups[ElemGroupIndex_phase].NElems;
        size_t NIP = mesh->ElementGroups[ElemGroupIndex_phase].BaseElem->ipcount;
        for (size_t i = 0; i < NElems; i++){
            for (size_t j = 0; j < NIP; j++){
                Energy_HistOld[i][j] = Energy_Hist[i][j];
                Plastic_StrainOld[i][j] = 1.0*Plastic_Strain[i][j];
                ViscDissOld[i][j] = ViscDiss[i][j];
                DamViscDissOld[i][j] = DamViscDiss[i][j];
                if (M->HistSize>0){
                    for (size_t k = 0; k < M->HistSize; k++){
                        MatHistOld[i][j][k] = MatHist[i][j][k];
                    }
                };
            }
        }
    }
}

void PhaseFieldDamage::ResetStep(){
    size_t NElems = mesh->ElementGroups[ElemGroupIndex_phase].NElems;
    size_t NIP = mesh->ElementGroups[ElemGroupIndex_phase].BaseElem->ipcount;
    for (size_t i = 0; i < NElems; i++){
        for (size_t j = 0; j < NIP; j++){
            Energy_Hist[i][j] = Energy_HistOld[i][j];
            Plastic_Strain[i][j] = 1.0*Plastic_StrainOld[i][j];
            ViscDiss[i][j] = ViscDissOld[i][j];
            DamViscDiss[i][j] = DamViscDissOld[i][j];
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
void PhaseFieldDamage::Setup(inputData& inputs){
    if (dim==3){
        DofNames_u.push_back("uz");
        g = {0, 0, -9.81};
    }

    //get indices for relevant element groups
    std::string ElemGroupName;
    inputs.GetRequired(ElemGroupName, {"Models", MyName, "ElementGroup_phase"});
    ElemGroupIndex_phase = mesh->GetElementGroupIdx(ElemGroupName);

    inputs.GetRequired(ElemGroupName, {"Models", MyName, "ElementGroup_u"});
    ElemGroupIndex_u = mesh->GetElementGroupIdx(ElemGroupName);

    //get relevant degree of freedmo steps and indices
    std::vector<size_t> dofSteps;
    dofs->getDofTypesSteps(DofNames_u, dofTypes_u, dofSteps);
    if (dofSteps[0] != dofSteps[1]) throw std::invalid_argument(ModelName+" requires "+ DofNames_u[0] + " and " + DofNames_u[1] + " to be in the same solver step\n");
    Step_u = dofSteps[0];

    dofs->getDofTypesSteps(DofNames_phase, dofType_phase, Step_phase);

    //material parameters
    inputs.GetRequired(SolidMatName, {"Models", MyName, "Material"});
    M = new IceMaterial(inputs, SolidMatName);

    inputs.GetRequired(viscNum, {"Models", MyName, "ViscousStab"});

    IncludeInertia = false;
    inputs.GetOptional(IncludeInertia, {"Models", MyName, "Intertia"});

    ConsistentDegradation = false;
    inputs.GetOptional(ConsistentDegradation, {"properties", "PhaseField", "ConsistentDegradation"});

    LumpedCap = false;
    inputs.GetOptional(LumpedCap, {"Models", MyName, "LumpedPFCapacity"});

    LumpedDrive = false;
    inputs.GetOptional(LumpedDrive, {"Models", MyName, "LumpedDriving"});

    bool HasGravity = false;
    inputs.GetOptional(HasGravity, {"Models", MyName, "Gravity"});
    if (HasGravity == false){
        g.setZero();
    }

    HasT = false;
    if (physics->dofspace->hasDofType(DofName_T)){
        HasT = true;
    }

    ViscousHeating = false;
    Step_T = 99;
    inputs.GetOptional(ViscousHeating, {"Models", MyName, "ViscousHeating"});

    if (HasT){
        dofs->getDofTypesSteps(DofName_T, dofType_T, Step_T);
        inputs.GetRequired(ElemGroupName, {"Models", MyName, "ElementGroup_T"});
        ElemGroupIndex_T = mesh->GetElementGroupIdx(ElemGroupName);
    } else {
        ElemGroupIndex_T = ElemGroupIndex_phase;
    }

    //damage functions
    PF_Util = new PhaseFieldUtility(inputs, SolidMatName, MyName);

    std::string DegFunctionName;
    inputs.GetRequired(DegFunctionName, {"properties", "PhaseField", "DamFunction"});
    if (DamageMethodNames.count(DegFunctionName)){
        DegFunction = DamageMethodNames[DegFunctionName];
    } else {
        throw std::invalid_argument("Damage function type "+DegFunctionName+" not defined,");
    }

    std::string DensityDegFunctionName;
    inputs.GetRequired(DensityDegFunctionName, {"properties", "PhaseField", "DensityDamage"});
    if (DamageMethodNames.count(DensityDegFunctionName)){
        DensityDegFunction = DamageMethodNames[DensityDegFunctionName];
    } else {
        throw std::invalid_argument("Damage function type "+DensityDegFunctionName+" not defined,");
    }

    std::string GravityDegFunctionName;
    inputs.GetRequired(GravityDegFunctionName, {"properties", "PhaseField", "DensityDamage"});
    if (DamageMethodNames.count(GravityDegFunctionName)){
        GravityDegFunction = DamageMethodNames[GravityDegFunctionName];
    } else {
        throw std::invalid_argument("Damage function type "+GravityDegFunctionName+" not defined,");
    }

    ViscDamageFact = 0.0;
    inputs.GetOptional(ViscDamageFact, {"properties", "PhaseField", "ViscoDrivingFactor"});
    if (ViscDamageFact>1.0e-12){
        std::string ViscoDegFunctionName;
        inputs.GetRequired(ViscoDegFunctionName, {"properties", "PhaseField", "ViscoDriving"});
        if (DamageMethodNames.count(ViscoDegFunctionName)){
            ViscoDegFunction = DamageMethodNames[ViscoDegFunctionName];
        } else {
            throw std::invalid_argument("Damage function type "+ViscoDegFunctionName+" not defined,");
        }
    }
    

    if (PF_Util->IrreversibleMethod==PF_Util->LagrangeMultiplier || PF_Util->IrreversibleMethod==PF_Util->AugLagrangeMultiplier){
        size_t Step_phase2;
        dofs->getDofTypesSteps(DofName_LagMult, dofType_lMult, Step_phase2);
        if (Step_phase != Step_phase2) throw std::invalid_argument(ModelName+" requires "+ DofNames_phase + " and " + DofName_LagMult + " to be in the same solver step\n");
    }

    if (PF_Util->IrreversibleMethod==PF_Util->AugLagrangeMultiplier2 || PF_Util->IrreversibleMethod==PF_Util->AugmentedLagrangeMult3){
        std::vector<size_t> Step_phase2;
        dofs->getDofTypesSteps(DofName_LagMult2, dofType_lMult2, Step_phase2);
        if (Step_phase != Step_phase2[0]) throw std::invalid_argument(ModelName+" requires "+ DofNames_phase + " and " + DofName_LagMult + " to be in the same solver step\n");
    }

    LCrack = 0.0;
    LCrackOld = -42;
}

/// @brief Assemble the tangent stiffness matrix and residual force vector
/// @param K Tangent matrix
/// @param f residual force vector
/// @param cons constraints applied by this model [none]
/// @param step input: Step for which to assemble the tangential matrix
void PhaseFieldDamage::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    double t = physics->time;
    double dt = physics->timeScheme->dt;
    // if (t>=0.0 && physics->timeScheme->dt<0.1){
    //     ACreep = 0;
    // }

    if (step == Step_u){    //assemble stiffness matrix for momentum balance
        Assemble_U(f, K);
    }

    if (step == Step_phase){    //phase field evolution equation
        Assemble_Phase(f, K);
    }

    if (step == Step_T && ViscousHeating){
        Assemble_T(f, K);
    }
}

void PhaseFieldDamage::Assemble_U(Vec &f, Mat &K){
    size_t ipcount = mesh->ElementGroups[ElemGroupIndex_phase].BaseElem->ipcount;    // number of integration points
    size_t nNodes_u = mesh->ElementGroups[ElemGroupIndex_u].NNodes_per_elem;         // number of nodes per displacement element
    size_t nNodes_Phase = mesh->ElementGroups[ElemGroupIndex_phase].NNodes_per_elem; // number of nodes per phasefield element
    size_t nNodes_T = mesh->ElementGroups[ElemGroupIndex_T].NNodes_per_elem; // number of nodes per phasefield element

    std::vector<double> w(ipcount); // integration weights
    std::vector<Eigen::RowVectorXd> Nu(ipcount); for (size_t ip = 0; ip < ipcount; ip++) Nu[ip].resize(nNodes_u); // displacement shape functions
    std::vector<Eigen::MatrixXd> Gu(ipcount); for (size_t ip = 0; ip < ipcount; ip++) Gu[ip].resize(dim, nNodes_u);        // displacement shape function gradients
    Eigen::MatrixXd B(6, dim * nNodes_u);    // stress to strain mapping matrix
    Eigen::MatrixXd N_uu(3, dim * nNodes_u); // state to displacement mapping matrix

    std::vector<Eigen::RowVectorXd> Nph(ipcount); for (size_t ip = 0; ip < ipcount; ip++) Nph[ip].resize(nNodes_Phase); // phasefield shape functions
    std::vector<Eigen::MatrixXd> Gph(ipcount); for (size_t ip = 0; ip < ipcount; ip++) Gph[ip].resize(dim, nNodes_Phase); // phasefield shape function gradients

    std::vector<size_t> El_u(nNodes_u), El_Phase(nNodes_Phase), El_T(nNodes_T); // nodes contained within this element

    std::vector<PetscInt> dofsU(dim * nNodes_u);     // displacement degree of freedom indices
    std::vector<PetscInt> dofsPhase(nNodes_Phase); // phase field degree of freedom indices


    Eigen::VectorXd U(dim * nNodes_u), dU(dim * nNodes_u), ddU(dim * nNodes_u), UOld(dim * nNodes_u); // current displacement, velocity, acceleration, and old displacement
    Eigen::VectorXd Phase(nNodes_Phase), PhaseOld(nNodes_Phase);                              // current phasefield parameter and phasefield time derivative
    Vector6d Stress_Dam, Stress_Undam;                                                 // damageable and undamageable stresses
    Vector6d Strain, StrainOld, Strain_Pl, Strain_Pl_Old;                              // total and plastic strain vectors
    Matrix6d D_Dam, D_UnDam, DPlast;                                                   // tangential matrices (damageable, undamageable, plastic)

    std::vector<Eigen::RowVectorXd> NT(ipcount); for (size_t ip = 0; ip < ipcount; ip++) NT[ip].resize(nNodes_T); // phasefield shape functions
    std::vector<Eigen::MatrixXd> GT(ipcount); for (size_t ip = 0; ip < ipcount; ip++) GT[ip].resize(dim, nNodes_T); // phasefield shape function gradients   
    std::vector<PetscInt> dofsT(nNodes_T);
    Eigen::VectorXd T(nNodes_T);

    double ddXdt2 = physics->timeScheme->ddu_dt;                       // state to acceleration derivative
    double dXdt = physics->timeScheme->du_dt;
    double dt= physics->timeScheme->dt;
    double Dam, dDamdPhi, DamI, dDamIdPhi, NotNeeded, DamG, dDamGdPhi; // damage variable and its derivatives

    Eigen::VectorXd F_U(dim * nNodes_u);                // element force vector
    Eigen::MatrixXd K_UU(dim * nNodes_u, dim * nNodes_u); // tangential matrix, displacement/displacement component

    Eigen::VectorXd WLumped(nNodes_u);

    for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex_u].NElems; el++){ // loop over all elements
        // get shape functions and gradients
        mesh->getShapeGrads(ElemGroupIndex_u, el, El_u, Nu, Gu, w);
        mesh->getShapeGrads(ElemGroupIndex_phase, el, El_Phase, Nph, Gph, w);

        // obtain degrees of freedom and state vectors
        dofs->getDofForNodes(El_u, dofTypes_u, dofsU);
        dofs->getDofForNodes(El_Phase, dofType_phase, dofsPhase);

        physics->StateVectors[Step_u].GetValues(dofsU, U);
        physics->StateVectorsOld[Step_u].GetValues(dofsU, UOld);
        physics->dStateVectors[Step_u].GetValues(dofsU, dU);
        physics->ddStateVectors[Step_u].GetValues(dofsU, ddU);
        physics->StateVectors[Step_phase].GetValues(dofsPhase, Phase);
        physics->StateVectorsOld[Step_phase].GetValues(dofsPhase, PhaseOld);

        if (HasT){
            mesh->getShapeGrads(ElemGroupIndex_T, el, El_T, NT, GT, w);
            dofs->getDofForNodes(El_T, dofType_T, dofsT);
            physics->StateVectors[Step_T].GetValues(dofsT, T);
        }

        // zero matrices before assembly
        F_U.setZero();
        K_UU.setZero();

        WLumped.setZero();

        // integration over element
        for (size_t ip = 0; ip < ipcount; ip++){
            M->SetTemperature(NT[ip]*T, 1.0);

            double phi = Nph[ip] * Phase; // phase field variable within integration point
            double phi_Old = Nph[ip] * PhaseOld;
            GetDamFunc(phi, Dam, dDamdPhi, NotNeeded);          // get damage function based on phasefield
            GetDamFuncInertia(phi, DamI, dDamIdPhi, NotNeeded); // get inertia damage function
            GetDamFuncGravity(phi, DamG, dDamGdPhi, NotNeeded);

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
                M->UpdatePlasticStrains(DPlast, Strain_Pl, Strain, Strain_Pl_Old, StrainOld, phi, phi_Old, dt, MatHist[el][ip], MatHistOld[el][ip], PlasticDiss, DamPlasticDiss);
            } else {
                std::vector<double> Unused, Unused2;
                M->UpdatePlasticStrains(DPlast, Strain_Pl, Strain, Strain_Pl_Old, StrainOld, phi, phi_Old, dt, Unused, Unused2, PlasticDiss, DamPlasticDiss);
            }
            Plastic_Strain[el][ip] = Strain_Pl;
            ViscDiss[el][ip] = ViscDissOld[el][ip] + PlasticDiss*dt;
            DamViscDiss[el][ip] = DamViscDissOld[el][ip] + DamPlasticDiss*dt;

            Strain = Strain - Strain_Pl;
            M->EnergySplit(D_Dam, D_UnDam, Stress_Dam, Stress_Undam, Strain);

            F_U += w[ip] * B.transpose() * (Dam * Stress_Dam + Stress_Undam);
            K_UU += w[ip] * B.transpose() * (Dam * D_Dam + D_UnDam) * DPlast * B;

            // gravity
            //F_U += -w[ip] * DamG * rho * N_uu.transpose() * g;

            //inertia
            if (IncludeInertia){
            //    F_U  += w[ip] * DamI * rho * N_uu.transpose()*N_uu*ddU;
            //    K_UU += w[ip] * DamI * rho * N_uu.transpose()*N_uu*ddXdt2;
            }

            // Lumped weights
            WLumped += w[ip] * Nu[ip];// * DamI;
        }
        for (size_t n = 0; n < nNodes_u; n++){
            int ny = n + nNodes_u;
            int nz = n + 2*nNodes_u;

            // inertia
            if (IncludeInertia){
                GetDamFuncInertia(Phase(n), DamI, dDamIdPhi, NotNeeded); // get inertia damage function

                F_U(n) += WLumped(n) * DamI * M->Density * ddU(n);
                K_UU(n, n) += WLumped(n) * DamI * M->Density * ddXdt2;

                F_U(ny) += WLumped(n) * DamI * M->Density * ddU(ny);
                K_UU(ny, ny) += WLumped(n) * DamI * M->Density * ddXdt2;

                if (dim==3){
                    F_U(nz) += WLumped(n) * DamI * M->Density * ddU(nz);
                    K_UU(nz, nz) += WLumped(n) * DamI * M->Density * ddXdt2;
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
            GetDamFuncGravity(Phase(n), DamG, dDamGdPhi, NotNeeded);
            F_U(n) += -WLumped(n) * DamG * M->Density * g(0);
            F_U(ny)+= -WLumped(n) * DamG * M->Density * g(1);
            if (dim==3){
                F_U(nz)+= -WLumped(n) * DamG * M->Density * g(2);
            }
        }

        // add to total vector and matrix
        VecAdd(f, dofsU, F_U);
        MatAdd(K, dofsU, dofsU, K_UU);
    }
}

void PhaseFieldDamage::Assemble_Phase(Vec &f, Mat &K){
    double t = physics->time;
    double dt = physics->timeScheme->dt;

    double unused, unused2;

    size_t ipcount = mesh->ElementGroups[ElemGroupIndex_phase].BaseElem->ipcount;   //number ofintegration points
    size_t nNodes_u = mesh->ElementGroups[ElemGroupIndex_u].NNodes_per_elem;        //number of nodes per displacement element
    size_t nNodes_Phase = mesh->ElementGroups[ElemGroupIndex_phase].NNodes_per_elem;// number of nodes per phasefield element
    size_t nNodes_T = mesh->ElementGroups[ElemGroupIndex_T].NNodes_per_elem;// number of nodes per phasefield element

    std::vector<double> w(ipcount); //integration weights

    std::vector<Eigen::RowVectorXd> Nu(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Nu[ip].resize(nNodes_u);  //displacement shape functions
    std::vector<Eigen::MatrixXd> Gu(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gu[ip].resize(dim, nNodes_u);  //displacement shape function gradients
    Eigen::MatrixXd B(6, dim*nNodes_u);   //stress to strain mapping matrix
    Eigen::MatrixXd N_uu(dim,dim*nNodes_u); //state to displacement mapping matrix

    std::vector<Eigen::RowVectorXd> Nph(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Nph[ip].resize(nNodes_Phase); //phasefield shape functions
    std::vector<Eigen::MatrixXd> Gph(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gph[ip].resize(dim, nNodes_Phase); //phasefield shape function gradients
    std::vector<Eigen::MatrixXd> G2ph(ipcount); for (size_t ip=0; ip<ipcount;  ip++) G2ph[ip].resize(3*(dim-1), nNodes_Phase);
    Eigen::VectorXd WLumped(nNodes_Phase), LumpedDriving(nNodes_Phase);

    std::vector<size_t> El_u(nNodes_u), El_Phase(nNodes_Phase), El_T(nNodes_T); //node indices for element

    std::vector<PetscInt> dofsU(dim*nNodes_u), dofsPhase(nNodes_Phase), dofsLag(nNodes_Phase), dofsLag2(nNodes_Phase); //degree of freedom indices for displacent and phase field

    Eigen::VectorXd U(dim*nNodes_u), dU(dim*nNodes_u); //current displacement, velocity
    Eigen::VectorXd Phase(nNodes_Phase), dPhase(nNodes_Phase), PhaseOld(nNodes_Phase), L(nNodes_Phase), L2(nNodes_Phase);  //current phasefield and change over time

    std::vector<Eigen::RowVectorXd> NT(ipcount); for (size_t ip = 0; ip < ipcount; ip++) NT[ip].resize(nNodes_T); // phasefield shape functions
    std::vector<Eigen::MatrixXd> GT(ipcount); for (size_t ip = 0; ip < ipcount; ip++) GT[ip].resize(dim, nNodes_T); // phasefield shape function gradients   
    std::vector<PetscInt> dofsT(nNodes_T);
    Eigen::VectorXd T(nNodes_T);

    double dXdt = physics->timeScheme->du_dt;   //state to velocity derivative

    double Dam, dDamdPhi, DamI, dDamIdPhi, ddDamdPhi2, ddDamIdPhi2, DamG, dDamGdPhi, ddDamGdPhi2; //damage functions and derivatives
    Vector6d Stress_Dam, Stress_Undam, Strain;   //damageable undamageable stresses, elastic strain
    Matrix6d D_Dam, D_UnDam; //damageable and undamageable material stiffness matrices

    Eigen::VectorXd Dist(nNodes_Phase), DistCap(nNodes_Phase);
    Eigen::MatrixXd dDistdPhi(nNodes_Phase,nNodes_Phase), dDistCap(nNodes_Phase, nNodes_Phase);

    Eigen::VectorXd F_Ph(nNodes_Phase); //phasefield residual vector
    Eigen::VectorXd F_L(nNodes_Phase); //Lagrange Multiplier residual vector
    Eigen::MatrixXd K_PhPh(nNodes_Phase, nNodes_Phase); //tangent matrix, phasefield/phasefield component
    Eigen::MatrixXd K_LPh(nNodes_Phase, nNodes_Phase); //tangent matrix, Lagrange/phasefield component
    Eigen::MatrixXd K_PhL(nNodes_Phase, nNodes_Phase); //tangent matrix, phasefield/Lagrange component
    Eigen::MatrixXd K_LL(nNodes_Phase, nNodes_Phase); //tangent matrix, Lagrange/Lagrange component
    Eigen::VectorXd F_L2(nNodes_Phase); //Lagrange Multiplier residual vector
    Eigen::MatrixXd K_L2Ph(nNodes_Phase, nNodes_Phase); //tangent matrix, Lagrange/phasefield component
    Eigen::MatrixXd K_PhL2(nNodes_Phase, nNodes_Phase); //tangent matrix, phasefield/Lagrange component
    Eigen::MatrixXd K_L2L2(nNodes_Phase, nNodes_Phase); //tangent matrix, Lagrange/Lagrange component

    LCrack = 0.0;
    for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex_u].NElems; el++){   //loop over all elements

        //get shape functions and gradients
        mesh->getShapeGrads(ElemGroupIndex_u, el, El_u, Nu, Gu, w);

        if (PF_Util->NeedsGrad2 == false){
            mesh->getShapeGrads(ElemGroupIndex_phase, el, El_Phase, Nph, Gph, w);
        } else {
            mesh->getShapeGrads(ElemGroupIndex_phase, el, El_Phase, Nph, Gph, G2ph, w);
        }

        //get degrees of freedom indices and state vectors
        dofs->getDofForNodes(El_u, dofTypes_u, dofsU);
        dofs->getDofForNodes(El_Phase, dofType_phase, dofsPhase);

        physics->StateVectors[Step_u].GetValues(dofsU, U);
        physics->dStateVectors[Step_u].GetValues(dofsU, dU);
        physics->StateVectors[Step_phase].GetValues(dofsPhase, Phase);
        physics->dStateVectors[Step_phase].GetValues(dofsPhase, dPhase);
        physics->StateVectorsOld[Step_phase].GetValues(dofsPhase, PhaseOld);

        if (HasT){
            mesh->getShapeGrads(ElemGroupIndex_T, el, El_T, NT, GT, w);
            physics->StateVectors[Step_T].GetValues(dofsT, T);
        }

        if (PF_Util->IrreversibleMethod==PF_Util->LagrangeMultiplier || PF_Util->IrreversibleMethod==PF_Util->AugLagrangeMultiplier){
            dofs->getDofForNodes(El_Phase, dofType_lMult, dofsLag);
            physics->StateVectors[Step_phase].GetValues(dofsLag, L);
        }
        if (PF_Util->IrreversibleMethod==PF_Util->AugLagrangeMultiplier2 || PF_Util->IrreversibleMethod==PF_Util->AugmentedLagrangeMult3){
            dofs->getDofForNodes(El_Phase, dofType_lMult2[0], dofsLag);
            dofs->getDofForNodes(El_Phase, dofType_lMult2[1], dofsLag2);
            physics->StateVectors[Step_phase].GetValues(dofsLag, L);
            physics->StateVectors[Step_phase].GetValues(dofsLag2, L2);
        }

        //zero before assembly
        F_Ph.setZero();
        K_PhPh.setZero();
        WLumped.setZero();
        LumpedDriving.setZero();
        if (PF_Util->IrreversibleMethod==PF_Util->LagrangeMultiplier || PF_Util->IrreversibleMethod==PF_Util->AugLagrangeMultiplier 
            || PF_Util->IrreversibleMethod==PF_Util->AugLagrangeMultiplier2 || PF_Util->IrreversibleMethod==PF_Util->AugmentedLagrangeMult3){
            F_L.setZero();
            K_LPh.setZero();
            K_PhL.setZero();
            K_LL.setZero();
            if (PF_Util->IrreversibleMethod==PF_Util->AugLagrangeMultiplier2 || PF_Util->IrreversibleMethod==PF_Util->AugmentedLagrangeMult3){
                F_L2.setZero();
                K_L2Ph.setZero();
                K_PhL2.setZero();
                K_L2L2.setZero();
            }
        }

        for (size_t ip = 0; ip < ipcount; ip++){ //integrate over element
            M->SetTemperature(NT[ip]*T, 1.0);

            double phi  = Nph[ip]*Phase;    //local phasefield
            GetDamFunc(phi, Dam, dDamdPhi, ddDamdPhi2);
            GetDamFuncInertia(phi, DamI, dDamIdPhi, ddDamIdPhi2);
            GetDamFuncGravity(phi, DamG, dDamGdPhi, ddDamGdPhi2);

            //special matrices
            ConstructBMat(B, Gu[ip], dim);

            N_uu.setZero();
            N_uu(0,Eigen::seq(0,nNodes_u-1)) = Nu[ip];
            N_uu(1,Eigen::seq(nNodes_u,2*nNodes_u-1)) = Nu[ip];
            if (dim==3){
                N_uu(2,Eigen::seq(2*nNodes_u,3*nNodes_u-1)) = Nu[ip];
            }

            //Energy split
            Strain = B*U - Plastic_Strain[el][ip];
            double Energy = 0;
            Energy = M->EnergySplit(D_Dam, D_UnDam, Stress_Dam, Stress_Undam, Strain);

            //allow for initialization period
            //if (t<0.0){
            //    Energy = E_min;
            //}
            //slow loading
        //    double tDecay = 20*0.25e-2;
        //    if (t<=tDecay){
        //       Energy = Energy/((1.0+50.0*((tDecay-t)/tDecay)));
        //    }

            //History Field
            if (PF_Util->IrreversibleMethod==PF_Util->HistParameter){
                if (Energy_HistOld[el][ip]>Energy){
                    Energy_Hist[el][ip] = Energy_HistOld[el][ip];
                    Energy = Energy_HistOld[el][ip];
                } else {
                    Energy_Hist[el][ip] = Energy;
                }
            } else {
                Energy_Hist[el][ip] = Energy;
            }

            //crack resistance
            double gamma = PF_Util->GetDistributionFunction(Dist, dDistdPhi, Nph[ip], Gph[ip], G2ph[ip], Phase);
            if (LumpedCap){
                PF_Util->GetCapacityOnly(DistCap, dDistCap, unused, unused2, Nph[ip], phi);
            } else {
                DistCap.setZero();
                dDistCap.setZero();
            }
            LCrack += gamma*w[ip];

            F_Ph   += w[ip]*M->Gc*(Dist-DistCap);
            K_PhPh += w[ip]*M->Gc*(dDistdPhi-dDistCap);

            //crack driving force (assuming staggered scheme!!!!)
            if (LumpedDrive==false){
                F_Ph   += w[ip]*dDamdPhi*Nph[ip].transpose()*Energy;
                K_PhPh += w[ip]*ddDamdPhi2*Nph[ip].transpose()*Nph[ip]*Energy;
            } else {
                LumpedDriving += w[ip]*Nph[ip]*Energy;
            }

            //plastic dissipation driving force
            if (ViscDamageFact>1.0e-9){
                double DamV, dDamV, ddDamV;
                double PlasticDiss = ViscDamageFact*DamViscDiss[el][ip];

                GetDamFuncVisc(phi, DamV, dDamV, ddDamV);

                F_Ph   += w[ip]*dDamV*Nph[ip].transpose()*PlasticDiss;
                K_PhPh += w[ip]*ddDamV*Nph[ip].transpose()*Nph[ip]*PlasticDiss;
            }

            //Inertial driving force
            if (IncludeInertia && ConsistentDegradation){
                F_Ph  += -w[ip]*0.5*M->Density*((N_uu*dU).transpose()*(N_uu*dU))*Nph[ip].transpose()*dDamIdPhi;
                K_PhPh+= -w[ip]*0.5*M->Density*((N_uu*dU).transpose()*(N_uu*dU))*Nph[ip].transpose()*Nph[ip]*ddDamIdPhi2;
            }

            //Gravitational Driving force
            if (ConsistentDegradation){
                F_Ph   += -w[ip]*M->Density*(g.transpose()*N_uu*U)*Nph[ip].transpose()*dDamGdPhi;
                K_PhPh += -w[ip]*M->Density*(g.transpose()*N_uu*U)*Nph[ip].transpose()*Nph[ip]*ddDamGdPhi2;
            }

            //if (IrreversibleMethod==AugmentedLagrangeMult3){
            //    AugLagrIrr3(Phase, PhaseOld, L, L2, F_L, F_L2, F_Ph, K_PhL, K_PhL2, K_LPh, K_L2Ph, K_PhPh, K_LL, K_L2L2, w[ip], Nph[ip], Gph[ip], G2ph[ip]);
            //}

            //lumped integration parameter
            WLumped += w[ip]*Nph[ip];
        }
        for (size_t n = 0; n < nNodes_Phase; n++){
            double LDistCap, dLDistCap;
            // Driving force
            if (LumpedDrive){
                GetDamFunc(Phase(n), Dam, dDamdPhi, ddDamdPhi2);
                F_Ph(n)     += LumpedDriving(n)*dDamdPhi;
                K_PhPh(n,n) += LumpedDriving(n)*ddDamdPhi2;
            }

            if (LumpedCap){
                PF_Util->GetCapacityOnly(DistCap, dDistCap, LDistCap, dLDistCap, Nph[0], Phase(n));
                F_Ph(n)     += WLumped(n)*M->Gc*LDistCap;
                K_PhPh(n,n) += WLumped(n)*M->Gc*dLDistCap;
            }

            //viscosity term
            F_Ph(n)    += WLumped(n)*M->pf_visc/dt*(Phase(n)-PhaseOld(n));
            K_PhPh(n,n)+= WLumped(n)*(M->pf_visc/dt+viscNum*M->Gc/M->pf_l);
        }
        PF_Util->EnforceIrreversible(Phase, PhaseOld, WLumped, L, L2, F_Ph, F_L, F_L2, K_PhPh, K_LL, K_L2L2, K_PhL, K_LPh, K_PhL2, K_L2Ph);

        //add to total force and tangential matrix
        VecAdd(f, dofsPhase, F_Ph);
        MatAdd(K, dofsPhase, dofsPhase, K_PhPh);
        if (PF_Util->IrreversibleMethod==PF_Util->LagrangeMultiplier || PF_Util->IrreversibleMethod==PF_Util->AugLagrangeMultiplier 
            || PF_Util->IrreversibleMethod==PF_Util->AugLagrangeMultiplier2 || PF_Util->IrreversibleMethod==PF_Util->AugmentedLagrangeMult3){
            VecAdd(f, dofsLag, F_L);
            MatAdd(K, dofsLag,   dofsPhase, K_LPh);
            MatAdd(K, dofsPhase, dofsLag,   K_PhL);
            MatAdd(K, dofsLag,   dofsLag,   K_LL);
            if (PF_Util->IrreversibleMethod==PF_Util->AugLagrangeMultiplier2 || PF_Util->IrreversibleMethod==PF_Util->AugmentedLagrangeMult3){
                VecAdd(f, dofsLag2, F_L2);
                MatAdd(K, dofsLag2,   dofsPhase, K_L2Ph);
                MatAdd(K, dofsPhase, dofsLag2,   K_PhL2);
                MatAdd(K, dofsLag2,   dofsLag2,   K_L2L2);
            }
        }
    }
}

void PhaseFieldDamage::Assemble_T(Vec &f, Mat &K){
    size_t ipcount = mesh->ElementGroups[ElemGroupIndex_T].BaseElem->ipcount;    // number of integration points
    size_t nNodes_T = mesh->ElementGroups[ElemGroupIndex_T].NNodes_per_elem; // number of nodes per phasefield element

    std::vector<double> w(ipcount); // integration weights
    std::vector<Eigen::RowVectorXd> N(ipcount); for (size_t ip = 0; ip < ipcount; ip++) N[ip].resize(nNodes_T); // displacement shape functions
    std::vector<Eigen::MatrixXd> G(ipcount); for (size_t ip = 0; ip < ipcount; ip++) G[ip].resize(dim, nNodes_T);        // displacement shape function gradients

    std::vector<size_t> El(nNodes_T); // nodes contained within this element
    std::vector<PetscInt> dofsT(nNodes_T);     // displacement degree of freedom indices

    double dt= physics->timeScheme->dt;

    Eigen::VectorXd F_T(nNodes_T);                // element force vector

    for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex_T].NElems; el++){ // loop over all elements
        // get shape functions and gradients
        mesh->getShapeGrads(ElemGroupIndex_T, el, El, N, G, w);

        // obtain degrees of freedom and state vectors
        dofs->getDofForNodes(El, dofType_T, dofsT);

        // zero matrices before assembly
        F_T.setZero();

        // integration over element
        for (size_t ip = 0; ip < ipcount; ip++){
            double HeatingRate = (ViscDiss[el][ip] - ViscDissOld[el][ip])/dt;
            //if (HeatingRate>0.0) std::cout << HeatingRate << "\n";
            F_T += -w[ip] * N[ip].transpose() * HeatingRate;
        }

        // add to total vector and matrix
        VecAdd(f, dofsT, F_T);
    }
}

/// @brief Gets the stress damage function and derivatives based on the set function type and input phasefield
/// @param phi input: phasefield parameter
/// @param D output: Damage
/// @param dD output: d(Damage)/d(phi)
/// @param ddD output: d2(Damage)/d(phi2)
void PhaseFieldDamage::GetDamFunc(double phi, double& D, double& dD, double& ddD){
    PF_Util->GenericDamageFunction(phi, DegFunction, D, dD, ddD);
}

void PhaseFieldDamage::GetDamFuncGravity(double phi, double& D, double& dD, double& ddD){
    PF_Util->GenericDamageFunction(phi, GravityDegFunction, D, dD, ddD);
}

/// @brief Gets the inertia damage function and derivatives based on the set function type and input phasefield
/// @param phi input: phasefield parameter
/// @param D output: Damage
/// @param dD output: d(Damage)/d(phi)
/// @param ddD output: d2(Damage)/d(phi2)
void PhaseFieldDamage::GetDamFuncInertia(double phi, double& D, double& dD, double& ddD){
    PF_Util->GenericDamageFunction(phi, DensityDegFunction, D, dD, ddD);
}

void PhaseFieldDamage::GetDamFuncVisc(double phi, double& D, double& dD, double& ddD){
    PF_Util->GenericDamageFunction(phi, ViscoDegFunction, D, dD, ddD);
}

/// @brief Save integration-point specific values to output files (for later post-processing)
/// @param SaveLoc input: location of the variables to save (ip/nodes)
/// @param DataName input: name of parameters to save (not necesarily coming from this model)
/// @param ElemGroup input: element group to save outputs for
/// @param Data in/outputs: object to save data into
/// @return output: has this saved any data?
bool PhaseFieldDamage::SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>>& Data){

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
    if (SaveLoc == "ip" && ElemGroup == ElemGroupIndex_u && (DataName == "sxx" || DataName == "syy" || DataName == "szz" || DataName == "sxy" || DataName == "syz" || DataName == "sxz" || DataName == "s1" || DataName == "s2" || DataName == "s3" || DataName == "dam")){
        SaveStressComponent(DataName, Data, true);
        return true;
    }

    if (SaveLoc == "ip" && ElemGroup == ElemGroupIndex_u && (DataName == "sxx0" || DataName == "syy0" || DataName == "szz0" || DataName == "sxy0" || DataName == "syz0" || DataName == "sxz0")){
        SaveStressComponent(DataName, Data, false);
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
            for (size_t ip = 0; ip < mesh->ElementGroups[ElemGroupIndex_phase].BaseElem->ipcount; ip++){
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

    //save phasefield energy histroy field
    if (SaveLoc == "ip" && ElemGroup == ElemGroupIndex_u && DataName == "phi_hist"){
        for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex_u].NElems; el++){
            for (size_t ip = 0; ip < mesh->ElementGroups[ElemGroupIndex_phase].BaseElem->ipcount; ip++){
                Data[el][ip] = Energy_Hist[el][ip];
            }
        }
        return true;
    }
    return false;
}


void PhaseFieldDamage::SaveStrainComponent(std::string DataName, std::vector<std::vector<double>> &Data){
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

    size_t ipcount = mesh->ElementGroups[ElemGroupIndex_phase].BaseElem->ipcount;
    size_t nNodes_u = mesh->ElementGroups[ElemGroupIndex_u].NNodes_per_elem;

    std::vector<double> w(ipcount);

    std::vector<Eigen::RowVectorXd> Nu(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Nu[ip].resize(nNodes_u);
    std::vector<Eigen::MatrixXd> Gu(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gu[ip].resize(dim, nNodes_u);
    Eigen::MatrixXd B(6, dim*nNodes_u);

    std::vector<size_t> El_u(nNodes_u);

    std::vector<PetscInt> dofsU(dim*nNodes_u);

    Eigen::VectorXd U(dim*nNodes_u);
    Vector6d Strain;

    double Dam, dDamdPhi, DamI, dDamIdPhi, NotNeeded;

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
void PhaseFieldDamage::SaveStressComponent(std::string DataName, std::vector<std::vector<double>>& Data, bool Damaged){
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

    size_t ipcount = mesh->ElementGroups[ElemGroupIndex_phase].BaseElem->ipcount;
    size_t nNodes_u = mesh->ElementGroups[ElemGroupIndex_u].NNodes_per_elem;
    size_t nNodes_Phase = mesh->ElementGroups[ElemGroupIndex_phase].NNodes_per_elem;

    std::vector<double> w(ipcount);

    std::vector<Eigen::RowVectorXd> Nu(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Nu[ip].resize(nNodes_u);
    std::vector<Eigen::MatrixXd> Gu(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gu[ip].resize(dim, nNodes_u);
    Eigen::MatrixXd B(6, dim*nNodes_u);

    std::vector<Eigen::RowVectorXd> Nph(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Nph[ip].resize(nNodes_Phase);
    std::vector<Eigen::MatrixXd> Gph(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gph[ip].resize(dim, nNodes_Phase);

    std::vector<size_t> El_u(nNodes_u), El_Phase(nNodes_Phase);

    std::vector<PetscInt> dofsU(dim*nNodes_u), dofsPhase(nNodes_Phase);

    Eigen::VectorXd U(dim*nNodes_u);
    Eigen::VectorXd Phase(nNodes_Phase);
    Vector6d Stress_Dam, Stress_Undam, Strain;
    Matrix6d D_Dam, D_UnDam;
    Eigen::Matrix3d Stress_Total;
    Eigen::Vector3d PriStresses;

    double Dam, dDamdPhi, DamI, dDamIdPhi, NotNeeded;

    for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex_u].NElems; el++){
        mesh->getShapeGrads(ElemGroupIndex_u, el, El_u, Nu, Gu, w);
        mesh->getShapeGrads(ElemGroupIndex_phase, el, El_Phase, Nph, Gph, w);

        dofs->getDofForNodes(El_u, dofTypes_u, dofsU);
        dofs->getDofForNodes(El_Phase, dofType_phase, dofsPhase);

        physics->StateVectors[Step_u].GetValues(dofsU, U);
        physics->StateVectors[Step_phase].GetValues(dofsPhase, Phase);

        for (size_t ip = 0; ip < ipcount; ip++){
            double phi  = Nph[ip]*Phase;
            GetDamFunc(phi, Dam, dDamdPhi, NotNeeded);
            GetDamFuncInertia(phi, DamI, dDamIdPhi, NotNeeded);

            //special matrices
            ConstructBMat(B, Gu[ip], dim);

            //stiffness
            Strain = B*U-Plastic_Strain[el][ip];
            M->EnergySplit(D_Dam, D_UnDam, Stress_Dam, Stress_Undam, Strain);

            if (DataName == "dam"){
                Data[el][ip] = Dam;
            } else {
                if (Pri==false){
                    if (Damaged){
                        Data[el][ip] = Dam*Stress_Dam(Comp) + Stress_Undam(Comp);
                    } else {
                        Data[el][ip] = Stress_Dam(Comp) + Stress_Undam(Comp);
                    }
                } else {
                    Stress_Total.setZero();
                    Stress_Total(0,0) = Dam*Stress_Dam(0) + Stress_Undam(0);
                    Stress_Total(1,1) = Dam*Stress_Dam(1) + Stress_Undam(1);
                    Stress_Total(2,2) = Dam*Stress_Dam(2) + Stress_Undam(2);
                    Stress_Total(0,1) = Dam*Stress_Dam(3) + Stress_Undam(3);
                    Stress_Total(1,0) = Dam*Stress_Dam(3) + Stress_Undam(3);
                    Stress_Total(1,2) = Dam*Stress_Dam(4) + Stress_Undam(4);
                    Stress_Total(2,1) = Dam*Stress_Dam(4) + Stress_Undam(4);
                    Stress_Total(0,2) = Dam*Stress_Dam(5) + Stress_Undam(5);
                    Stress_Total(2,0) = Dam*Stress_Dam(5) + Stress_Undam(5);

                    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(Stress_Total);
                    PriStresses = eigensolver.eigenvalues();
                    Data[el][ip] = PriStresses(Comp);
                }
            }
        }
    }
}

size_t PhaseFieldDamage::hasTimeData(std::vector<std::string>& DataNames){
    size_t nData = 2;
    DataNames.resize(nData);
    DataNames[0] = MyName+"/L_crack";
    DataNames[1] = MyName+"/dL_crack";
    return nData;
}

void PhaseFieldDamage::GetTimeData(std::vector<double>& DataValues){
    double LCrackRed = LCrack;
    MPI_Allreduce(MPI_IN_PLACE, &LCrackRed, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

    size_t nData = 2;
    DataValues.resize(0);
    DataValues.push_back(LCrackRed);
    if (LCrackOld<0.0){
        LCrackOld = LCrack;
    }
    DataValues.push_back(LCrackRed-LCrackOld);
}