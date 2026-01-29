#include "SymFluidInterface.h"
#include "../../../Physics/physics.h"
#include <Eigen/Dense>

/// @brief Registers the model into the model creator
void Register_SymFluidInterfaceModel(){
    ModelNames.push_back("SymFluidInterfaceModel");
    ModelCreators.push_back(New_SymFluidInterfaceModel);
}

/// @brief Creates a new symmetric fluid interface model
/// @param My_Physics input: physics object
/// @param MyNameIn input: name of this model, from where to use input properties
/// @return pointer to the newly created model
BaseModel* New_SymFluidInterfaceModel(Physics& My_Physics, std::string MyNameIn){
    return new SymFluidInterfaceModel(My_Physics, MyNameIn);
}

/// @brief Initializer, forwards to baseModel
/// @param My_Physics //input: physics object, reference copied into model
/// @param MyName   //name of this model, from where to use input properties
SymFluidInterfaceModel::SymFluidInterfaceModel(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "SymFluidInterfaceModel";
}   

SymFluidInterfaceModel::~SymFluidInterfaceModel(){

}

/// @brief Initialize this model from scratch
/// @param inputs input: JSON object representing the input file
void SymFluidInterfaceModel::init(inputData& inputs){
    Setup(inputs);

    std::vector<size_t> UniqueNodes_u = mesh->ElementGroups[ElemGroupIndex_u].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes_u, dofTypes_u);

    std::vector<size_t> UniqueNodes_p = mesh->ElementGroups[ElemGroupIndex_p].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes_p, dofType_p);

    size_t nNodes = mesh->ElementGroups[ElemGroupIndex_u].BaseElem->NodeCount;
    hMax.resize(mesh->ElementGroups[ElemGroupIndex_u].NElems);
    hMaxOld.resize(mesh->ElementGroups[ElemGroupIndex_u].NElems);
    Broken.resize(mesh->ElementGroups[ElemGroupIndex_u].NElems);
    BrokenOld.resize(mesh->ElementGroups[ElemGroupIndex_u].NElems);
    for (size_t i = 0; i < mesh->ElementGroups[ElemGroupIndex_u].NElems; i++){
        hMax[i].resize(nNodes);
        hMaxOld[i].resize(nNodes);
        Broken[i] = false;
        BrokenOld[i] = false;
        for (size_t j = 0; j < nNodes; j++){
            hMax[i][j] = 0.0;
            hMaxOld[i][j] = 0.0;
        }
    }

    double h0, R0, RBase;
    inputs.GetRequired(h0, {"Models", MyName, "h0"});
    inputs.GetOptional(R0, {"Models", MyName, "r0"});
    inputs.GetOptional(RBase, {"Models", MyName, "rBase"});

    std::vector<size_t> Nodes(nNodes);
    std::vector<double> coordsX(nNodes), coordsY(nNodes), coordsZ(nNodes);
    double maxX, maxY, maxZ;
    maxX = 0.0;
    maxY = 0.0;
    maxZ = 0.0;
    
    for (size_t i = 0; i < mesh->ElementGroups[ElemGroupIndex_u].NElems; i++){
        mesh->GetNodesForElem(Nodes, ElemGroupIndex_u, i);
        mesh->GetCoordsForNodes(coordsX, coordsY, coordsZ, Nodes);
        if (maxX<*std::max_element(coordsX.begin(), coordsX.end())) maxX = *std::max_element(coordsX.begin(), coordsX.end());
        if (maxY<*std::max_element(coordsY.begin(), coordsY.end())) maxY = *std::max_element(coordsY.begin(), coordsY.end());
        if (maxZ<*std::max_element(coordsZ.begin(), coordsZ.end())) maxZ = *std::max_element(coordsZ.begin(), coordsZ.end());
    }

    MPI_Allreduce(MPI_IN_PLACE, &maxX, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &maxY, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &maxZ, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Barrier(PETSC_COMM_WORLD);
    X = maxY;
    Z = maxZ;

    TopCentre_u = 0;
    BotCentre_u = 0;
    size_t nNodesu = mesh->ElementGroups[ElemGroupIndex_u].BaseElem->NodeCount;
    for (size_t i = 0; i < mesh->ElementGroups[ElemGroupIndex_u].NElems; i++){
        mesh->GetNodesForElem(Nodes, ElemGroupIndex_u, i);
        mesh->GetCoordsForNodes(coordsX, coordsY, coordsZ, Nodes);

        double MeanZ = sum(coordsZ) / nNodes;
        double MeanY = sum(coordsY) / nNodes;
        double MeanX = sum(coordsX) / nNodes;
        if (MeanZ > Z-h0){
            Broken[i] = true;
            BrokenOld[i] = true;
            for (size_t j = 0; j < nNodes; j++){
                hMax[i][j] = 1e6;//u0;
                hMaxOld[i][j] = 1e6;//u0;
            }
        } else if (MeanY*MeanY+(MeanZ-(Z-h0))*(MeanZ-(Z-h0)) < R0*R0){
            Broken[i] = true;
            BrokenOld[i] = true;
            for (size_t j = 0; j < nNodes; j++){
                hMax[i][j] = 1e6;//u0;
                hMaxOld[i][j] = 1e6;//u0;
            }
        } else if ((MeanY-maxY)*(MeanY-maxY)+MeanZ*MeanZ < RBase*RBase){
            Broken[i] = true;
            BrokenOld[i] = true;
            for (size_t j = 0; j < nNodes; j++){
                hMax[i][j] = 1e6;//u0;
                hMaxOld[i][j] = 1e6;//u0;
            }
        } else {
            Broken[i] = false;
            BrokenOld[i] = false;
        }

        for (size_t j = 0; j < nNodesu; j++){
            if (std::abs(coordsZ[j]-Z)<1.0e-9 && std::abs(coordsY[j])<1.0e-9 && coordsX[j] < 1.0e-9){
                TopCentre_u = Nodes[j];
            }
            if (std::abs(coordsZ[j])<1.0e-9 && std::abs(coordsY[j]-X)<1.0e-9 && coordsX[j] < 1.0e-9){
                BotCentre_u = Nodes[j];
            }
        }
    }
    // check if any core has a top node
    size_t TopReduced = TopCentre_u;
    size_t BotReduced = BotCentre_u;
    MPI_Allreduce(MPI_IN_PLACE, &TopReduced, 1, MPIU_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &BotReduced, 1, MPIU_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
    if (TopReduced == 0) throw std::invalid_argument("No top centre u node found for SymFluidInterfaceModel "+MyName);
    if (BotReduced == 0) throw std::invalid_argument("No bottom centre u node found for SymFluidInterfaceModel "+MyName);

    size_t nNodesp = mesh->ElementGroups[ElemGroupIndex_p].BaseElem->NodeCount;
    std::vector<size_t> NodesP(nNodesp);
    std::vector<double> coordsXP(nNodesp), coordsYP(nNodesp), coordsZP(nNodesp);
    TopCentre_p = 0;
    BotCentre_p = 0;
    for (size_t i = 0; i < mesh->ElementGroups[ElemGroupIndex_p].NElems; i++){
        mesh->GetNodesForElem(NodesP, ElemGroupIndex_p, i);
        mesh->GetCoordsForNodes(coordsXP, coordsYP, coordsZP, NodesP);

        for (size_t j = 0; j < nNodesp; j++){
            if (std::abs(coordsZP[j]-Z)<1.0e-9 && coordsYP[j] < R0 && coordsXP[j] < 1.0e-9){
                TopNodes.push_back(NodesP[j]);
            }
            if (std::abs(coordsZP[j])<1.0e-9 && coordsXP[j] < 1.0e-9){
                BotNodes.push_back(NodesP[j]);
            }
            if (std::abs(coordsZP[j]-Z)<1.0e-9 && std::abs(coordsYP[j])<1.0e-9 && coordsXP[j] < 1.0e-9){
                TopCentre_p = NodesP[j];
            }
            if (std::abs(coordsZP[j])<1.0e-9 && std::abs(coordsYP[j]-X)<1.0e-9 && coordsXP[j] < 1.0e-9){
                BotCentre_p = NodesP[j];
            }
        }
    }
    // check if any core has a top node
    TopReduced = TopCentre_p;
    BotReduced = BotCentre_p;
    MPI_Allreduce(MPI_IN_PLACE, &TopReduced, 1, MPIU_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &BotReduced, 1, MPIU_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
    if (TopReduced == 0) throw std::invalid_argument("No top centre p node found for SymFluidInterfaceModel "+MyName);
    if (BotReduced == 0) throw std::invalid_argument("No bottom centre p node found for SymFluidInterfaceModel "+MyName);

    //Ensure Unique top nodes
    std::sort(TopNodes.begin(), TopNodes.end());
    TopNodes.erase(std::unique(TopNodes.begin(), TopNodes.end()), TopNodes.end());

    //Ensure Unique bottom nodes
    std::sort(BotNodes.begin(), BotNodes.end());
    BotNodes.erase(std::unique(BotNodes.begin(), BotNodes.end()), BotNodes.end());

}

/// @brief Loads model from a previous restart savefile
/// @param inputs input: JSON object representing the input file
/// @param data  input: restart data from where to load variables
void SymFluidInterfaceModel::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);
    Commit(CommitTypes::TIMEDEP_COMMIT_TYPE);
}

/// @brief Save model to a restart file
/// @param data input: pointer to object in which to save restart data
void SymFluidInterfaceModel::save(SaveDataFile& data){

}

/// @brief Commits time dependent history variables (Energy history and plastic strains)
/// @param CommitType
void SymFluidInterfaceModel::Commit(int CommitType){
    if (CommitType==CommitTypes::TIMEDEP_COMMIT_TYPE){
        size_t NElems = mesh->ElementGroups[ElemGroupIndex_u].NElems;
        size_t NNodes = mesh->ElementGroups[ElemGroupIndex_u].BaseElem->NodeCount;
        for (size_t i = 0; i < NElems; i++){
            BrokenOld[i] = Broken[i];
            for (size_t j = 0; j < NNodes; j++){
                hMaxOld[i][j] = 1.0*hMax[i][j];
                if (Broken[i]){
                    //hMaxOld[i][j] = std::max(hMaxOld[i][j], u0);
                }
            }
        }
    }
}

void SymFluidInterfaceModel::ResetStep(){
    size_t NElems = mesh->ElementGroups[ElemGroupIndex_u].NElems;
    size_t NNodes = mesh->ElementGroups[ElemGroupIndex_u].BaseElem->NodeCount;
    for (size_t i = 0; i < NElems; i++){
        for (size_t j = 0; j < NNodes; j++){
            hMax[i][j] = 1.0*hMaxOld[i][j];
        }
    }
}

/// @brief Performs set-up of this model based on input file
/// @param inputs input: Input json data file
void SymFluidInterfaceModel::Setup(inputData& inputs){
    if (dim==3){
        DofNames_u.push_back("uz");
        g = {0, 0, -9.81};
    }

    //get indices for relevant element groups
    std::string ElemGroupName;
    inputs.GetRequired(ElemGroupName, {"Models", MyName, "ElementGroup_u"});
    ElemGroupIndex_u = mesh->GetElementGroupIdx(ElemGroupName);

    inputs.GetRequired(ElemGroupName, {"Models", MyName, "ElementGroup_p"});
    ElemGroupIndex_p = mesh->GetElementGroupIdx(ElemGroupName);

    if (ElemGroupIndex_u != ElemGroupIndex_p) {
        LumpedCap = false;
    }

    //get relevant degree of freedmo steps and indices
    std::vector<size_t> dofSteps;
    dofs->getDofTypesSteps(DofNames_u, dofTypes_u, dofSteps);
    if (dofSteps[0] != dofSteps[1]) throw std::invalid_argument(ModelName+" requires "+ DofNames_u[0] + " and " + DofNames_u[1] + " to be in the same solver step\n");
    Step_u = dofSteps[0];

    dofs->getDofTypesSteps(DofName_p, dofType_p, Step_p);
    if (Step_u != Step_p) throw std::invalid_argument(ModelName+" requires "+ DofNames_u[0] + " and " + DofName_p + " to be in the same solver step\n");

    //material parameters
    inputs.GetRequired(SolidMatName, {"Models", MyName, "Solid"});
    S = new SolidMaterial(inputs, SolidMatName);

    inputs.GetRequired(FluidMatName, {"Models", MyName, "Fluid"});
    F = new FluidMaterial(inputs, FluidMatName);

    inputs.GetRequired(FlowModel, {"properties", FluidMatName,"FractureFlow", "FlowType"});
    inputs.GetRequired(kDummy, {"properties", SolidMatName, "FractureProperties", "kDummy"});
    inputs.GetRequired(Gc, {"properties", SolidMatName, "FractureProperties", "Gc"});
    inputs.GetRequired(ft, {"properties", SolidMatName, "FractureProperties", "ft"});

    hOver = 0.0;
    TopBC = inputs.GetOptional(hOver, {"Models", MyName, "hOver"});

    Q_Top = 0.0;
    inputs.GetOptional(Q_Top, {"Models", MyName, "Q_Top"});

    hBase = 0.0;
    Q_Bot = 0.0;
    BotBC = inputs.GetOptional(hBase, {"Models", MyName, "hBase"});
    inputs.GetOptional(Q_Bot, {"Models", MyName, "QBase"});

    u0 = ft/kDummy;

    //saturation model
    std::string satModelString = "Saturated";
    inputs.GetOptional(satModelString, {"properties", FluidMatName, "FractureFlow", "SaturationModel"});
    if (SaturationModels.find(satModelString) == SaturationModels.end()){
        throw std::invalid_argument("Saturation model " + satModelString + " not recognized. Available models are Exponent, ExponentSmooth, Tanh, Saturated and Linear");
    } else {
        satModel = SaturationModels[satModelString];
        inputs.GetOptional(sat_exp_pref, {"properties", FluidMatName, "FractureFlow", "ReferencePressure"});
    }

    //Fracture criterion
    std::string fracCritString = "SingleIP";
    inputs.GetOptional(fracCritString, {"Models", MyName, "FractureCriterion"});
    if (FracCriteriaNames.find(fracCritString) == FracCriteriaNames.end()){
        throw std::invalid_argument("Fracture criterion " + fracCritString + " not recognized. Available criteria are SingleIP and AllIP");
    } else {
        fracCrit = FracCriteriaNames[fracCritString];
    }
}

void SymFluidInterfaceModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    double t = physics->time;
    double dt = physics->timeScheme->dt;

    if (step == Step_u){    //assemble stiffness matrix for momentum balance
        Assemble_UP(f, K);
    }

}

/// @brief Assembles the stiffness matrix and force vector for the fluid interface model
/// @param f output: force vector
/// @param K output: stiffness matrix
void SymFluidInterfaceModel::Assemble_UP(Vec &f, Mat &K){
    size_t ipcount = mesh->ElementGroups[ElemGroupIndex_u].BaseElem->ipcount;    // number of integration points
    size_t nNodes_u = mesh->ElementGroups[ElemGroupIndex_u].NNodes_per_elem;         // number of nodes per displacement element
    size_t nNodes_p = mesh->ElementGroups[ElemGroupIndex_p].NNodes_per_elem;         // number of nodes per pressure element

    std::vector<double> w(ipcount); // integration weights
    std::vector<Eigen::RowVectorXd> Nu(ipcount); for (size_t ip = 0; ip < ipcount; ip++) Nu[ip].resize(nNodes_u); // displacement shape functions
    std::vector<Eigen::MatrixXd> Gu(ipcount); for (size_t ip = 0; ip < ipcount; ip++) Gu[ip].resize(dim-1, nNodes_u);        // displacement shape function gradients
    Eigen::MatrixXd N_uu(3, dim * nNodes_u); // state to displacement mapping matrix

    std::vector<Eigen::RowVectorXd> Np(ipcount); for (size_t ip = 0; ip < ipcount; ip++) Np[ip].resize(nNodes_p); // pressure shape functions
    std::vector<Eigen::MatrixXd> Gp(ipcount); for (size_t ip = 0; ip < ipcount; ip++) Gp[ip].resize(dim-1, nNodes_p);        // pressure shape function gradients

    std::vector<size_t> El_u(nNodes_u), El_p(nNodes_p); // element node indices
    std::vector<PetscInt> dofsU(dim * nNodes_u), dofsP(nNodes_p); // dof indices for displacements and pressure

    Eigen::VectorXd U(nNodes_u*dim), UOld(nNodes_u*dim), P(nNodes_p), POld(nNodes_p); // displacement and pressure vectors
    Eigen::VectorXd WLumped(nNodes_u), WLumpedHOld(nNodes_p);

    Eigen::VectorXd Xc(nNodes_p), Yc(nNodes_p), Zc(nNodes_p); // coordinates of displacement nodes
    Eigen::MatrixXd Coords(nNodes_p, dim); // coordinates matrix for tangent calculation

    double ddXdt2 = physics->timeScheme->ddu_dt;   // state to acceleration derivative
    double dXdt = physics->timeScheme->du_dt;
    double dt= physics->timeScheme->dt;
    double t = physics->time;

    Eigen::VectorXd f_u(nNodes_u*dim), f_p(nNodes_p); // force vectors for displacements and pressure
    Eigen::MatrixXd K_uu(nNodes_u*dim, nNodes_u*dim);
    Eigen::MatrixXd K_pp(nNodes_p, nNodes_p); // stiffness matrices for displacements and pressure
    Eigen::MatrixXd K_up(nNodes_u*dim, nNodes_p); // coupling stiffness matrix
    Eigen::MatrixXd K_pu(nNodes_p, nNodes_u*dim); // coupling stiffness matrix

    Eigen::Vector3d fInt;
    Eigen::Matrix3d df_du;
    Eigen::Vector2d GradP, dq_dSw;
    Eigen::Vector2d q, dq_dh;
    Eigen::Matrix2d dq_dp;
    double Sw, dSw, ddSw;
    double SwOld, dSwOld, ddSwOld;

    Eigen::VectorXd SurfNorm(dim); SurfNorm.setZero(); SurfNorm(0) = 1.0; // surface normal vector

    // loop over all elements in the element group
    LCrack = Z;
    LMeanCrack = 0.0;
    h_Bot = -1.0e99;
    h_Top = -1.0e99;
    p_Bot = -1.0e99;
    p_Top = -1.0e99;
    for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex_u].NElems; el++){
        mesh->getShapeGrads(ElemGroupIndex_u, el, El_u, Nu, Gu, w);
        dofs->getDofForNodes(El_u, dofTypes_u, dofsU);

        mesh->getShapeGrads(ElemGroupIndex_p, el, El_p, Np, Gp, w);
        dofs->getDofForNodes(El_p, dofType_p, dofsP);
        mesh->GetCoordsForNodes(Xc, Yc, Zc, El_p);
        
        // Build coordinates matrix for tangent vector calculation
        if (dim == 3) {
            for (size_t i = 0; i < nNodes_p; i++) {
                Coords(i, 0) = Xc(i);
                Coords(i, 1) = Yc(i);
                Coords(i, 2) = Zc(i);
            }
        }

        physics->StateVectors[Step_u].GetValues(dofsU, U);
        physics->StateVectorsOld[Step_u].GetValues(dofsU, UOld);
        physics->StateVectors[Step_p].GetValues(dofsP, P);
        physics->StateVectorsOld[Step_p].GetValues(dofsP, POld);

        f_u.setZero();
        f_p.setZero();
        K_uu.setZero();
        K_pp.setZero();
        K_up.setZero();
        K_pu.setZero();
        WLumped.setZero();
        WLumpedHOld.setZero();

        if (fracCrit == AllIP){
            Broken[el] = true; // Assume broken, check if any ip is unbroken
        } else {
            Broken[el] = false; // Assume unbroken, check if any ip is broken
        }

        for (size_t ip = 0; ip < ipcount; ip++){
            N_uu.setZero();
            N_uu(0, Eigen::seq(0, nNodes_u - 1)) = Nu[ip];
            N_uu(1, Eigen::seq(nNodes_u, 2 * nNodes_u - 1)) = Nu[ip];
            if (dim==3){
                N_uu(2, Eigen::seq(2*nNodes_u, 3 * nNodes_u - 1)) = Nu[ip];
            }

            double uJump = SurfNorm.transpose() * (N_uu * U);
            double uJumpOld = SurfNorm.transpose() * N_uu * UOld;
            double h = 2.0*uJump;
            double hOld = 2.0*uJumpOld;
            double dh = 2.0;
            double hOldLim = std::max(hOld, 1e-3);

            // Check propagation
            if (BrokenOld[el]){
                Broken[el] = true; // If old was broken, we are still broken
            } else {
                if (fracCrit == SingleIP){
                    if (uJump > u0 && t>0.0){
                        Broken[el] = true; // If jump exceeds threshold, we are broken
                    }
                } else {
                    if (uJump <= u0 || t<=0.0){
                        Broken[el] = false; // If any ip is unbroken, we are unbroken
                    }
                }
            }

            //Fluid Pressure
            double pWall = Np[ip]*P;
            double pWallOld = Np[ip]*POld;
            SwOld = GetSaturation(pWallOld, dSw, ddSw);
            Sw = GetSaturation(pWall, dSw, ddSw);
            if (BrokenOld[el]){
                f_u  += -w[ip] * N_uu.transpose() * SurfNorm * Sw * pWall;
                K_up += -w[ip] * N_uu.transpose() * SurfNorm * Np[ip]*(Sw + dSw*pWall);

                LMeanCrack += w[ip];
                if (Np[ip]*Zc < LCrack){
                    LCrack = Np[ip]*Zc;
                }
            }
            if (FlowModel == "HydroStatic"){
                    // Penalty enforce P
                    double z_ip = Np[ip]*Zc;
                    double pw;
                    pw = F->Density*9.81*(-(z_ip-Z) + hOver);

                    f_p  += 0.5*w[ip] * kDummy * Np[ip].transpose()*(Np[ip]*P-pw);
                    K_pp += 0.5*w[ip] * kDummy * Np[ip].transpose() * Np[ip];
            } else {
                if (BrokenOld[el] == true){ 
                    if (LumpedCap == false){
                        // // Capacity
                        f_p  += 0.5*w[ip]*Np[ip].transpose()*(h - hOld)*Sw/dt;
                        K_pu += 0.5*w[ip]*Np[ip].transpose()*SurfNorm.transpose()*Sw *dh*N_uu / dt;
                        K_pp += 0.5*w[ip]*Np[ip].transpose()*(h - hOld)*dSw/dt*Np[ip];

                        if (h<1e-3){
                            h = 1.0e-3;
                            dh = 0.0;
                        }

                        //saturation changes
                        // f_p  += 0.5*w[ip]*Np[ip].transpose()*h*dSw*Np[ip]*(P - POld)/dt;
                        // K_pp += 0.5*w[ip]*Np[ip].transpose()*h*(ddSw*Np[ip]*(P - POld)+dSw)/dt*Np[ip];
                        // K_pu += 0.5*w[ip]*Np[ip].transpose()*dh*SurfNorm.transpose()*dSw*(Np[ip]*(P - POld))*N_uu/dt;


                        // if (hOld<0.0){
                        //     hOld = 0.0;
                        // }

                        // f_p  += 0.5*w[ip]*Np[ip].transpose()*(h*Sw - hOld*SwOld)/dt;
                        // K_pp += 0.5*w[ip]*Np[ip].transpose()*h*dSw/dt*Np[ip];
                        // K_pu += 0.5*w[ip]*Np[ip].transpose()*dh*Sw/dt*SurfNorm.transpose()*N_uu;
                    }

                    // flow
                    GradP = Gp[ip] * P;
                    
                    // Compute tangent vectors at this integration point
                    Eigen::MatrixXd dXdxi = Gp[ip] * Coords; // (dim-1) x dim matrix
                    Eigen::Vector3d tangent1 = dXdxi.row(0); // first tangent vector (d/dxi1)
                    Eigen::Vector3d tangent2 = dXdxi.row(1); // second tangent vector (d/dxi2)
                    tangent1.normalize();
                    tangent2.normalize();

                    // Project gravity onto the tangent plane
                    Eigen::Vector3d gVec(g[0], g[1], g[2]);
                    Eigen::Vector2d gLocal;
                    gLocal(0) = gVec.dot(tangent1);
                    gLocal(1) = gVec.dot(tangent2);
                    
                    // Add gravity contribution to pressure gradient
                    GetFluidFlux(GradP, hOldLim, Sw, q, dq_dp, dq_dh, dq_dSw, gLocal);
                    f_p  += -0.5*w[ip]    * Gp[ip].transpose() * q;
                    K_pp += -0.5*w[ip]    * Gp[ip].transpose() * (dq_dp * Gp[ip] + dq_dSw*dSw * Np[ip]);
                    //K_pu += -0.5*dh*w[ip] * Gp[ip].transpose() * dq_dh * SurfNorm.transpose() * N_uu;
                }
            }

            WLumped += w[ip]*Nu[ip];
            WLumpedHOld += w[ip]*Np[ip]*hOldLim;
        }

        if ((FlowModel != "HydroStatic")){
            // Compressibility
            for (size_t n = 0; n < nNodes_p; n++){
                double stab = 1;
                if (t < 0.0){
                    stab = 1.0e6;
                } else {
                    stab = 1.0e6;
                }
                f_p(n)    += stab*0.5*WLumpedHOld[n] / F->Bulk*(P(n) - POld(n)) / dt;
                K_pp(n,n) += stab*0.5*WLumpedHOld[n] / F->Bulk / dt;
            }

            // Saturation changes
            if (LumpedCap == false){
                for (size_t n = 0; n < nNodes_p; n++){
                    Sw = GetSaturation(P(n), dSw, ddSw);
                    f_p(n)    += 0.5*WLumpedHOld[n]*(dSw*(P(n)-POld(n)))/dt;
                    K_pp(n,n) += 0.5*WLumpedHOld[n]*(ddSw*(P(n)-POld(n))+dSw)/dt;
                }
            }

        }

        //pressure capacity
        if ((FlowModel != "HydroStatic") && (LumpedCap == true) ){
            for (size_t n = 0; n < nNodes_u; n++){
                int ny = n + nNodes_u;
                int nz = n + 2*nNodes_u;

                Eigen::Vector3d UNode, UNodeOld; 
                UNode[0] = U[n]; UNodeOld[0] = UOld[n];
                UNode[1] = U[ny]; UNodeOld[1] = UOld[ny];
                UNode[2] = U[nz]; UNodeOld[2] = UOld[nz];

                double uJump = SurfNorm.transpose()*UNode;
                double uJumpOld = SurfNorm.transpose()*UNodeOld;
                double h = 2.0*uJump;
                double hOld = 2.0*uJumpOld;
                double dh = 2.0;

                Sw = GetSaturation(P(n), dSw, ddSw);

                // Capacity
                if (BrokenOld[el] == true){
                    // if (t < 0.0){
                    //     // Initial period, set h to 1e-3 to allow notch to open
                    //     h = 1e-3;
                    //     hOld = 1e-3;
                    //     dh = 0.0;
                    // }

                    //height changes
                    f_p(n)     += 0.5*WLumped[n]*Sw*(h-hOld)/dt;

                    K_pu(n,n)  += 0.5*WLumped[n]*Sw*SurfNorm[0]*dh/dt;
                    K_pu(n,ny) += 0.5*WLumped[n]*Sw*SurfNorm[1]*dh/dt;
                    K_pu(n,nz) += 0.5*WLumped[n]*Sw*SurfNorm[2]*dh/dt;
                    K_pp(n,n)  += 0.5*WLumped[n]*(h-hOld)*dSw/dt;


                    if (h < 1e-3) {
                        h = 1e-3;
                        dh = 0.0;
                    }
                    if (hOld < 1e-3){
                        hOld = 1e-3;
                    }
                    //saturation changes
                    f_p(n)     += 0.5*WLumped[n]*h*dSw*(P(n)-POld(n))/dt;
                    K_pp(n,n)  += 0.5*WLumped[n]*h*(ddSw*(P(n)-POld(n))+dSw)/dt;
                    K_pu(n,n)  += 0.5*WLumped[n]*dh*SurfNorm[0]*dSw*(P(n)-POld(n))/dt;
                    K_pu(n,ny) += 0.5*WLumped[n]*dh*SurfNorm[1]*dSw*(P(n)-POld(n))/dt;
                    K_pu(n,nz) += 0.5*WLumped[n]*dh*SurfNorm[2]*dSw*(P(n)-POld(n))/dt;
                }
            }
        }


        // Cohesive zone model
        for (size_t n = 0; n < nNodes_u; n++){
            fInt.setZero();
            df_du.setZero();

            int ny = n + nNodes_u;
            int nz = n + 2*nNodes_u;

            Eigen::Vector3d UNode; 
            UNode[0] = U[n];
            UNode[1] = U[ny];
            UNode[2] = U[nz];

            double uJump = 2.0*SurfNorm.transpose()*UNode;
            double uMaxOld = hMaxOld[el][n];
            hMax[el][n] = std::max(uJump, uMaxOld);

            if (BrokenOld[el]){ //BrokenOld[el]
                if (uJump < 0.0){ // No Penetration
                    fInt = kDummy * uJump * SurfNorm;
                    df_du = kDummy  * SurfNorm* SurfNorm.transpose();
                } else if (uJump < u0 && uMaxOld < u0){ //never been fractured, elastic
                    fInt = kDummy * uJump * SurfNorm;
                    df_du = kDummy  * SurfNorm* SurfNorm.transpose();
                } else if (uJump >= uMaxOld){ 
                    // Loading, exponential softening through CZM
                //     fInt = ft * std::exp(- (uJump - u0) * ft / Gc) * SurfNorm;
                //     df_du = -ft*ft/Gc * std::exp(- (uJump - u0) * ft / Gc) *SurfNorm* SurfNorm.transpose();

                    // Loading, softening through Linear CZM
                    double uc = 2.0 * Gc / ft + u0;  // Critical separation
                    if (uJump < uc) {
                        double traction = ft * (1.0 - (uJump - u0) / (uc - u0));
                        fInt = traction * SurfNorm;
                        df_du = -ft / (uc - u0) * SurfNorm * SurfNorm.transpose();
                    } else {
                        fInt.setZero();
                        df_du.setZero();
                    }

                } else {
                    // Unloading below threshold, use linear from last max
                    // double fMax = ft * std::exp(- (uMaxOld - u0) * ft / Gc);
                    // fInt = fMax * (uJump / uMaxOld) * SurfNorm;
                    // df_du = fMax/uMaxOld  * SurfNorm* SurfNorm.transpose();

                    // Unloading below threshold, use linear from last max - Linear CZM
                    double uc = 2.0 * Gc / ft + u0;  // Critical separation
                    double fMax;
                    if (uMaxOld < uc)
                        fMax = ft * (1.0 - (uMaxOld - u0) / (uc - u0));
                    else
                        fMax = 0.0;
                    fInt = fMax * (uJump / uMaxOld) * SurfNorm;
                    df_du = fMax / uMaxOld * SurfNorm * SurfNorm.transpose();

                }
            } else {
                fInt = kDummy * uJump * SurfNorm;
                df_du = kDummy  * SurfNorm* SurfNorm.transpose();
            }

            // check for NaN's
            if (std::isnan(fInt[0]) || std::isnan(fInt[1]) || std::isnan(fInt[2])){
                std::cout << "============================================" << std::endl;
                std::cout << "ft" << ft << ", Gc: " << Gc << "u0: " << u0 << std::endl;
                std::cout << "uJump: " << uJump << ", uMaxOld: " << uMaxOld << ", fInt: " << fInt.transpose() << std::endl;
                std::cout << "df_du: " << std::endl << df_du << std::endl;
                std::cout << "SurfNorm: " << SurfNorm.transpose() << std::endl;
                std::cout << "============================================" << std::endl;

                throw std::runtime_error("NaN detected in cohesive force calculation in SymFluidInterfaceModel "+MyName);
            }

            f_u[n]  += WLumped[n]*fInt[0];  
            f_u[ny] += WLumped[n]*fInt[1];  
            f_u[nz] += WLumped[n]*fInt[2];
            K_uu(n, n)   += WLumped[n] * 2.0*df_du(0, 0); 
            K_uu(n, ny)  += WLumped[n] * 2.0*df_du(0, 1); 
            K_uu(n, nz)  += WLumped[n] * 2.0*df_du(0, 2);
            K_uu(ny, n)  += WLumped[n] * 2.0*df_du(1, 0);
            K_uu(ny, ny) += WLumped[n] * 2.0*df_du(1, 1);
            K_uu(ny, nz) += WLumped[n] * 2.0*df_du(1, 2);
            K_uu(nz, n)  += WLumped[n] * 2.0*df_du(2, 0);
            K_uu(nz, ny) += WLumped[n] * 2.0*df_du(2, 1);
            K_uu(nz, nz) += WLumped[n] * 2.0*df_du(2, 2);
        }


        VecAdd(f, dofsU, f_u);
        VecAdd(f, dofsP, f_p);
        MatAdd(K, dofsU, dofsU, K_uu);
        MatAdd(K, dofsP, dofsP, K_pp);
        MatAdd(K, dofsU, dofsP, K_up);
        MatAdd(K, dofsP, dofsU, K_pu);
    }

    // Top Boundary condition
    if (TopBC && TopNodes.size() > 0){
        size_t nTopNodes = TopNodes.size();
        PetscInt dofTop;
        double PTop;
        double pTarget = hOver * F->Density * 9.81; // target pressure at top nodes
        for (size_t i = 0; i < nTopNodes; i++){
            dofs->getDofForNodes(TopNodes[i], dofType_p, dofTop);
            PTop = physics->StateVectors[Step_p].GetValue(dofTop);

            double fTop = -1000*kDummy * (PTop - pTarget);
            double kTop = -1000*kDummy;
            VecAdd(f, dofTop, fTop);
            MatAdd(K, dofTop, dofTop, kTop);
        }
    }

    if (TopCentre_p > 0 && t>=0.0){
        PetscInt dofTop, uDofTop;
        dofs->getDofForNodes(TopCentre_p, dofType_p, dofTop);
        dofs->getDofForNodes(TopCentre_u, dofTypes_u[0], uDofTop);
        p_Top = physics->StateVectors[Step_p].GetValue(dofTop);
        h_Top = 2.0*physics->StateVectors[Step_u].GetValue(uDofTop);
        double fTop = -Q_Top;
        VecAdd(f, dofTop, fTop);
    }

    //Bottom Boundary condition
    if (BotBC && BotNodes.size() > 0){
        size_t nBotNodes = BotNodes.size();
        PetscInt dofBot, dofUBot;
        double PBot;
        double UZBot;
        double pTarget;
        for (size_t i = 0; i < nBotNodes; i++){
            dofs->getDofForNodes(BotNodes[i], dofType_p, dofBot);
            //dofs->getDofForNodes(BotNodes[i], dofTypes_u[2], dofUBot);
            PBot = physics->StateVectors[Step_p].GetValue(dofBot);
            //UZBot = physics->StateVectors[Step_u].GetValue(dofUBot);
            UZBot = 0.0;
            pTarget = (hBase-UZBot) * F->Density * 9.81; // target pressure at bottom nodes
            double fBot = -1000*kDummy * (PBot - pTarget);
            double kBot = -1000*kDummy;
            double kUBot = -1000*kDummy * F->Density * 9.81;
            VecAdd(f, dofBot, fBot);
            MatAdd(K, dofBot, dofBot, kBot);
            //MatAdd(K, dofBot, dofUBot, kUBot);
        }
    }

    if (BotCentre_p > 0 && t>=0.0){
        PetscInt dofBot, udofBot;
        dofs->getDofForNodes(BotCentre_p, dofType_p, dofBot);
        dofs->getDofForNodes(BotCentre_u, dofTypes_u[0], udofBot);
        p_Bot = physics->StateVectors[Step_p].GetValue(dofBot);
        h_Bot = 2.0*physics->StateVectors[Step_u].GetValue(udofBot);
        double fBot = -Q_Bot;
        VecAdd(f, dofBot, fBot);
    }

}

void SymFluidInterfaceModel::GetFluidFlux(Eigen::Vector2d gradP, double h, double Sw, Eigen::Vector2d&q, Eigen::Matrix2d &dq_dGradP, Eigen::Vector2d &dq_dh, Eigen::Vector2d &dq_dSw, Eigen::Vector2d& gLocal){
    double Sw_eps = 1e-2; // to avoid zero flow at zero saturation
    if (FlowModel == "Laminar"){
        double expFact = 2;
        //cubic Law flow
        q = -std::pow((1.0-Sw_eps)*Sw+Sw_eps,expFact)*h*h*h/12.0/F->visc * gradP;
        dq_dGradP = -std::pow((1.0-Sw_eps)*Sw+Sw_eps,expFact)*h*h*h/12.0/F->visc*Eigen::Matrix2d::Identity();
        dq_dh = -3*std::pow((1.0-Sw_eps)*Sw+Sw_eps,expFact)*h*h/12.0/F->visc * gradP;
        dq_dSw = -h*h*h*expFact*(1.0-Sw_eps)*std::pow((1.0-Sw_eps)*Sw+Sw_eps,expFact-1)/12.0/F->visc * gradP;

        q += std::pow(Sw,expFact)*h*h*h/12.0/F->visc * F->Density * gLocal;
        dq_dGradP += Eigen::Matrix2d::Zero();
        dq_dh += std::pow(Sw,expFact)*h*h/4.0/F->visc * F->Density * gLocal;
        dq_dSw += h*h*h*expFact*std::pow(Sw,expFact-1)/12.0/F->visc * F->Density * gLocal;
    } else if (FlowModel == "Turbulent"){
        //Friction factor turbulent flow
        double kwall = 1e-2; // roughness
        double f0 = 0.143; // friction factor
        double pNorm = gradP.norm() + 1.0e-10;
        double preFac = 2.0*std::pow(F->Density, -0.5) * std::pow(kwall, -0.16666) * std::pow(f0, -0.5);
        double hPow = std::pow(h, 5.0/3.0);
        double invSqrtPNorm =std::pow(pNorm, -0.5);
        double outerCoeff = -0.5 * invSqrtPNorm;

        q         = -Sw*preFac*         hPow * invSqrtPNorm * gradP;
        dq_dGradP = -Sw*preFac*         hPow * (invSqrtPNorm * Eigen::Matrix2d::Identity() + outerCoeff * (gradP * gradP.transpose()));
        dq_dh     = -Sw*preFac*5.0/3.0* std::pow(h, 2.0/3.0)*invSqrtPNorm * gradP;
        dq_dSw    = -   preFac*         hPow * invSqrtPNorm * gradP;
    }
}

bool SymFluidInterfaceModel::SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>>& Data){
    if (SaveLoc=="Nodes" && DataName == "hCrack" && ElemGroup == ElemGroupIndex_u){
        size_t nNodes_u = mesh->ElementGroups[ElemGroupIndex_u].NNodes_per_elem;
        size_t NElems = mesh->ElementGroups[ElemGroupIndex_u].NElems;
        size_t NExport = mesh->ElementGroups[ElemGroupIndex_u].BaseElem->NExport.size();

        std::vector<Eigen::RowVectorXd> NExportShape(NExport);
        for (size_t i = 0; i < NExport; i++) NExportShape[i].resize(nNodes_u);

        std::vector<size_t> El_u(nNodes_u);
        std::vector<PetscInt> dofsU(dim*nNodes_u);

        Eigen::MatrixXd N_uu(3,dim*nNodes_u);
        Eigen::VectorXd U(dim*nNodes_u);

        Data.resize(NElems);
        for (size_t i = 0; i < NElems; i++){
            mesh->getExportShape(ElemGroupIndex_u, i, El_u, NExportShape);

            mesh->GetNodesForElem(El_u, ElemGroupIndex_u, i);
            dofs->getDofForNodes(El_u, dofTypes_u, dofsU);
            physics->StateVectors[Step_u].GetValues(dofsU, U);

            Data[i].resize(NExport); // Only save the four corner nodes
            for (size_t j = 0; j < NExport; j++){
                N_uu.setZero();
                N_uu(0,Eigen::seq(0,nNodes_u-1)) = NExportShape[j];
                N_uu(1,Eigen::seq(nNodes_u,2*nNodes_u-1)) = NExportShape[j];
                if (dim==3){
                    N_uu(2,Eigen::seq(2*nNodes_u,3*nNodes_u-1)) = NExportShape[j];
                }

                Eigen::Vector3d disp = N_uu*U;
                // hCrack is twice the normal jump (assuming normal is x-direction)
                double hCrack = 2.0 * disp(0);
                Data[i][j] = hCrack;
            }
        }
        return true;
    }
    if (SaveLoc=="Nodes" && DataName == "Sw" && ElemGroup == ElemGroupIndex_p){
        size_t nNodes_p = mesh->ElementGroups[ElemGroupIndex_p].NNodes_per_elem;
        size_t NElems = mesh->ElementGroups[ElemGroupIndex_p].NElems;
        size_t NExport = mesh->ElementGroups[ElemGroupIndex_p].BaseElem->NExport.size();

        std::vector<Eigen::RowVectorXd> NExportShape(NExport);
        for (size_t i = 0; i < NExport; i++) NExportShape[i].resize(nNodes_p);

        std::vector<size_t> El_p(nNodes_p);
        std::vector<PetscInt> dofsP(nNodes_p);
        Eigen::VectorXd P(nNodes_p);

        Data.resize(NElems);
        for (size_t i = 0; i < NElems; i++){
            mesh->getExportShape(ElemGroupIndex_p, i, El_p, NExportShape);

            mesh->GetNodesForElem(El_p, ElemGroupIndex_p, i);
            dofs->getDofForNodes(El_p, dofType_p, dofsP);
            physics->StateVectors[Step_p].GetValues(dofsP, P);

            Data[i].resize(NExport); // Only save the four corner nodes
            for (size_t j = 0; j < NExport; j++){
                double p = NExportShape[j]*P;
                double Unused1, Unused2;
                Data[i][j] = GetSaturation(p, Unused1, Unused2);
            }
        }
        return true;
    }
    if (SaveLoc=="Nodes" && DataName == "pOver" && ElemGroup == ElemGroupIndex_p){
        size_t nNodes_p = mesh->ElementGroups[ElemGroupIndex_p].NNodes_per_elem;
        size_t NElems = mesh->ElementGroups[ElemGroupIndex_p].NElems;
        size_t NExport = mesh->ElementGroups[ElemGroupIndex_p].BaseElem->NExport.size();
        Eigen::VectorXd Xc(nNodes_p), Yc(nNodes_p), Zc(nNodes_p); // coordinates of displacement nodes

        std::vector<Eigen::RowVectorXd> NExportShape(NExport);
        for (size_t i = 0; i < NExport; i++) NExportShape[i].resize(nNodes_p);

        std::vector<size_t> El_p(nNodes_p);
        std::vector<PetscInt> dofsP(nNodes_p);
        Eigen::VectorXd P(nNodes_p);

        Data.resize(NElems);
        for (size_t i = 0; i < NElems; i++){
            mesh->getExportShape(ElemGroupIndex_p, i, El_p, NExportShape);

            mesh->GetNodesForElem(El_p, ElemGroupIndex_p, i);
            mesh->GetCoordsForNodes(Xc, Yc, Zc, El_p);
            dofs->getDofForNodes(El_p, dofType_p, dofsP);
            physics->StateVectors[Step_p].GetValues(dofsP, P);

            Data[i].resize(NExport); // Only save the four corner nodes
            for (size_t j = 0; j < NExport; j++){
                double p = NExportShape[j]*P;
                double z = NExportShape[j]*Zc;

                double Unused1, Unused2;
                double pOver = p-9.81*S->Density*(Z-z);
                Data[i][j] = std::max(0.0, pOver);
            }
        }
        return true;
    }
    if (SaveLoc=="Nodes" && DataName == "broken" && ElemGroup == ElemGroupIndex_u){
        size_t nNodes_u = mesh->ElementGroups[ElemGroupIndex_u].NNodes_per_elem;
        size_t NElems = mesh->ElementGroups[ElemGroupIndex_u].NElems;
        size_t NExport = mesh->ElementGroups[ElemGroupIndex_u].BaseElem->NExport.size();

        Data.resize(NElems);
        for (size_t i = 0; i < NElems; i++){
            Data[i].resize(NExport);
            for (size_t j = 0; j < NExport; j++){
                Data[i][j] = (Broken[i]+BrokenOld[i]) / 2.0; // Average of current and old state
            }
        }

        return true;
    }
    if (SaveLoc=="ip" && (DataName == "qx" || DataName == "qy") && ElemGroup == ElemGroupIndex_p){
        size_t nNodes_p = mesh->ElementGroups[ElemGroupIndex_p].NNodes_per_elem;
        size_t nNodes_u = mesh->ElementGroups[ElemGroupIndex_u].NNodes_per_elem;
        size_t NElems = mesh->ElementGroups[ElemGroupIndex_p].NElems;
        size_t ipcount = mesh->ElementGroups[ElemGroupIndex_p].BaseElem->ipcount;

        std::vector<Eigen::RowVectorXd> Np(ipcount); for (size_t ip = 0; ip < ipcount; ip++) Np[ip].resize(nNodes_p); // pressure shape functions
        std::vector<Eigen::MatrixXd> Gp(ipcount); for (size_t ip = 0; ip < ipcount; ip++) Gp[ip].resize(dim-1, nNodes_p);        // pressure shape function gradients
        std::vector<Eigen::RowVectorXd> Nu(ipcount); for (size_t ip = 0; ip < ipcount; ip++) Nu[ip].resize(nNodes_u); // displacement shape functions
        std::vector<Eigen::MatrixXd> Gu(ipcount); for (size_t ip = 0; ip < ipcount; ip++) Gu[ip].resize(dim-1, nNodes_u);        // displacement shape function gradients
        std::vector<double> w(ipcount);

        std::vector<size_t> El_p(nNodes_p), El_u(nNodes_u);
        std::vector<PetscInt> dofsP(nNodes_p), dofsU(nNodes_u);
        Eigen::VectorXd P(nNodes_p), U(nNodes_u);
        Eigen::VectorXd Xc(nNodes_p), Yc(nNodes_p), Zc(nNodes_p);
        Eigen::MatrixXd Coords(nNodes_p, dim);

        Data.resize(NElems);
        for (size_t el = 0; el < NElems; el++){
            Data[el].resize(ipcount);
            mesh->getShapeGrads(ElemGroupIndex_u, el, El_u, Nu, Gu, w);
            mesh->getShapeGrads(ElemGroupIndex_p, el, El_p, Np, Gp, w);
            mesh->GetCoordsForNodes(Xc, Yc, Zc, El_p);
            
            // Build coordinates matrix for tangent vector calculation
            if (dim == 3) {
                for (size_t i = 0; i < nNodes_p; i++) {
                    Coords(i, 0) = Xc(i);
                    Coords(i, 1) = Yc(i);
                    Coords(i, 2) = Zc(i);
                }
            }
            
            dofs->getDofForNodes(El_p, dofType_p, dofsP);
            dofs->getDofForNodes(El_u, dofTypes_u[0], dofsU);
            physics->StateVectors[Step_u].GetValues(dofsU, U);
            physics->StateVectors[Step_p].GetValues(dofsP, P);

            for (size_t ip = 0; ip < ipcount; ip++){
                double h = 2.0 * Nu[ip] * U; // crack opening height
                double p = Np[ip] * P;            // fluid pressure
                double Sw, dSw, ddSw;
                Sw = GetSaturation(p, dSw, ddSw);
                if (h < 1.0e-6) h = 1.0e-6;

                Eigen::Vector2d GradP = Gp[ip] * P;
                
                // Compute tangent vectors at this integration point
                Eigen::MatrixXd dXdxi = Gp[ip] * Coords; // (dim-1) x dim matrix
                Eigen::Vector3d tangent1 = dXdxi.row(0); // first tangent vector (d/dxi1)
                Eigen::Vector3d tangent2 = dXdxi.row(1); // second tangent vector (d/dxi2)
                tangent1.normalize();
                tangent2.normalize();
                
                // Project gravity onto the tangent plane
                Eigen::Vector3d gVec(g[0], g[1], g[2]);
                Eigen::Vector2d gLocal;
                gLocal(0) = gVec.dot(tangent1);
                gLocal(1) = gVec.dot(tangent2);
                
                // Add gravity contribution to pressure gradient
                Eigen::Vector2d q, dq_dh, dq_dSw;
                Eigen::Matrix2d dq_dGradP;
                if (Broken[el]==false){
                    q.setZero();
                } else {
                    GetFluidFlux(GradP, h, Sw, q, dq_dGradP, dq_dh, dq_dSw, gLocal);
                    if (std::isnan(q[0]) || std::isnan(q[1])){
                        std::cout << "NaN flux detected at element " << el << ", ip " << ip << ", h=" << h << ", p=" << p << std::endl;
                        q.setZero();
                    }
                }

                if (DataName == "qx"){
                    Data[el][ip] = q[0];
                } else if (DataName == "qy"){
                    Data[el][ip] = q[1];
                }
            }
        }
        return true;
    }
    return false; // Not implemented yet
}

size_t SymFluidInterfaceModel::hasTimeData(std::vector<std::string>& DataNames){
    size_t nData = 8;
    DataNames.resize(nData);
    DataNames[0] = MyName+"/L_crack";
    DataNames[1] = MyName+"/LMean_crack";
    DataNames[2] = MyName+"/h_Top";
    DataNames[3] = MyName+"/h_Bot";
    DataNames[4] = MyName+"/p_Top";
    DataNames[5] = MyName+"/p_Bot";
    DataNames[6] = MyName+"/peff_Top";
    DataNames[7] = MyName+"/peff_Bot";
    return nData;
}

void SymFluidInterfaceModel::GetTimeData(std::vector<double>& DataValues){
    double LCrackRed = LCrack;
    double LMeanCrackRed = LMeanCrack;
    double p_TopRed = p_Top;
    double p_BotRed = p_Bot;
    double h_TopRed = h_Top;
    double h_BotRed = h_Bot;
    MPI_Allreduce(MPI_IN_PLACE, &LCrackRed, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &LMeanCrackRed, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &p_TopRed, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &p_BotRed, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &h_TopRed, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &h_BotRed, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);

    size_t nData = 8;
    DataValues.resize(0);
    DataValues.push_back(Z-LCrackRed);
    DataValues.push_back(LMeanCrackRed/X);
    DataValues.push_back(h_TopRed);
    DataValues.push_back(h_BotRed);
    DataValues.push_back(p_TopRed);
    DataValues.push_back(p_BotRed);

    double unused1, unused2;
    double Sw = GetSaturation(p_TopRed, unused1, unused2);
    double peff_Top = Sw*p_TopRed;
    Sw = GetSaturation(p_BotRed, unused1, unused2);
    double peff_Bot = Sw*p_BotRed;
    DataValues.push_back(peff_Top);
    DataValues.push_back(peff_Bot);
}

double SymFluidInterfaceModel::GetSaturation(double p, double& dSw, double& ddSw){
    double Sw;
    double S_0 = 0.0;

    switch (satModel){
        case FullSaturation:{ 
            Sw = 1.0;
            dSw = 0.0;
            ddSw = 0.0;
        } break;
        case ExponentSmoothSaturation:{
            double Offset = 0.0;
            double x = p / sat_exp_pref + Offset;
            Sw   = S_0 + (1.0-S_0)*sigmoid(x);
            dSw  = (1.0-S_0)*(1.0 / sat_exp_pref) * sigmoid_derivative(x);
            ddSw =-(1.0-S_0)*(1.0 / (sat_exp_pref * sat_exp_pref)) * sigmoid_second_derivative(x);

        } break;
        case TanhSaturation:{ //note: functionally the same as ExponentSmooth, but this one is used in other literature
            Sw   = S_0 + (1.0-S_0)*(0.5+0.5*tanh(p/sat_exp_pref));
            dSw  = (1.0-S_0)*0.5/sat_exp_pref*std::pow(cosh(p/sat_exp_pref), -2);
            ddSw =-(1.0-S_0)/sat_exp_pref/sat_exp_pref*std::pow(cosh(p/sat_exp_pref), -3)*sinh(p/sat_exp_pref);
        } break;
        case ExponentSaturation:{
            if (p<0){
                Sw = S_0 + (1.0-S_0)*std::exp(p/sat_exp_pref);
                dSw = (1.0-S_0)*1.0/sat_exp_pref*std::exp(p/sat_exp_pref);
                ddSw = (1.0-S_0)*1.0/sat_exp_pref/sat_exp_pref*std::exp(p/sat_exp_pref);
            } else {
                Sw = 1.0;
                dSw = 0.0;
                ddSw = 0.0;
            }
        } break;
        case LinearSaturation:{
            if (p<-sat_exp_pref){
                Sw = S_0;
                dSw = 0.0;
                ddSw = 0.0;
            } else if (p<0){
                Sw = S_0 + (1.0-S_0)*(p+sat_exp_pref)/sat_exp_pref;
                dSw = (1.0-S_0)/sat_exp_pref;
                ddSw = 0.0;
            } else {
                Sw = 1.0;
                dSw = 0.0;
                ddSw = 0.0;
            }
        } break;
        default:{
            throw std::invalid_argument("Saturation function not defined in PoroElasticRelations.cpp,");
        }
    }

    if (std::isnan(Sw) || std::isnan(dSw) || std::isnan(ddSw)){
        std::stringstream  ErrString;
        ErrString << "Saturnation function is NaN in PoroElasticRelations.cpp,\n" << "Sw=" << Sw << ", dSw=" << dSw << ", ddSw=" << ddSw << ", pw=" << p << "\n";
        std::cout << ErrString.str();
        throw std::invalid_argument(ErrString.str());
    }

    return Sw;
}
