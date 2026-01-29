#include "LargeDefModel.h"
#include "../../Physics/physics.h"
#include <Eigen/Dense>



/// NOT YET TESTED DO NOT USE IN CURRENT STATE


/// @brief Registers the model into the model creator
void Register_LargeDeformationModel(){
    ModelNames.push_back("LargeDeformationModel");
    ModelCreators.push_back(New_LargeDeformationModel);
}

/// @brief Creates a new general solid model
/// @param My_Physics input: physics object
/// @param MyNameIn input: name of this model, from where to use input properties
/// @return pointer to the newly created model
BaseModel* New_LargeDeformationModel(Physics& My_Physics, std::string MyNameIn){
    // Throw an error as this model is not yet tested
    std::invalid_argument("LargeDeformationModel is not yet tested and should not be used in its current state.\n");

    return new LargeDeformationModel(My_Physics, MyNameIn);
}

/// @brief Initializer, forwards to baseModel
/// @param My_Physics //input: physics object, reference copied into model
/// @param MyName   //name of this model, from where to use input properties
LargeDeformationModel::LargeDeformationModel(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "LargeDeformationModel";
}

LargeDeformationModel::~LargeDeformationModel(){

}

/// @brief Initialize this model from scratch
/// @param inputs input: JSON object representing the input file
void LargeDeformationModel::init(inputData& inputs){
    Setup(inputs);

    std::vector<size_t> UniqueNodes_u = mesh->ElementGroups[ElemGroupIndex_u].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes_u, dofTypes_u);

    size_t NElems = mesh->ElementGroups[ElemGroupIndex_u].NElems;
    size_t NIP = mesh->ElementGroups[ElemGroupIndex_u].BaseElem->ipcount;
    // No plasticity: nothing to initialize here for large-strain elasticity
}

/// @brief Loads model from a previous restart savefile
/// @param inputs input: JSON object representing the input file
/// @param data  input: restart data from where to load variables
void LargeDeformationModel::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);
    //data.Load("Energy", Energy_Hist);
    // No plasticity: nothing to load
}

/// @brief Save model to a restart file
/// @param data input: pointer to object in which to save restart data
void LargeDeformationModel::save(SaveDataFile& data){
    //data.Save("Energy", Energy_Hist);
    // No plasticity: nothing to save
}

/// @brief Commits time dependent history variables (Energy history and plastic strains)
/// @param CommitType
void LargeDeformationModel::Commit(int CommitType){
    // No plasticity: nothing to commit
}

void LargeDeformationModel::ResetStep(){
    // No plasticity: nothing to reset
}

/// @brief Performs set-up of this model based on input file
/// @param inputs input: Input json data file
void LargeDeformationModel::Setup(inputData& inputs){
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

void LargeDeformationModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    double t = physics->time;
    double dt = physics->timeScheme->dt;

    if (step == Step_u){    //assemble stiffness matrix for momentum balance
        Assemble_U(f, K);
    }

}

void LargeDeformationModel::Assemble_U(Vec &f, Mat &K){
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
    // Large-strain: will use F, C, E, and 2nd Piola-Kirchhoff stress

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

            // Compute deformation gradient F using displacement gradients
            // Vectorized computation of F
            Eigen::MatrixXd Umat = Eigen::Map<const Eigen::MatrixXd>(U.data(), dim, nNodes_u);
            Eigen::MatrixXd F = Eigen::MatrixXd::Identity(dim, dim) + Umat * Gu[ip].transpose();

            double J = F.determinant();
            Eigen::MatrixXd C = F.transpose() * F;
            Eigen::MatrixXd Cinv = C.inverse();

            // Material parameters
            double mu = M->Shear; // Shear modulus
            double lambda = M->Lame; // Lame's first parameter

            // 2nd Piola-Kirchhoff stress for compressible Neo-Hookean
            // S = mu*(I - Cinv) + lambda*log(J)*Cinv
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(dim, dim);
            Eigen::MatrixXd S = mu * (I - Cinv) + lambda * std::log(J) * Cinv;

            // Convert S to Voigt notation (for 3D: xx, yy, zz, xy, yz, xz)
            Vector6d S_voigt;
            S_voigt.setZero();
            S_voigt(0) = S(0,0);
            S_voigt(1) = S(1,1);
            S_voigt(2) = S(2,2);
            S_voigt(3) = S(0,1); // xy
            S_voigt(4) = S(1,2); // yz
            S_voigt(5) = S(0,2); // xz

            // Internal force: Fint = B_G^T * S_voigt (total Lagrangian)
            // Construct large-strain B_G matrix (mapping nodal displacements to Green-Lagrange strain)
            // Fully nonlinear B_G (large-strain B-matrix) including all F terms, for block/component-major ordering
            // U = [U_x1, U_x2, ..., U_y1, U_y2, ..., U_z1, U_z2, ...]
            Eigen::MatrixXd B_G = Eigen::MatrixXd::Zero(6, dim * nNodes_u);
            if (dim == 2) {
                // Voigt: [E11, E22, -, E12, -, -]
                // E11: sum_j F(0,j) * Gu[ip](j,:)
                B_G.row(0).segment(0, nNodes_u) = (F.row(0) * Gu[ip]).transpose();
                // E22: sum_j F(1,j) * Gu[ip](j,:)
                B_G.row(1).segment(nNodes_u, nNodes_u) = (F.row(1) * Gu[ip]).transpose();
                // E12 (symmetric): 0.5 * [sum_j F(0,j) * Gu[ip](j,:) for Uy + sum_j F(1,j) * Gu[ip](j,:) for Ux]
                B_G.row(3).segment(0, nNodes_u) = 0.5 * (F.row(1) * Gu[ip]).transpose();
                B_G.row(3).segment(nNodes_u, nNodes_u) = 0.5 * (F.row(0) * Gu[ip]).transpose();
            } else if (dim == 3) {
                // Voigt: [E11, E22, E33, E12, E23, E13]
                // E11: sum_j F(0,j) * Gu[ip](j,:)
                B_G.row(0).segment(0 * nNodes_u, nNodes_u) = (F.row(0) * Gu[ip]).transpose();
                // E22: sum_j F(1,j) * Gu[ip](j,:)
                B_G.row(1).segment(1 * nNodes_u, nNodes_u) = (F.row(1) * Gu[ip]).transpose();
                // E33: sum_j F(2,j) * Gu[ip](j,:)
                B_G.row(2).segment(2 * nNodes_u, nNodes_u) = (F.row(2) * Gu[ip]).transpose();
                // E12 (symmetric): 0.5 * [sum_j F(1,j) * Gu[ip](j,:) for Ux + sum_j F(0,j) * Gu[ip](j,:) for Uy]
                B_G.row(3).segment(0 * nNodes_u, nNodes_u) = 0.5 * (F.row(1) * Gu[ip]).transpose();
                B_G.row(3).segment(1 * nNodes_u, nNodes_u) = 0.5 * (F.row(0) * Gu[ip]).transpose();
                // E23 (symmetric): 0.5 * [sum_j F(2,j) * Gu[ip](j,:) for Uy + sum_j F(1,j) * Gu[ip](j,:) for Uz]
                B_G.row(4).segment(1 * nNodes_u, nNodes_u) = 0.5 * (F.row(2) * Gu[ip]).transpose();
                B_G.row(4).segment(2 * nNodes_u, nNodes_u) = 0.5 * (F.row(1) * Gu[ip]).transpose();
                // E13 (symmetric): 0.5 * [sum_j F(2,j) * Gu[ip](j,:) for Ux + sum_j F(0,j) * Gu[ip](j,:) for Uz]
                B_G.row(5).segment(0 * nNodes_u, nNodes_u) = 0.5 * (F.row(2) * Gu[ip]).transpose();
                B_G.row(5).segment(2 * nNodes_u, nNodes_u) = 0.5 * (F.row(0) * Gu[ip]).transpose();
            }
            F_U += w[ip] * B_G.transpose() * S_voigt;


            // Full consistent tangent for compressible Neo-Hookean (material part)
            // 4th order tensor in Voigt notation (6x6 for 3D)
            // Vectorized construction of Cmat (material tangent) in Voigt notation
            Matrix6d Cmat;
            Cmat.setZero();
            // Voigt index mapping arrays
            constexpr int vmap[6][2] = { {0,0}, {1,1}, {2,2}, {0,1}, {1,2}, {0,2} };
            Eigen::Matrix3d sqrt2 = Eigen::Matrix3d::Ones();
            sqrt2(0,1) = sqrt2(1,0) = sqrt2(1,2) = sqrt2(2,1) = sqrt2(0,2) = sqrt2(2,0) = std::sqrt(2.0);
            for (int I = 0; I < 6; ++I) {
                int i = vmap[I][0], j = vmap[I][1];
                for (int J = 0; J < 6; ++J) {
                    int k = vmap[J][0], l = vmap[J][1];
                    double val = lambda * Cinv(i, j) * Cinv(k, l)
                        + mu * 0.5 * (Cinv(i, k) * Cinv(j, l) + Cinv(i, l) * Cinv(j, k));
                    double scale = 1.0;
                    if (I >= 3) scale *= std::sqrt(2.0);
                    if (J >= 3) scale *= std::sqrt(2.0);
                    Cmat(I, J) = val * scale;
                }
            }

            // Geometric stiffness (initial stress part)
            // Kgeo = sum over i,j of (Gu[ip](i,:)^T * S(i,j) * Gu[ip](j,:))
            Eigen::MatrixXd Kgeo = Eigen::MatrixXd::Zero(dim * nNodes_u, dim * nNodes_u);
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j) {
                    Kgeo.block(i * nNodes_u, j * nNodes_u, nNodes_u, nNodes_u) += S(i, j) * Gu[ip].row(i).transpose() * Gu[ip].row(j);
                }
            }

            // Material tangent: Kmat = B_G^T * Cmat * B_G (total Lagrangian)
            Eigen::MatrixXd Kmat = B_G.transpose() * Cmat * B_G;

            // Total tangent
            K_UU += w[ip] * (Kmat + Kgeo);

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
bool LargeDeformationModel::SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>>& Data){

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

    // No plasticity: skip saving plastic strains

    //save total strains
    if (SaveLoc == "ip" && ElemGroup == ElemGroupIndex_u && (DataName == "exx" || DataName == "eyy" || DataName == "ezz" || DataName == "exy" || DataName == "eyz" || DataName == "exz")){
        SaveStrainComponent(DataName, Data);
        return true;
    }

    return false;
}


void LargeDeformationModel::SaveStrainComponent(std::string DataName, std::vector<std::vector<double>> &Data){
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

            // For large-strain, output Green-Lagrange strain E = 0.5*(C - I)
            // Vectorized computation of F
            Eigen::MatrixXd Umat = Eigen::Map<const Eigen::MatrixXd>(U.data(), dim, nNodes_u);
            Eigen::MatrixXd F = Eigen::MatrixXd::Identity(dim, dim) + Umat * Gu[ip].transpose();
            Eigen::MatrixXd C = F.transpose() * F;
            Eigen::MatrixXd E = 0.5 * (C - Eigen::MatrixXd::Identity(dim, dim));
            // Output requested component in Voigt order
            if (Comp == 0) Data[el][ip] = E(0,0);
            if (Comp == 1) Data[el][ip] = E(1,1);
            if (Comp == 2) Data[el][ip] = E(2,2);
            if (Comp == 3) Data[el][ip] = E(0,1);
            if (Comp == 4) Data[el][ip] = E(1,2);
            if (Comp == 5) Data[el][ip] = E(0,2);
        }
    }
}


/// @brief Saves the stress components to a file
/// @param DataName
/// @param Data
/// @param Damaged
void LargeDeformationModel::SaveStressComponent(std::string DataName, std::vector<std::vector<double>>& Data, bool Damaged){
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
    Eigen::Matrix3d Stress_Total;
    Eigen::Vector3d PriStresses;

    for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex_u].NElems; el++){
        mesh->getShapeGrads(ElemGroupIndex_u, el, El_u, Nu, Gu, w);
        dofs->getDofForNodes(El_u, dofTypes_u, dofsU);
        physics->StateVectors[Step_u].GetValues(dofsU, U);

        for (size_t ip = 0; ip < ipcount; ip++){

            //special matrices
            ConstructBMat(B, Gu[ip], dim);

            // Large-strain: output Neo-Hookean 2nd Piola-Kirchhoff stress
            // Vectorized computation of F
            Eigen::MatrixXd Umat = Eigen::Map<const Eigen::MatrixXd>(U.data(), dim, nNodes_u);
            Eigen::MatrixXd F = Eigen::MatrixXd::Identity(dim, dim) + Umat * Gu[ip].transpose();
            double J = F.determinant();
            Eigen::MatrixXd C = F.transpose() * F;
            Eigen::MatrixXd Cinv = C.inverse();
            double mu = M->Shear;
            double lambda = M->Lame;
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(dim, dim);
            Eigen::MatrixXd S = mu * (I - Cinv) + lambda * std::log(J) * Cinv;
            // Output requested component in Voigt order
            if (Pri==false){
                if (Comp == 0) Data[el][ip] = S(0,0);
                if (Comp == 1) Data[el][ip] = S(1,1);
                if (Comp == 2) Data[el][ip] = S(2,2);
                if (Comp == 3) Data[el][ip] = S(0,1);
                if (Comp == 4) Data[el][ip] = S(1,2);
                if (Comp == 5) Data[el][ip] = S(0,2);
            } else {
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(S);
                PriStresses = eigensolver.eigenvalues();
                Data[el][ip] = PriStresses(Comp);
            }

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
