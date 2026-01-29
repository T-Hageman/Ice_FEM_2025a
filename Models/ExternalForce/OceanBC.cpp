#include "OceanBC.h"
#include "../../Physics/physics.h"

// Register the model
void Register_OceanBCModel() {
    ModelNames.push_back("OceanBC");
    ModelCreators.push_back(New_OceanBCModel);
}

// Create a new instance of the model
BaseModel* New_OceanBCModel(Physics& My_Physics, std::string MyNameIn) {
    return new OceanBCModel(My_Physics, MyNameIn);
}

// Constructor
OceanBCModel::OceanBCModel(Physics& My_Physics, std::string MyName) : BaseModel(My_Physics, MyName) {
    ModelName = "OceanBC";
}

// Destructor
OceanBCModel::~OceanBCModel() {

}

// Load model state
void OceanBCModel::load(inputData& inputs, SaveDataFile& data) {

}

// Save model state
void OceanBCModel::save(SaveDataFile& data) {

}

// Initialize the model
void OceanBCModel::init(inputData& inputs) {
    // Perform setup based on inputs
    Setup(inputs);

    std::vector<size_t> UniqueNodes = mesh->ElementGroups[ElementGroupIndex_u].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes, dofTypes_U);

    std::vector<size_t> UniqueNodes2 = mesh->ElementGroups[ElementGroupIndex_Phase].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes2, dofType_phase);

    if (HasPW){
        std::vector<size_t> UniqueNodes3 = mesh->ElementGroups[ElementGroupIndex_p].GetUniqueNodes();
        dofs->AddDofs(UniqueNodes3, dofType_p);
    }
}

// Setup the model
void OceanBCModel::Setup(inputData& inputs) {

    HasPW = dofs->hasDofType(DofNames_p);

    std::string ElementGroupName;
    inputs.GetRequired(ElementGroupName, {"Models", MyName, "ElementGroup_u"});
    ElementGroupIndex_u = mesh->GetElementGroupIdx(ElementGroupName);

    inputs.GetRequired(ElementGroupName, {"Models", MyName, "ElementGroup_phase"});
    ElementGroupIndex_Phase = mesh->GetElementGroupIdx(ElementGroupName);

    if (HasPW){
        inputs.GetRequired(ElementGroupName, {"Models", MyName, "ElementGroup_p"});
        ElementGroupIndex_p = mesh->GetElementGroupIdx(ElementGroupName);
        dofs->getDofTypesSteps(DofNames_p, dofType_p, dofStep_p);
    }

    // Initialize DOF types
    std::vector<size_t> dofSteps;
    dofs->getDofTypesSteps(DofNames_u, dofTypes_U, dofSteps);
    dofStep_U = dofSteps[0];

    dofs->getDofTypesSteps(Dofname_PF, dofType_phase, dofStep_Phase);
    
    inputs.GetRequired(hw, {"Models", MyName, "hw"});
    
    inputs.GetRequired(FluidMatName, {"Models", MyName, "Fluid"});
    Fluid = new FluidMaterial(inputs, FluidMatName);

    inputs.GetRequired(pwAboveWater, {"Models", MyName, "pwAboveWater"});
}

// Assemble the force vector and stiffness matrix
void OceanBCModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step) {
    if (step == dofStep_U) {
        Assemble_U(K, f);
    }
    if (HasPW && step == dofStep_p) {
        Assemble_P(K, f);
    }
}

// Assemble the displacement part of the tangent stiffness matrix and residual force vector
void OceanBCModel::Assemble_U(Mat& K, Vec& f){
    size_t ipcount = mesh->ElementGroups[ElementGroupIndex_u].BaseElem->ipcount;   //number of integration points
    size_t nNodes_u = mesh->ElementGroups[ElementGroupIndex_u].NNodes_per_elem;//number of nodes per pressure element
    size_t nNodes_phase = mesh->ElementGroups[ElementGroupIndex_Phase].NNodes_per_elem;//number of nodes per phasefield element

    std::vector<double> w(ipcount); //integration weights
    std::vector<Eigen::RowVectorXd> Nu(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Nu[ip].resize(nNodes_u); //pressure shape functions
    std::vector<Eigen::MatrixXd> Gu(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gu[ip].resize(1, nNodes_u); //pressure shape function gradients
    std::vector<size_t> El_u(nNodes_u); //nodes contained within this element
    std::vector<PetscInt> dofsU(2*nNodes_u);  //pressure field degree of freedom indices

    std::vector<Eigen::RowVectorXd> Nphase(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Nphase[ip].resize(nNodes_phase); //phasefield shape functions
    std::vector<Eigen::MatrixXd> Gphase(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gphase[ip].resize(1, nNodes_phase); //pressure shape function gradients
    std::vector<size_t> El_phase(nNodes_phase); //nodes contained within this element
    std::vector<PetscInt> dofsPhase(nNodes_phase);  //phase field degree of freedom indices
    Eigen::VectorXd Phase(nNodes_phase);  //current phasefield parameter and phasefield time derivative

    Eigen::MatrixXd coordsNodes(nNodes_u, 2);

    Eigen::VectorXd F_U(2*nNodes_u);    //element force vector
    Eigen::VectorXd WLumped(nNodes_u);

    for (size_t el = 0; el < mesh->ElementGroups[ElementGroupIndex_u].NElems; el++){   //loop over all elements
        //get shape functions and gradients
        mesh->getShapeGrads(ElementGroupIndex_u, el, El_u, Nu, Gu, w);
        mesh->getShapeGrads(ElementGroupIndex_Phase, el, El_phase, Nphase, Gphase, w);
        mesh->GetCoordsForNodes(coordsNodes, El_u);

        //obtain degrees of freedom and state vectors
        dofs->getDofForNodes(El_u, dofTypes_U, dofsU);
        dofs->getDofForNodes(El_phase, dofType_phase, dofsPhase);
        physics->StateVectors[dofStep_Phase].GetValues(dofsPhase, Phase);

        F_U.setZero();
        WLumped.setZero();

        for (size_t ip = 0; ip < ipcount; ip++){
            WLumped += w[ip] * (Nphase[ip]*Phase)* Nu[ip];
        }
        for (size_t n=0; n < nNodes_u; n++){
            double hIP = hw - coordsNodes(n, 1);
            double pTarget;
            if (hIP < 0.0){
                pTarget = 0.0;
            } else {
                pTarget = g*Fluid->Density * hIP;
            }

            //only apply to X component
            F_U(n)   += WLumped[n] * pTarget;
        }

        VecAdd(f, dofsU, F_U);
    }
}

void OceanBCModel::Assemble_P(Mat& K, Vec& f){
    size_t ipcount = mesh->ElementGroups[ElementGroupIndex_p].BaseElem->ipcount;   //number of integration points
    size_t nNodes_p = mesh->ElementGroups[ElementGroupIndex_p].NNodes_per_elem;//number of nodes per pressure element

    std::vector<double> w(ipcount); //integration weights
    std::vector<Eigen::RowVectorXd> Np(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Np[ip].resize(nNodes_p); //pressure shape functions
    std::vector<Eigen::MatrixXd> Gp(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gp[ip].resize(1, nNodes_p); //pressure shape function gradients
    std::vector<size_t> El_p(nNodes_p); //nodes contained within this element
    std::vector<PetscInt> dofsP(nNodes_p);  //pressure field degree of freedom indices
    Eigen::VectorXd P(nNodes_p);

    Eigen::MatrixXd coordsNodes(nNodes_p, 2);

    Eigen::VectorXd F_P(nNodes_p);    //element force vector
    Eigen::MatrixXd K_PP(nNodes_p, nNodes_p);   //tangential matrix
    Eigen::VectorXd WLumped(nNodes_p);

    for (size_t el = 0; el < mesh->ElementGroups[ElementGroupIndex_p].NElems; el++){   //loop over all elements
        //get shape functions and gradients
        mesh->getShapeGrads(ElementGroupIndex_p, el, El_p, Np, Gp, w);
        mesh->GetCoordsForNodes(coordsNodes, El_p);

        //obtain degrees of freedom and state vectors
        dofs->getDofForNodes(El_p, dofType_p, dofsP);
        physics->StateVectors[dofStep_p].GetValues(dofsP, P);

        F_P.setZero();
        K_PP.setZero();
        WLumped.setZero();

        for (size_t ip = 0; ip < ipcount; ip++){
            WLumped += w[ip] * Np[ip];
        }
        for (size_t n=0; n < nNodes_p; n++){
            double hIP = hw - coordsNodes(n, 1);
            double pTarget;
            if (hIP < 0.0){
                pTarget = pwAboveWater;
            } else {
                pTarget = g*Fluid->Density * hIP;
            }

            F_P(n)   += kDummy * WLumped[n] * (P(n) - pTarget);
            K_PP(n,n)+= kDummy * WLumped[n];
        }

        VecAdd(f, dofsP, F_P);
        MatAdd(K, dofsP, dofsP, K_PP);
    }
}