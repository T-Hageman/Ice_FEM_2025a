#include "FloatingBC.h"
#include "../../Physics/physics.h"

// Register the model
void Register_FloatingBCModel() {
    ModelNames.push_back("FloatingBC");
    ModelCreators.push_back(New_FloatingBCModel);
}

// Create a new instance of the model
BaseModel* New_FloatingBCModel(Physics& My_Physics, std::string MyNameIn) {
    return new FloatingBCModel(My_Physics, MyNameIn);
}

// Constructor
FloatingBCModel::FloatingBCModel(Physics& My_Physics, std::string MyName) : BaseModel(My_Physics, MyName) {
    ModelName = "FloatingBC";
}

// Destructor
FloatingBCModel::~FloatingBCModel() {

}

// Load model state
void FloatingBCModel::load(inputData& inputs, SaveDataFile& data) {

}

// Save model state
void FloatingBCModel::save(SaveDataFile& data) {

}

// Initialize the model
void FloatingBCModel::init(inputData& inputs) {
    // Perform setup based on inputs
    Setup(inputs);

    std::vector<size_t> UniqueNodes = mesh->ElementGroups[ElementGroupIndex_u].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes, dofTypes_U);
}

// Setup the model
void FloatingBCModel::Setup(inputData& inputs) {
    std::string ElementGroupName;
    inputs.GetRequired(ElementGroupName, {"Models", MyName, "ElementGroup_u"});
    ElementGroupIndex_u = mesh->GetElementGroupIdx(ElementGroupName);

    // Initialize DOF types
    if (mesh->dim == 3){
        DofName_u = "uz";
    }
    std::vector<size_t> dofSteps;
    dofs->getDofTypesSteps(DofName_u, dofTypes_U, dofStep_U);

    inputs.GetRequired(hw, {"Models", MyName, "hw"});
    
    inputs.GetRequired(FluidMatName, {"Models", MyName, "Fluid"});
    Fluid = new FluidMaterial(inputs, FluidMatName);
}

// Assemble the force vector and stiffness matrix
void FloatingBCModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step) {
    if (step == dofStep_U) {
        Assemble_U(K, f);
    }
}

// Assemble the displacement part of the tangent stiffness matrix and residual force vector
void FloatingBCModel::Assemble_U(Mat& K, Vec& f){
    size_t ipcount = mesh->ElementGroups[ElementGroupIndex_u].BaseElem->ipcount;   //number of integration points
    size_t nNodes_u = mesh->ElementGroups[ElementGroupIndex_u].NNodes_per_elem;//number of nodes per pressure element

    std::vector<double> w(ipcount); //integration weights
    std::vector<Eigen::RowVectorXd> Nu(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Nu[ip].resize(nNodes_u); //pressure shape functions
    std::vector<Eigen::MatrixXd> Gu(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gu[ip].resize(1, nNodes_u); //pressure shape function gradients
    std::vector<size_t> El_u(nNodes_u); //nodes contained within this element
    std::vector<PetscInt> dofsU(nNodes_u);  //pressure field degree of freedom indices
    Eigen::VectorXd UZ(nNodes_u); //displacement field values at element nodes

    Eigen::MatrixXd coordsNodes(nNodes_u, mesh->dim);

    Eigen::VectorXd F_U(nNodes_u);    //element force vector
    Eigen::MatrixXd K_UU(nNodes_u, nNodes_u); //element stiffness matrix
    Eigen::VectorXd WLumped(nNodes_u);

    for (size_t el = 0; el < mesh->ElementGroups[ElementGroupIndex_u].NElems; el++){   //loop over all elements
        //get shape functions and gradients
        mesh->getShapeGrads(ElementGroupIndex_u, el, El_u, Nu, Gu, w);
        mesh->GetCoordsForNodes(coordsNodes, El_u);

        //obtain degrees of freedom and state vectors
        dofs->getDofForNodes(El_u, dofTypes_U, dofsU);
        physics->StateVectors[dofStep_U].GetValues(dofsU, UZ);

        F_U.setZero();
        K_UU.setZero();
        WLumped.setZero();

        for (size_t ip = 0; ip < ipcount; ip++){
            WLumped += w[ip] * Nu[ip];
        }
        for (size_t n=0; n < nNodes_u; n++){
            double hIP = coordsNodes(n, mesh->dim-1) + UZ(n) - hw; //height of the node above water, accounting for displacement
            double pTarget;
            if (hIP > 0.0){
                pTarget = 0.0;
            } else {
                pTarget = g*Fluid->Density * hIP;

                //only apply to Z component
                F_U(n)    += WLumped[n] * pTarget;
                K_UU(n,n) += WLumped[n] * g*Fluid->Density;
            }
        }

        VecAdd(f, dofsU, F_U);
        MatAdd(K, dofsU, dofsU, K_UU);
    }
}