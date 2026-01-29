#include "LinearElasticModel.h"
#include "../../Physics/physics.h"

void Register_LinearElasticModel(){
    ModelNames.push_back("LinearElastic");
    ModelCreators.push_back(New_LinearElasticModel);
}

BaseModel* New_LinearElasticModel(Physics& My_Physics, std::string MyNameIn){
    return new LinearElasticModel(My_Physics, MyNameIn);
}

LinearElasticModel::LinearElasticModel(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "LinearElastic";
};

LinearElasticModel::~LinearElasticModel(){

};

void LinearElasticModel::init(inputData& inputs){
    Setup(inputs);

    std::vector<size_t> UniqueNodes = mesh->ElementGroups[ElemGroupIndex].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes, dofTypes);
}

void LinearElasticModel::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);
}

void LinearElasticModel::save(SaveDataFile& data){

}

void LinearElasticModel::Setup(inputData& inputs){
    //std::cout << ModelName + " Init\n";
    std::string Groupname; inputs.GetRequired(Groupname, {"Models",MyName,"ElementGroup"});
    ElemGroupIndex = mesh->GetElementGroupIdx(Groupname);
    
    // Get and Add dofs
    std::vector<size_t> dofSteps;
    dofs->getDofTypesSteps(DofNames, dofTypes, dofSteps);
    if (dofSteps[0] != dofSteps[1]) throw std::invalid_argument(ModelName+" requires "+ DofNames[0] + " and " + DofNames[1] + " to be in the same solver step\n");
    MyStep = dofSteps[0];

    //material properties
    inputs.GetRequired(Young, {"properties","solid","Young"});
    inputs.GetRequired(Poisson, {"properties","solid","Poisson"});

    D.setZero();
    double preFac = Young/((1.0+Poisson)*(1-2*Poisson));
    D(0,0) = preFac*(1.0-Poisson);
    D(0,1) = preFac*Poisson;
    D(0,2) = preFac*Poisson;
    D(1,0) = preFac*Poisson;
    D(1,1) = preFac*(1.0-Poisson);
    D(1,2) = preFac*Poisson;
    D(2,0) = preFac*Poisson;
    D(2,1) = preFac*Poisson;
    D(2,2) = preFac*(1.0-Poisson);
    D(3,3) = 0.5*preFac*(1-2.0*Poisson);
}

void LinearElasticModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    if (step == MyStep){
        size_t ipcount = mesh->ElementGroups[ElemGroupIndex].BaseElem->ipcount;
        size_t nNodes = mesh->ElementGroups[ElemGroupIndex].NNodes_per_elem;
        std::vector<size_t> Nodes(nNodes);

        Eigen::MatrixXd K_el(2*nNodes, 2*nNodes); 
        Eigen::VectorXd F_el(2*nNodes);
        
        std::vector<Eigen::RowVectorXd> N(ipcount); for (size_t ip=0; ip<ipcount;  ip++) N[ip].resize(nNodes); 
        std::vector<Eigen::MatrixXd> G(ipcount); for (size_t ip=0; ip<ipcount;  ip++) G[ip].resize(2, nNodes); 
        std::vector<double> w(ipcount);
        Eigen::MatrixXd B(4, 2*nNodes);
        Eigen::VectorXd XY(2*nNodes);
        std::vector<size_t> dofsX(nNodes), dofsY(nNodes);
        std::vector<PetscInt> dofsXY(2*nNodes);

        for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex].NElems; el++){
            mesh->getShapeGrads(ElemGroupIndex, el, Nodes, N, G, w);

            dofs->getDofForNodes(Nodes, dofTypes[0], dofsX);
            dofs->getDofForNodes(Nodes, dofTypes[1], dofsY);
            for (size_t i = 0; i < nNodes; i++){
                dofsXY[i] = dofsX[i];
                dofsXY[i+nNodes] = dofsY[i];
            }
            physics->StateVectors[MyStep].GetValues(dofsXY, XY);

            K_el.setZero(); F_el.setZero();
            for (size_t ip = 0; ip < ipcount; ip++){
                getB(G[ip], B);
                F_el += w[ip] * B.transpose()*D*B*XY;
                K_el += w[ip] * B.transpose()*D*B;
            }
            VecAdd(f, dofsXY, F_el);
            MatAdd(K, dofsXY, dofsXY, K_el);
        } 
    }
}

void LinearElasticModel::getB(Eigen::MatrixXd &G, Eigen::MatrixXd &B){
    B.setZero();
    size_t nNodes = G.cols();
    for (size_t i = 0; i < nNodes; i++){
		//dx
		B(0, i) = G(0, i);
		B(3, i) = G(1, i);

		//dy
		B(1, i + nNodes) = G(1, i);
		B(3, i + nNodes) = G(0, i);        
    }
}

bool LinearElasticModel::SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>>& Data){
    if (SaveLoc == "ip" && ElemGroup == ElemGroupIndex && (DataName == "sxx" || DataName == "syy" || DataName == "szz" || DataName == "sxy")){
        SaveStressComponent(DataName, Data);
        return true;
    }
    return false;
}

void LinearElasticModel::SaveStressComponent(std::string DataName, std::vector<std::vector<double>>& Data){
    size_t Comp;
    if (DataName == "sxx"){
        Comp = 0;
    } 
    if (DataName == "syy"){
        Comp = 1;
    } 
    if (DataName == "szz"){
        Comp = 2;
    } 
    if (DataName == "sxy"){
        Comp = 3;
    } 

    size_t ipcount = mesh->ElementGroups[ElemGroupIndex].BaseElem->ipcount;
    size_t nNodes = mesh->ElementGroups[ElemGroupIndex].NNodes_per_elem;
    std::vector<Eigen::RowVectorXd> N(ipcount); for (size_t ip=0; ip<ipcount;  ip++) N[ip].resize(nNodes); 
    std::vector<Eigen::MatrixXd> G(ipcount); for (size_t ip=0; ip<ipcount;  ip++) G[ip].resize(2, nNodes); 
    std::vector<double> w(ipcount);
    Eigen::MatrixXd B(4, 2*nNodes);
    Eigen::MatrixXd coordsNodes(nNodes, 2);
    Eigen::VectorXd XY(2*nNodes);
    std::vector<size_t> dofsX(nNodes), dofsY(nNodes);
    std::vector<PetscInt> dofsXY(2*nNodes);
    Eigen::Vector4d Stresses;

    size_t i_el = 0;
    for (std::vector<size_t> El : mesh->ElementGroups[ElemGroupIndex].Elems){
        mesh->GetCoordsForNodes(coordsNodes, El);
        mesh->ElementGroups[ElemGroupIndex].BaseElem->getShapeGrads(N, G, w, coordsNodes);

        dofs->getDofForNodes(El, dofTypes[0], dofsX);
        dofs->getDofForNodes(El, dofTypes[1], dofsY);
        for (size_t i = 0; i < nNodes; i++){
            dofsXY[i] = dofsX[i];
            dofsXY[i+nNodes] = dofsY[i];
        }
        physics->StateVectors[MyStep].GetValues(dofsXY, XY);
        for (size_t ip = 0; ip < ipcount; ip++){
            getB(G[ip], B);
            Stresses = D*B*XY;
            Data[i_el][ip] = Stresses(Comp);
        }
        i_el += 1;
    }
}