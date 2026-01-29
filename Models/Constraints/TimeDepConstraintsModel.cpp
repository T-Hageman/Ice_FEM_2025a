#include "TimeDepConstraintsModel.h"
#include "../../Physics/physics.h"

void Register_TimeDepConstraintsModel(){
    ModelNames.push_back("TimeDepConstraints");
    ModelCreators.push_back(New_TimeDepConstraintsModel);
}

BaseModel* New_TimeDepConstraintsModel(Physics& My_Physics, std::string MyNameIn){
    return new TimeDepConstraintsModel(My_Physics, MyNameIn);
}

TimeDepConstraintsModel::TimeDepConstraintsModel(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "TimeDepConstraints";
};

TimeDepConstraintsModel::~TimeDepConstraintsModel(){

};

void TimeDepConstraintsModel::init(inputData& inputs){
    Setup(inputs);
    dofs->AddDofs(mesh->NodeGroups[NodeGroupIndex].Nodes, dofTypes);
}

void TimeDepConstraintsModel::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);

}

void TimeDepConstraintsModel::save(SaveDataFile& data){
    
}

void TimeDepConstraintsModel::Setup(inputData& inputs){
    std::string NodeGroupName; inputs.GetRequired(NodeGroupName, {"Models", MyName.c_str(), "NodeGroup"});
    NodeGroupIndex = physics->mesh->GetNodeGroupIdx(NodeGroupName);
    
    // Get dofs
    inputs.GetRequired(dofNames, {"Models",MyName,"dofs"});
    inputs.GetRequired(conVals, {"Models",MyName,"values"});
    inputs.GetRequired(dconVals, {"Models",MyName,"dvalues"});

    dofs->getDofTypesSteps(dofNames, dofTypes, dofSteps);
}

void TimeDepConstraintsModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    for (size_t i = 0; i < dofNames.size(); i++){
        if (dofSteps[i] == step){
            double time = physics->time;
            double dt = physics->timeScheme->dt;
            double value = conVals[i] + dconVals[i]*time;

            std::vector<size_t> dofIDs(mesh->NodeGroups[NodeGroupIndex].NNodes);
            dofs->getDofForNodes(mesh->NodeGroups[NodeGroupIndex].Nodes, dofTypes[i], dofIDs);
            cons->AddConstraint(step, dofIDs, value);
        }
    }
}