#include "BasalFriction.h"
#include "../../Physics/physics.h"

void Register_BasalFrictionModel(){
    ModelNames.push_back("BasalFriction");
    ModelCreators.push_back(New_BasalFrictionModel);
}

BaseModel* New_BasalFrictionModel(Physics& My_Physics, std::string MyNameIn){
    return new BasalFrictionModel(My_Physics, MyNameIn);
}

BasalFrictionModel::BasalFrictionModel(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "BasalFriction";
};

BasalFrictionModel::~BasalFrictionModel(){

};

void BasalFrictionModel::init(inputData& inputs){
    Setup(inputs);

    std::vector<size_t> UniqueNodes = mesh->ElementGroups[ElementGroupIndex_u].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes, dofTypes);

    if (Damage){
        std::vector<size_t> UniqueNodes2 = mesh->ElementGroups[ElementGroupIndex_phase].GetUniqueNodes();
        dofs->AddDofs(UniqueNodes2, dofType_phase);     
    }
}

void BasalFrictionModel::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);
}

void BasalFrictionModel::Setup(inputData& inputs){
    std::string EGroupName; inputs.GetRequired(EGroupName, {"Models",MyName,"ElementGroup_u"});
    ElementGroupIndex_u = mesh->GetElementGroupIdx(EGroupName);
    
    // Get and Add dofs
    std::vector<size_t> Steps_U(2);
    dofs->getDofTypesSteps(DofNames_u, dofTypes, Steps_U);
    dofSteps = Steps_U[0];

    //other parameters
    inputs.GetRequired(friction, {"Models",MyName,"friction"});
    inputs.GetRequired(k, {"Models",MyName,"k"});

    Lumped = false;
    inputs.GetOptional(Lumped, {"Models",MyName,"Lumped"});

    if (Damage){
        inputs.GetRequired(EGroupName, {"Models",MyName,"ElementGroup_phase"});
        ElementGroupIndex_phase = mesh->GetElementGroupIdx(EGroupName);
        dofs->getDofTypesSteps(DofNames_phase, dofType_phase, dofStep_phase);
    }

    basalStrength = 1.0e20;
    inputs.GetOptional(basalStrength, {"Models",MyName,"BasalStrength"});

    u0 = 5e-5;
    inputs.GetOptional(u0, {"Models",MyName,"u0"});
}

void BasalFrictionModel::save(SaveDataFile& data){

}

void BasalFrictionModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    if (dofSteps == step){
        size_t ipcount = mesh->ElementGroups[ElementGroupIndex_u].BaseElem->ipcount;
        size_t nNodes = mesh->ElementGroups[ElementGroupIndex_u].NNodes_per_elem;
        size_t nNodesPhase = mesh->ElementGroups[ElementGroupIndex_phase].NNodes_per_elem;

        std::vector<size_t> Nodes(nNodes), NodesPhase(nNodesPhase);
        Eigen::MatrixXd coordsNodes(nNodes, 2);
        std::vector<Eigen::RowVectorXd> N(ipcount); for (size_t ip=0; ip<ipcount;  ip++) N[ip].resize(nNodes); 
        std::vector<Eigen::MatrixXd> G(ipcount); for (size_t ip=0; ip<ipcount;  ip++) G[ip].resize(1, nNodes); 

        std::vector<Eigen::RowVectorXd> Nph(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Nph[ip].resize(nNodesPhase); 
        std::vector<Eigen::MatrixXd> Gph(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gph[ip].resize(1, nNodesPhase); 

        std::vector<double> w(ipcount);
        Eigen::MatrixXd Nu(2, 2*nNodes); 

        Eigen::VectorXd dU(2*nNodes), U(2*nNodes), UOld(2*nNodes), Phase(nNodesPhase);
        Eigen::VectorXd WLumped(nNodes);

        Eigen::VectorXd F_el(2*nNodes);
        Eigen::MatrixXd K_el(2*nNodes, 2*nNodes);
        std::vector<PetscInt> dofsU(2*nNodes), dofsPhase(nNodesPhase);
        double dXdt = physics->timeScheme->du_dt;

        for (size_t el = 0; el < mesh->ElementGroups[ElementGroupIndex_u].NElems; el++){
            mesh->getShapeGrads(ElementGroupIndex_u, el, Nodes, N, G, w);
            mesh->GetCoordsForNodes(coordsNodes, Nodes);

            dofs->getDofForNodes(Nodes, dofTypes, dofsU);
            physics->StateVectors[step].GetValues(dofsU, U);
            physics->StateVectorsOld[step].GetValues(dofsU, UOld);
            physics->dStateVectors[step].GetValues(dofsU, dU);

            if (Damage){
                mesh->getShapeGrads(ElementGroupIndex_phase, el, NodesPhase, Nph, Gph, w);
                dofs->getDofForNodes(NodesPhase, dofType_phase, dofsPhase);
                physics->dStateVectors[dofStep_phase].GetValues(dofsPhase, Phase);
            }

            F_el.setZero(); 
            K_el.setZero();
            WLumped.setZero();
            for (size_t ip = 0; ip < ipcount; ip++){
                double d = 1.0;
                if (Damage){
                    double phi = Nph[ip]*Phase;
                    if (phi<0.8){
                        d = 1.0;//std::pow(1.0-std::min(0.0, std::max(1.0, phi)), 2);
                    } else {
                        d = 0.0;
                    }
                }
                if (Lumped == false){
                    Nu.setZero();
                    Nu(0,Eigen::seq(0,nNodes-1)) = N[ip];
                    Nu(1,Eigen::seq(nNodes,2*nNodes-1)) = N[ip];
                    
                    Eigen::Vector2d norm, tangent;
                    norm[0] = 0.0; norm[1] = 1.0;
                    tangent[0] = 1.0; tangent[1] = 0.0;

                    // No-Pen
                    double Fy = k * (norm.transpose()*Nu*U-0.0);
                    F_el += w[ip] * Nu.transpose() * norm * Fy ;
                    K_el += w[ip] * k * Nu.transpose() * norm * norm.transpose() * Nu;

                    //Friction
                    if (Fy>0.0){
                        double Fx = friction * Fy * d * tangent.transpose()*Nu*dU;

                        F_el += w[ip] * Nu.transpose() * tangent * Fx;
                        K_el += w[ip] * friction * d * Nu.transpose() * tangent * (Fy*tangent.transpose() * Nu + (tangent.transpose() * Nu * dU) * k * norm.transpose()*Nu);
                    }
                } else {
                    WLumped += w[ip] * d * N[ip];
                }
            }
            if (Lumped == true){
                for (size_t n = 0; n < nNodes; n++){
                    size_t ny = n + nNodes;

                    Nu.setZero();
                    Nu(0,n) = 1.0;
                    Nu(1,ny)= 1.0;

                    Eigen::Vector2d norm, tangent;
                    norm[0] = 0.0; norm[1] = 1.0;
                    tangent[0] = 1.0; tangent[1] = 0.0;

                    // No-Pen
                    double Fy = k * (0.0-norm.transpose()*Nu*U);
                    double dFy = k;
                    if (-Fy>basalStrength){
                        Fy = -basalStrength;
                        dFy = 0.0;
                    }

                    F_el += -WLumped[n] * Nu.transpose() * norm * Fy ;
                    K_el += -WLumped[n] * -dFy * Nu.transpose() * norm * norm.transpose() * Nu;

                    //Friction
                    if (true){
                        double FyOld = k *  (0.0-norm.transpose()*Nu*UOld);
                        double Velocity = tangent.transpose()*Nu*dU;
                        if (FyOld>0.0){
                            double Fx, dFx_dVel;
                            if(false){
                                Fx= -friction * FyOld * Velocity;
                                dFx_dVel = -friction * FyOld;
                            } else {
                                Fx = -friction*FyOld/(1+std::abs(Velocity)/u0)*Velocity;
                                dFx_dVel = -friction*FyOld*(1-std::abs(Velocity)/u0/(1+std::abs(Velocity)/u0))/(1+std::abs(Velocity)/u0);
                            }

                            F_el += -WLumped[n] * Nu.transpose() * tangent * Fx;
                            K_el += -WLumped[n] * Nu.transpose() * tangent * dFx_dVel * (tangent.transpose() * Nu * dXdt );
                        }
                    } else {
                        if (Fy>0.0){
                            double Fx = -friction * Fy * tangent.transpose()*Nu*dU;

                            F_el += -WLumped[n] * Nu.transpose() * tangent * Fx;
                            K_el += -WLumped[n] * -friction * Nu.transpose() * tangent * (Fy*tangent.transpose() * Nu * dXdt + (tangent.transpose() * Nu * dU) * -k * norm.transpose()*Nu);
                        }
                    }
                }
            }
            
            VecAdd(f, dofsU, F_el);
            MatAdd(K, dofsU, dofsU, K_el);
        }
    }
}
