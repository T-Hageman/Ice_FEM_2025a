#include "TriAxial.h"
#include "../../Physics/physics.h"

void Register_TriAxialModel(){
    ModelNames.push_back("TriAxial");
    ModelCreators.push_back(New_TriAxialModel);
}

BaseModel* New_TriAxialModel(Physics& My_Physics, std::string MyNameIn){
    return new TriAxialModel(My_Physics, MyNameIn);
}

TriAxialModel::TriAxialModel(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "TriAxialModel";
};

TriAxialModel::~TriAxialModel(){

};

void TriAxialModel::init(inputData& inputs){
    Setup(inputs);

    std::vector<size_t> UniqueNodes = mesh->ElementGroups[ElementGroupIndex_Side].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes, dofTypes);

    std::vector<size_t> UniqueNodes2 = mesh->ElementGroups[ElementGroupIndex_Top].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes2, dofTypes);
}

void TriAxialModel::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);
}

void TriAxialModel::Setup(inputData& inputs){
    std::string EGroupName; 
    
    inputs.GetRequired(EGroupName, {"Models",MyName,"SideGroup"});
    ElementGroupIndex_Side = mesh->GetElementGroupIdx(EGroupName);

    inputs.GetRequired(EGroupName, {"Models",MyName,"TopGroup"});
    ElementGroupIndex_Top = mesh->GetElementGroupIdx(EGroupName);
    
    // Get and Add dofs
    if (dim==2){
        DofNames_u.resize(2); //remove uz
    }

    std::vector<size_t> dofSteps(dim);
    dofs->getDofTypesSteps(DofNames_u, dofTypes, dofSteps);
    dofStep = dofSteps[0];

    //other parameters
    inputs.GetRequired(kdummy, {"Models",MyName,"dummy"});
    inputs.GetRequired(p0,{"Models",MyName,"p0"});
    inputs.GetRequired(uRate,{"Models",MyName,"dU"});
    inputs.GetRequired(t_postFail, {"Models",MyName,"t_postFail"});

    inputs.GetRequired(DamageForces,{"Models",MyName,"DamagedForces"});
    if (DamageForces){
        dofs->getDofTypesSteps(Dofname_PF, dofType_phase, dofStep_Phase);
        std::string SolidMatName;
        inputs.GetRequired(SolidMatName, {"Models", MyName, "Material"});
        PF_Util = new PhaseFieldUtility(inputs, SolidMatName);
    }

    uz_Sum = 0.0;
    F_Ext = 0.0;
    Area = 1.0;
    u0 = -1.0e8;

    tbreak = 0.0;
    MinLoadPassed = false;
    maxLoad = -1e9;
}

void TriAxialModel::save(SaveDataFile& data){

}

void TriAxialModel::Commit(int CommitType){
    if (CommitType==CommitTypes::TIMEDEP_COMMIT_TYPE){
        MPI_Allreduce(MPI_IN_PLACE, &uz_Sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &Area, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &F_Ext, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

        uz_SumOld = uz_Sum;
        F_ExtOld = F_Ext;

        maxLoad = std::max(maxLoad, F_Ext/Area*sgn(uRate));

        if (sgn(uRate)<0){//loaded downwards
            if (MinLoadPassed == false && (F_Ext/Area)*sgn(uRate)>std::max(1.0e3, std::abs(p0)) && physics->time>0.0){
                MinLoadPassed = true;
                Logs.PrintEvery("MinLoadPassed\n",2);
            }

            double MinLoad = std::max(0.05*std::abs(p0), 0.05*maxLoad);
            if (MinLoadPassed && abs(F_Ext/Area)<MinLoad && tbreak<1e-8){ //broken, remove external load from sides
                tbreak = physics->time; 
                physics->MaxTime = tbreak + t_postFail;
                Logs.PrintEvery("Broken\n",2);
            }
        } else { //loaded upwards
            double MinLoad = std::max(0.0, 0.25*maxLoad);
            if (MinLoadPassed && abs(F_Ext/Area)<MinLoad && tbreak<1e-8){ //broken, remove external load from sides
                tbreak = physics->time; 
                physics->MaxTime = tbreak + t_postFail;
                Logs.PrintEvery("Broken\n",2);
            }
            if (MinLoadPassed == false && (F_Ext/Area)*sgn(uRate)>0.0e3 && physics->time>0.0){
                MinLoadPassed = true;
                Logs.PrintEvery("MinLoadPassed\n",2);
            }  

        }

        if (physics->time<0.0){
            u0 = uz_Sum/Area;
        }
    }
}

void TriAxialModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    //side pressure
    double dt = physics->timeScheme->dt;
    double t = physics->time + dt;

    // MPI_Allreduce(MPI_IN_PLACE, &uz_Sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    // MPI_Allreduce(MPI_IN_PLACE, &Area, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    // MPI_Allreduce(MPI_IN_PLACE, &F_Ext, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

    // double MinLoad = std::max(std::abs(p0), 0.25*maxLoad);
    // if (MinLoadPassed && abs(F_Ext/Area)<MinLoad && tbreak<1e-8){ //broken, remove external load from sides
    //     tbreak = physics->time; 
    //     physics->MaxTime = tbreak + t_postFail;
    //     Logs.PrintEvery("Broken\n",2);
    // }

    if (tbreak>0.0){
        p0 = 0.0;
    }

    // Logs.PrintSingle("p_ax: "+std::to_string(F_Ext/Area*1e-6)+" MPa \n p_rad: "+std::to_string(p0*1e-6)+" MPa\n",2);

    if (dofStep == step){
        size_t ipcount = mesh->ElementGroups[ElementGroupIndex_Side].BaseElem->ipcount;
        size_t nNodes = mesh->ElementGroups[ElementGroupIndex_Side].NNodes_per_elem;
        Eigen::VectorXd F_el(dim*nNodes);
        std::vector<size_t> Nodes(nNodes);
        Eigen::MatrixXd coordsNodes(nNodes, dim), coordsIP(ipcount, dim);
        std::vector<Eigen::RowVectorXd> N(ipcount); for (size_t ip=0; ip<ipcount;  ip++) N[ip].resize(nNodes); 
        std::vector<Eigen::MatrixXd> G(ipcount); for (size_t ip=0; ip<ipcount;  ip++) G[ip].resize(dim-1, nNodes); 
        std::vector<Eigen::VectorXd> normals(ipcount); for (size_t ip=0; ip<ipcount;  ip++) normals[ip].resize(dim); 
        Eigen::MatrixXd N_uu(dim, dim * nNodes); // state to displacement mapping matrix
        std::vector<double> w(ipcount);
        Eigen::VectorXd PF(nNodes), dU(dim*nNodes);
        Eigen::MatrixXd K_el(dim*nNodes, dim*nNodes);
        Eigen::RowVectorXd dF_du(dim*nNodes);
        std::vector<PetscInt> dofsU(dim*nNodes), dofsPF(nNodes);

        double dXdt = physics->timeScheme->du_dt;
        double Force, targetX;

        for (size_t el = 0; el < mesh->ElementGroups[ElementGroupIndex_Side].NElems; el++){
            mesh->getShapeGrads(ElementGroupIndex_Side, el, Nodes, N, G, w);
            mesh->GetCoordsForNodes(coordsNodes, Nodes);
            mesh->getCoordsIP(ElementGroupIndex_Side, el, coordsIP, coordsNodes);
            mesh->getNormals(ElementGroupIndex_Side, el, normals);

            dofs->getDofForNodes(Nodes, dofTypes, dofsU);
            physics->dStateVectors[step].GetValues(dofsU, dU);
            if (DamageForces){
                dofs->getDofForNodes(Nodes, dofType_phase, dofsPF);
                physics->StateVectors[dofStep_Phase].GetValues(dofsPF, PF);
            }

            F_el.setZero(); 
            K_el.setZero();

            for (size_t ip = 0; ip < ipcount; ip++){
                double Dam = 1.0;
                if (DamageForces){
                    double phi = N[ip]*PF;
                    double unused2, unused3;
                    PF_Util->GenericDamageFunction(phi, DamageMethods::QuadraticDamageNR, Dam, unused2, unused3);
                }

                N_uu.setZero();
                N_uu(0, Eigen::seq(0, nNodes - 1)) = N[ip];
                N_uu(1, Eigen::seq(nNodes, 2 * nNodes - 1)) = N[ip];
                if (dim==3){
                    N_uu(2, Eigen::seq(2*nNodes, 3 * nNodes - 1)) = N[ip];
                }

                double f_ext;
                if (false){
                    double dUnorm = normals[ip].transpose()*N_uu*dU;
                    double k = 1.0e8;

                    f_ext = p0/(1+k*dUnorm*dUnorm);
                    dF_du = -2*k*p0*dUnorm/(1+k*dUnorm*dUnorm)/(1+k*dUnorm*dUnorm)*normals[ip].transpose()*N_uu*dXdt;
                } else {
                    f_ext = p0;
                    dF_du.setZero();
                }


                F_el += Dam*w[ip]*N_uu.transpose()*normals[ip]*f_ext;
                K_el += Dam*w[ip]*(N_uu.transpose()*normals[ip])*dF_du;
            }
            
            VecAdd(f, dofsU, F_el);
            MatAdd(K, dofsU, dofsU, K_el);
        }
    }   
    //top displacement/pressure
    if (dofStep == step){
        uz_Sum = 0.0;
        F_Ext = 0.0;
        Area = 0.0;
        size_t ipcount = mesh->ElementGroups[ElementGroupIndex_Top].BaseElem->ipcount;
        size_t nNodes = mesh->ElementGroups[ElementGroupIndex_Top].NNodes_per_elem;
        Eigen::VectorXd F_el(dim*nNodes);
        std::vector<size_t> Nodes(nNodes);
        Eigen::MatrixXd coordsNodes(nNodes, dim), coordsIP(ipcount, dim);
        std::vector<Eigen::RowVectorXd> N(ipcount); for (size_t ip=0; ip<ipcount;  ip++) N[ip].resize(nNodes); 
        std::vector<Eigen::MatrixXd> G(ipcount); for (size_t ip=0; ip<ipcount;  ip++) G[ip].resize(dim-1, nNodes); 
        std::vector<Eigen::VectorXd> normals(ipcount); for (size_t ip=0; ip<ipcount;  ip++) normals[ip].resize(dim); 
        Eigen::MatrixXd N_uu(dim, dim * nNodes); // state to displacement mapping matrix
        std::vector<double> w(ipcount);
        Eigen::VectorXd U(dim*nNodes);
        Eigen::MatrixXd K_el(dim*nNodes, dim*nNodes);
        std::vector<PetscInt> dofsU(dim*nNodes);
        Eigen::VectorXd disps(dim);

        double Force, targetX;
        for (size_t el = 0; el < mesh->ElementGroups[ElementGroupIndex_Top].NElems; el++){
            mesh->getShapeGrads(ElementGroupIndex_Top, el, Nodes, N, G, w);
            //mesh->GetCoordsForNodes(coordsNodes, Nodes);
            //mesh->getCoordsIP(ElementGroupIndex_Top, el, coordsIP, coordsNodes);
            mesh->getNormals(ElementGroupIndex_Top, el, normals);

            dofs->getDofForNodes(Nodes, dofTypes, dofsU);
            physics->StateVectors[step].GetValues(dofsU, U);

            F_el.setZero(); 
            K_el.setZero();

            for (size_t ip = 0; ip < ipcount; ip++){
                N_uu.setZero();
                N_uu(0, Eigen::seq(0, nNodes - 1)) = N[ip];
                N_uu(1, Eigen::seq(nNodes, 2 * nNodes - 1)) = N[ip];
                if (dim==3){
                    N_uu(2, Eigen::seq(2*nNodes, 3 * nNodes - 1)) = N[ip];
                }

                if (dim==2){
                    normals[ip](0)=0.0; normals[ip](1)=1.0;
                } else {
                    normals[ip](0)=0.0; normals[ip](1)=0.0; normals[ip](2)=1.0; 
                }
                if (t<0.0){
                    F_el += w[ip]*N_uu.transpose()*normals[ip]*p0;

                    F_Ext+= -w[ip]*p0;
                } else {
                     double uTarget = u0 + t*uRate;
                     double pExt = kdummy*(normals[ip].transpose()*N_uu*U - uTarget);

                     F_el += w[ip]*N_uu.transpose()*normals[ip]*pExt;
                     K_el += w[ip]*kdummy*N_uu.transpose()*normals[ip]*normals[ip].transpose()*N_uu;

                     F_Ext+= -w[ip]*pExt;
                }
                Area += w[ip];

                disps = N_uu*U;
                uz_Sum += w[ip]*disps(dim-1);
            }
            VecAdd(f, dofsU, F_el);
            MatAdd(K, dofsU, dofsU, K_el);
        }
    }   
}

size_t TriAxialModel::hasTimeData(std::vector<std::string>& DataNames){
    size_t nData = 2;
    DataNames.resize(nData);
    DataNames[0] = MyName+"/p_Ext";
    DataNames[1] = MyName+"/U_Ext";
    
    return nData;
}

void TriAxialModel::GetTimeData(std::vector<double>& DataValues){
    size_t nData = 2;
    DataValues.resize(nData);

    std::vector<double> mpiData(3); mpiData[0] = Area; mpiData[1] = F_Ext; mpiData[2] = uz_Sum;

    MPI_Allreduce(MPI_IN_PLACE, mpiData.data(), 3, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

    DataValues[0] = mpiData[1]/mpiData[0];
    DataValues[1] = mpiData[2]/mpiData[0];
}