#include "Initializations.h"
#include <algorithm>

/// @brief Sets initial conditions based on input parameter "nonlinSolver.Initialization.Type"
/// @param inputs Input properties object
/// @param MyPhysics pointer to physics object
void Initialize(inputData& inputs, Physics& MyPhysics){
    std::string InitType; inputs.GetRequired(InitType,{"nonlinSolver","Initialization","Type"});

    if (InitType == "Zero") Initialize_Zero(MyPhysics);
    if (InitType == "PartSolid") Initialize_PartSolid(inputs, MyPhysics);
    if (InitType == "ElectroChemistry") Initialize_ElectroChemistry(inputs, MyPhysics);
    if (InitType == "Glacier") Initialize_Glacier(inputs, MyPhysics);
    if (InitType == "DamageChallenge") Initialize_DamageChallenge(inputs, MyPhysics);
    if (InitType == "Defect") Initialize_Defect(inputs, MyPhysics);
    if (InitType == "InterfaceCrack") Initialize_InterfaceCrack(inputs, MyPhysics);


    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        MyPhysics.StateVectors[i].AssemblyStart();
        MyPhysics.dStateVectors[i].AssemblyStart();
        MyPhysics.ddStateVectors[i].AssemblyStart();
    }
    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        MyPhysics.StateVectors[i].AssemblyEnd();
        MyPhysics.dStateVectors[i].AssemblyEnd();
        MyPhysics.ddStateVectors[i].AssemblyEnd();
        MyPhysics.StateVectors[i].SyncStart(INSERT_VALUES);
        MyPhysics.dStateVectors[i].SyncStart(INSERT_VALUES);
        MyPhysics.ddStateVectors[i].SyncStart(INSERT_VALUES);
    } 
    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        MyPhysics.StateVectors[i].SyncEnd(INSERT_VALUES);
        MyPhysics.dStateVectors[i].SyncEnd(INSERT_VALUES);
        MyPhysics.ddStateVectors[i].SyncEnd(INSERT_VALUES);
    }
}

/// @brief Initializer to set properties to zero (called from Initialize, do not call directly)
/// @param MyPhysics pointer to physics object
void Initialize_Zero(Physics& MyPhysics){
    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        VecZeroEntries(MyPhysics.StateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.dStateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.ddStateVectors[i].DataVector);
    }
}

/// @brief Initializer to set phasefield in part of the domain (called from Initialize, do not call directly)
/// @param inputs Input properties object
/// @param MyPhysics pointer to physics object
void Initialize_PartSolid(inputData& inputs, Physics& MyPhysics){
    std::string phaseGroup;
    inputs.GetRequired(phaseGroup, {"nonlinSolver","Initialization","PhaseGroup"});
    size_t nDofs = MyPhysics.dofspace->nDofs;
    std::vector<std::string> GroupNames(nDofs);

    double Lx; inputs.GetRequired(Lx, {"nonlinSolver","Initialization","Lx"});
    double Ly; inputs.GetRequired(Ly, {"nonlinSolver","Initialization","Ly"});

    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        VecZeroEntries(MyPhysics.StateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.dStateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.ddStateVectors[i].DataVector);
    }
    size_t phaseStep, phaseType;
    size_t phaseGroupIdx = MyPhysics.mesh->GetElementGroupIdx(phaseGroup);
    MyPhysics.dofspace->getDofTypesSteps("phase", phaseType, phaseStep);
    std::vector<size_t> Nodes = MyPhysics.mesh->ElementGroups[phaseGroupIdx].GetUniqueNodes();
    std::vector<double> NodeValues(Nodes.size());
    std::vector<PetscInt> NodeDofs(Nodes.size());
    MyPhysics.dofspace->getDofForNodes(Nodes, phaseType, NodeDofs);

    for (size_t i = 0; i < Nodes.size(); i++){
        double x = MyPhysics.mesh->Xcoords.GetValue(Nodes[i]);
        double y = MyPhysics.mesh->Ycoords.GetValue(Nodes[i]);

    //    double dx = Lx-x;
    //    double dy = Ly-y;
    //    double l = 0.1;
    //    NodeValues[i] = std::max(0.0,std::min(1.0,dx/l))*std::max(0.0,std::min(1.0,dy/l));

        if (x<Lx && y<Ly){
            NodeValues[i] = 1.0;
        } else {
            NodeValues[i] = 0.0;
        }
    }
    
    MyPhysics.StateVectors[phaseStep].Set(NodeDofs, NodeValues, INSERT_VALUES);
}

void Initialize_InterfaceCrack(inputData& inputs, Physics& MyPhysics){
    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        VecZeroEntries(MyPhysics.StateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.dStateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.ddStateVectors[i].DataVector);
    }

    std::string CrackGroup; inputs.GetRequired(CrackGroup,{"nonlinSolver","Initialization","pGroup"});
    size_t CrackGroupIdx = MyPhysics.mesh->GetNodeGroupIdx(CrackGroup);
    std::vector<size_t> Nodes = MyPhysics.mesh->NodeGroups[CrackGroupIdx].Nodes;
    size_t DofTypes, dofSteps;
    MyPhysics.dofspace->getDofTypesSteps("p", DofTypes, dofSteps);

    double p0 = 0.0; inputs.GetRequired(p0, {"properties","Water","FractureFlow","ReferencePressure"});
    p0 = -4.0*p0;

    std::vector<double> NodeValues(Nodes.size());
    std::vector<PetscInt> NodeDofs(Nodes.size());
    MyPhysics.dofspace->getDofForNodes(Nodes, DofTypes, NodeDofs);
    
    for (size_t j = 0; j < NodeDofs.size(); j++){
        NodeValues[j] = p0;
    }
    MyPhysics.StateVectors[dofSteps].Set(NodeDofs, NodeValues, INSERT_VALUES);

    
}

void Initialize_ElectroChemistry(inputData& inputs, Physics& MyPhysics){
    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        VecZeroEntries(MyPhysics.StateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.dStateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.ddStateVectors[i].DataVector);
    }


    bool ChemPotBased; inputs.GetRequired(ChemPotBased, {"properties","Species","UseChemPotentials"});
    size_t nSpecies;
    std::vector<std::string> Species; inputs.GetRequired(Species, {"properties","Species","Types"});
    std::vector<double> C0, z;
    int i_from_neutrality = -1;
    double Z_Filler = 0;


    nSpecies = Species.size();
    C0.resize(nSpecies);
    z.resize(nSpecies);
    for (size_t i = 0; i < nSpecies; i++){
        inputs.GetRequired(z[i], {"properties","Species",Species[i],"z"});
        std::string C0Type = inputs.GetType({"properties","Species",Species[i].c_str(),"C0"});
        if (C0Type == "Number"){
            inputs.GetRequired(C0[i], {"properties","Species",Species[i],"C0"});
            if (ChemPotBased && C0[i]<std::exp(-30.0)){
                C0[i] = std::exp(-30.0);
            }
            Z_Filler += z[i]*C0[i];
        } else {
            if (i_from_neutrality<0){
                i_from_neutrality = i;
            } else {
                throw std::invalid_argument("Multiple species are designated as initialized from electroneutrality\n");
            }
        }
    }
    if (i_from_neutrality>=0){
        C0[i_from_neutrality] = -Z_Filler/z[i_from_neutrality];
        if (C0[i_from_neutrality]<1.0e-20){
            throw std::invalid_argument("Initial concentrations do not follow from electroneutrality\n");
        }
    }

    std::string InteriorGroup; inputs.GetRequired(InteriorGroup,{"nonlinSolver","Initialization","InteriorGroup"});
    size_t InteriorGroupIdx = MyPhysics.mesh->GetNodeGroupIdx(InteriorGroup);

    std::vector<size_t> DofTypes(nSpecies), dofSteps(nSpecies);
    MyPhysics.dofspace->getDofTypesSteps(Species, DofTypes, dofSteps);
    std::vector<size_t> Nodes = MyPhysics.mesh->NodeGroups[InteriorGroupIdx].Nodes;
    std::vector<double> NodeValues(Nodes.size());
    std::vector<PetscInt> NodeDofs(Nodes.size());

    for (size_t i = 0; i < nSpecies; i++){
        MyPhysics.dofspace->getDofForNodes(Nodes, DofTypes[i], NodeDofs);
        double cCon;
        if (ChemPotBased){
            cCon = std::log(C0[i]);
        } else {
            cCon = C0[i];
        }
        
        for (size_t j = 0; j < Nodes.size(); j++){
            NodeValues[j] = cCon;
        }
        MyPhysics.StateVectors[dofSteps[i]].Set(NodeDofs, NodeValues, INSERT_VALUES);
    }
}


void Initialize_Glacier(inputData& inputs, Physics& MyPhysics){
    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        VecZeroEntries(MyPhysics.StateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.dStateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.ddStateVectors[i].DataVector);
    }

    double HFirn; inputs.GetRequired(HFirn, {"nonlinSolver","Initialization","HFirn"});

    std::string InteriorGroup; 
    size_t DofType, dofStep;
    size_t InteriorGroupIdx ;
    double HMax, HMin; 
    VecMax(MyPhysics.mesh->Ycoords.DataVector, NULL, &HMax);
    VecMin(MyPhysics.mesh->Ycoords.DataVector, NULL, &HMin);

    //pressure
    std::string pName = "pw";
    if (MyPhysics.dofspace->hasDofType(pName)){
        double p0; inputs.GetRequired(p0, {"nonlinSolver","Initialization","p0"});
        inputs.GetRequired(InteriorGroup,{"nonlinSolver","Initialization","InteriorGroup_p"});
        InteriorGroupIdx = MyPhysics.mesh->GetElementGroupIdx(InteriorGroup);

        MyPhysics.dofspace->getDofTypesSteps(pName, DofType, dofStep);
        std::vector<size_t> Nodes = MyPhysics.mesh->ElementGroups[InteriorGroupIdx].GetUniqueNodes();

        std::vector<double> NodeValues(Nodes.size());
        std::vector<PetscInt> NodeDofs(Nodes.size());
        MyPhysics.dofspace->getDofForNodes(Nodes, DofType, NodeDofs);
        for (size_t j = 0; j < Nodes.size(); j++){
            NodeValues[j] = p0;
        }
        MyPhysics.StateVectors[dofStep].Set(NodeDofs, NodeValues, INSERT_VALUES);
    }

    //temperature
    std::string TName = "T";
    double TSurf, TBase;
    inputs.GetRequired(TSurf, {"nonlinSolver","Initialization","TSurf"});
    inputs.GetRequired(TBase, {"nonlinSolver","Initialization","TBase"});

    inputs.GetRequired(InteriorGroup,{"nonlinSolver","Initialization","InteriorGroup_T"});
    InteriorGroupIdx = MyPhysics.mesh->GetElementGroupIdx(InteriorGroup);
    MyPhysics.dofspace->getDofTypesSteps(TName, DofType, dofStep);
    std::vector<size_t> Nodes_T = MyPhysics.mesh->ElementGroups[InteriorGroupIdx].GetUniqueNodes();

    std::vector<double> NodeValues_T(Nodes_T.size());
    std::vector<PetscInt> NodeDofs_T(Nodes_T.size());  
    MyPhysics.dofspace->getDofForNodes(Nodes_T, DofType, NodeDofs_T);

    for (size_t j = 0; j < Nodes_T.size(); j++){
        Eigen::VectorXd coords(2);
        MyPhysics.mesh->GetCoordsForNodes(coords, Nodes_T[j]);
        if (coords(1)>HMax-HFirn){
            NodeValues_T[j] = TBase + (TSurf-TBase)*(coords(1)-(HMax-HFirn))/HFirn;
        } else {
            NodeValues_T[j] = TBase;
        }
    }
    MyPhysics.StateVectors[dofStep].Set(NodeDofs_T, NodeValues_T, INSERT_VALUES);


    //Porosity
    std::string PorName = "poros";
    bool VariablePoros = MyPhysics.dofspace->hasDofType(PorName);
    if (VariablePoros){
        double poros_ice, poros_firn;
        inputs.GetRequired(poros_ice, {"nonlinSolver","Initialization","poros_ice"});
        inputs.GetRequired(poros_firn, {"nonlinSolver","Initialization","poros_firn"});

        MyPhysics.dofspace->getDofTypesSteps(PorName, DofType, dofStep);

        std::vector<size_t> Nodes_poros = MyPhysics.mesh->ElementGroups[InteriorGroupIdx].GetUniqueNodes();
        std::vector<double> NodeValues_poros(Nodes_poros.size());
        std::vector<PetscInt> NodeDofs_poros(Nodes_poros.size());  
        MyPhysics.dofspace->getDofForNodes(Nodes_poros, DofType, NodeDofs_poros);

        for (size_t j = 0; j < Nodes_poros.size(); j++){
            Eigen::VectorXd coords(2);
            MyPhysics.mesh->GetCoordsForNodes(coords, Nodes_poros[j]);
            if (coords(1)>HMax-HFirn){
                NodeValues_poros[j] = poros_ice+(poros_firn-poros_ice)*(coords(1)-(HMax-HFirn))/HFirn;
            } else {
                NodeValues_poros[j] = poros_ice;
            }
        }
        MyPhysics.StateVectors[dofStep].Set(NodeDofs_poros, NodeValues_poros, INSERT_VALUES);

    }

    //Phasefield
    double LFrac, ell;
    inputs.GetRequired(LFrac, {"nonlinSolver","Initialization","InitDepth"});
    inputs.GetRequired(ell, {"properties","Ice","FractureProperties","l"});
    std::string PhaseName = "phase";
    MyPhysics.dofspace->getDofTypesSteps(PhaseName, DofType, dofStep);

    inputs.GetRequired(InteriorGroup,{"nonlinSolver","Initialization","InteriorGroup_phase"});
    InteriorGroupIdx = MyPhysics.mesh->GetElementGroupIdx(InteriorGroup);

    std::vector<size_t> Nodes_phase = MyPhysics.mesh->ElementGroups[InteriorGroupIdx].GetUniqueNodes();
    std::vector<double> NodeValues_phase(Nodes_phase.size());
    std::vector<PetscInt> NodeDofs_phase(Nodes_phase.size());  
    MyPhysics.dofspace->getDofForNodes(Nodes_phase, DofType, NodeDofs_phase);

    for (size_t j = 0; j < Nodes_phase.size(); j++){
        Eigen::VectorXd coords(2);
        MyPhysics.mesh->GetCoordsForNodes(coords, Nodes_phase[j]);
        NodeValues_phase[j] = 0.0;
        if (coords(1)>HMax-LFrac){
            double dx = coords(0);
            //NodeValues_phase[j] = std::exp(-(std::abs(dx))/ell);
            if (std::abs(dx)<ell) NodeValues_phase[j] = 1.0;
        }
    }
    MyPhysics.StateVectors[dofStep].Set(NodeDofs_phase, NodeValues_phase, INSERT_VALUES);
}

void Initialize_DamageChallenge(inputData& inputs, Physics& MyPhysics){
    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        VecZeroEntries(MyPhysics.StateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.dStateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.ddStateVectors[i].DataVector);
    }

    size_t DofType, dofStep;
 
    //Phasefield
    double ell;
    std::vector<double> Notch_x(2), Notch_y(2);
    inputs.GetRequired(Notch_x, {"nonlinSolver","Initialization","Notch_x"});
    inputs.GetRequired(Notch_y, {"nonlinSolver","Initialization","Notch_y"});
    inputs.GetRequired(ell, {"properties","PhaseField","l"});

    std::string PhaseName = "phase";
    MyPhysics.dofspace->getDofTypesSteps(PhaseName, DofType, dofStep);

    std::string InteriorGroup; inputs.GetRequired(InteriorGroup,{"nonlinSolver","Initialization","InteriorGroup"});
    size_t InteriorGroupIdx = MyPhysics.mesh->GetElementGroupIdx(InteriorGroup);

    std::vector<size_t> Nodes_phase = MyPhysics.mesh->ElementGroups[InteriorGroupIdx].GetUniqueNodes();
    std::vector<double> NodeValues_phase(Nodes_phase.size());
    std::vector<PetscInt> NodeDofs_phase(Nodes_phase.size());  
    MyPhysics.dofspace->getDofForNodes(Nodes_phase, DofType, NodeDofs_phase);

    for (size_t j = 0; j < Nodes_phase.size(); j++){
        Eigen::VectorXd coords(2);
        MyPhysics.mesh->GetCoordsForNodes(coords, Nodes_phase[j]);
        NodeValues_phase[j] = 0.0;

        if (coords(0)>Notch_x[0] && coords(0)<Notch_x[1] && coords(1)>Notch_y[0] && coords(1)<Notch_y[1]){
            NodeValues_phase[j] = 1.0;
        }
    }
    MyPhysics.StateVectors[dofStep].Set(NodeDofs_phase, NodeValues_phase, INSERT_VALUES);
}

void Initialize_Defect(inputData& inputs, Physics& MyPhysics){
    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        VecZeroEntries(MyPhysics.StateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.dStateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.ddStateVectors[i].DataVector);
    }

    // phasefield defect
    std::string PhaseName = "phase";
    size_t DofType, dofStep;
    MyPhysics.dofspace->getDofTypesSteps(PhaseName, DofType, dofStep);

    std::string InteriorGroup; inputs.GetRequired(InteriorGroup,{"nonlinSolver","Initialization","InteriorGroup"});
    size_t InteriorGroupIdx = MyPhysics.mesh->GetElementGroupIdx(InteriorGroup);

    std::vector<size_t> Nodes_phase = MyPhysics.mesh->ElementGroups[InteriorGroupIdx].GetUniqueNodes();
    std::vector<double> NodeValues_phase(Nodes_phase.size());
    std::vector<PetscInt> NodeDofs_phase(Nodes_phase.size());  
    MyPhysics.dofspace->getDofForNodes(Nodes_phase, DofType, NodeDofs_phase);

    std::string MatName;
    inputs.GetRequired(MatName,{"nonlinSolver","Initialization","Material"});
    double ell;
    double Z; inputs.GetRequired(Z,{"nonlinSolver","Initialization","Z"});
    double Y; inputs.GetRequired(Y,{"nonlinSolver","Initialization","Y"});
    double R; inputs.GetRequired(R,{"nonlinSolver","Initialization","R"});
    inputs.GetRequired(ell, {"properties",MatName,"FractureProperties","l"});
    for (size_t j = 0; j < Nodes_phase.size(); j++){
        if (MyPhysics.mesh->dim==3){
            Eigen::VectorXd coords(3);
            MyPhysics.mesh->GetCoordsForNodes(coords, Nodes_phase[j]);
            NodeValues_phase[j] = 0.0;

            if ((coords(0)+Y)*(coords(0)+Y)+coords(1)*coords(1)+(coords(2)-Z)*(coords(2)-Z)<R*R*ell*ell){
                NodeValues_phase[j] = 1.0;
            }
        } else {
            Eigen::VectorXd coords(2);
            MyPhysics.mesh->GetCoordsForNodes(coords, Nodes_phase[j]);
            NodeValues_phase[j] = 0.0;

            if ((coords(0)+Y)*(coords(0)+Y)+(coords(1)-Z)*(coords(1)-Z)<R*R*ell*ell){
                NodeValues_phase[j] = 1.0;
            }
        }
    }
    MyPhysics.StateVectors[dofStep].Set(NodeDofs_phase, NodeValues_phase, INSERT_VALUES);

    //Temperature
    std::string TName = "T";
    if (MyPhysics.dofspace->hasDofType(TName)){

        MyPhysics.dofspace->getDofTypesSteps(TName, DofType, dofStep);

        std::vector<size_t> Nodes_T = MyPhysics.mesh->ElementGroups[InteriorGroupIdx].GetUniqueNodes();
        std::vector<double> NodeValues_T(Nodes_T.size());
        std::vector<PetscInt> NodeDofs_T(Nodes_T.size());  
        MyPhysics.dofspace->getDofForNodes(Nodes_T, DofType, NodeDofs_T);

        double T0 = 0.0; inputs.GetOptional(T0,{"nonlinSolver","Initialization","T0"});
        for (size_t j = 0; j < Nodes_T.size(); j++){
            NodeValues_T[j] = T0;
        }
        MyPhysics.StateVectors[dofStep].Set(NodeDofs_T, NodeValues_T, INSERT_VALUES);
    }
}