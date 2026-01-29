#include "DofSpace.h"
#include "../utility/utility.h"

#include <petsc.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <highfive/H5File.hpp>

DofSpace::DofSpace(inputData& inputs) {
    PetscCallThrow(MPI_Comm_size(PETSC_COMM_WORLD,&size));
    PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));

    inputs.GetRequired(DofNames, {"Dofs","DofNames"});
    inputs.GetRequired(DofStep, {"Dofs","DofStep"});
    nDofs = DofNames.size();
    if (nDofs!=DofStep.size()){
        throw std::invalid_argument("Please specify dof steps for all dofs");
    }

    maxSteps = *std::max_element(DofStep.begin(), DofStep.end())+1;
    DofNumbering.resize(maxSteps);
    DofsToAdd.resize(nDofs);
    GhostDofs.resize(maxSteps);
    LocalDofNumRange.resize(maxSteps);
    TotalDofRange.resize(maxSteps);

    inputs.GetRequired(FilenamePrefix,{"Outputs","SaveFolder"});
}

DofSpace::~DofSpace() {

}

/// @brief Saves dofspace to restartable data file
/// @param data 
void DofSpace::save(SaveDataFile& data){
    data.Save("LocalDofNumRange", LocalDofNumRange);
    data.Save("TotalDofRange", TotalDofRange);
    data.Save("DofNumbering", DofNumbering);
    data.Save("GhostDofs", GhostDofs);
}

/// @brief Load dofspace from restartable data file
/// @param data 
void DofSpace::load(SaveDataFile& data){
    data.Load("LocalDofNumRange", LocalDofNumRange);
    data.Load("TotalDofRange", TotalDofRange);
    data.Load("DofNumbering", DofNumbering);
    data.Load("GhostDofs", GhostDofs);
}

/// @brief For a provided dof name, return the index and step in which they are resolved
/// @param names input: Strings to identify the type of degree of freedom
/// @param DofTypes_out output: Index of the degree of freedom type
/// @param DofSteps_out output: Step in which degree of freedom is resolved
void DofSpace::getDofTypesSteps(std::vector<std::string> names, std::vector<size_t> &DofTypes_out, std::vector<size_t> &DofSteps_out){
    DofTypes_out.resize(names.size());
    DofSteps_out.resize(names.size());

    for (size_t i = 0; i < names.size(); i++){
        getDofTypesSteps(names[i], DofTypes_out[i], DofSteps_out[i]);
    }
    return;
}

/// @brief For a provided dof name, return the index and step in which they are resolved
/// @param name input: Strings to identify the type of degree of freedom
/// @param DofType_out output: Index of the degree of freedom type
/// @param DofStep_out output: Step in which degree of freedom is resolved
void DofSpace::getDofTypesSteps(std::string name, size_t &DofType_out, size_t &DofStep_out){
    bool HasDOF = hasDofType(name, DofType_out, DofStep_out);
    if (HasDOF == false){
        throw std::invalid_argument("Dof of type \"" + name + "\" is not defined");
    }  
}

/// @brief Check if degree of freedom exists, if it does, provide info about it
/// @param name input: Strings to identify the type of degree of freedom
/// @param DofType_out output: Index of the degree of freedom type (only defined if output = true)
/// @param DofStep_out output: Step in which degree of freedom is resolved (only defined if output = true)
/// @return does this dof exist, yes/no
bool DofSpace::hasDofType(std::string name, size_t &DofType_out, size_t &DofStep_out){
    auto it = std::find(DofNames.begin(), DofNames.end(), name);
    if (it != DofNames.end()) {
        DofType_out = it - DofNames.begin();
        DofStep_out = DofStep[DofType_out];
        return true;
    } else {
        return false;
    }    
}

/// @brief Check if degree of freedom exists
/// @param name input: Strings to identify the type of degree of freedom
/// @return does this dof exist, yes/no
bool DofSpace::hasDofType(std::string name){
    size_t dummy1, dummy2;
    return hasDofType(name, dummy1, dummy2);
}

/// @brief Selects degrees of freedom to be added to nodes
/// @param Nodes Nodes to which dofs should be added
/// @param DofInd Index of the dof to add (applies ALL DOF types to ALL nodes)
void DofSpace::AddDofs(std::vector<size_t> &Nodes, std::vector<size_t> &DofInd){
    if (DofInd.empty()) {
        // Empty DOF index array is invalid
        throw std::runtime_error("DofSpace::AddDofs: Empty DOF index array");
    }
    if (Nodes.empty()) {
        // Empty nodes array is valid but nothing to do - return early
        return;
    }
    
    // Apply each DOF type to all nodes
    for (size_t idof = 0; idof < DofInd.size(); idof++){
        if (DofInd[idof] >= nDofs) {
            throw std::runtime_error("DofSpace::AddDofs: DOF index " + std::to_string(DofInd[idof]) + 
                                   " (position " + std::to_string(idof) + ") exceeds maximum " + std::to_string(nDofs-1));
        }
        AddDofs(Nodes, DofInd[idof]);
    }
}

/// @brief Selects degrees of freedom to be added to nodes
/// @param Nodes Nodes to which dofs should be added
/// @param DofInd Index of the dof to add
void DofSpace::AddDofs(std::vector<size_t> &Nodes, size_t DofInd){
    if (DofInd >= nDofs) {
        throw std::runtime_error("DofSpace::AddDofs: DOF index " + std::to_string(DofInd) + 
                               " exceeds maximum " + std::to_string(nDofs-1));
    }
    if (DofInd >= DofsToAdd.size()) {
        throw std::runtime_error("DofSpace::AddDofs: DOF index " + std::to_string(DofInd) + 
                               " exceeds DofsToAdd array size " + std::to_string(DofsToAdd.size()));
    }
    
    // If nodes array is empty, nothing to do - return early
    if (Nodes.empty()) {
        return;
    }
    
    size_t initialSize = DofsToAdd[DofInd].size();
    DofsToAdd[DofInd].insert(DofsToAdd[DofInd].end(), Nodes.begin(), Nodes.end());
    std::sort(DofsToAdd[DofInd].begin(), DofsToAdd[DofInd].end());
    auto last = std::unique(DofsToAdd[DofInd].begin(), DofsToAdd[DofInd].end());
    DofsToAdd[DofInd].erase(last, DofsToAdd[DofInd].end());
    
    // Validate we didn't lose any valid nodes (only check if we had nodes to add)
    if (DofsToAdd[DofInd].size() < initialSize) {
        throw std::runtime_error("DofSpace::AddDofs: Node list size decreased from " + std::to_string(initialSize) + 
                               " to " + std::to_string(DofsToAdd[DofInd].size()) + " for DOF " + std::to_string(DofInd));
    }
}

/// @brief Selects degrees of freedom to be added to nodes
/// @param Node Nodes to which dofs should be added
/// @param DofInd Index of the dof to add
void DofSpace::AddDofs(size_t Node, size_t DofInd){
    if (DofInd >= nDofs) {
        throw std::runtime_error("DofSpace::AddDofs: DOF index " + std::to_string(DofInd) + 
                               " exceeds maximum " + std::to_string(nDofs-1));
    }
    if (DofInd >= DofsToAdd.size()) {
        throw std::runtime_error("DofSpace::AddDofs: DOF index " + std::to_string(DofInd) + 
                               " exceeds DofsToAdd array size " + std::to_string(DofsToAdd.size()));
    }
    
    DofsToAdd[DofInd].insert(DofsToAdd[DofInd].end(), Node);
    std::sort(DofsToAdd[DofInd].begin(), DofsToAdd[DofInd].end());
    auto last = std::unique(DofsToAdd[DofInd].begin(), DofsToAdd[DofInd].end());
    DofsToAdd[DofInd].erase(last, DofsToAdd[DofInd].end());
}

/// @brief Flushes degrees of freedom to a permanent structure, and syncs across nodes
/// @param mesh pointer to mesh (to get owner ranges)
void DofSpace::SyncDofs(Mesh* mesh){
    if (!mesh) {
        throw std::invalid_argument("DofSpace::SyncDofs: mesh pointer is null");
    }
    
    // Get Node Owner Ranges
    std::vector<size_t> NodeRanges(size+1);
    const PetscInt *NodeRangesPETSC;
    PetscCallThrow(VecGetOwnershipRanges(mesh->Xcoords.DataVector, &NodeRangesPETSC));
    for (int i = 0; i < size+1; i++){
        NodeRanges[i] = NodeRangesPETSC[i];
    }
    
    // Validate node ranges
    if (NodeRanges[0] != 0) {
        throw std::runtime_error("DofSpace::SyncDofs: First node range should start at 0, got " + std::to_string(NodeRanges[0]));
    }
    for (int i = 1; i < size+1; i++) {
        if (NodeRanges[i] <= NodeRanges[i-1]) {
            throw std::runtime_error("DofSpace::SyncDofs: Node ranges not monotonically increasing at index " + std::to_string(i) + 
                                   " (" + std::to_string(NodeRanges[i-1]) + " -> " + std::to_string(NodeRanges[i]) + ")");
        }
    }
    
    // Validate DofsToAdd before processing
    for (size_t dofInd = 0; dofInd < nDofs; dofInd++) {
        for (size_t nd = 0; nd < DofsToAdd[dofInd].size(); nd++) {
            size_t nodeId = DofsToAdd[dofInd][nd];
            if (nodeId >= NodeRanges[size]) {
                throw std::runtime_error("DofSpace::SyncDofs: Node ID " + std::to_string(nodeId) + 
                                       " for DOF type " + std::to_string(dofInd) + 
                                       " exceeds maximum node range " + std::to_string(NodeRanges[size]-1));
            }
        }
    }

    //send requested dofs to other cores to verify they are all added (doing this first to not have to fix dof numbering later on)
    {
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> RequestNum(size,size); RequestNum.setZero();
        Eigen::Matrix<std::vector<size_t>, Eigen::Dynamic, Eigen::Dynamic> RequestNodes(size,size), RequestDofs(size,size);
        
        size_t totalRequests = 0;
        for (size_t dofInd = 0; dofInd < nDofs; dofInd++){
            if (dofInd >= DofStep.size()) {
                throw std::runtime_error("DofSpace::SyncDofs: DOF index " + std::to_string(dofInd) + " exceeds DofStep size " + std::to_string(DofStep.size()));
            }
            
            for (size_t nd = 0; nd < DofsToAdd[dofInd].size(); nd++){
                int c_to_add_to = -1;
                for (int c = 0; c < size; c++){
                    if (DofsToAdd[dofInd][nd] >= NodeRanges[c] && DofsToAdd[dofInd][nd] < NodeRanges[c+1]){
                        c_to_add_to = c;
                        break;
                    }
                }
                if (c_to_add_to == -1) {
                    throw std::runtime_error("DofSpace::SyncDofs: Node " + std::to_string(DofsToAdd[dofInd][nd]) + 
                                           " (DOF type " + std::to_string(dofInd) + ") not found in any core range [0-" + std::to_string(NodeRanges[size]-1) + "]");
                }
                if (c_to_add_to < 0 || c_to_add_to >= size) {
                    throw std::runtime_error("DofSpace::SyncDofs: Invalid core assignment " + std::to_string(c_to_add_to) + 
                                           " for node " + std::to_string(DofsToAdd[dofInd][nd]) + " (valid range: 0-" + std::to_string(size-1) + ")");
                }
                
                // CRITICAL FIX: Always add nodes to their owner core, regardless of who requests them
                // This ensures that each node gets a DOF assignment on the core that owns it
                if (c_to_add_to == rank) {
                    // Node belongs to this core - ensure it's in our DOF list (already done above)
                } else {
                    // Node belongs to another core - request that core to add it
                    RequestNum(rank,c_to_add_to) += 1;
                    RequestNodes(rank,c_to_add_to).push_back(DofsToAdd[dofInd][nd]);
                    RequestDofs(rank,c_to_add_to).push_back(dofInd);
                    totalRequests++;
                }
            }
        }
        
        // Validate request matrix bounds
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (RequestNum(i,j) > 1000000) { // Sanity check for unreasonable request counts
                    throw std::runtime_error("DofSpace::SyncDofs: Excessive request count " + std::to_string(RequestNum(i,j)) + 
                                           " from core " + std::to_string(i) + " to core " + std::to_string(j));
                }
                if (RequestNum(i,j) != RequestNodes(i,j).size() || RequestNum(i,j) != RequestDofs(i,j).size()) {
                    throw std::runtime_error("DofSpace::SyncDofs: Request count mismatch for cores " + std::to_string(i) + 
                                           "->" + std::to_string(j) + ": count=" + std::to_string(RequestNum(i,j)) + 
                                           ", nodes=" + std::to_string(RequestNodes(i,j).size()) + 
                                           ", dofs=" + std::to_string(RequestDofs(i,j).size()));
                }
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, RequestNum.data(), size*size, my_MPI_SIZE_T, MPI_SUM, PETSC_COMM_WORLD);
        
        // Validate AllReduce results
        size_t totalRequestsGlobal = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                totalRequestsGlobal += RequestNum(i,j);
                if (RequestNum(i,j) > 1000000) {
                    throw std::runtime_error("DofSpace::SyncDofs: Post-AllReduce excessive request count " + std::to_string(RequestNum(i,j)) + 
                                           " from core " + std::to_string(i) + " to core " + std::to_string(j));
                }
            }
        }

        for (size_t c = 0; c < size; c++){
            MPI_Status mySendRequests;
            if (c==rank){ //receiving and adding
                for (size_t c2 = 0; c2 < size; c2++){
                    if (c2!=rank){
                        size_t nToCheck = RequestNum(c2, rank);
                        if (nToCheck > 0) {
                            if (nToCheck > 1000000) {
                                throw std::runtime_error("DofSpace::SyncDofs: Receiving excessive request count " + std::to_string(nToCheck) + 
                                                       " from core " + std::to_string(c2));
                            }
                            
                            std::vector<size_t> RecNodes(nToCheck), RecDofs(nToCheck);
                            // Use unique tags to avoid collision: sender_id * 1000 + receiver_id * 10 + message_type
                            MPI_Recv(RecNodes.data(), nToCheck, my_MPI_SIZE_T, c2, c2*1000+rank*10+1, PETSC_COMM_WORLD, &mySendRequests);
                            MPI_Recv(RecDofs.data(), nToCheck, my_MPI_SIZE_T, c2, c2*1000+rank*10+2, PETSC_COMM_WORLD, &mySendRequests);

                            // Validate received data
                            for (size_t i = 0; i < RecNodes.size(); i++){
                                if (RecNodes[i] >= NodeRanges[size]) {
                                    throw std::runtime_error("DofSpace::SyncDofs: Received invalid node ID " + std::to_string(RecNodes[i]) + 
                                                           " from core " + std::to_string(c2) + " (max valid: " + std::to_string(NodeRanges[size]-1) + ")");
                                }
                                if (RecDofs[i] >= nDofs) {
                                    throw std::runtime_error("DofSpace::SyncDofs: Received invalid DOF type " + std::to_string(RecDofs[i]) + 
                                                           " from core " + std::to_string(c2) + " (max valid: " + std::to_string(nDofs-1) + ")");
                                }
                                
                                // Verify node belongs to this core
                                if (RecNodes[i] < NodeRanges[rank] || RecNodes[i] >= NodeRanges[rank+1]) {
                                    throw std::runtime_error("DofSpace::SyncDofs: Received node " + std::to_string(RecNodes[i]) + 
                                                           " from core " + std::to_string(c2) + " but it doesn't belong to core " + std::to_string(rank) + 
                                                           " (range: " + std::to_string(NodeRanges[rank]) + "-" + std::to_string(NodeRanges[rank+1]-1) + ")");
                                }
                                
                                AddDofs(RecNodes[i], RecDofs[i]);
                            }
                        }
                    }
                }
            } else { //sending
                size_t n_to_send = RequestNum(rank, c);
                if (n_to_send > 0) {
                    if (n_to_send != RequestNodes(rank, c).size() || n_to_send != RequestDofs(rank, c).size()) {
                        throw std::runtime_error("DofSpace::SyncDofs: Send size mismatch to core " + std::to_string(c) + 
                                               ": expected " + std::to_string(n_to_send) + 
                                               ", got nodes=" + std::to_string(RequestNodes(rank, c).size()) + 
                                               ", dofs=" + std::to_string(RequestDofs(rank, c).size()));
                    }
                    
                    MPI_Send(RequestNodes(rank, c).data(), n_to_send, my_MPI_SIZE_T, c, rank*1000+c*10+1, PETSC_COMM_WORLD);
                    MPI_Send(RequestDofs(rank, c).data(), RequestNum(rank, c), my_MPI_SIZE_T, c, rank*1000+c*10+2, PETSC_COMM_WORLD);
                }
            }
        }
    }

    //set up dof space
    for (size_t CurdofStep = 0; CurdofStep < maxSteps; CurdofStep++){
        if (CurdofStep >= DofNumbering.size() || CurdofStep >= GhostDofs.size() || 
            CurdofStep >= LocalDofNumRange.size() || CurdofStep >= TotalDofRange.size()) {
            throw std::runtime_error("DofSpace::SyncDofs: Current step " + std::to_string(CurdofStep) + 
                                   " exceeds allocated arrays (max: " + std::to_string(maxSteps-1) + ")");
        }
        
        //sort Dofs to Add into "buckets"
        std::vector<std::vector<std::vector<size_t>>> DofsToAddPerCore; DofsToAddPerCore.resize(size);   //(Core;dof -> Nodenumbers)
        for (int i = 0; i < size; i++){
            DofsToAddPerCore[i].resize(nDofs);
        }
        
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> RequestNum(size,size); RequestNum.setZero();
        Eigen::Matrix<std::vector<size_t>, Eigen::Dynamic, Eigen::Dynamic> RequestNodes(size,size), RequestDofs(size,size);
        for (int i = 0; i < size; i++){
            for (int j = 0; j < size; j++){
                RequestNodes(i,j).resize(0);
            } 
        }
        
        // Count DOFs for current step and validate
        size_t dofsForThisStep = 0;
        for (size_t dofInd = 0; dofInd < nDofs; dofInd++){
            if (DofStep[dofInd] == CurdofStep){
                dofsForThisStep++;
                
                for (size_t nd = 0; nd < DofsToAdd[dofInd].size(); nd++){
                    int c_to_add_to = -1;
                    for (int c = 0; c < size; c++){
                        if (DofsToAdd[dofInd][nd] >= NodeRanges[c] && DofsToAdd[dofInd][nd] < NodeRanges[c+1]){
                            c_to_add_to = c;
                            break;
                        }
                    }
                    if (c_to_add_to == -1) {
                        throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Node " + 
                                               std::to_string(DofsToAdd[dofInd][nd]) + " not found in any core range");
                    }
                    if (c_to_add_to < 0 || c_to_add_to >= size) {
                        throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Invalid core " + 
                                               std::to_string(c_to_add_to) + " for node " + std::to_string(DofsToAdd[dofInd][nd]));
                    }
                    
                    DofsToAddPerCore[c_to_add_to][dofInd].insert(DofsToAddPerCore[c_to_add_to][dofInd].end(), DofsToAdd[dofInd][nd]);
                    if (c_to_add_to!=rank){
                        RequestNum(rank,c_to_add_to) += 1;
                        RequestNodes(rank,c_to_add_to).push_back(DofsToAdd[dofInd][nd]);
                        RequestDofs(rank,c_to_add_to).push_back(dofInd);
                    }
                }
            }
        }
        
        if (dofsForThisStep == 0) {
            // This is normal - some steps might not have DOFs, just continue
            continue;
        }

        std::vector<size_t> TotalDofs(size); for (int i = 0; i < size; i++) TotalDofs[i] = 0;
        for (int i = 0; i < size; i++){
            for (size_t j = 0; j < nDofs; j++){
                TotalDofs[i] += DofsToAddPerCore[i][j].size();
                
                // Remove duplicates within each core's DOF list
                std::sort(DofsToAddPerCore[i][j].begin(), DofsToAddPerCore[i][j].end());
                auto last = std::unique(DofsToAddPerCore[i][j].begin(), DofsToAddPerCore[i][j].end());
                if (last != DofsToAddPerCore[i][j].end()) {
                    size_t duplicates = DofsToAddPerCore[i][j].end() - last;
                    DofsToAddPerCore[i][j].erase(last, DofsToAddPerCore[i][j].end());
                    TotalDofs[i] -= duplicates; // Adjust count after removing duplicates
                }
            }
            
            if (TotalDofs[i] > 10000000) { // Sanity check
                throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Core " + std::to_string(i) + 
                                       " has excessive DOF count: " + std::to_string(TotalDofs[i]));
            }
        }

        // sync dofs per core and request table
        std::vector<size_t> OwnedDofs(size); for (int i = 0; i < size; i++) OwnedDofs[i] = 0;
        OwnedDofs[rank] = TotalDofs[rank];
        
        // Validate before AllReduce
        if (rank < 0 || rank >= size) {
            throw std::runtime_error("DofSpace::SyncDofs: Invalid rank " + std::to_string(rank) + " (size: " + std::to_string(size) + ")");
        }
        
        MPI_Allreduce(MPI_IN_PLACE, OwnedDofs.data(), size, my_MPI_SIZE_T, MPI_SUM, PETSC_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, RequestNum.data(), size*size, my_MPI_SIZE_T, MPI_SUM, PETSC_COMM_WORLD);
        
        // Validate AllReduce results
        size_t totalDofCount = 0;
        for (int i = 0; i < size; i++) {
            if (OwnedDofs[i] > 10000000) {
                throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Core " + std::to_string(i) + 
                                       " reports excessive owned DOFs: " + std::to_string(OwnedDofs[i]));
            }
            totalDofCount += OwnedDofs[i];
        }
        
        if (totalDofCount == 0) {
            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: No DOFs found across all cores");
        }
        
        PetscCallThrow(PetscBarrier(NULL)); 

        //set up local dof numberings
        DofNumbering[CurdofStep].resize(nDofs);

        LocalDofNumRange[CurdofStep].resize(2);  
        LocalDofNumRange[CurdofStep][0] = 0; 
        LocalDofNumRange[CurdofStep][1] = 0;
        TotalDofRange[CurdofStep] = 0;
        
        // Calculate ranges with validation
        for (int i = 0; i < rank; i++){
            LocalDofNumRange[CurdofStep][0] += OwnedDofs[i];
        }
        for (int i = 0; i < rank+1; i++){
            LocalDofNumRange[CurdofStep][1] += OwnedDofs[i];
        }
        for (int i = 0; i < size; i++){
            TotalDofRange[CurdofStep] += OwnedDofs[i];
        }
        
        // Validate ranges
        if (LocalDofNumRange[CurdofStep][0] > LocalDofNumRange[CurdofStep][1]) {
            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Invalid local range [" + 
                                   std::to_string(LocalDofNumRange[CurdofStep][0]) + ", " + std::to_string(LocalDofNumRange[CurdofStep][1]) + "]");
        }
        
        size_t expectedLocalCount = LocalDofNumRange[CurdofStep][1] - LocalDofNumRange[CurdofStep][0];
        if (expectedLocalCount != OwnedDofs[rank]) {
            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Local DOF count mismatch: " + 
                                   "range gives " + std::to_string(expectedLocalCount) + 
                                   " but owned count is " + std::to_string(OwnedDofs[rank]));
        }
        
        size_t DofCounter = LocalDofNumRange[CurdofStep][0];
        size_t assignedDofs = 0;
        
        for (size_t nDof = 0; nDof < nDofs; nDof++){
            if (DofStep[nDof] == CurdofStep) {
                for (size_t dof = 0; dof < DofsToAddPerCore[rank][nDof].size(); dof++){
                    size_t nodeId = DofsToAddPerCore[rank][nDof][dof];
                    
                    // Validate node ownership
                    if (nodeId < NodeRanges[rank] || nodeId >= NodeRanges[rank+1]) {
                        throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Trying to assign DOF " + 
                                               std::to_string(DofCounter) + " to node " + std::to_string(nodeId) + 
                                               " which doesn't belong to core " + std::to_string(rank) + 
                                               " (valid range: " + std::to_string(NodeRanges[rank]) + "-" + std::to_string(NodeRanges[rank+1]-1) + ")");
                    }
                    
                    // Check for duplicate assignments
                    auto insertResult = DofNumbering[CurdofStep][nDof].insert({nodeId, DofCounter});
                    if (!insertResult.second) {
                        throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Duplicate DOF assignment for node " + 
                                               std::to_string(nodeId) + " DOF type " + std::to_string(nDof) + 
                                               ": existing " + std::to_string(insertResult.first->second) + 
                                               ", new " + std::to_string(DofCounter));
                    }
                    
                    DofCounter += 1;
                    assignedDofs++;
                }
            }
        }
        
        // Validate final DOF assignment count
        if (assignedDofs != OwnedDofs[rank]) {
            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: DOF assignment count mismatch: " + 
                                   "assigned " + std::to_string(assignedDofs) + 
                                   " but expected " + std::to_string(OwnedDofs[rank]));
        }
        
        if (DofCounter != LocalDofNumRange[CurdofStep][1]) {
            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Final DOF counter mismatch: " + 
                                   "counter at " + std::to_string(DofCounter) + 
                                   " but expected " + std::to_string(LocalDofNumRange[CurdofStep][1]));
        }
        
        if (size>1){// Get ghosted dof numbering
            //send out requests
            std::vector<MPI_Request> mySendRequests; mySendRequests.reserve(2*size*size+10); 
            std::vector<size_t> mySendIDs; mySendIDs.reserve(2*size*size+10); 
            std::vector<size_t> myReceiveIDs; myReceiveIDs.reserve(2*size*size+10); 
            
            // Validate request sizes before sending
            size_t totalRequestsToSend = 0;
            for (int RequestTo = 0; RequestTo < size; RequestTo++){
                if (RequestNum(rank, RequestTo)>0){
                    if (RequestNum(rank, RequestTo) != RequestNodes(rank,RequestTo).size() || 
                        RequestNum(rank, RequestTo) != RequestDofs(rank,RequestTo).size()) {
                        throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Request size mismatch to core " + 
                                               std::to_string(RequestTo) + ": count=" + std::to_string(RequestNum(rank, RequestTo)) + 
                                               ", nodes=" + std::to_string(RequestNodes(rank,RequestTo).size()) + 
                                               ", dofs=" + std::to_string(RequestDofs(rank,RequestTo).size()));
                    }
                    totalRequestsToSend += RequestNum(rank, RequestTo);
                    
                    // Validate nodes and DOF types in requests
                    for (size_t i = 0; i < RequestNodes(rank,RequestTo).size(); i++) {
                        if (RequestNodes(rank,RequestTo)[i] >= NodeRanges[size]) {
                            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Invalid node ID " + 
                                                   std::to_string(RequestNodes(rank,RequestTo)[i]) + " in request to core " + std::to_string(RequestTo));
                        }
                        if (RequestDofs(rank,RequestTo)[i] >= nDofs) {
                            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Invalid DOF type " + 
                                                   std::to_string(RequestDofs(rank,RequestTo)[i]) + " in request to core " + std::to_string(RequestTo));
                        }
                        if (DofStep[RequestDofs(rank,RequestTo)[i]] != CurdofStep) {
                            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: DOF type " + 
                                                   std::to_string(RequestDofs(rank,RequestTo)[i]) + " belongs to step " + 
                                                   std::to_string(DofStep[RequestDofs(rank,RequestTo)[i]]) + " not current step");
                        }
                    }
                    
                    mySendRequests.push_back(0); mySendIDs.push_back(rank*1000+RequestTo*10+3);
                    MPI_Isend(RequestNodes(rank,RequestTo).data(), RequestNodes(rank,RequestTo).size(), my_MPI_SIZE_T, RequestTo, mySendIDs.back(), PETSC_COMM_WORLD, &mySendRequests.back());
                    mySendRequests.push_back(0); mySendIDs.push_back(rank*1000+RequestTo*10+4);
                    MPI_Isend(RequestDofs(rank,RequestTo).data(), RequestDofs(rank,RequestTo).size(), my_MPI_SIZE_T, RequestTo, mySendIDs.back(), PETSC_COMM_WORLD, &mySendRequests.back());
                }
            }
            //receive requests
            size_t totalRequestsToReceive = 0;
            for (int RequestFrom = 0; RequestFrom < size; RequestFrom++){
                if (RequestNum(RequestFrom, rank)>0){
                    size_t expectedSize = RequestNum(RequestFrom, rank);
                    if (expectedSize > 1000000) {
                        throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Excessive request size " + 
                                               std::to_string(expectedSize) + " from core " + std::to_string(RequestFrom));
                    }
                    
                    totalRequestsToReceive += expectedSize;
                    RequestNodes(RequestFrom,rank).resize(expectedSize);
                    RequestDofs(RequestFrom,rank).resize(expectedSize);

                    mySendRequests.push_back(0);  myReceiveIDs.push_back(RequestFrom*1000+rank*10+3);
                    MPI_Irecv(RequestNodes(RequestFrom,rank).data(), expectedSize, my_MPI_SIZE_T, RequestFrom, myReceiveIDs.back(), PETSC_COMM_WORLD, &mySendRequests.back());
                    mySendRequests.push_back(0);  myReceiveIDs.push_back(RequestFrom*1000+rank*10+4);
                    MPI_Irecv(RequestDofs(RequestFrom,rank).data(), expectedSize, my_MPI_SIZE_T, RequestFrom, myReceiveIDs.back(), PETSC_COMM_WORLD, &mySendRequests.back());
                }
            }
            
            if (mySendRequests.size() > 10000) {
                throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Excessive MPI request count: " + 
                                       std::to_string(mySendRequests.size()));
            }
            
            std::vector<MPI_Status> myStatus(mySendRequests.size());
            MPI_Waitall(mySendRequests.size(), mySendRequests.data(), myStatus.data());
            
            // Validate received requests
            for (int RequestFrom = 0; RequestFrom < size; RequestFrom++){
                if (RequestNum(RequestFrom, rank)>0){
                    for (size_t i = 0; i < RequestNodes(RequestFrom,rank).size(); i++) {
                        size_t nodeId = RequestNodes(RequestFrom,rank)[i];
                        size_t dofType = RequestDofs(RequestFrom,rank)[i];
                        
                        // Validate received node belongs to this core
                        if (nodeId < NodeRanges[rank] || nodeId >= NodeRanges[rank+1]) {
                            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Received request for node " + 
                                                   std::to_string(nodeId) + " from core " + std::to_string(RequestFrom) + 
                                                   " but node doesn't belong to core " + std::to_string(rank) + 
                                                   " (range: " + std::to_string(NodeRanges[rank]) + "-" + std::to_string(NodeRanges[rank+1]-1) + ")");
                        }
                        
                        if (dofType >= nDofs) {
                            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Received invalid DOF type " + 
                                                   std::to_string(dofType) + " from core " + std::to_string(RequestFrom));
                        }
                        
                        if (DofStep[dofType] != CurdofStep) {
                            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Received DOF type " + 
                                                   std::to_string(dofType) + " from core " + std::to_string(RequestFrom) + 
                                                   " but it belongs to step " + std::to_string(DofStep[dofType]));
                        }
                    }
                }
            }
            
            PetscCallThrow(PetscBarrier(NULL)); 
            //answer requests
            mySendRequests.resize(0);   mySendRequests.reserve(2*size*size+10); 
            mySendIDs.resize(0);        mySendIDs.reserve(2*size*size+10); 
            myReceiveIDs.resize(0);     myReceiveIDs.reserve(2*size*size+10); 
            std::vector<std::vector<size_t>> DofsToSend;
            std::vector<size_t> RequestVec;
            
            for (int RequestFrom = 0; RequestFrom < size; RequestFrom++){
                if (RequestNum(RequestFrom, rank)>0){
                    RequestVec.push_back(RequestFrom);
                    mySendRequests.push_back(0); mySendIDs.push_back(rank*1000+RequestFrom*10+5);
                    DofsToSend.resize(DofsToSend.size()+1); 
                    DofsToSend[DofsToSend.size()-1].resize(RequestNum(RequestFrom, rank));
                    
                    try {
                        getDofForNodesSeries(RequestNodes(RequestFrom,rank), RequestDofs(RequestFrom,rank), DofsToSend[DofsToSend.size()-1]);
                    } catch (const std::exception& e) {
                        throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Failed to get DOFs for request from core " + 
                                               std::to_string(RequestFrom) + ": " + e.what());
                    }
                    
                    // Validate DOF numbers before sending
                    for (size_t i = 0; i < DofsToSend[DofsToSend.size()-1].size(); i++) {
                        size_t dofNum = DofsToSend[DofsToSend.size()-1][i];
                        if (dofNum < LocalDofNumRange[CurdofStep][0] || dofNum >= LocalDofNumRange[CurdofStep][1]) {
                            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: DOF number " + 
                                                   std::to_string(dofNum) + " for core " + std::to_string(RequestFrom) + 
                                                   " is outside local range [" + std::to_string(LocalDofNumRange[CurdofStep][0]) + 
                                                   ", " + std::to_string(LocalDofNumRange[CurdofStep][1]) + ")");
                        }
                    }
                }
            }
            
            if (RequestVec.size() != DofsToSend.size()) {
                throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Request vector size mismatch: " + 
                                       std::to_string(RequestVec.size()) + " vs " + std::to_string(DofsToSend.size()));
            }
            
            for (size_t i = 0; i < mySendIDs.size(); i++){
                if (i >= RequestVec.size() || i >= DofsToSend.size()) {
                    throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Index " + std::to_string(i) + 
                                           " out of bounds for send arrays");
                }
                MPI_Isend(DofsToSend[i].data(), DofsToSend[i].size(), my_MPI_SIZE_T, RequestVec[i], mySendIDs[i], PETSC_COMM_WORLD, &mySendRequests[i]);
            }

            //process requests
            std::vector<std::vector<size_t>> DofsReceived(size);
            for (int RequestTo = 0; RequestTo < size; RequestTo++){
                if (RequestNum(rank, RequestTo)>0){
                    size_t expectedSize = RequestNum(rank, RequestTo);
                    mySendRequests.push_back(0); myReceiveIDs.push_back(RequestTo*1000+rank*10+5);
                    DofsReceived[RequestTo].resize(expectedSize);
                    MPI_Irecv(DofsReceived[RequestTo].data(), expectedSize, my_MPI_SIZE_T, RequestTo, myReceiveIDs.back(), PETSC_COMM_WORLD, &mySendRequests.back());
                }
            }
            
            myStatus.resize(mySendRequests.size());
            MPI_Waitall(mySendRequests.size(), mySendRequests.data(), myStatus.data());
            PetscCallThrow(PetscBarrier(NULL)); 

            size_t totalGhostDofs = 0;
            for (int c = 0; c < size; c++){
                if (RequestNum(rank, c)>0){
                    if (DofsReceived[c].size() != RequestNum(rank, c)) {
                        throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Received DOF count mismatch from core " + 
                                               std::to_string(c) + ": expected " + std::to_string(RequestNum(rank, c)) + 
                                               ", got " + std::to_string(DofsReceived[c].size()));
                    }
                    
                    for (size_t idof = 0; idof < DofsReceived[c].size(); idof++){
                        size_t dofType = RequestDofs(rank,c)[idof];
                        size_t correctStep = DofStep[dofType];  // Use the correct step for this DOF type
                        size_t nodeId = RequestNodes(rank,c)[idof];
                        size_t receivedDofNum = DofsReceived[c][idof];
                        
                        // Verify we're in the correct step for this DOF type
                        if (correctStep != CurdofStep) {
                            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Ghost DOF type " + 
                                                   std::to_string(dofType) + " belongs to step " + std::to_string(correctStep) + 
                                                   " but we're processing step " + std::to_string(CurdofStep));
                        }
                        
                        // Validate received DOF number is in valid global range
                        if (receivedDofNum >= TotalDofRange[correctStep]) {
                            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Received invalid DOF number " + 
                                                   std::to_string(receivedDofNum) + " from core " + std::to_string(c) + 
                                                   " (max valid: " + std::to_string(TotalDofRange[correctStep]-1) + ")");
                        }
                        
                        // Validate DOF number is NOT in our local range (should be ghost)
                        if (receivedDofNum >= LocalDofNumRange[correctStep][0] && receivedDofNum < LocalDofNumRange[correctStep][1]) {
                            throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Received DOF number " + 
                                                   std::to_string(receivedDofNum) + " from core " + std::to_string(c) + 
                                                   " is in our local range [" + std::to_string(LocalDofNumRange[correctStep][0]) + 
                                                   ", " + std::to_string(LocalDofNumRange[correctStep][1]) + ") - should be ghost");
                        }
                        
                        // Check for duplicate ghost DOF assignments
                        auto insertResult = DofNumbering[correctStep][dofType].insert({nodeId, receivedDofNum});
                        if (!insertResult.second) {
                            if (insertResult.first->second != receivedDofNum) {
                                throw std::runtime_error("DofSpace::SyncDofs[Step " + std::to_string(CurdofStep) + "]: Conflicting ghost DOF assignment for node " + 
                                                       std::to_string(nodeId) + " DOF type " + std::to_string(dofType) + 
                                                       ": existing " + std::to_string(insertResult.first->second) + 
                                                       ", received " + std::to_string(receivedDofNum) + " from core " + std::to_string(c));
                            }
                            // If same value, it's ok (duplicate request)
                        } else {
                            // New ghost DOF added
                            GhostDofs[correctStep].push_back(receivedDofNum);
                            totalGhostDofs++;
                        }
                    }
                }
            }
        }
        
        // Remove duplicates from GhostDofs after collecting all (for current step only)
        std::sort(GhostDofs[CurdofStep].begin(), GhostDofs[CurdofStep].end());
        auto last = std::unique(GhostDofs[CurdofStep].begin(), GhostDofs[CurdofStep].end());
        GhostDofs[CurdofStep].erase(last, GhostDofs[CurdofStep].end());

        printStats(CurdofStep);
    }
}

/// @brief Prints statistics regarding degrees of freedom
/// @param curStep Step to print info for
void DofSpace::printStats(size_t curStep){
    size_t GhostTotal = GhostDofs[curStep].size();
    size_t LocalTotal = LocalDofNumRange[curStep][1]-LocalDofNumRange[curStep][0];
    // combine results from all processes
    MPI_Allreduce(MPI_IN_PLACE, &GhostTotal, 1, my_MPI_SIZE_T, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &LocalTotal, 1, my_MPI_SIZE_T, MPI_SUM, PETSC_COMM_WORLD);

    std::string Message = "Degrees of Freedom Added: ";
    for (size_t i = 0; i < nDofs; i++){
        if (DofStep[i]==curStep){
            Message += DofNames[i] + ", ";
        }
    }
    Message += "Total amount of DOFs: "+std::to_string(TotalDofRange[curStep]);
    Logs.PrintSingle(Message+"\n",1);
    Logs.PrintSingle("Total Local dofs:"+std::to_string(LocalTotal)+"\n",2);
    Logs.PrintSingle("Total Ghost dofs:"+std::to_string(GhostTotal)+"\n",2);
    PetscCallThrow(PetscBarrier(NULL)); 
    Message = "Local Dofs: "+std::to_string(LocalDofNumRange[curStep][1]-LocalDofNumRange[curStep][0])+", Ghost Dofs: "+std::to_string(GhostDofs[curStep].size())+"\n";
    Logs.PrintEvery(Message,2);
    PetscCallThrow(PetscBarrier(NULL));
}

void DofSpace::ExportDofSpace(){
    // Create file
    std::string FileName = FilenamePrefix+"/DofSpace.hdf5";
    
    FileAccessProps fapl;
    fapl.add(MPIOFileAccess{PETSC_COMM_WORLD, MPI_INFO_NULL});

    File file(FileName, File::ReadWrite | File::Create | File::Truncate, fapl);

    // Loop through each step and save the data
    for (size_t step = 0; step < maxSteps; step++){
        std::string stepName = "Step_" + std::to_string(step);
        file.createGroup(stepName);
        
        // Save the DofNumbering for this step
        for (size_t dofInd = 0; dofInd < nDofs; dofInd++){
            std::string DofGroupName = stepName + "/" + DofNames[dofInd];
            file.createGroup(DofGroupName);

            std::map<size_t, size_t> m = DofNumbering[step][dofInd];
            size_t dofNums = m.size();
            std::vector<size_t> nodeNumbers;
            std::vector<size_t> dofNumbers;
            
            for(std::map<size_t,size_t>::iterator it = m.begin(); it != m.end(); ++it) {
                nodeNumbers.push_back(it->first);
                dofNumbers.push_back(it->second);
            }

            std::vector<size_t> TotalDofNums(size);
            MPI_Allgather(&dofNums, 1, my_MPI_SIZE_T, TotalDofNums.data(), 1, my_MPI_SIZE_T, PETSC_COMM_WORLD);
            for (size_t i = 0; i < size; i++){
                DataSet dataset = file.createDataSet<size_t>(DofGroupName + "/Core_" + std::to_string(i), DataSpace({TotalDofNums[i],2}));
                if(i == rank){
                    dataset.select({0,0},{dofNums,1} ).write(nodeNumbers);
                    dataset.select({0,1},{dofNums,1} ).write(dofNumbers);
                }
            }
        }
    }

}

/// @brief Obtains degree of freedom indices for the provided node
/// @param Node input: node to obtain the degree of freedom index for (global numbering)
/// @param DofInd input: Index of degree of freedom type to obtain
/// @param Dof_out output: degree of freedom numbering
void DofSpace::getDofForNodes(size_t &Node, size_t DofInd, PetscInt &Dof_out){
    if (DofInd >= nDofs) {
        throw std::runtime_error("DofSpace::getDofForNodes: DOF index " + std::to_string(DofInd) + 
                               " exceeds maximum " + std::to_string(nDofs-1));
    }
    if (DofStep[DofInd] >= maxSteps) {
        throw std::runtime_error("DofSpace::getDofForNodes: DOF step " + std::to_string(DofStep[DofInd]) + 
                               " for DOF " + std::to_string(DofInd) + " exceeds maximum " + std::to_string(maxSteps-1));
    }
    
    auto it = DofNumbering[DofStep[DofInd]][DofInd].find(Node);
    if (it == DofNumbering[DofStep[DofInd]][DofInd].end()) {
        throw std::runtime_error("DofSpace::getDofForNodes: Node " + std::to_string(Node) + 
                               " not found for DOF type " + std::to_string(DofInd) + 
                               " (step " + std::to_string(DofStep[DofInd]) + ")");
    }
    
    Dof_out = it->second;
}



/// @brief Obtains degrees of freedom indices for the provided nodes
/// @param Nodes input: nodes to obtain the degree of freedom index for (global numbering)
/// @param DofInd input: Index of degree of freedom type to obtain
/// @param Dofs_out output: degree of freedom numbering (needs to be pre-sized)
void DofSpace::getDofForNodes(std::vector<size_t> &Nodes, size_t DofInd, std::vector<size_t> &Dofs_out){
    if (DofInd >= nDofs) {
        throw std::runtime_error("DofSpace::getDofForNodes: DOF index " + std::to_string(DofInd) + 
                               " exceeds maximum " + std::to_string(nDofs-1));
    }
    if (Dofs_out.size() != Nodes.size()) {
        throw std::runtime_error("DofSpace::getDofForNodes: Output array size " + std::to_string(Dofs_out.size()) + 
                               " doesn't match input nodes size " + std::to_string(Nodes.size()));
    }
    if (DofStep[DofInd] >= maxSteps) {
        throw std::runtime_error("DofSpace::getDofForNodes: DOF step " + std::to_string(DofStep[DofInd]) + 
                               " for DOF " + std::to_string(DofInd) + " exceeds maximum " + std::to_string(maxSteps-1));
    }
    
    for (size_t i = 0; i < Nodes.size(); i++){
        auto it = DofNumbering[DofStep[DofInd]][DofInd].find(Nodes[i]);
        if (it == DofNumbering[DofStep[DofInd]][DofInd].end()) {
            // Provide detailed diagnostic information
            size_t step = DofStep[DofInd];
            std::string diagnostics = "\nDiagnostic info:\n";
            diagnostics += "  - Total nodes in DOF numbering for this type: " + std::to_string(DofNumbering[step][DofInd].size()) + "\n";
            diagnostics += "  - DOF name: " + (DofInd < DofNames.size() ? DofNames[DofInd] : "unknown") + "\n";
            diagnostics += "  - Current MPI rank: " + std::to_string(rank) + "/" + std::to_string(size) + "\n";
            diagnostics += "  - Local DOF range for this step: [" + 
                          std::to_string(LocalDofNumRange[step][0]) + ", " + std::to_string(LocalDofNumRange[step][1]) + ")\n";
            
            // Check node ownership ranges (if available from SyncDofs context)
            // This will help determine if the node should be owned by this processor or is a ghost node
            diagnostics += "  - Requested node: " + std::to_string(Nodes[i]) + "\n";
            
            // Show a sample of available nodes (first 10 and last 10)
            if (!DofNumbering[step][DofInd].empty()) {
                diagnostics += "  - Available nodes range: [" + 
                              std::to_string(DofNumbering[step][DofInd].begin()->first) + ", " + 
                              std::to_string(DofNumbering[step][DofInd].rbegin()->first) + "]\n";
                diagnostics += "  - Sample available nodes (first 10): ";
                int count = 0;
                for (auto& pair : DofNumbering[step][DofInd]) {
                    if (count >= 10) {
                        diagnostics += "...";
                        break;
                    }
                    diagnostics += std::to_string(pair.first) + " ";
                    count++;
                }
                diagnostics += "\n";
                
                // Show last few nodes if we have many nodes
                if (DofNumbering[step][DofInd].size() > 10) {
                    diagnostics += "  - Sample available nodes (last 5): ";
                    auto it_rev = DofNumbering[step][DofInd].rbegin();
                    for (int i = 0; i < 5 && it_rev != DofNumbering[step][DofInd].rend(); ++i, ++it_rev) {
                        diagnostics += std::to_string(it_rev->first) + " ";
                    }
                    diagnostics += "\n";
                }
            }
            
            // Check if this could be a ghost node issue
            size_t requestedNode = Nodes[i];
            bool isInRange = false;
            if (!DofNumbering[step][DofInd].empty()) {
                size_t minNode = DofNumbering[step][DofInd].begin()->first;
                size_t maxNode = DofNumbering[step][DofInd].rbegin()->first;
                if (requestedNode < minNode) {
                    diagnostics += "  - Issue: Requested node " + std::to_string(requestedNode) + 
                                  " is BELOW the range of owned nodes [" + std::to_string(minNode) + 
                                  ", " + std::to_string(maxNode) + "] - likely a ghost node from lower rank\n";
                } else if (requestedNode > maxNode) {
                    diagnostics += "  - Issue: Requested node " + std::to_string(requestedNode) + 
                                  " is ABOVE the range of owned nodes [" + std::to_string(minNode) + 
                                  ", " + std::to_string(maxNode) + "] - likely a ghost node from higher rank\n";
                } else {
                    diagnostics += "  - Issue: Requested node " + std::to_string(requestedNode) + 
                                  " is WITHIN the range [" + std::to_string(minNode) + 
                                  ", " + std::to_string(maxNode) + "] but missing - possible DOF assignment issue\n";
                }
            }
            
            throw std::runtime_error("DofSpace::getDofForNodes: Node " + std::to_string(Nodes[i]) + 
                                   " (index " + std::to_string(i) + ") not found for DOF type " + std::to_string(DofInd) + 
                                   " (step " + std::to_string(DofStep[DofInd]) + ")" + diagnostics);
        }
        Dofs_out[i] = it->second;
    }
}

/// @brief Obtains degrees of freedom indices for the provided nodes
/// @param Nodes input: nodes to obtain the degree of freedom index for (global numbering)
/// @param DofInd input: Index of degree of freedom type to obtain
/// @param Dofs_out output: degree of freedom numbering (needs to be pre-sized)
void DofSpace::getDofForNodes(std::vector<size_t> &Nodes, size_t DofInd, std::vector<PetscInt> &Dofs_out){
    if (DofInd >= nDofs) {
        throw std::runtime_error("DofSpace::getDofForNodes: DOF index " + std::to_string(DofInd) + 
                               " exceeds maximum " + std::to_string(nDofs-1));
    }
    if (Dofs_out.size() != Nodes.size()) {
        throw std::runtime_error("DofSpace::getDofForNodes: Output array size " + std::to_string(Dofs_out.size()) + 
                               " doesn't match input nodes size " + std::to_string(Nodes.size()));
    }
    if (DofStep[DofInd] >= maxSteps) {
        throw std::runtime_error("DofSpace::getDofForNodes: DOF step " + std::to_string(DofStep[DofInd]) + 
                               " for DOF " + std::to_string(DofInd) + " exceeds maximum " + std::to_string(maxSteps-1));
    }
    
    for (size_t i = 0; i < Nodes.size(); i++){
        auto it = DofNumbering[DofStep[DofInd]][DofInd].find(Nodes[i]);
        if (it == DofNumbering[DofStep[DofInd]][DofInd].end()) {
            throw std::runtime_error("DofSpace::getDofForNodes: Node " + std::to_string(Nodes[i]) + 
                                   " (index " + std::to_string(i) + ") not found for DOF type " + std::to_string(DofInd) + 
                                   " (step " + std::to_string(DofStep[DofInd]) + ")");
        }
        Dofs_out[i] = it->second;
    }
}

/// @brief Obtains degrees of freedom indices for the provided nodes
/// @param Nodes input: nodes to obtain the degree of freedom index for (global numbering)
/// @param DofInd input: Index of degree of freedom type to obtain
/// @param Dofs_out output: degree of freedom numbering (needs to be pre-sized)
void DofSpace::getDofForNodes(std::vector<size_t> &Nodes, std::vector<size_t> DofInd, std::vector<PetscInt> &Dofs_out){
    if (DofInd.empty()) {
        throw std::runtime_error("DofSpace::getDofForNodes: Empty DOF index array");
    }
    if (Nodes.empty()) {
        throw std::runtime_error("DofSpace::getDofForNodes: Empty nodes array");
    }
    if (Dofs_out.size() != Nodes.size() * DofInd.size()) {
        throw std::runtime_error("DofSpace::getDofForNodes: Output array size " + std::to_string(Dofs_out.size()) + 
                               " doesn't match expected size " + std::to_string(Nodes.size() * DofInd.size()) + 
                               " (nodes: " + std::to_string(Nodes.size()) + ", DOF types: " + std::to_string(DofInd.size()) + ")");
    }
    
    for (size_t j = 0; j < DofInd.size(); j++){
        if (DofInd[j] >= nDofs) {
            throw std::runtime_error("DofSpace::getDofForNodes: DOF index " + std::to_string(DofInd[j]) + 
                                   " (position " + std::to_string(j) + ") exceeds maximum " + std::to_string(nDofs-1));
        }
        if (DofStep[DofInd[j]] >= maxSteps) {
            throw std::runtime_error("DofSpace::getDofForNodes: DOF step " + std::to_string(DofStep[DofInd[j]]) + 
                                   " for DOF " + std::to_string(DofInd[j]) + " exceeds maximum " + std::to_string(maxSteps-1));
        }
        
        for (size_t i = 0; i < Nodes.size(); i++){
            size_t outputIndex = i + j * Nodes.size();
            auto it = DofNumbering[DofStep[DofInd[j]]][DofInd[j]].find(Nodes[i]);
            if (it == DofNumbering[DofStep[DofInd[j]]][DofInd[j]].end()) {
                throw std::runtime_error("DofSpace::getDofForNodes: Node " + std::to_string(Nodes[i]) + 
                                       " (index " + std::to_string(i) + ") not found for DOF type " + std::to_string(DofInd[j]) + 
                                       " (step " + std::to_string(DofStep[DofInd[j]]) + ")");
            }
            Dofs_out[outputIndex] = it->second;
        }
    }
}

/// @brief Obtains degrees of freedom indices for the provided combination of nodes and dof index
/// @param Nodes input: nodes to obtain the degree of freedom index for (global numbering)
/// @param DofInd input: Index of degree of freedom type to obtain
/// @param Dofs_out output: degree of freedom numbering (needs to be pre-sized)
void DofSpace::getDofForNodesSeries(std::vector<size_t> &Nodes, std::vector<size_t> &DofInd, std::vector<size_t> &Dofs_out){
    if (Nodes.size() != DofInd.size()) {
        throw std::runtime_error("DofSpace::getDofForNodesSeries: Nodes array size " + std::to_string(Nodes.size()) + 
                               " doesn't match DOF index array size " + std::to_string(DofInd.size()));
    }
    if (Dofs_out.size() != Nodes.size()) {
        throw std::runtime_error("DofSpace::getDofForNodesSeries: Output array size " + std::to_string(Dofs_out.size()) + 
                               " doesn't match input array size " + std::to_string(Nodes.size()));
    }
    
    for (size_t i = 0; i < Nodes.size(); i++){
        if (DofInd[i] >= nDofs) {
            throw std::runtime_error("DofSpace::getDofForNodesSeries: DOF index " + std::to_string(DofInd[i]) + 
                                   " (position " + std::to_string(i) + ") exceeds maximum " + std::to_string(nDofs-1));
        }
        if (DofStep[DofInd[i]] >= maxSteps) {
            throw std::runtime_error("DofSpace::getDofForNodesSeries: DOF step " + std::to_string(DofStep[DofInd[i]]) + 
                                   " for DOF " + std::to_string(DofInd[i]) + " (position " + std::to_string(i) + 
                                   ") exceeds maximum " + std::to_string(maxSteps-1));
        }
        
        auto it = DofNumbering[DofStep[DofInd[i]]][DofInd[i]].find(Nodes[i]);
        if (it == DofNumbering[DofStep[DofInd[i]]][DofInd[i]].end()) {
            throw std::runtime_error("DofSpace::getDofForNodesSeries: Node " + std::to_string(Nodes[i]) + 
                                   " (position " + std::to_string(i) + ") not found for DOF type " + std::to_string(DofInd[i]) + 
                                   " (step " + std::to_string(DofStep[DofInd[i]]) + ")");
        }
        Dofs_out[i] = it->second;
    }
}

/// @brief Obtain ghost numbering schemes
/// @param Step input: staggered step
/// @return ghost degree of freedom indices
std::vector<PetscInt> DofSpace::getGhostDofNumbers(size_t Step){
    if (Step >= maxSteps) {
        throw std::runtime_error("DofSpace::getGhostDofNumbers: Step " + std::to_string(Step) + 
                               " exceeds maximum " + std::to_string(maxSteps-1));
    }
    if (Step >= GhostDofs.size()) {
        throw std::runtime_error("DofSpace::getGhostDofNumbers: Step " + std::to_string(Step) + 
                               " exceeds GhostDofs array size " + std::to_string(GhostDofs.size()));
    }
    
    return GhostDofs[Step];
}