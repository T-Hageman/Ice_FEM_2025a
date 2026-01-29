#ifndef VARIABLEPOROSITY_H
#define VARIABLEPOROSITY_H

#include <vector>
#include <iostream>
#include <unordered_map>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"
#include "PoroElasticRelations.h"
#include "../Materials/ThermalMaterial.h"

class VariablePorosity: public BaseModel{
    public:
        std::string DofNames_p = "pw";            
        std::string DofNames_T = "T";   
        std::string DofNames_Por = "poros";           


        VariablePorosity(Physics& My_Physics, std::string MyName);
        ~VariablePorosity();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);

        bool SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>>& Data);

    protected:

    private:
        bool Fluid = false;

        ThermalMaterial* S;
        ThermalMaterial* F;

        size_t ElemGroupIndex_p, ElemGroupIndex_T, ElemGroupIndex_Por;  //index of the element groups
        size_t dofType_p, dofType_T, dofType_Por;                           // dof index for phasefield
        size_t Step_p, Step_T, Step_Por;                      // steps in which displacements and phasefield are resolved (NOT THE SAME STEP)

        std::string SolidMatName;   //name of solid material
        std::string FluidMatName;

        PoroElasticRelations* PoroRelations;

        bool PhaseChange;
};

void Register_VariablePorosity();
BaseModel* New_VariablePorosity(Physics& My_Physics, std::string MyNameIn);

#endif

