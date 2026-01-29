#ifndef THERMAL_H
#define THERMAL_H

#include <vector>
#include <iostream>
#include <unordered_map>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"
#include "../Materials/ThermalMaterial.h"

class ThermalModel: public BaseModel{
    public:
        std::string DofNames_T = "T";                   //name of wetting pressure degree of freedom 

        ThermalModel(Physics& My_Physics, std::string MyName);
        ~ThermalModel();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);

        bool SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>>& Data);

    protected:

    private:
        size_t ElemGroupIndex_T;  //index of the element groups
        size_t dofType_T;                          // dof index for phasefield
        size_t Step_T;                      // steps in which displacements and phasefield are resolved (NOT THE SAME STEP)

        std::string SolidMatName;   //name of solid material
        ThermalMaterial* M;

};

void Register_ThermalModel();
BaseModel* New_ThermalModel(Physics& My_Physics, std::string MyNameIn);

#endif

