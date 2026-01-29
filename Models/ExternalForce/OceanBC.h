#ifndef OCEANBCMODEL_H
#define OCEANBCMODEL_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"
#include "../Fracture/PhaseField/PhaseFieldUtil.h"
#include "../Materials/FluidMaterial.h"

class OceanBCModel: public BaseModel {
    public:
        std::string Name;
        std::vector<std::string> DofNames_u = {"ux","uy"};      //names for displacement degrees of freedom
        std::string Dofname_PF = "phase";
        std::string DofNames_p = "pw";

        PhaseFieldUtility* PF_Util;

        OceanBCModel(Physics& My_Physics, std::string MyName);
        ~OceanBCModel();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);

    protected:

    private:
        void Assemble_U(Mat& K, Vec& f);
        void Assemble_P(Mat& K, Vec& f);

        size_t ElementGroupIndex_u, ElementGroupIndex_Phase, ElementGroupIndex_p;  //index of the element groups
        size_t dofStep_U, dofStep_Phase, dofType_phase, dofStep_p, dofType_p; 
        std::vector<size_t> dofTypes_U;

        double hw;
        double kDummy = 1.0e9;
        double g=9.81;

        bool HasPW;
        double pwAboveWater;

        std::string FluidMatName;
        FluidMaterial* Fluid;
};

void Register_OceanBCModel();
BaseModel* New_OceanBCModel(Physics& My_Physics, std::string MyNameIn);

#endif

