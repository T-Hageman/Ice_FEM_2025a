#ifndef BASALFRICTIONMODEL_H
#define BASALFRICTIONMODEL_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"

class BasalFrictionModel: public BaseModel {
    public:
        std::string Name;
        std::vector<std::string> DofNames_u = {"ux", "uy"};
        std::string DofNames_phase = "phase";

        BasalFrictionModel(Physics& My_Physics, std::string MyName);
        ~BasalFrictionModel();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);

    protected:

    private:
        size_t ElementGroupIndex_u, ElementGroupIndex_phase;
        size_t dofSteps, dofStep_phase;
        std::vector<size_t> dofTypes;
        size_t dofType_phase;

        double friction, k, basalStrength, u0;
        bool Lumped;

        bool Damage = true;
};

void Register_BasalFrictionModel();
BaseModel* New_BasalFrictionModel(Physics& My_Physics, std::string MyNameIn);

#endif

