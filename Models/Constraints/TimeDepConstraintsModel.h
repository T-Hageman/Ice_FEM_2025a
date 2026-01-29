#ifndef TDEPCONSTRAINTSMODEL_H
#define TDEPCONSTRAINTSMODEL_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"

class TimeDepConstraintsModel: public BaseModel {
    public:
        std::string Name;

        TimeDepConstraintsModel(Physics& My_Physics, std::string MyName);
        ~TimeDepConstraintsModel();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);
    protected:

    private:
        size_t NodeGroupIndex;
        std::vector<std::string> dofNames;
        std::vector<size_t> dofSteps;
        std::vector<size_t> dofTypes;
        std::vector<double> conVals;
        std::vector<double> dconVals;

};

void Register_TimeDepConstraintsModel();
BaseModel* New_TimeDepConstraintsModel(Physics& My_Physics, std::string MyNameIn);

#endif

