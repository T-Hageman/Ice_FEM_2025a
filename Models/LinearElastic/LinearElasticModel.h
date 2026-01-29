#ifndef LINEARELASTIC_H
#define LINEARELASTIC_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"

class LinearElasticModel: public BaseModel{
    public:
        std::vector<std::string> DofNames = {"dx","dy"};

        LinearElasticModel(Physics& My_Physics, std::string MyName);
        ~LinearElasticModel();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);
        bool SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>>& Data);
        
    protected:

    private:
        void getB(Eigen::MatrixXd &G, Eigen::MatrixXd &B);
        void SaveStressComponent(std::string DataName, std::vector<std::vector<double>>& Data);

        size_t ElemGroupIndex;
        std::vector<size_t> dofTypes;
        size_t MyStep;

        double Young;
        double Poisson;
        Eigen::Matrix4d D;
};

void Register_LinearElasticModel();
BaseModel* New_LinearElasticModel(Physics& My_Physics, std::string MyNameIn);

#endif

