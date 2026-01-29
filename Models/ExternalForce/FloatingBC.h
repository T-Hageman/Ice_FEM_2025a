#ifndef FLOATINGBCMODEL_H
#define FLOATINGBCMODEL_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"
#include "../Materials/FluidMaterial.h"

class FloatingBCModel: public BaseModel {
    public:
        std::string Name;
        std::string DofName_u = "uy";      //names for displacement degrees of freedom

        FloatingBCModel(Physics& My_Physics, std::string MyName);
        ~FloatingBCModel();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);

    protected:

    private:
        void Assemble_U(Mat& K, Vec& f);

        size_t ElementGroupIndex_u;  //index of the element groups
        size_t dofStep_U; 
        size_t dofTypes_U;

        double hw;
        double g=9.81;

        std::string FluidMatName;
        FluidMaterial* Fluid;
};

void Register_FloatingBCModel();
BaseModel* New_FloatingBCModel(Physics& My_Physics, std::string MyNameIn);

#endif

