#ifndef BASEMATERIAL_H
#define BASEMATERIAL_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../../InputsOutputs/inputData.h"
#include "../../mesh/mesh.h"
#include "../../Physics/DofSpace.h"
#include "../../Physics/Constrainer.h"
#include "../../InputsOutputs/SaveData.h"
#include "../../utility/utility.h"



class BaseMaterial{
    public:
        BaseMaterial(inputData& inputs, std::string SolidMatName);
        ~BaseMaterial();

        std::string MaterialType, MaterialName;

        double Density;

        double R = 8.31446261815324;
    protected:

    private:

};



#endif