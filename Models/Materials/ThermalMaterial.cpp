#include "ThermalMaterial.h"

#include "../../Physics/physics.h"

ThermalMaterial::ThermalMaterial(inputData &inputs, std::string MatName): BaseMaterial(inputs, MatName){
    inputs.GetRequired(k, {"properties", MatName, "ThermalProperties", "ThermalConductivity"});
    inputs.GetRequired(cp, {"properties", MatName, "ThermalProperties", "HeatCapacity"});

    LatentHeat = 0.0;
    inputs.GetOptional(LatentHeat, {"properties", MatName, "ThermalProperties", "LatentHeat"});
}

ThermalMaterial::~ThermalMaterial(){

}

