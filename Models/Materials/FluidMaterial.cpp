#include "FluidMaterial.h"

#include "../../Physics/physics.h"

FluidMaterial::FluidMaterial(inputData &inputs, std::string MatName): BaseMaterial(inputs, MatName){
    inputs.GetRequired(visc, {"properties", MatName, "FluidProperties", "Viscosity"});
    Bulk = 1e12;
    inputs.GetOptional(Bulk, {"properties", MatName, "ThermalProperties", "Bulk"});
}

FluidMaterial::~FluidMaterial(){

}

