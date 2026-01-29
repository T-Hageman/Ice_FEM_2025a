#include "BaseMaterial.h"

#include "../../Physics/physics.h"

BaseMaterial::BaseMaterial(inputData &inputs, std::string MatName){
    MaterialType = "BaseModel";
    MaterialName = MatName;

    inputs.GetRequired(Density, {"properties", MatName, "SolidProperties", "Density"});
}

BaseMaterial::~BaseMaterial(){

}

