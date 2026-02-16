#include "ModelRegister.h"
#include <iostream>

std::vector<std::string> ModelNames; //vector of available model names
std::vector<std::function<BaseModel*(Physics&, std::string)>> ModelCreators; //pointer to model creation functions


/// @brief Registers the available physics models
void RegisterModels(){
    //Basic models
    Register_BaseModel();
    Register_LinearElasticModel();
    Register_ConstraintsModel();
    Register_ExternalForceModel();
    Register_WeakExternalModel();
    Register_WeakAreaModel();
    Register_BasalFrictionModel();
    Register_AreaConstraintsModel();
    Register_GeneralSolidModel();
    Register_TimeDepConstraintsModel();

    //Eulerian reference frame
    // Register_TwoPhaseGravity();
    // Register_TwoPhaseAllenCahn();
    // Register_TwoPhaseCahnHilliard();
    // Register_ConstitutiveLinearElastic();
    // Register_TwoPhaseMassBalance();
    // Register_TwoPhaseMomentum();
    // Register_TwoPhaseViscoElastic();

    //fracture mechanics
    Register_PhaseFieldDamage();
    Register_CossPhaseField();
    Register_TriAxialModel();
    Register_SymFluidInterfaceModel();

    //Thermals
    Register_ThermalModel();

    //poroelasticity
    // Register_PoroElasticFlow();
    // Register_ThermalPorous();
    Register_VariablePorosity();
    // Register_SurfaceFlow();
    // Register_FractureFlow();
    // Register_SurfaceFlux();
    Register_OceanBCModel();
    Register_FloatingBCModel();
}


/// @brief Creates the physics models
/// @param physics reference to physics object
/// @param ModelName characteristic model to be created
/// @param MyName unique identifier to name model
/// @return pointer to model (remember to destroy when not needed)
BaseModel* CreateModel(Physics& physics, std::string ModelName, std::string MyName){
    BaseModel* Model;

    size_t it = 0;
    while (true){
        if (ModelName == ModelNames[it]){
            Model = ModelCreators[it](physics, MyName);
            return Model;
        }
        it+= 1;
        if (it == ModelNames.size()){
            std::string validTypes = "";
            for (size_t i = 0; i < ModelNames.size(); i++){
                validTypes.append(ModelNames[i]);
                validTypes.append(", ");
            }
            throw std::invalid_argument("Model type "+ModelName+" not defined, valid options are: " + validTypes);
            return 0;
        }
    } 
}