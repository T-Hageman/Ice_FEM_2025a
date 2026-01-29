#ifndef MODELTYPEREGISTER_H
#define MODELTYPEREGISTER_H

#include "../InputsOutputs/inputData.h"
#include <string>

#include "BaseModel/BaseModel.h"
#include "Constraints/ConstraintsModel.h"
#include "LinearElastic/LinearElasticModel.h"
#include "ExternalForce/ExternalForce.h"
#include "ExternalForce/WeakExternal.h"
#include "ExternalForce/WeakArea.h"
#include "ExternalForce/BasalFriction.h"
// #include "Eulerian/TwoPhase/AllenCahn.h"
// #include "Eulerian/TwoPhase/CahnHilliard.h"
// #include "Eulerian/TwoPhase/Gravity.h"
// #include "Eulerian/TwoPhase/MassBalance.h"
// #include "Eulerian/TwoPhase/MomentumBalance.h"
// #include "Eulerian/TwoPhase/ConstitutiveLinearElastic.h"
// #include "Eulerian/TwoPhase/ConstitutiveViscoElastic.h"
#include "Fracture/PhaseField/PhaseField.h"
// #include "PoroElasticity/PoroElasticFlow.h"
// #include "PoroElasticity/ThermalPorous.h"
// #include "PoroElasticity/VariablePorosity.h"
// #include "PoroElasticity/SurfaceFlow.h"
// #include "PoroElasticity/FractureFlow.h"
// #include "PoroElasticity/SurfaceFlux.h"
#include "ExternalForce/TriAxial.h"
#include "Thermal/ThermalModel.h"
#include "Fracture/PhaseField/CossPhaseField.h"
#include "Constraints/AreaConstraint.h"
#include "ExternalForce/OceanBC.h"
#include "LinearElastic/GeneralSolidModel.h"
#include "Fracture/InterfaceElements/SymFluidInterface.h"
#include "ExternalForce/FloatingBC.h"
#include "LinearElastic/LargeDefModel.h"
#include "Constraints/TimeDepConstraintsModel.h"


extern std::vector<std::string> ModelNames;
extern std::vector<std::function<BaseModel*(Physics&, std::string)>> ModelCreators;

void RegisterModels();
BaseModel* CreateModel(Physics& physics, std::string ModelName, std::string MyName);

#endif