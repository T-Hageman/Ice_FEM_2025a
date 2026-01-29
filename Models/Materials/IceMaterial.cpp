#include "IceMaterial.h"

#include "../../Physics/physics.h"

IceMaterial::IceMaterial(inputData &inputs, std::string MatName): 
    ThermalMaterial(inputs, MatName), CosseratMaterial(inputs, MatName), PhaseFieldMaterial(inputs, MatName), SolidMaterial(inputs, MatName), BaseMaterial(inputs, MatName){

    TDependent = false;
    if (inputs.HasKey({"properties", MatName, "TemperatureDependent"})){
        TDependent = true;
        inputs.GetOptional(TDependent, {"properties", MatName, "TemperatureDependent","TemperatureDependent"});
    }

    TOffset = 0.0;
    inputs.GetOptional(TOffset, {"properties", MatName, "TemperatureDependent", "TOffset"});
    inputs.GetRequired(TRef, {"properties", MatName, "TemperatureDependent", "TRef"});

    A0 = ACreep;
    inputs.GetRequired(QCreep, {"properties", MatName, "TemperatureDependent", "QCreep"});

    inputs.GetRequired(ft0, {"properties", MatName, "TemperatureDependent", "ft0"});
    inputs.GetRequired(ftdeg, {"properties", MatName, "TemperatureDependent", "ftdeg"});

    if (StressSplit==SplittingMethods::ICE){
        Ice_eRef0 = Ice_eRef;
        Ice_cReduce0 = Ice_cReduce;
        inputs.GetRequired(Ice_cReduceFact, {"properties", MatName, "TemperatureDependent", "Ice_cReduceFact"});
        inputs.GetRequired(Ice_eReduceFact, {"properties", MatName, "TemperatureDependent", "Ice_eReduceFact"});
    }

    sy0 = sy;
    inputs.GetRequired(sydeg, {"properties", MatName, "TemperatureDependent", "sydeg"});


    D0 = 1.0*D;
    D_coss0 = 1.0*D_coss;
    G_coss0 = G_coss;
    Young0 = Young;


    SetTemperature(TRef, 0.0);

}

IceMaterial::~IceMaterial(){

}

void IceMaterial::SetTemperature(double T, double poros){
    double dam = 1-poros;
    double TAbs = T+TOffset;
    if (TAbs>273.15){
        TAbs = 273.15;
    }

    //elasticity
    Young = Young0 * dam;
    Shear  = Young/(2.0*(1+Poisson));
    Lame   = Young*Poisson/((1.0+Poisson)*(1.0-2.0*Poisson));
    Bulk = Lame + 2.0/3.0*Shear;
    G_coss = G_coss0 * dam;

    D = dam*D0;
    D_coss = dam*D_coss0;

    //creep
    double ExpTerm = -QCreep/R*(1/TAbs-1/TRef);
    ACreep = A0*std::exp(ExpTerm);

    //Fracture energy
    double ft = dam*(ft0 - ftdeg*(TAbs-TRef));
    Gc = 2.0 * ft*ft/Young*pf_l;

    //phasefield
    if (StressSplit==SplittingMethods::ICE){
        Ice_eRef = Ice_eRef0 + Ice_eReduceFact*std::sqrt(std::max(1e-3,-(TAbs-TRef)));
        Ice_cReduce = Ice_cReduce0 - (TAbs-TRef)*Ice_cReduceFact;
    }

    //plasticity
    sy = dam*(sy0 - sydeg*(TAbs-TRef));


}
