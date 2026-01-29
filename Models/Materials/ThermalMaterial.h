#ifndef THERMALMATERIAL_H
#define THERMALMATERIAL_H

#include "BaseMaterial.h"
#include "../Fracture/PhaseField/PhaseFieldUtil.h"

class ThermalMaterial: virtual public BaseMaterial{
    public:
        ThermalMaterial(inputData& inputs, std::string SolidMatName);
        ~ThermalMaterial();

        double k;
        double cp;
        double LatentHeat;

    protected:

    private:

};


#endif