#ifndef FLUIDMATERIAL_H
#define FLUIDMATERIAL_H

#include "BaseMaterial.h"
#include "../Fracture/PhaseField/PhaseFieldUtil.h"

class FluidMaterial: virtual public BaseMaterial{
    public:
        FluidMaterial(inputData& inputs, std::string SolidMatName);
        ~FluidMaterial();

        double visc;
        double Bulk;

    protected:

    private:

};


#endif