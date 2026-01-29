#ifndef ICEMATERIAL_H
#define ICEMATERIAL_H

#include "PhaseFieldMaterial.h"
#include "CosseratMaterial.h"
#include "ThermalMaterial.h"

class IceMaterial: virtual public CosseratMaterial, ThermalMaterial{
    public:
        IceMaterial(inputData& inputs, std::string SolidMatName);
        ~IceMaterial();

        void SetTemperature(double T, double poros);

        bool TDependent;
        double TRef, TOffset;
        double A0, QCreep;
        double ft0, ftdeg;
        double sy0, sydeg;
        double G_coss0;
        double Young0;

        Matrix18d D_coss0;
        Matrix6d D0;

        double Ice_eRef0, Ice_cReduce0, Ice_cReduceFact, Ice_eReduceFact;

    protected:

    private:

};
#endif