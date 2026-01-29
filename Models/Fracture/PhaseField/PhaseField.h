#ifndef PHASEFIELDFRACTURE_H
#define PHASEFIELDFRACTURE_H

#include <vector>
#include <iostream>
#include <unordered_map>
#include <Eigen/Dense>

#include "../../BaseModel/BaseModel.h"
#include "PhaseFieldUtil.h"
#include "../../Materials/IceMaterial.h"

/// @brief Model to solve phase field fracture combined with a dynamic momentum balance.
/**
Required input parameters, within the "Models" node:

\code{.json}
"PhaseField1":{
    "Name": "Fracture/PhaseField",
    "ElementGroup_u": "internal",
    "ElementGroup_phase": "internal",
    "Material": "Ice",
    "ViscousStab": 0.0
},
\endcode

and within the parameters node:

\code{.json}
"PhaseField":{
        "_Comment": "Properties for fracture phase-field: length scale l [m], interfacial energy Gc [J/m^2]",
        "DamFunction": "quadratic",
        "DensityDamage": "none",
        "Split": "SpecStrains",
        "l": 2.5,
        "k0": 1e-6,
        "Viscosity": 0.0
    },
"Ice":{
        "_Comment": "Properties describing the metal: intact material density Density [kg/m^3]",
        "Density": 910,
        "Young": 5.0e9,
        "Poisson": 0.35,
        "Gc": 500.0,
        "Rheology":{
            "Type": "ViscoElastic",
            "A": 5.0e-24,
            "n": 3.0
        }
\endcode

*/
class PhaseFieldDamage: public BaseModel{
    public:
        std::vector<std::string> DofNames_u = {"ux","uy"};      //names for displacement degrees of freedom
        std::string DofNames_phase = "phase";                   //name of phase field degree of freedom
        std::string DofName_LagMult = "lMult_pf";
        std::string DofName_T = "T";
        std::vector<std::string> DofName_LagMult2 = {"lMult_pf","lMult_pf2"};

        PhaseFieldDamage(Physics& My_Physics, std::string MyName);
        ~PhaseFieldDamage();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);
        void ResetStep();

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat &K, Vec &f, Constrainer *cons, size_t step);
        void Commit(int CommitType);

        bool SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>> &Data);
        void SaveStressComponent(std::string DataName, std::vector<std::vector<double>> &Data, bool Damaged);
        void SaveStrainComponent(std::string DataName, std::vector<std::vector<double>> &Data);

        size_t hasTimeData(std::vector<std::string> &DataNames);
        void GetTimeData(std::vector<double>& DataValues);
    protected:

    private:
        IceMaterial* M;
        PhaseFieldUtility* PF_Util;

        size_t ElemGroupIndex_phase, ElemGroupIndex_u, ElemGroupIndex_T;  //index of the element groups
        std::vector<size_t> dofTypes_u;                 // dof indices for displacements
        size_t dofType_phase, dofType_lMult, dofType_T;                           // dof index for phasefield
        std::vector<size_t> dofType_lMult2;                 // dof indices for displacements
        size_t Step_u, Step_phase, Step_T;                      // steps in which displacements and phasefield are resolved (NOT THE SAME STEP)

        std::string SolidMatName;   //name of solid material

        double viscNum;     //artificial additional viscosity
        double k0;          //residual stiffness
        Eigen::Vector3d g={0, -9.81, 0};   //gravity vector
        DamageMethods DegFunction, ViscoDegFunction;        //degradation function used for stresses.
        DamageMethods DensityDegFunction; //density degradation function,
        DamageMethods GravityDegFunction; //density degradation function,
        double ViscDamageFact;

        bool ConsistentEnergy;
        bool ConsistentDegradation;
        bool IncludeInertia;
        bool LumpedCap, LumpedDrive;
        bool ViscousHeating, HasT;

        void Assemble_U(Vec &f, Mat &K);
        void Assemble_Phase(Vec &f, Mat &K);
        void Assemble_T(Vec &f, Mat &K);

        // void AugLagrIrr3(Eigen::VectorXd &Phase, Eigen::VectorXd &PhaseOld, Eigen::VectorXd &L, Eigen::VectorXd &L2, 
        //                 Eigen::VectorXd &F_L, Eigen::VectorXd &F_L2, Eigen::VectorXd &F_Ph, 
        //                 Eigen::MatrixXd &K_PhL, Eigen::MatrixXd &K_PhL2, Eigen::MatrixXd &K_LPh, Eigen::MatrixXd &K_L2Ph, Eigen::MatrixXd &K_PhPh, Eigen::MatrixXd &K_LL, Eigen::MatrixXd &K_L2L2, 
        //                 double w_ip, Eigen::RowVectorXd &Nph_ip, Eigen::MatrixXd &Gph_ip, Eigen::MatrixXd &G2ph_ip);


        std::vector<std::vector<double>> Energy_Hist, Energy_HistOld;            //driving energy history field
        std::vector<std::vector<Vector6d>> Plastic_Strain, Plastic_StrainOld;    //plastic strain history
        std::vector<std::vector<double>> ViscDiss, ViscDissOld, DamViscDiss, DamViscDissOld;
        std::vector<std::vector<std::vector<double>>> MatHist, MatHistOld;

        void GetDamFunc(double phi, double& D, double& dD, double& ddD);
        void GetDamFuncInertia(double phi, double& D, double& dD, double& ddD);
        void GetDamFuncGravity(double phi, double& D, double& dD, double& ddD);
        void GetDamFuncVisc(double phi, double& D, double& dD, double& ddD);

        double LCrack, LCrackOld;
};

void Register_PhaseFieldDamage();
BaseModel* New_PhaseFieldDamage(Physics& My_Physics, std::string MyNameIn);

#endif

