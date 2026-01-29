#ifndef COSSPHASEFIELDFRACTURE_H
#define COSSPHASEFIELDFRACTURE_H

#include <vector>
#include <iostream>
#include <unordered_map>
#include <Eigen/Dense>

#include "../../BaseModel/BaseModel.h"
#include "PhaseFieldUtil.h"
#include "../../Materials/IceMaterial.h"

class CossPhaseField: public BaseModel{
    public:
        std::vector<std::string> DofNames_u = {"ux","uy","oz"};      //names for displacement degrees of freedom
        std::string DofNames_phase = "phase";                   //name of phase field degree of freedom
        std::string DofName_LagMult = "lMult_pf";
        std::string DofName_T = "T";
        std::string DofName_poros = "poros";
        std::vector<std::string> DofName_LagMult2 = {"lMult_pf","lMult_pf2"};

        CossPhaseField(Physics& My_Physics, std::string MyName);
        ~CossPhaseField();

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

        size_t ElemGroupIndex_phase, ElemGroupIndex_u, ElemGroupIndex_T, ElemGroupIndex_poros;  //index of the element groups
        std::vector<size_t> dofTypes_u;                 // dof indices for displacements
        size_t dofType_phase, dofType_lMult, dofType_T, dofType_poros;                           // dof index for phasefield
        std::vector<size_t> dofType_lMult2;                 // dof indices for displacements
        size_t Step_u, Step_phase, Step_T, Step_poros;                      // steps in which displacements and phasefield are resolved (NOT THE SAME STEP)

        std::string SolidMatName;   //name of solid material
        size_t nDisp = 3; 

        double viscNum;     //artificial additional viscosity
        double k0;          //residual stiffness
        Vector6d g={0, -9.81, 0, 0, 0, 0};   //gravity vector
        DamageMethods DegFunction, ViscoDegFunction;        //degradation function used for stresses.
        DamageMethods DensityDegFunction; //density degradation function,
        DamageMethods GravityDegFunction; //density degradation function,
        double ViscDamageFact;

        bool VariablePorosity;
        bool ConsistentEnergy;
        bool ConsistentDegradation;
        bool IncludeInertia;
        bool ViscousHeating, HasT;

        void Assemble_U(Vec &f, Mat &K);
        void Assemble_Phase(Vec &f, Mat &K);
        void Assemble_T(Vec &f, Mat &K);

        void ConstructBMatCoss(Eigen::MatrixXd& B, Eigen::RowVectorXd& N, Eigen::MatrixXd& G);
        void ConstructNMatCoss(Eigen::MatrixXd& Nu, Eigen::RowVectorXd& N);

        std::vector<std::vector<double>> Energy_Hist, Energy_HistOld;            //driving energy history field
        std::vector<std::vector<Vector18d>> Plastic_Strain, Plastic_StrainOld;    //plastic strain history
        std::vector<std::vector<double>> ViscDiss, ViscDissOld, DamViscDiss, DamViscDissOld;
        std::vector<std::vector<std::vector<double>>> MatHist, MatHistOld;

        void GetDamFunc(double phi, double& D, double& dD, double& ddD);
        void GetDamFuncInertia(double phi, double& D, double& dD, double& ddD);
        void GetDamFuncGravity(double phi, double& D, double& dD, double& ddD);
        void GetDamFuncVisc(double phi, double& D, double& dD, double& ddD);

        double LCrack, LCrackOld;
};

void Register_CossPhaseField();
BaseModel* New_CossPhaseField(Physics& My_Physics, std::string MyNameIn);

#endif

