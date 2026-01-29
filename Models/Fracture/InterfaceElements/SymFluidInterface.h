#ifndef SYMFLUIDINTERFACE_H
#define SYMFLUIDINTERFACE_H

#include <vector>
#include <iostream>
#include <unordered_map>
#include <Eigen/Dense>

#include "../../BaseModel/BaseModel.h"
#include "../../Materials/SolidMaterial.h"
#include "../../Materials/FluidMaterial.h"

class SymFluidInterfaceModel: public BaseModel{
    public:
        std::vector<std::string> DofNames_u = {"ux","uy"};      //names for displacement degrees of freedom
        std::string DofName_p = "p";           //names for pressure degrees of freedom

        SymFluidInterfaceModel(Physics& My_Physics, std::string MyName);
        ~SymFluidInterfaceModel();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);
        void ResetStep();

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat &K, Vec &f, Constrainer *cons, size_t step);
        void Commit(int CommitType);

        bool SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>> &Data);
        size_t hasTimeData(std::vector<std::string>& DataNames);
        void GetTimeData(std::vector<double>& DataValues);
    protected:

    private:
        void Assemble_UP(Vec &f, Mat &K);

        double GetSaturation(double p, double& dSw, double& ddSw);
        void GetFluidFlux(Eigen::Vector2d gradP, double h, double Sw, Eigen::Vector2d&q, Eigen::Matrix2d &dq_dGradP, Eigen::Vector2d &dq_dh, Eigen::Vector2d &dq_dSw, Eigen::Vector2d &gVec);
        enum SaturationModel{ExponentSaturation, ExponentSmoothSaturation, TanhSaturation, FullSaturation, LinearSaturation};
        std::unordered_map<std::string, SaturationModel> SaturationModels = {
            {"Exponent", ExponentSaturation},
            {"ExponentSmooth", ExponentSmoothSaturation},
            {"Tanh", TanhSaturation},
            {"Saturated", FullSaturation},
            {"Linear", LinearSaturation}
        };
        double sat_exp_pref = 1.0e4, sat_exp_diffexp = 1.0;
        SaturationModel satModel;

        bool LumpedCap = true;

        enum FractureCriteria{SingleIP, AllIP};
        std::unordered_map<std::string, FractureCriteria> FracCriteriaNames = {
            {"SingleIP", SingleIP},
            {"AllIP", AllIP}
        };

        FractureCriteria fracCrit = SingleIP;

        std::vector<size_t> TopNodes, BotNodes;
        size_t TopCentre_p, BotCentre_p, TopCentre_u, BotCentre_u;

        std::string FlowModel;
        std::string SolidMatName, FluidMatName;   //name of solid material
        SolidMaterial* S;
        FluidMaterial* F;

        double kDummy, Gc, ft, u0;
        double hOver, hBase, Z, X;
        bool TopBC = false, BotBC = false;
        double Q_Top, Q_Bot;

        double LCrack, LMeanCrack;
        double h_Top, h_Bot;
        double p_Top, p_Bot;

        size_t ElemGroupIndex_u, ElemGroupIndex_p;  //index of the element groups
        std::vector<size_t> dofTypes_u;                 // dof indices for displacements
        size_t dofType_p;
        size_t Step_u, Step_p;                      // steps in which displacements and phasefield are resolved (NOT THE SAME STEP)

        Eigen::Vector3d g={0, -9.81, 0};   //gravity vector
        std::vector<std::vector<double>> hMax, hMaxOld;
        std::vector<bool> Broken, BrokenOld;
};

void Register_SymFluidInterfaceModel();
BaseModel* New_SymFluidInterfaceModel(Physics& My_Physics, std::string MyNameIn);

#endif
