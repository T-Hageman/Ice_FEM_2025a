#ifndef PHASEFIELD_UTIL_H
#define PHASEFIELD_UTIL_H

#include <vector>
#include <iostream>
#include <unordered_map>
#include <Eigen/Dense>

#include "../../BaseModel/BaseModel.h"

//options for phase-field damage formulation
enum DamageMethods{NoDamage, Step2, Step5, Step8, LinearDamage, QuadraticDamage, CubicDamage, LinearDamageNR, QuadraticDamageNR, CubicDamageNR, CZMDamage, PorderDamage};
inline std::unordered_map<std::string, DamageMethods> DamageMethodNames = {
    {"None", NoDamage},
    {"Linear", LinearDamage},
    {"LinearNoRes", LinearDamageNR},
    {"Quadratic", QuadraticDamage},
    {"QuadraticNoRes", QuadraticDamageNR},
    {"Cubic", CubicDamage},
    {"CubicNoRes", CubicDamageNR},
    {"Porder", PorderDamage},
    {"CZM", CZMDamage},
    {"Step2", Step2},
    {"Step5", Step5},
    {"Step8", Step8}
};

//Distribution function
enum DistributionMethods{LinearDistribution, CZMDistribution, Linear2Distribution, QuadraticDistribution, HO_LinearDistribution, HO_Linear2Distribution, HO_CZMDistribution, HO_QuadraticDistribution, LocalDistribution, CustomDistribution};
inline std::unordered_map<std::string, DistributionMethods> DistributionMethodNames = {
    {"Linear", LinearDistribution},
    {"CZM", CZMDistribution},
    {"Quadratic", QuadraticDistribution},
    {"HO_Linear",HO_LinearDistribution},
    {"HO_Linear2",HO_Linear2Distribution},
    {"HO_CZM", HO_CZMDistribution},
    {"HO_Quadratic", HO_QuadraticDistribution},
    {"Local", LocalDistribution},
    {"Linear2", Linear2Distribution},
    {"Custom", CustomDistribution}
};


class PhaseFieldUtility{
    public:
        PhaseFieldUtility(inputData& inputs, std::string SolidMatName, std::string ModelName);
        PhaseFieldUtility(inputData& inputs, std::string SolidMatName);
        ~PhaseFieldUtility();

        //options for enforcing irreversible phase-field damage
        enum IrreversibilityMethods{NoIrreversibleEnforced, HistParameter, Penalty, LagrangeMultiplier, AugLagrangeMultiplier, AugLagrangeMultiplier2, AugmentedLagrangeMult3};
        std::unordered_map<std::string, IrreversibilityMethods> IrreversibilityMethodsNames = {
            {"None", NoIrreversibleEnforced},
            {"Hist", HistParameter},
            {"Penalty", Penalty},
            {"LagrangeMult", LagrangeMultiplier},
            {"AugmentedLagrangeMult", AugLagrangeMultiplier},
            {"AugmentedLagrangeMult2", AugLagrangeMultiplier2},
            {"AugmentedLagrangeMult3", AugmentedLagrangeMult3}
        };

        void EnforceIrreversible(Eigen::VectorXd& Phase, Eigen::VectorXd& PhaseOld, Eigen::VectorXd& WLumped, Eigen::VectorXd& L, Eigen::VectorXd& L2,
            Eigen::VectorXd& F_Ph, Eigen::VectorXd& F_L, Eigen::VectorXd& F_L2, 
            Eigen::MatrixXd& K_PhPh, Eigen::MatrixXd& K_LL, Eigen::MatrixXd& K_L2L2, Eigen::MatrixXd& K_PhL, Eigen::MatrixXd& K_LPh, Eigen::MatrixXd& K_PhL2, Eigen::MatrixXd& K_L2Ph);

        void AugLagrIrr(Eigen::VectorXd &F_L, size_t n, Eigen::VectorXd &WLumped, Eigen::VectorXd &L, Eigen::MatrixXd &K_LL, Eigen::VectorXd &Phase, Eigen::VectorXd &PhaseOld, std::vector<double> &PhaseLim, Eigen::VectorXd &F_Ph, Eigen::MatrixXd &K_PhPh, Eigen::MatrixXd &K_PhL, Eigen::MatrixXd &K_LPh);
        void PenaltyIrr(Eigen::VectorXd &Phase, size_t n, std::vector<double> &PhaseLim, Eigen::VectorXd &F_Ph, Eigen::VectorXd &WLumped, Eigen::MatrixXd &K_PhPh);
        void LagrMultIrr(Eigen::VectorXd &F_Ph, size_t n, Eigen::VectorXd &WLumped, Eigen::VectorXd &L, Eigen::MatrixXd &K_LL, Eigen::MatrixXd &K_PhL, Eigen::VectorXd &Phase, Eigen::VectorXd &F_L, std::vector<double> &PhaseLim, Eigen::MatrixXd &K_LPh, Eigen::VectorXd &PhaseOld);
        void AugLagrIrr2(Eigen::VectorXd &F_L, Eigen::VectorXd &F_L2, size_t n, Eigen::VectorXd &WLumped, Eigen::VectorXd &L, Eigen::VectorXd &L2,
            Eigen::MatrixXd &K_LL, Eigen::VectorXd &Phase, Eigen::VectorXd &PhaseOld, std::vector<double> &PhaseLim, Eigen::VectorXd &F_Ph, Eigen::MatrixXd &K_PhPh, Eigen::MatrixXd &K_PhL, Eigen::MatrixXd &K_LPh,
            Eigen::MatrixXd &K_L2L2, Eigen::MatrixXd &K_PhL2, Eigen::MatrixXd &K_L2Ph);

        IrreversibilityMethods IrreversibleMethod;
        double LagFact = 1.0e0;
        double PenaltyFactor = 1.0e6;

        void GenericDamageFunction(double phi, DamageMethods damfunc, double& D, double& dD, double& ddD);
        double GetDistributionFunction(Eigen::VectorXd& Dist, Eigen::MatrixXd& dDist, Eigen::RowVectorXd& N, Eigen::MatrixXd& G, Eigen::VectorXd& Phase);
        double GetDistributionFunction(Eigen::VectorXd& Dist, Eigen::MatrixXd& dDist, Eigen::RowVectorXd& N, Eigen::MatrixXd& G, Eigen::MatrixXd& G2, Eigen::VectorXd& Phase);

        void GetCapacityOnly(Eigen::VectorXd& Dist, Eigen::MatrixXd& dDist, double& DistLumped, double& dDistLumped, Eigen::RowVectorXd& N, double phi);

        bool NeedsGrad2;

    private:
        double l;
        double k0;
        double Gc;
        DistributionMethods DistFunctionMethod;

        //degradation parameters
        double p;
        size_t dim;

        //costum specific
        std::vector<double> Coefs; // w, a_1, a_2, b, c
};

#endif