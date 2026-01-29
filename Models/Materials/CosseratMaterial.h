#ifndef COSSERATMATERIAL_H
#define COSSERATMATERIAL_H

#include "PhaseFieldMaterial.h"

class CosseratMaterial: virtual public PhaseFieldMaterial{
    public:
    //xx/yy/zz/xy/yx/yz/zy/xz/zx  0-8: strain, 9-17: curv

        CosseratMaterial(inputData& inputs, std::string SolidMatName);
        ~CosseratMaterial();

        Matrix18d D_coss;
        Vector18d ic;
        Vector9d icl;
        Matrix18d J2Mat_Coss;
        Matrix9d SkewMat;
        Eigen::Matrix<double, 9, 18> MacroStrains, MicroStrains;

        double G_coss, l_coss;
        double Hardening_e_ref;

        using PhaseFieldMaterial::UpdatePlasticStrains;
        void UpdatePlasticStrains(Matrix18d& DMat, Vector18d& Strain_Pl, Vector18d& Strain, Vector18d& Strain_Pl_Old, Vector18d& StrainOld, 
                                  double phi, double phiOld, double dt, 
                                  std::vector<double>& hist, std::vector<double>& histOld, double& DissRate, double& DamDissRate);

        using PhaseFieldMaterial::EnergySplit;
        double EnergySplit(Matrix18d& D_Dam, Matrix18d& D_UnDam, Vector18d& Stress_Dam, Vector18d& Stress_Undam, Vector18d& Strain);
    protected:

    private:
        double SpecStressSplit(Matrix18d& D_Dam, Matrix18d& D_UnDam, Vector18d& Stress_Dam, Vector18d& Stress_Undam, Vector18d& Strain);
        double SpecStressSplitConsistent(Matrix18d& D_Dam, Matrix18d& D_UnDam, Vector18d& Stress_Dam, Vector18d& Stress_Undam, Vector18d& Strain);
        double DruckerPragerSplit(Matrix18d& D_Dam, Matrix18d& D_UnDam, Vector18d& Stress_Dam, Vector18d& Stress_Undam, Vector18d& Strain);
        double StarConvexSplit(Matrix18d& D_Dam, Matrix18d& D_UnDam, Vector18d& Stress_Dam, Vector18d& Stress_Undam, Vector18d& Strain);

        void return_mapping_VEP(Matrix18d& DMat, Vector18d& Strain_Pl, Vector18d& Strain, Vector18d& Strain_Pl_Old, double dt, double Damage, 
                                    std::vector<double>& hist, std::vector<double>& histOld, double& DissRate, double& DamDissRate);
        void return_mapping_VEP_BackUp(Matrix20d& K, Vector20d& sol, Vector18d& Strain, Vector18d& Strain_Pl_Old, double dt, double Damage, 
                                    std::vector<double>& histOld);
        void return_mapping_VEP_getKF(Matrix20d& K, Vector20d& f,
                                      Vector20d& sol, Vector18d& strain_total,
                                      Vector18d& strain_vep_old, double dt, double Damage, std::vector<double>& histOld); 
        void return_mapping_VEP_updateStrain(Vector18d& strain_total, Vector18d& Strain_New, Vector18d& strain_vep_old, 
                                               Vector20d& sol, double dt, double Damage, 
                                               std::vector<double>& hist, std::vector<double>& histOld, 
                                               double& DissRate, double& DamDissRate);
        double YieldFunction(Vector18d& Stress, Vector18d& TrialStress, Vector18d& dF, Matrix18d& ddF);
        double YieldFunctionVE(Vector18d& Stress, Vector18d& dF, Matrix18d& ddF);
        void PotentialFunction(Vector18d& Stress, Vector18d& TrialStress, Vector18d& dG, Matrix18d& ddG);

        double sc=1e6;
        double scale1 = 1.0e3;
        double scale2 = -1.0e-6;
        double scale3 = -1.0e-6;

        double eps = 1.0e-12;
};


#endif