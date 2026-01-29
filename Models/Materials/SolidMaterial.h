#ifndef SOLIDMATERIAL_H
#define SOLIDMATERIAL_H

#include "BaseMaterial.h"
#include "../Fracture/PhaseField/PhaseFieldUtil.h"

class SolidMaterial: virtual public BaseMaterial{
    public:
        SolidMaterial(inputData& inputs, std::string SolidMatName);
        ~SolidMaterial();

        // options for elastic energy splitting schemes
        enum SolidMaterialModel{VE, VE_expl, dVE, dVE_expl, VEP, dVEP, LE};
        std::unordered_map<std::string, SolidMaterialModel> SolidMaterialNames = {
            {"LinearElastic", LE},
            {"ViscoElastic", VE},
            {"ViscoElasticExpl", VE_expl},
            {"DamViscoElastic", dVE},
            {"DamViscoElasticExpl", dVE_expl},
            {"ViscoElasticPlastic", VEP},
            {"DamViscoElasticPlastic", dVEP}
        };

        double Young;
        double Poisson;
        double Bulk;
        double Lame;
        double Shear;
        double Damping;

        uint HistSize;

        SolidMaterialModel Rheology;
        double ACreep;
        double nCreep;
        double PhiViscLim;
        PhaseFieldUtility* PF_Util;
        DamageMethods DegFunction;

        double sy, yield_alpha, pot_alpha, hardening, pl_visc;

        Matrix6d D;  //Linear-elastic stiffness matrix
        Vector6d i;  // trace operator, [1 1 1 0]^T
        Matrix6d P;
        Matrix6d VoightCorrect, VoightDoubleDot;
        Matrix6d J2Mat_stress, J2Mat_strain;

        void UpdatePlasticStrains(Matrix6d& DMat, Vector6d& Strain_Pl, Vector6d& Strain, Vector6d& Strain_Pl_Old, Vector6d& StrainOld, 
                                  double phi, double phiOld, double dt, 
                                  std::vector<double>& hist, std::vector<double>& histOld, double& DissRate, double& DamDissRate);
    protected:

    private:
        double return_mapping_Viscous(Matrix6d& DMat, Vector6d& Strain_Pl, Vector6d& Strain, Vector6d& Strain_Pl_Old, double phi, double dt, 
                                      std::vector<double>& hist, std::vector<double>& histOld);
        void return_mapping_getKF(  Matrix7d& K, Vector7d& f,
                                    Vector7d& sol, Vector6d& strain_total,
                                    Vector6d& strain_ve_old,
                                    double Damage, double dt);

        void return_mapping_VEP(Matrix6d& DMat, Vector6d& Strain_Pl, Vector6d& Strain, Vector6d& Strain_Pl_Old, double dt, double Damage, 
                                    std::vector<double>& hist, std::vector<double>& histOld, double& DissRate, double& DamDissRate);
        void return_mapping_VEP_getKF(Eigen::Matrix<double, 8, 8>& K, Eigen::Matrix<double, 8, 1>& f,
                                    Eigen::Matrix<double, 8, 1>& sol, Vector6d& strain_total,
                                    Vector6d& strain_vep_old, double dt, double Damage, std::vector<double>& hist, std::vector<double>& histOld); 
        void return_mapping_VEP_updateStrain(Vector6d& strain_total, Vector6d& Strain_New, Vector6d& strain_vep_old, 
                                               Eigen::Matrix<double, 8, 1>& sol, double dt, double Damage, 
                                               std::vector<double>& hist, std::vector<double>& histOld, 
                                               double& DissRate, double& DamDissRate);
        double YieldFunction(Vector6d& Stress, Vector6d& TrialStress, Vector6d& dF, Matrix6d& ddF);
        void PotentialFunction(Vector6d& Stress, Vector6d& TrialStress, Vector6d& dG, Matrix6d& ddG);
        double YieldFunctionVE(Vector6d& Stress, Vector6d& dF, Matrix6d& ddF);
};


#endif