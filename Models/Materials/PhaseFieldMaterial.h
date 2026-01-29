#ifndef PHASEFIELDMATERIAL_H
#define PHASEFIELDMATERIAL_H

#include "SolidMaterial.h"

/// @brief Implements phase field splitting schemes
/**

Splitting schemes included via setting:
\code{.json}
    "properties":{
        "PhaseField":{
            "Split": "SpecStress"
            }
        }
\endcode

where split is either:
- "None": Do not perform any split, defining the elastic energy as 
        \f[ \Psi_\text{el}^+=\frac{1}{2}K\text{tr}(\varepsilon)^2 + \mu \text{tr}(\varepsilon^2) \qquad \Psi_\text{el}^-= 0 \f]
- "VolStrains": Split based on volumetric and deviatoric strains as \cite Miehe2010 
        \f[ \Psi_\text{el}^+=\frac{1}{2}K\left<\text{tr}(\varepsilon)^2\right>^+ + \mu \text{tr}(\varepsilon^2) \qquad \Psi_\text{el}^-= \frac{1}{2}K\left<\text{tr}(\varepsilon)^2\right>^- \f]
- "SpecStrains": Split based on principal strains as \cite Zhang2022 
        \f[ \Psi_\text{el}^+=\frac{1}{2}\lambda \left(\left< \varepsilon_1 + \varepsilon_2 + \varepsilon_3\right>^+\right)^2 + \mu \sum_i \left(\left<\varepsilon_i\right>^+\right)^2 \qquad \Psi_\text{el}^-= \frac{1}{2}\lambda \left(\left< \varepsilon_1 + \varepsilon_2 + \varepsilon_3\right>^-\right)^2 + \mu \sum_i \left(\left<\varepsilon_i\right>^-\right)^2 \f]
- "VolStress": Split based on volumetric and deviatoric stresses as \cite Zhang2022 
        \f[ \Psi_\text{el}^\pm = \frac{1}{2}\sigma^{\pm T} \varepsilon \qquad \sigma^+ = \frac{1}{3}\left<\text{tr}(\sigma)\right>^+\boldsymbol{I}+(\sigma - \frac{1}{3}\text{tr}(\sigma)\boldsymbol{I}) \qquad \sigma^- = \frac{1}{3}\left<\text{tr}(\sigma)\right>^-\boldsymbol{I} \f]
- "SpecStress": Split based on principal stresses as \cite Miehe2010 
        \f[  \Psi_\text{el}^\pm = \frac{1}{2}\sigma^{\pm T} \varepsilon \qquad \text{where} \sigma^+ = \left<\text{eig}(\sigma)\right>^+  \qquad \sigma^- = \left<\text{eig}(\sigma)\right>^-      \f]
- "DruckerPrager": Drucker-Prager like split, following \cite Navidtehrani2022
        \f[ \Psi_\text{el}^+=\begin{cases} \frac{1}{2}K I_1^2 + 2\mu J_2 \qquad \text{if} \qquad -6B\sqrt{J_2}<I_1 \\ 0 \qquad \text{if} \qquad -6\mu \sqrt{J_2}<3BKI_1 \\ \frac{1}{18B^2K+2\mu}\left( -3BKI_1+2\mu\sqrt{J_2}\right)^2\qquad \text{otherwise} \end{cases} \f]
        \f[ \Psi_\text{el}^-=\begin{cases} 0 \qquad \text{if} \qquad -6B\sqrt{J_2}<I_1 \\ \frac{1}{2}K I_1^2 + 2\mu J_2 \qquad \text{if} \qquad -6\mu \sqrt{J_2}<3BKI_1 \\ \frac{1}{18B^2K+2\mu}\left( I_1+6B\sqrt{J_2}\right)^2\qquad \text{otherwise} \end{cases} \f]
        where \f$ I_1 \f$ and \f$ J_2 \f$ are the first invariant of the principal strain, and second invariant of the deviatoric strain

*/
class PhaseFieldMaterial: virtual public SolidMaterial{
    public:
        PhaseFieldMaterial(inputData& inputs, std::string SolidMatName);
        ~PhaseFieldMaterial();

        // options for elastic energy splitting schemes
        enum SplittingMethods{NoSplit, VolStrains, SpecStrains, VolStress, SpecStress, DruckerPrager, ICE, StarConvex};
        std::unordered_map<std::string, SplittingMethods> SplittingMethodNames = {
            {"None", NoSplit},
            {"VolStrains", VolStrains},
            {"SpecStrains", SpecStrains},
            {"VolStress", VolStress},
            {"SpecStress", SpecStress},
            {"DruckerPrager", DruckerPrager},
            {"ICE", ICE},
            {"StarConvex", StarConvex}
        };

        double pf_l, Gc, pf_visc;

        double EnergySplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain);
        double B_dp;
        double Ice_eRef, Ice_cReduce;
        double StarConvexGamma;

        SplittingMethods StressSplit;        //stress split for fracture driving energy: "none", "VolStrains" or "SpecStrains" (all components included, Amor Vol/Dev split, Miehe Spectral split)
        bool ConsistentEnergy;
    protected:

    private:
        void GetEigenValues2D(Vector6d& s_xxyyzzxy, std::vector<double>& s, std::vector<Vector6d>& ds, std::vector<Matrix6d>& dds);
        void GetEigenVectors2D(Vector6d& s_xxyyzzxy, std::vector<Eigen::Vector3d>& EigVec);

        double VolStressSplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain);
        double SpecStressSplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain);
        double VolStrainSplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain);
        double SpecStrainSplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain);
        double DruckerPragerSplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain);
        double IceSplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain);
        double StarConvexSplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain);

};


#endif