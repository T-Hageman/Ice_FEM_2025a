#include "SolidMaterial.h"

#include "../../Physics/physics.h"

SolidMaterial::SolidMaterial(inputData &inputs, std::string MatName): BaseMaterial(inputs, MatName){
    inputs.GetRequired(Young, {"properties", MatName, "SolidProperties", "Young"});
    inputs.GetRequired(Poisson, {"properties", MatName, "SolidProperties", "Poisson"});

    Shear  = Young/(2.0*(1+Poisson));
    Lame   = Young*Poisson/((1.0+Poisson)*(1.0-2.0*Poisson));
    Bulk = Lame + 2.0/3.0*Shear;

    HistSize = 0;

    Damping = 0.0;
    inputs.GetOptional(Damping, {"properties", MatName, "SolidProperties", "Damping"});

    std::string RheologyName = "LinearElastic";
    inputs.GetOptional(RheologyName, {"properties", MatName, "SolidProperties", "Rheology"});
    if (SolidMaterialNames.count(RheologyName)){
        Rheology = SolidMaterialNames[RheologyName];
    } else {
        throw std::invalid_argument("Solid material model type "+RheologyName+" not defined,");
    }

    if (Rheology == SolidMaterialModel::VE || Rheology == SolidMaterialModel::dVE || 
        Rheology == SolidMaterialModel::VE_expl || Rheology == SolidMaterialModel::dVE_expl ||
        Rheology == SolidMaterial::VEP || Rheology == SolidMaterial::dVEP){
        HistSize = 1;

        inputs.GetRequired(ACreep, {"properties", MatName, "SolidProperties", "A"});
        inputs.GetRequired(nCreep, {"properties", MatName, "SolidProperties", "n"});
        
        PhiViscLim = 1.0e3;
        inputs.GetOptional(PhiViscLim, {"properties", MatName, "SolidProperties", "PhiViscLim"});

        if (Rheology == SolidMaterialModel::dVE || Rheology == SolidMaterialModel::dVE_expl || Rheology == SolidMaterialModel::dVEP){
            PF_Util = new PhaseFieldUtility(inputs, MatName);

            std::string DegFunctionName;
            inputs.GetRequired(DegFunctionName, {"properties", "PhaseField", "DamFunction"});
            if (DamageMethodNames.count(DegFunctionName)){
                DegFunction = DamageMethodNames[DegFunctionName];
            } else {
                throw std::invalid_argument("Damage function type "+DegFunctionName+" not defined,");
            }
        }
    }

    if (Rheology==SolidMaterialModel::VEP || Rheology==SolidMaterialModel::dVEP){
        HistSize = 4;

        inputs.GetRequired(sy, {"properties", MatName, "SolidProperties", "sy"});
        yield_alpha = 0.0;
        inputs.GetOptional(yield_alpha, {"properties", MatName, "SolidProperties", "F_deg"});
        pot_alpha = yield_alpha;
        inputs.GetOptional(pot_alpha, {"properties", MatName, "SolidProperties", "G_deg"});
        hardening = 0.0;
        inputs.GetOptional(hardening, {"properties", MatName, "SolidProperties", "hardening"});

        yield_alpha = M_PI/180.0 * yield_alpha; yield_alpha = 6.0*std::sin(yield_alpha)/(3.0-std::sin(yield_alpha));
        pot_alpha = M_PI/180.0 * pot_alpha; pot_alpha = 6.0*std::sin(pot_alpha)/(3.0-std::sin(pot_alpha));
        //sy = sy * std::cos(yield_alpha)/(1.0-1.0/3.0*std::sin(yield_alpha));

        pl_visc = 0.0;
        inputs.GetOptional(pl_visc, {"properties", MatName, "SolidProperties", "pl_visc"});
    }

    //utility matrices
    i.setZero(); i(0) = 1.0; i(1) = 1.0; i(2) = 1.0;

    VoightCorrect.setZero();
    VoightCorrect(0,0) = 1.0;
    VoightCorrect(1,1) = 1.0;
    VoightCorrect(2,2) = 1.0;
    VoightCorrect(3,3) = 0.5;
    VoightCorrect(4,4) = 0.5;
    VoightCorrect(5,5) = 0.5;

    VoightDoubleDot.setZero();
    VoightDoubleDot(0,0) = 1.0;
    VoightDoubleDot(1,1) = 1.0;
    VoightDoubleDot(2,2) = 1.0;
    VoightDoubleDot(3,3) = 2.0;
    VoightDoubleDot(4,4) = 2.0;
    VoightDoubleDot(5,5) = 2.0;

    P = Matrix6d::Identity() - 1.0/3.0*i*i.transpose();

    D.setZero();
    D = Lame*i*i.transpose()+2*Shear*Matrix6d::Identity(); 
    D(3,3) = Shear;
    D(4,4) = Shear;
    D(5,5) = Shear;

    J2Mat_stress.setZero();
    J2Mat_stress = Matrix6d::Identity()-1.0/3.0*i*i.transpose();
    J2Mat_stress(3,3) = 2.0;
    J2Mat_stress(4,4) = 2.0;
    J2Mat_stress(5,5) = 2.0;

    J2Mat_strain.setZero();
    J2Mat_strain = Matrix6d::Identity()-1.0/3.0*i*i.transpose();
    J2Mat_strain(3,3) = 2.0/4.0; //correcting for voight notation in strains
    J2Mat_strain(4,4) = 2.0/4.0;
    J2Mat_strain(5,5) = 2.0/4.0;
}

SolidMaterial::~SolidMaterial(){

}




/// @brief updates the plastic strains, based on total stresses and history
/// @param DMat output: Tangent matrix
/// @param Strain_Pl output/input: Current plastic strain (output, and used as initial guess)
/// @param Strain input: total strains
/// @param Strain_Pl_Old input history parameter: plastic strains at the end of last time increment
/// @param StrainOld input history: total strains at the end of the last time increment
/// @param phi input: current damage
void SolidMaterial::UpdatePlasticStrains(Matrix6d& DMat, Vector6d& Strain_Pl, Vector6d& Strain, Vector6d& Strain_Pl_Old, Vector6d& StrainOld,
                         double phi, double phiOld, double dt, std::vector<double>& hist, std::vector<double>& histOld,
                         double& DissRate, double& DamDissRate){
    DissRate = 0.0;
    DamDissRate = 0.0;
    double damage = 1.0;
    double Unused1, Unused2;

    switch (Rheology){
    case SolidMaterialModel::LE:
            //no plastic strains, no need to update
            DMat = Matrix6d::Identity();
            Strain_Pl = 1.0*Strain_Pl_Old;
            DissRate = 0.0;
        break;
    case SolidMaterial::VE:
            if (phi>PhiViscLim){
                //no additional plastic strains, no need to update
                DMat = Matrix6d::Identity();
                Strain_Pl = 1.0*Strain_Pl_Old;
            } else {
                DissRate = return_mapping_Viscous(DMat, Strain_Pl, Strain, Strain_Pl_Old, 1.0, dt, hist, histOld);
            }
        break;
    case SolidMaterial::dVE:
            if (phi>PhiViscLim){
                //no additional plastic strains, no need to update
                DMat = Matrix6d::Identity();
                Strain_Pl = 1.0*Strain_Pl_Old;
            } else {
                PF_Util->GenericDamageFunction(phi, DegFunction, damage, Unused1, Unused2);
                DissRate = return_mapping_Viscous(DMat, Strain_Pl, Strain, Strain_Pl_Old, damage, dt, hist, histOld);
            }
        break;
    case SolidMaterial::VE_expl:
            if (phiOld>PhiViscLim){
                //no additional plastic strains, no need to update
                DMat = Matrix6d::Identity();
                Strain_Pl = 1.0*Strain_Pl_Old;
            } else {
                DissRate = return_mapping_Viscous(DMat, Strain_Pl, StrainOld, Strain_Pl_Old, 1.0, dt, hist, histOld);
                DMat = Matrix6d::Identity();
            }
        break;
    case SolidMaterial::dVE_expl:
            if (phiOld>PhiViscLim){
                //no additional plastic strains, no need to update
                DMat = Matrix6d::Identity();
                Strain_Pl = 1.0*Strain_Pl_Old;
            } else {
                PF_Util->GenericDamageFunction(phiOld, DegFunction, damage, Unused1, Unused2);
                DissRate = return_mapping_Viscous(DMat, Strain_Pl, StrainOld, Strain_Pl_Old, damage, dt, hist, histOld);
                DMat = Matrix6d::Identity();
            }
        break;
    case SolidMaterial::VEP:
        if (phiOld>PhiViscLim){
            //no additional plastic strains, no need to update
            DMat = Matrix6d::Identity();
            Strain_Pl = 1.0*Strain_Pl_Old;
        } else {
            //PF_Util->GenericDamageFunction(phiOld, DegFunction, damage, Unused1, Unused2);
            return_mapping_VEP(DMat, Strain_Pl, Strain, Strain_Pl_Old, dt, damage, hist, histOld, DissRate, DamDissRate);
        }
        break;
    default:
        throw std::invalid_argument("Rheology type not defined in SolidMaterial.cpp,");
        break;
    }
}

/// @brief Performs the return mapping required to obtain plastic strains using a visco-elastic rheology
/// @param DMat output: tangent matrix
/// @param Strain_Pl //in/output: plastic strains (and initial guess for these)
/// @param Strain //input: total strains
/// @param Strain_Pl_Old //input (history): Plastic strains at the end of the previous time increment
/// @param phi input:
double SolidMaterial::return_mapping_Viscous(Matrix6d& DMat, Vector6d& Strain_Pl, Vector6d& Strain, Vector6d& Strain_Pl_Old, double Damage, double dt, std::vector<double>& hist, std::vector<double>& histOld){
    double DissipatedRate = 0.0;

    Matrix7d K, Kinv;
    Vector7d f, fOld;
    Vector7d sol, dsol;
    double e, e0, lineSearch, ls1, ls2;
    uint it;

    sol(Eigen::seq(0,5)) = Strain-Strain_Pl;
    sol(6) = hist[0];

    return_mapping_getKF(K, f, sol, Strain, Strain_Pl_Old, Damage, dt);
    e0 = f.transpose()*f;
    e0 = std::max(e0,1e-12);

    bool conv = false;
    it = 0;
    while (conv == false){
        it += 1;
        if (false){
            dsol = K.fullPivLu().solve(f); //This does not seem to work? always returns zero
            assert(f.isApprox(K*dsol)); 
            dsol *= -1.0;
        } else {
            dsol = -K.inverse()*f;
        }
        sol += dsol;

        fOld = f;
        return_mapping_getKF(K, f, sol, Strain, Strain_Pl_Old, Damage, dt);
        ls1 = f.transpose()*dsol;
        ls2 = fOld.transpose()*dsol;
        if ((ls1-ls2)==0.0){
            lineSearch = 1.0;
        } else {
            lineSearch = -(ls2)/(abs(ls1-ls2)+1e-12)*sgn(ls1-ls2);
		    lineSearch = std::min(std::max(lineSearch, 0.2), 1.0);
        }
        sol = sol-(1.0-lineSearch)*dsol;

        return_mapping_getKF(K, f, sol, Strain, Strain_Pl_Old, Damage, dt);
        e = f.transpose()*f;
        if (e<1.0e-16 || e/e0<1.0e-9){
            conv = true;
        }
        if (it>250){
            std::string msg = "Viscous return mapping did not provide a solution within 250 iterations: "+std::to_string(e)+"/"+std::to_string(e0)+"\n";
            Logs.PrintEvery(msg,2);
            conv = true;
        }
    }

    Vector6d Stress = D*sol(Eigen::seq(0,5));
    Vector6d dF; 
    Matrix6d unused;
    YieldFunctionVE(Stress, dF, unused);

    Strain_Pl = Strain_Pl_Old + sol(6)*dF;

    Kinv = K.inverse();
    DMat = Kinv(Eigen::seq(0,5),Eigen::seq(0,5));

    //dissipation rate: \sigma^T*\dot{\varepsilon_p}
    DissipatedRate = Damage*Stress.transpose()*(Strain_Pl-Strain_Pl_Old);
    DissipatedRate = std::abs(DissipatedRate)/dt;

    hist[0] = sol(6);
    return DissipatedRate;
}

/// @brief Part of visco-elastic return mapping, provides the residual and tangent
/// @param K output: tangent matrix
/// @param f output: residual vector
/// @param sol input: current solution state
/// @param strain_total input: total strains
/// @param strain_ve_old input: plastic strains at previous time increment
/// @param phi input:
void SolidMaterial::return_mapping_getKF(  Matrix7d& K, Vector7d& f,
                                              Vector7d& sol, Vector6d& strain_total,
                                              Vector6d& strain_ve_old,
                                              double Damage, double dt){
    Vector6d Stress = D*sol(Eigen::seq(0,5));

    double F_VE;
    Vector6d dF_VE;
    Matrix6d ddF_VE;

    F_VE = YieldFunctionVE(Stress, dF_VE, ddF_VE);
    ddF_VE.setZero(); 

    f(Eigen::seq(0,5)) = sol(Eigen::seq(0,5)) - strain_total + strain_ve_old + sol(6)*dF_VE;
    K(Eigen::seq(0,5),Eigen::seq(0,5)) = Matrix6d::Identity() + Damage*sol(6)*ddF_VE*D; 
    K(Eigen::seq(0,5),6) = dF_VE;

    double sc = 1.0e6;
    double PreFac = ACreep*dt*std::pow(sc, nCreep);
    f(6) = PreFac*std::pow(F_VE/sc, nCreep) - sol(6);
    K(6,6) = - 1.0;
    K(6,Eigen::seq(0,5)) = nCreep*PreFac/sc*std::pow(F_VE/sc, nCreep-1.0)*( Damage*dF_VE.transpose()*D );

}

void SolidMaterial::return_mapping_VEP(Matrix6d& DMat, Vector6d& Strain_Pl, Vector6d& Strain, Vector6d& Strain_Pl_Old, double dt, double Damage, 
                                        std::vector<double>& hist, std::vector<double>& histOld, double& DissRate, double& DamDissRate){
    Eigen::Matrix<double, 8, 8> K, Kinv;
    Eigen::Matrix<double, 8, 1> f, fOld;
    Eigen::Matrix<double, 8, 1> sol, dsol;
    double e, e0, lineSearch, ls1, ls2;
    uint it;

    sol(Eigen::seq(0,5)) = Strain-Strain_Pl; //stress
    sol(6) = hist[1]; // dLambda_visc
    sol(7) = hist[2]; // dLambda_p

    return_mapping_VEP_getKF(K, f, sol, Strain, Strain_Pl_Old, dt, Damage, hist, histOld);
    e0 = f.transpose()*f;
    e0 = std::max(e0,1e-12);

    bool conv = false;
    it = 0;
    while (conv == false){
        it += 1;
        if (false){
            dsol = K.fullPivLu().solve(f); //This does not seem to work? always returns zero
            assert(f.isApprox(K*dsol)); 
            dsol *= -1.0;
        } else {
            dsol = -K.inverse()*f;
        }
        sol += dsol;

        fOld = f;
        return_mapping_VEP_getKF(K, f, sol, Strain, Strain_Pl_Old, dt, Damage, hist, histOld);
        ls1 = f.transpose()*dsol;
        ls2 = fOld.transpose()*dsol;
        lineSearch = -(ls2)/(abs(ls1-ls2)+1e-12)*sgn(ls1-ls2);
		lineSearch = std::min(std::max(lineSearch, 0.2), 1.0);
        sol = sol-(1.0-lineSearch)*dsol;

        return_mapping_VEP_getKF(K, f, sol, Strain, Strain_Pl_Old, dt, Damage, hist, histOld);

        e = f.transpose()*f;
        if (e<1.0e-16 || e/e0<1.0e-8){
            conv = true;
        }
        if (it>250){
            std::string msg = "ViscoPlastic return mapping did not provide a solution within 250 iterations: "+std::to_string(e)+"/"+std::to_string(e0)+"\n";
            Logs.PrintEvery(msg,2);
            conv = true;
        }
    }

    return_mapping_VEP_updateStrain(Strain, Strain_Pl, Strain_Pl_Old, sol, dt, Damage, hist, histOld, DissRate, DamDissRate);

    Kinv = K.inverse();
    DMat = Kinv(Eigen::seq(0,5),Eigen::seq(0,5));
}

void SolidMaterial::return_mapping_VEP_getKF(Eigen::Matrix<double, 8, 8>& K, Eigen::Matrix<double, 8, 1>& f,
                            Eigen::Matrix<double, 8, 1>& sol, Vector6d& strain_total,
                            Vector6d& strain_vep_old, double dt, double Damage, std::vector<double>& hist, std::vector<double>& histOld){
    f.setZero();
    K.setZero();

    Vector6d Stress = Damage*D*sol(Eigen::seq(0,5));          
    Vector6d SDev = P*Stress;            
    Vector6d StressTrial = Damage*D*(strain_total-strain_vep_old-sol(6)*SDev);

    double F, FTrial;
    Vector6d dF, dG;
    Matrix6d ddF, ddG;

    FTrial = YieldFunction(StressTrial, StressTrial, dF, ddF);
    F = YieldFunction(Stress, StressTrial, dF, ddF);
    PotentialFunction(Stress, StressTrial, dG, ddG);

    double h = histOld[0] + hardening*std::sqrt(sol[7]*sol[7]);
    double dh = hardening*sgn(sol[7]);
    if (h<-0.9*sy){
        h = -0.9*sy;
        dh = -0.0;
    }

    // 0 = s - D(strain - strain_old - dL_v s_dev - dLp df/ds)
    // 0 = e_el - strain + strain_old + dL_v s_dev + dLp df/ds
    f(Eigen::seq(0,5)) = sol(Eigen::seq(0,5)) - strain_total + strain_vep_old + sol(6)*SDev + sol(7)*dG;
    K(Eigen::seq(0,5),Eigen::seq(0,5)) = Matrix6d::Identity() + Damage*sol(6)*P*D + Damage*sol(7)*ddG*D;
    K(Eigen::seq(0,5),6) = SDev;
    K(Eigen::seq(0,5),7) = dG;

    // 0 = dL_v - A*dt*(|sDev|)^(n-1)
    f(6) = sol(6) - ACreep*dt*std::pow(SDev.transpose()*SDev, (nCreep-1.0)/2.0);
    K(6,6) = 1.0;
    K(6,Eigen::seq(0,5)) = - ACreep*dt*std::pow(SDev.transpose()*SDev, (nCreep-3.0)/2.0 ) * (nCreep-1.0)/2.0 * 2.0 * Damage*SDev.transpose()*P*D;

    // 0 = F;
    if (FTrial-histOld[0]<0.0){
        f(7) = -(pl_visc+1.0e0)/dt*sol(7);
        K(7,7) = -(pl_visc+1.0e0)/dt;
    } else {
       f(7) = F - h - pl_visc/dt*std::sqrt(sol(7)*sol(7)+1.0e-12);
       K(7, Eigen::seq(0,5)) = Damage*dF.transpose()*D;
       K(7,7) = -pl_visc/dt*sgn(sol(7)) - dh;
    }
    //std::cout << f << "\n" << K << "\n\n\n";
}

void SolidMaterial::return_mapping_VEP_updateStrain(Vector6d& strain_total, Vector6d& Strain_New, Vector6d& strain_vep_old, 
                                               Eigen::Matrix<double, 8, 1>& sol, double dt, double Damage, 
                                               std::vector<double>& hist, std::vector<double>& histOld, double& DissRate, double& DamDissRate){
    Vector6d Stress = Damage*D*sol(Eigen::seq(0,5));   
    Vector6d SDev = P*Stress;                           
    Vector6d StressTrial = Damage*D*(strain_total-strain_vep_old-sol(6)*SDev);

    Vector6d dG;
    Matrix6d ddG;
    Vector6d Strain_VE, Strain_Pl;

    PotentialFunction(Stress, StressTrial, dG, ddG);

    double DissV = 0.0, DissP = 0.0;
    Strain_New = strain_vep_old;

    // Creep
    Strain_VE = sol(6)*SDev;
    Strain_New += Strain_VE;

    //Plasticity
    Strain_Pl = sol(7)*dG;
    Strain_New += Strain_Pl;

    //thermal dissipation
    for (size_t i = 0; i < 6; i++){
        DissV += 1.0/dt*std::abs(Strain_VE[i]*Stress[i]);
        DissP += 1.0/dt*std::abs(Strain_Pl[i]*Stress[i]);
    }

    //hardening
    hist[0] = histOld[0] + hardening*std::sqrt(sol(7)*sol(7));
    hist[1] = sol(6);
    hist[2] = sol(7);

    DissRate = DissV+DissP;
    DamDissRate = DissP;
}

double SolidMaterial::YieldFunction(Vector6d& Stress, Vector6d& TrialStress, Vector6d& dF, Matrix6d& ddF){
    double F;
    
    double I1 = i.transpose()*Stress;
    Vector6d dI1ds = i;

    double J2 = 0.5*Stress.transpose()*J2Mat_stress*Stress;
    Vector6d dJ2 = J2Mat_stress*Stress;
    Matrix6d ddJ2 = J2Mat_stress;
    
    bool Apex = false;
    if (yield_alpha>0.0){
        double I1Apex = sy*3.0/yield_alpha;
        if (i.transpose()*TrialStress>I1Apex) Apex = true;
    }

    if (Apex==false){
        F = std::sqrt(3.0*J2+1.0e-6)+yield_alpha/3.0*I1-sy;
        dF  = 0.5*std::sqrt(3.0)/std::sqrt(J2+1.0e-6) * dJ2 + yield_alpha/3.0*dI1ds;
        ddF = 0.5*std::sqrt(3.0)/std::sqrt(J2+1.0e-6)*ddJ2
            -0.25*std::sqrt(3.0)*std::pow(J2+1.0e-6,-1.5)*dJ2*dJ2.transpose();
    } else {
        F = yield_alpha/3.0*I1-sy;
        dF = yield_alpha/3.0*dI1ds;
        ddF.setZero();
    }

    return F;
}

double SolidMaterial::YieldFunctionVE(Vector6d& Stress, Vector6d& dF, Matrix6d& ddF){
    double F;
    double sc = 1.0e6;
    double eps = 1.0e-12;

    Vector6d Stress_Scaled = Stress/sc;
    double J2 = 0.5*Stress_Scaled.transpose()*J2Mat_stress*Stress_Scaled;
    Vector6d dJ2 = J2Mat_stress*Stress_Scaled;
    Matrix6d ddJ2 = J2Mat_stress;
    
    F = sc*std::sqrt(2.0*J2);
    dF  = 0.5*std::sqrt(2.0)/std::sqrt(J2+eps) * dJ2;
    ddF = 0.5/sc*std::sqrt(3.0)/std::sqrt(J2+eps)*(ddJ2-0.5/(J2+eps)*dJ2*dJ2.transpose());

    return F;
}

void SolidMaterial::PotentialFunction(Vector6d& Stress, Vector6d& TrialStress, Vector6d& dG, Matrix6d& ddG){
    double I1 = i.transpose()*Stress;
    Vector6d dI1ds = i;

    double J2 = 0.5*Stress.transpose()*J2Mat_stress*Stress;
    Vector6d dJ2 = J2Mat_stress*Stress;
    Matrix6d ddJ2 = J2Mat_stress;
    
    bool Apex = false;
    if (yield_alpha>0.0){
        double I1Apex = sy*3.0/yield_alpha;
        if (i.transpose()*TrialStress>I1Apex) Apex = true;
    }

    if (Apex==false){
        dG  = 0.5*std::sqrt(3.0)/std::sqrt(J2+1.0e-6) * dJ2 + pot_alpha/3.0*dI1ds;
        ddG = 0.5*std::sqrt(3.0)/std::sqrt(J2+1.0e-6)*ddJ2
            -0.25*std::sqrt(3.0)*std::pow(J2+1.0e-6,-1.5)*dJ2*dJ2.transpose();
    } else {
        Vector6d dStress = TrialStress - 3.0*sy/yield_alpha*i;
        Vector6d dStrain = D.inverse()*dStress;

        dG = dStrain;
        ddG.setZero();
    }

    dG.normalize();
}