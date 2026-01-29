#include "PhaseFieldUtil.h"


PhaseFieldUtility::PhaseFieldUtility(inputData& inputs, std::string SolidMatName, std::string ModelName){
    inputs.GetRequired(l, {"properties", SolidMatName, "FractureProperties", "l"});
    inputs.GetRequired(Gc, {"properties", SolidMatName, "FractureProperties", "Gc"});
    inputs.GetRequired(k0, {"properties", "PhaseField", "k0"});
    inputs.GetRequired(dim, {"mesh", "dim"});

    //distribution function
    std::string DistributionName;
    DistFunctionMethod = QuadraticDistribution;
    bool distdefined = inputs.GetOptional(DistributionName, {"properties", "PhaseField", "DistributionFunction"});
    if (distdefined){
        DistFunctionMethod = DistributionMethodNames[DistributionName];
    }

    if (DistFunctionMethod==HO_LinearDistribution || DistFunctionMethod==HO_Linear2Distribution || DistFunctionMethod==HO_QuadraticDistribution || DistFunctionMethod==HO_CZMDistribution || DistFunctionMethod==CustomDistribution){
        NeedsGrad2 = true;
    } else {
        NeedsGrad2 = false;
    }

    p=2;
    inputs.GetOptional(p, {"properties", "PhaseField", "p"});

    if (DistFunctionMethod==CustomDistribution){
        Coefs.resize(5);
        inputs.GetRequired(Coefs, {"properties", "PhaseField", "Coefs"});
    }

    //irreversibility enforcement
    std::string IrreversibleMethodName = "Hist";
    inputs.GetOptional(IrreversibleMethodName, {"Models", ModelName, "Irreversible"});
    if (IrreversibilityMethodsNames.count(IrreversibleMethodName)){
        IrreversibleMethod = IrreversibilityMethodsNames[IrreversibleMethodName];
    } else {
        std::invalid_argument(ModelName+" has no valid method to enforce irreversible cracks \n");
    }
}

PhaseFieldUtility::PhaseFieldUtility(inputData& inputs, std::string SolidMatName){
    inputs.GetRequired(l, {"properties", SolidMatName, "FractureProperties", "l"});
    inputs.GetRequired(Gc, {"properties", SolidMatName, "FractureProperties", "Gc"});
    inputs.GetRequired(k0, {"properties", "PhaseField", "k0"});
    inputs.GetRequired(dim, {"mesh", "dim"});

    //distribution function
    std::string DistributionName;
    DistFunctionMethod = QuadraticDistribution;
    bool distdefined = inputs.GetOptional(DistributionName, {"properties", "PhaseField", "DistributionFunction"});
    if (distdefined){
        DistFunctionMethod = DistributionMethodNames[DistributionName];
    }

    if (DistFunctionMethod==HO_LinearDistribution || DistFunctionMethod==HO_Linear2Distribution || DistFunctionMethod==HO_QuadraticDistribution || DistFunctionMethod==HO_CZMDistribution || DistFunctionMethod==CustomDistribution){
        NeedsGrad2 = true;
    } else {
        NeedsGrad2 = false;
    }

    p=2;
    inputs.GetOptional(p, {"properties", "PhaseField", "p"});

    if (DistFunctionMethod==CustomDistribution){
        Coefs.resize(5);
        inputs.GetRequired(Coefs, {"properties", "PhaseField", "Coefs"});
    } 

    IrreversibleMethod = IrreversibilityMethods::AugLagrangeMultiplier;
}

PhaseFieldUtility::~PhaseFieldUtility(){

}

void PhaseFieldUtility::GenericDamageFunction(double phi, DamageMethods damfunc, double& D, double& dD, double& ddD){
    double phiLim = phi;
    if (phiLim<0.0){
        phiLim = 0.0;
    }
    if (phiLim>1.0){
        phiLim = 1.0;
    }

    switch (damfunc){
        case NoDamage:{
            D   = 1.0;
            dD  = 0.0;
            ddD = 0.0;
        } break;
        case LinearDamage:{
            D   =      (1.0-k0)*(1.0-phiLim)+k0;
            dD  = -1.0*(1.0-k0);
            ddD =  0.0;
        } break;
        case LinearDamageNR:{
            D   =  (1.0-phiLim);
            dD  = -1.0;
            ddD =  0.0;
        } break;
        case QuadraticDamage:{
            D   =      (1.0-k0)*std::pow((1.0-phiLim),2)+k0;
            dD  = -2.0*(1.0-k0)*(1.0-phi);
            ddD =  2.0*(1.0-k0);
        } break;
        case QuadraticDamageNR:{
            D   =      std::pow((1.0-phiLim),2);
            dD  = -2.0*(1.0-phi);
            ddD =  2.0;
        } break;
        case CubicDamage:{
            D   =      (1.0-k0)*std::pow((1.0-phiLim),3)+k0;
            dD  = -3.0*(1.0-k0)*(1.0-phi)*(1.0-phi);
            ddD =  6.0*(1.0-k0)*(1.0-phi);
        } break;
        case CubicDamageNR:{
            D   =      std::pow((1.0-phiLim),3);
            dD  = -3.0*(1.0-phi)*(1.0-phi);
            ddD =  6.0*(1.0-phi);
        } break;
        case PorderDamage:{
            D   =    (1.0-k0)*std::pow((1.0-phiLim),p)+k0;
            dD  = -p*(1.0-k0)*std::pow((1.0-phi),p-1);
            if (p>1){
                ddD =  p*(p-1)*(1.0-k0)*std::pow((1.0-phi),p-2);
            } else {
                ddD = 0.0;
            }
        } break;
        case CZMDamage:{
            //D   =      (1.0-k0)* ((1-phi)*(1-phi)/((1-phi)*(1-phi)+a1*phi*(1.0+a2*phi+a2*a3*phi*phi))) + k0;
            //double c1 = a1*(1.0-phi)*(a2*phi*(a3*(phi-3)*phi-2.0)-phi-1.0);
            //double c2 = std::max(1.0e-9, a1*phi*(a2*phi*(a3*phi+1.0)+1.0)+(phi-1.0)*(phi-1.0));
            //dD  = (1.0-k0)* c1/c2/c2;

            //double c3 = 2.0*a1*(a1*(a2*phi*(a3*((2.0-3.0*phi)*phi+3.0)*phi+a2*phi*(a3*phi*(a3*((phi-6.0)*phi+6.0)*phi-6*phi+8.0)-2.0*phi+3.0)-phi*phi+3.0)+1.0)-(phi-1.0)*(phi-1.0)*(a2*(3.0*a3*phi+2.0*phi+1.0)+phi+2.0));
            //ddD =  (1-k0)*c3/c2/c2/c2;

            // double s = 1.05;
            // double fac = (s-1.0)/s;

            // D = s*( 1.0-std::pow(fac, (1.0-phi)*(1.0-phi)) );
            // dD = 2.0*s*(1.0-phi)*std::pow(fac, (1.0-phi)*(1.0-phi))*std::log(fac);
            // ddD = -2.0*s*std::pow(fac, (1.0-phi)*(1.0-phi))*std::log(fac)*(1.0+2.0*(1.0-phi)*(1.0-phi)*std::log(fac));

            D = (1.0-k0)*(2.0*std::pow(1-phi,2)-1.0*std::pow(1-phi,3)) + k0;
            dD= (1.0-k0)*(3.0*(1-phi)*(1-phi) - 4.0*(1-phi));
            ddD = (1.0-k0)*(-6.0*(1-phi) + 4.0);

        } break;
        case Step2:{
            if (phi<0.2){
                D = 1.0;
            } else {
                D = 0.0;
            }
            dD = 0.0;
            ddD = 0.0;
        } break;
        case Step5:{
            if (phi<0.5){
                D = 1.0;
            } else {
                D = 0.0;
            }
            dD = 0.0;
            ddD = 0.0;
        } break;
        case Step8:{
            if (phi<0.8){
                D = 1.0;
            } else {
                D = 0.0;
            }
            dD = 0.0;
            ddD = 0.0;
        } break;
        default:{
            throw std::invalid_argument("Damage function type "+std::to_string(damfunc)+" not defined in phasefieldUtil.cpp,");
        }
    }
}

double PhaseFieldUtility::GetDistributionFunction(Eigen::VectorXd& Dist, Eigen::MatrixXd& dDist, Eigen::RowVectorXd& N, Eigen::MatrixXd& G, Eigen::VectorXd& Phase){
    double gamma;
    double phi = N*Phase;
    Eigen::VectorXd dG(dim); dG = G*Phase;
    switch (DistFunctionMethod){
        case QuadraticDistribution:{
            gamma = 1.0/(2.0*l)*(phi*phi             + l*l*dG.transpose()*dG);
            Dist  = 1.0/(2.0*l)*(2.0*N.transpose()*phi + 2.0*l*l*G.transpose()*dG);
            dDist = 1.0/(2.0*l)*(2.0*N.transpose()*N   + 2.0*l*l*G.transpose()*G);
        } break;
        case CZMDistribution:{
            gamma = 1.0/(3.141592*l)*(2.0*phi-phi*phi             + l*l*dG.transpose()*dG);
            Dist  = 1.0/(3.141592*l)*(2.0*N.transpose()*(1.0-phi) + 2.0*l*l*G.transpose()*dG);
            dDist = 1.0/(3.141592*l)*(-2.0*N.transpose()*N        + 2.0*l*l*G.transpose()*G);
        } break;
        case LinearDistribution:{
            gamma = 3.0/(8.0*l)*(phi + l*l*dG.transpose()*dG);
            Dist  = 3.0/(8.0*l)*(N.transpose() + 2.0*l*l*G.transpose()*dG);
            dDist = 3.0/(8.0*l)*(                2.0*l*l*G.transpose()*G);
        } break;
        case Linear2Distribution:{
            gamma = 1.0/(4.0*l)*(2.0*phi + l*l*dG.transpose()*dG);
            Dist  = 1.0/(4.0*l)*(2.0*N.transpose() + 2.0*l*l*G.transpose()*dG);
            dDist = 1.0/(4.0*l)*(                2.0*l*l*G.transpose()*G);
        } break;
        case LocalDistribution:{
            gamma = 1.0/(l)*(phi*phi);
            Dist  = 1.0/(l)*(2*N.transpose()*phi);
            dDist = 1.0/(l)*(2*N.transpose()*N);
        } break;
        default:{
            throw std::invalid_argument("Distribution function type "+std::to_string(DistFunctionMethod)+" not defined in phasefield.cpp,");
        }
    }

    return gamma;
}

double PhaseFieldUtility::GetDistributionFunction(Eigen::VectorXd& Dist, Eigen::MatrixXd& dDist, Eigen::RowVectorXd& N, Eigen::MatrixXd& G, Eigen::MatrixXd& G2, Eigen::VectorXd& Phase){
    double gamma;
    double phi = N*Phase;
    Eigen::VectorXd dG(dim); dG = G*Phase;

    size_t dimG2 = (dim-1)*3;
    Eigen::VectorXd IM(dimG2); 
    Eigen::MatrixXd DoubleCon(dimG2,dimG2); 
    if (dim==2){
        IM.setZero(); IM(0) = 1.0; IM(1) = 1.0;
        DoubleCon.setZero(); DoubleCon(0,0) = 1.0; DoubleCon(1,1) = 1.0; DoubleCon(2,2) = 2.0;
    } else {
        IM.setZero(); IM(0) = 1.0; IM(1) = 1.0; IM(2) = 1.0;
        DoubleCon.setZero(); 
        DoubleCon(0,0) = 1.0; DoubleCon(1,1) = 1.0; DoubleCon(2,2) = 1.0;
        DoubleCon(3,3) = 2.0; DoubleCon(4,4) = 2.0; DoubleCon(5,5) = 2.0;
    }

    Eigen::VectorXd ddG(dimG2);
    ddG = G2*Phase;
    switch (DistFunctionMethod){
        case HO_QuadraticDistribution:{
            gamma = 1.0/(2.0*l)*(phi*phi             + 0.5*l*l*dG.transpose()*dG + l*l*l*l/16.0*ddG.transpose()*DoubleCon*ddG);
            Dist  = 1.0/(2.0*l)*(2.0*N.transpose()*phi + l*l*G.transpose()*dG    + l*l*l*l/8.0*G2.transpose()*DoubleCon*ddG);
            dDist = 1.0/(2.0*l)*(2.0*N.transpose()*N + l*l*G.transpose()*G       + l*l*l*l/8.0*G2.transpose()*DoubleCon*G2);
        } break;
        case HO_CZMDistribution:{
            gamma = 1.0/(3.020*l)*(2.0*phi-phi*phi             + 0.5*l*l*dG.transpose()*dG + l*l*l*l/16.0*ddG.transpose()*DoubleCon*ddG);
            Dist  = 1.0/(3.020*l)*(2.0*N.transpose()*(1.0-phi) + l*l*G.transpose()*dG       + l*l*l*l/8.0*G2.transpose()*DoubleCon*ddG);
            dDist = 1.0/(3.020*l)*(-2.0*N.transpose()*N        + l*l*G.transpose()*G       + l*l*l*l/8.0*G2.transpose()*DoubleCon*G2);
        } break;
        case HO_LinearDistribution:{
            gamma = 1.0/(19.2*l)*(9.0*phi            + l*l*l*l*ddG.transpose()*DoubleCon*ddG);
            Dist  = 1.0/(19.2*l)*(9.0*N.transpose()  + 2.0*l*l*l*l*G2.transpose()*DoubleCon*ddG);
            dDist = 1.0/(19.2*l)*(                   + 2.0*l*l*l*l*G2.transpose()*DoubleCon*G2);
        } break;
        case HO_Linear2Distribution:{
            // gamma = 3.0/(8.0*l)*(phi             + 0.5*l*l*dG.transpose()*dG + l*l*l*l/16.0*ddG.transpose()*DoubleCon*ddG);
            // Dist  = 3.0/(8.0*l)*(N.transpose()   + l*l*G.transpose()*dG      + l*l*l*l/8.0*G2.transpose()*DoubleCon*ddG);
            // dDist = 3.0/(8.0*l)*(                  l*l*G.transpose()*G       + l*l*l*l/8.0*G2.transpose()*DoubleCon*G2);

            gamma = 1.0/(4.0*l)*(2.0*phi             + 0.5*l*l*dG.transpose()*dG + l*l*l*l/16.0*ddG.transpose()*DoubleCon*ddG);
            Dist  = 1.0/(4.0*l)*(2.0*N.transpose()   + l*l*G.transpose()*dG      + l*l*l*l/8.0*G2.transpose()*DoubleCon*ddG);
            dDist = 1.0/(4.0*l)*(                  l*l*G.transpose()*G       + l*l*l*l/8.0*G2.transpose()*DoubleCon*G2);

        } break;
        case CustomDistribution:{
            gamma = 1.0/(Coefs[0]*l)*(Coefs[1]*phi          +    Coefs[2]*phi*phi           +     Coefs[3]*l*l*dG.transpose()*dG +     Coefs[4]*l*l*l*l*ddG.transpose()*DoubleCon*ddG);
            Dist  = 1.0/(Coefs[0]*l)*(Coefs[1]*N.transpose()+2.0*Coefs[2]*N.transpose()*phi + 2.0*Coefs[3]*l*l*G.transpose()*dG  + 2.0*Coefs[4]*l*l*l*l*G2.transpose()*DoubleCon*ddG);
            dDist = 1.0/(Coefs[0]*l)*(                       2.0*Coefs[2]*N.transpose()*N   + 2.0*Coefs[3]*l*l*G.transpose()*G   + 2.0*Coefs[4]*l*l*l*l*G2.transpose()*DoubleCon*G2);
        } break;
        default:{
            gamma = GetDistributionFunction(Dist, dDist, N, G, Phase);
        }
    }

    return gamma;
}

void PhaseFieldUtility::GetCapacityOnly(Eigen::VectorXd& Dist, Eigen::MatrixXd& dDist, double& DistLumped, double& dDistLumped, Eigen::RowVectorXd& N, double phi){
    switch (DistFunctionMethod){
        case QuadraticDistribution:{
            Dist  = 1.0/(2.0*l)*(2.0*N.transpose()*phi );
            dDist = 1.0/(2.0*l)*(2.0*N.transpose()*N   );
            DistLumped  = 1.0/l*phi;
            dDistLumped = 1.0/l;
        } break;
        case CZMDistribution:{
            Dist  = 1.0/(3.141592*l)*(2.0*N.transpose()*(1.0-phi) );
            dDist = 1.0/(3.141592*l)*(-2.0*N.transpose()*N        );
            DistLumped  = 2.0/(3.141592*l)*(1.0-phi);
            dDistLumped = -2.0/(3.141592*l);
        } break;
        case LinearDistribution:{
            Dist  = 3.0/(8.0*l)*(N.transpose() );
            dDist.setZero();
            DistLumped  = 3.0/(8.0*l);
            dDistLumped = 0.0;
        } break;
        case Linear2Distribution:{
            Dist  = 1.0/(4.0*l)*(2.0*N.transpose() );
            dDist.setZero();
            DistLumped  = 1.0/(4.0*l);
            dDistLumped = 0.0;
        } break;
        case LocalDistribution:{
            Dist  = 1.0/(l)*(2*N.transpose()*phi);
            dDist = 1.0/(l)*(2*N.transpose()*N);
            DistLumped  = 2.0/(l)*phi;
            dDistLumped = 2.0/l;
        } break;
        case HO_QuadraticDistribution:{
            Dist  = 1.0/(2.0*l)*(2.0*N.transpose()*phi );
            dDist = 1.0/(2.0*l)*(2.0*N.transpose()*N   );
            DistLumped  = 1.0/(l)*phi;
            dDistLumped = 1.0/l;
        } break;
        case HO_CZMDistribution:{
            Dist  = 1.0/(3.020*l)*(2.0*N.transpose()*(1.0-phi));
            dDist = 1.0/(3.020*l)*(-2.0*N.transpose()*N       );
            DistLumped  = 2.0/(3.020*l)*(1-phi);
            dDistLumped = -2.0/(3.020*l);
        } break;
        case HO_LinearDistribution:{
            Dist  = 1.0/(19.2*l)*(9.0*N.transpose());
            dDist.setZero();
            DistLumped  = 9.0/(19.2*l);
            dDistLumped = 0.0;
        } break;
        case HO_Linear2Distribution:{
            Dist  = 3.0/(8.0*l)*N.transpose();
            dDist.setZero();
            DistLumped  = 3.0/(8.0*l);
            dDistLumped = 0.0;
        } break;
        case CustomDistribution:{
            Dist  = 1.0/(Coefs[0]*l)*(Coefs[1]*N.transpose()+2.0*Coefs[2]*N.transpose()*phi);
            dDist = 1.0/(Coefs[0]*l)*(                       2.0*Coefs[2]*N.transpose()*N  );
            DistLumped  = 1.0/(Coefs[0]*l)*(Coefs[1]+2.0*Coefs[2]*phi);
            dDistLumped = 1.0/(Coefs[0]*l)*          2.0*Coefs[2];
        } break;
        default:{
            throw std::invalid_argument("Distribution function type "+std::to_string(DistFunctionMethod)+" not defined in phasefield.cpp,");
        }
    }
}

void PhaseFieldUtility::EnforceIrreversible(Eigen::VectorXd& Phase, Eigen::VectorXd& PhaseOld, Eigen::VectorXd& WLumped, Eigen::VectorXd& L, Eigen::VectorXd& L2,
    Eigen::VectorXd& F_Ph, Eigen::VectorXd& F_L, Eigen::VectorXd& F_L2, 
    Eigen::MatrixXd& K_PhPh, Eigen::MatrixXd& K_LL, Eigen::MatrixXd& K_L2L2, Eigen::MatrixXd& K_PhL, Eigen::MatrixXd& K_LPh, Eigen::MatrixXd& K_PhL2, Eigen::MatrixXd& K_L2Ph){

    for (size_t n = 0; n < Phase.size(); n++){
        std::vector<double> PhaseLim(2); PhaseLim[0] = std::min(1.0-1.0e-3, std::max(0.0, PhaseOld(n))); PhaseLim[1] = 1.0+1e-1;

        if (IrreversibleMethod==Penalty){
            PenaltyIrr(Phase, n, PhaseLim, F_Ph, WLumped, K_PhPh);
        }
        if (IrreversibleMethod==LagrangeMultiplier){
            LagrMultIrr(F_Ph, n, WLumped, L, K_LL, K_PhL, Phase, F_L, PhaseLim, K_LPh, PhaseOld);
        }
        if (IrreversibleMethod==AugLagrangeMultiplier){
            AugLagrIrr(F_L, n, WLumped, L, K_LL, Phase, PhaseOld, PhaseLim, F_Ph, K_PhPh, K_PhL, K_LPh);
        }
        if (IrreversibleMethod==AugLagrangeMultiplier2){
            AugLagrIrr2(F_L, F_L2, n, WLumped, L, L2, K_LL, Phase, PhaseOld, PhaseLim, F_Ph, K_PhPh, K_PhL, K_LPh, K_L2L2, K_PhL2, K_L2Ph);
        }
    }
}

void PhaseFieldUtility::PenaltyIrr(Eigen::VectorXd &Phase, size_t n, std::vector<double> &PhaseLim, Eigen::VectorXd &F_Ph, Eigen::VectorXd &WLumped, Eigen::MatrixXd &K_PhPh){
    if (Phase(n) < PhaseLim[0]){
        F_Ph(n) += WLumped(n) * Gc * PenaltyFactor * (Phase(n) - PhaseLim[0]);
        K_PhPh(n, n) += WLumped(n) * Gc * PenaltyFactor;
    }
    if (Phase(n) > PhaseLim[1]){
        F_Ph(n) += WLumped(n) * Gc * PenaltyFactor * (Phase(n) - PhaseLim[1]);
        K_PhPh(n, n) += WLumped(n) * Gc * PenaltyFactor;
    }
}

void PhaseFieldUtility::LagrMultIrr(Eigen::VectorXd &F_Ph, size_t n, Eigen::VectorXd &WLumped, Eigen::VectorXd &L, Eigen::MatrixXd &K_LL, Eigen::MatrixXd &K_PhL, Eigen::VectorXd &Phase, Eigen::VectorXd &F_L, std::vector<double> &PhaseLim, Eigen::MatrixXd &K_LPh, Eigen::VectorXd &PhaseOld){
    F_Ph(n) += WLumped(n) * L(n);
    K_PhL(n, n) += WLumped(n);
    if (Phase(n) >= 1.0){
        // E = L*(Phase-1)
        F_L(n) += WLumped(n) * (Phase(n) - PhaseLim[1]);
        K_LPh(n, n) += WLumped(n);
    } else if (Phase(n) < PhaseOld(n)){
        // E = L*<Phase-PhaseOld>
        F_L(n) += WLumped(n) * (Phase(n) - PhaseLim[0]);
        K_LPh(n, n) += WLumped(n);
    } else {
        F_L(n) += WLumped(n) * L(n);
        K_LL(n, n) += WLumped(n);
    }
}

void PhaseFieldUtility::AugLagrIrr(Eigen::VectorXd &F_L, size_t n, Eigen::VectorXd &WLumped, Eigen::VectorXd &L, Eigen::MatrixXd &K_LL, Eigen::VectorXd &Phase, Eigen::VectorXd &PhaseOld, std::vector<double> &PhaseLim, Eigen::VectorXd &F_Ph, Eigen::MatrixXd &K_PhPh, Eigen::MatrixXd &K_PhL, Eigen::MatrixXd &K_LPh){
    double preFac, preFac2, dpreFac;

    double scaleFact =  Gc / l / LagFact;

    F_L(n) += -WLumped(n) * scaleFact * L(n) * (1.0+1.0e-10);
    K_LL(n, n) += -WLumped(n) * scaleFact * (1.0+1.0e-10);

    preFac  = L(n) + LagFact * (Phase(n) - PhaseLim[0]);
    preFac2 = L(n) + LagFact * (Phase(n) - PhaseLim[1]);
    dpreFac = LagFact;

    F_Ph(n)      += WLumped(n) * scaleFact * L(n)  * dpreFac;
    K_PhL(n, n)  += WLumped(n) * scaleFact * dpreFac;

    //if (Phase(n)<PhaseLim[1]){
        //E = 1/(2gamma)*(L+gamma*phiDot/dXdt)^2-1/2gamma * L^2   //see e.g. Geelen et al, A phase-field formulation for dynamic cohesive fracture
        //E = 1/2/gamma*preFac^2

        if (preFac <= 0){
            F_Ph(n)      += WLumped(n) * scaleFact * (preFac-L(n))  * dpreFac;
            K_PhPh(n, n) += WLumped(n) * scaleFact * dpreFac * dpreFac;

            F_L(n)      += WLumped(n) * scaleFact * preFac;
            K_LL(n, n)  += WLumped(n) * scaleFact;
            K_LPh(n, n) += WLumped(n) * scaleFact * dpreFac;
        }
    //} else {
        // if (preFac2 >= 0){
        //     F_Ph(n)      += WLumped(n) * scaleFact * (preFac2-L(n)) * dpreFac;
        //     K_PhPh(n, n) += WLumped(n) * scaleFact * dpreFac * dpreFac;
        //     //K_PhL(n, n) += WLumped(n) * Gc / l / LagFact * dpreFac;

        //     F_L(n)      += WLumped(n) * scaleFact * preFac2;
        //     K_LL(n, n)  += WLumped(n) * scaleFact;
        //     K_LPh(n, n) += WLumped(n) * scaleFact * dpreFac;
        // }
    //}
}


void PhaseFieldUtility::AugLagrIrr2(Eigen::VectorXd &F_L, Eigen::VectorXd &F_L2, size_t n, Eigen::VectorXd &WLumped, Eigen::VectorXd &L, Eigen::VectorXd &L2,
            Eigen::MatrixXd &K_LL, Eigen::VectorXd &Phase, Eigen::VectorXd &PhaseOld, std::vector<double> &PhaseLim, Eigen::VectorXd &F_Ph, Eigen::MatrixXd &K_PhPh, Eigen::MatrixXd &K_PhL, Eigen::MatrixXd &K_LPh,
            Eigen::MatrixXd &K_L2L2, Eigen::MatrixXd &K_PhL2, Eigen::MatrixXd &K_L2Ph){
    double preFac, preFac2, dpreFac;
    dpreFac = LagFact;

    F_L(n) += -WLumped(n) * Gc / l / LagFact * L(n) * (1.0-1.0e-8);
    K_LL(n, n) += -WLumped(n) * Gc / l / LagFact * (1.0-1.0e-8);
    F_L2(n) += -WLumped(n) * Gc / l / LagFact * L2(n) * (1.0-1.0e-8);
    K_L2L2(n, n) += -WLumped(n) * Gc / l / LagFact * (1.0-1.0e-8);

    F_Ph(n)      += WLumped(n) * Gc / l / LagFact * L(n)  * dpreFac;
    K_PhL(n, n)  += WLumped(n) * Gc / l / LagFact * dpreFac;
    F_Ph(n)         += WLumped(n) * Gc / l / LagFact * L2(n) * dpreFac;
    K_PhL2(n, n)    += WLumped(n) * Gc / l / LagFact * dpreFac;

    // lower limit
    preFac = L(n) + LagFact * (Phase(n) - PhaseLim[0]);
    if (preFac <= 0){
        F_Ph(n)      += WLumped(n) * Gc / l / LagFact * (preFac-L(n))  * dpreFac;
        K_PhPh(n, n) += WLumped(n) * Gc / l / LagFact * dpreFac * dpreFac;

        F_L(n)      += WLumped(n) * Gc / l / LagFact * preFac;
        K_LL(n, n)  += WLumped(n) * Gc / l / LagFact;
        K_LPh(n, n) += WLumped(n) * Gc / l / LagFact * dpreFac;
    }// else {
        //upper limit
        preFac2 = L2(n) + LagFact * (Phase(n) - PhaseLim[1]);
        if (preFac2 >= 0){
            F_Ph(n)         += WLumped(n) * Gc / l / LagFact * (preFac2-L2(n)) * dpreFac;
            K_PhPh(n, n)    += WLumped(n) * Gc / l / LagFact * dpreFac * dpreFac;

            F_L2(n)         += WLumped(n) * Gc / l / LagFact * preFac2;
            K_L2L2(n, n)    += WLumped(n) * Gc / l / LagFact;
            K_L2Ph(n, n)    += WLumped(n) * Gc / l / LagFact * dpreFac;
        }
    //}
}

// void PhaseFieldDamage::AugLagrIrr3(Eigen::VectorXd &Phase, Eigen::VectorXd &PhaseOld, Eigen::VectorXd &L, Eigen::VectorXd &L2, 
//                             Eigen::VectorXd &F_L, Eigen::VectorXd &F_L2, Eigen::VectorXd &F_Ph, 
//                             Eigen::MatrixXd &K_PhL, Eigen::MatrixXd &K_PhL2, Eigen::MatrixXd &K_LPh, Eigen::MatrixXd &K_L2Ph, Eigen::MatrixXd &K_PhPh, Eigen::MatrixXd &K_LL, Eigen::MatrixXd &K_L2L2, 
//                             double w_ip, Eigen::RowVectorXd &Nph_ip, Eigen::MatrixXd &Gph_ip, Eigen::MatrixXd &G2ph_ip){
//     uint nNodes = Nph_ip.size();

//     Eigen::VectorXd dDist(nNodes);        
//     Eigen::MatrixXd ddDist(nNodes, nNodes);      

//     Eigen::RowVectorXd dPreFac_dL(nNodes), dPreFac_dPhi(nNodes);
//     Eigen::MatrixXd ddPrefac_ddPhi(nNodes, nNodes);
//     double preFac, prefix, distOld, distNew;

//     prefix = M->Gc/M->pf_l/LagFact;

//     distOld = PF_Util->GetDistributionFunction(dDist, ddDist, Nph_ip, Gph_ip, G2ph_ip, PhaseOld);
//     distNew = PF_Util->GetDistributionFunction(dDist, ddDist, Nph_ip, Gph_ip, G2ph_ip, Phase);

//     // lower limit, dist>distOld,   E = 1/(2LagFact) * <L + LagFact*(dist-distOld) >^2 - 1/(2LagFact)*L^2
//     preFac = Nph_ip*L + LagFact * (distNew - distOld);
//     dPreFac_dL = Nph_ip;
//     dPreFac_dPhi = LagFact * dDist.transpose();
//     ddPrefac_ddPhi = LagFact * ddDist.transpose();

//     F_L  += -w_ip * prefix * Nph_ip.transpose()*Nph_ip*L;
//     K_LL += -w_ip * prefix * Nph_ip.transpose()*Nph_ip;

//     F_Ph    += w_ip * prefix * dPreFac_dPhi.transpose() * Nph_ip*L;
//     K_PhL   += w_ip * prefix * dPreFac_dPhi.transpose() * dPreFac_dL;
//     K_PhPh  += w_ip * prefix * ddPrefac_ddPhi * (Nph_ip*L);

//     if (preFac<=0.0){
//         F_Ph    += w_ip * prefix * dPreFac_dPhi.transpose() * (preFac-Nph_ip*L);
//         K_PhPh  += w_ip * prefix * (dPreFac_dPhi.transpose()*dPreFac_dPhi + (preFac-Nph_ip*L)*ddPrefac_ddPhi);
//         //K_PhL   += w_ip * prefix * dPreFac_dPhi.transpose() * dPreFac_dL;

//         F_L     += w_ip * prefix * dPreFac_dL.transpose() * preFac;
//         K_LL    += w_ip * prefix * dPreFac_dL.transpose() * dPreFac_dL;
//         K_LPh   += w_ip * prefix * dPreFac_dL.transpose() * dPreFac_dPhi;
//     }

//     //upper limit, phi=1
//     preFac = Nph_ip*L2 + LagFact * (Nph_ip*Phase - 1.001);
//     dPreFac_dL = Nph_ip;
//     dPreFac_dPhi = LagFact * Nph_ip;
//     ddPrefac_ddPhi.setZero(); 

//     F_L2   += -w_ip * prefix * Nph_ip.transpose()*Nph_ip*L2;
//     K_L2L2 += -w_ip * prefix * Nph_ip.transpose()*Nph_ip;

//     F_Ph    += w_ip * prefix * dPreFac_dPhi.transpose() * Nph_ip*L2;
//     K_PhL2   += w_ip * prefix * dPreFac_dPhi.transpose() * dPreFac_dL;
//     K_PhPh  += w_ip * prefix * ddPrefac_ddPhi * (Nph_ip*L2);

//     if (preFac>=0.0){
//         F_Ph      += w_ip * prefix * dPreFac_dPhi.transpose() * (preFac-Nph_ip*L2);
//         K_PhPh    += w_ip * prefix * (dPreFac_dPhi.transpose()*dPreFac_dPhi + (preFac-Nph_ip*L2)*ddPrefac_ddPhi);
//         //K_PhL2    += w_ip * prefix * dPreFac_dPhi.transpose() * dPreFac_dL;

//         F_L2      += w_ip * prefix * dPreFac_dL.transpose() * preFac;
//         K_L2L2    += w_ip * prefix * dPreFac_dL.transpose() * dPreFac_dL;
//         K_L2Ph    += w_ip * prefix * dPreFac_dL.transpose() * dPreFac_dPhi;
//     }

// }