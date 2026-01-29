#include "PhaseFieldMaterial.h"

PhaseFieldMaterial::PhaseFieldMaterial(inputData &inputs, std::string MatName): SolidMaterial(inputs, MatName), BaseMaterial(inputs, MatName){

    inputs.GetRequired(pf_l, {"properties", MatName, "FractureProperties", "l"});
    inputs.GetRequired(Gc, {"properties", MatName, "FractureProperties", "Gc"});

    pf_visc = 0.0;
    inputs.GetOptional(pf_visc, {"properties", MatName, "FractureProperties", "Viscosity"});

    //stress splitting scheme
    std::string SplitName;
    inputs.GetRequired(SplitName, {"properties", MatName, "FractureProperties", "Split"});
    if (SplittingMethodNames.count(SplitName)){
        StressSplit = SplittingMethodNames[SplitName];
    } else {
        throw std::invalid_argument("Stress split function type "+SplitName+" not defined,");
    }

    ConsistentEnergy = true;
    inputs.GetOptional(ConsistentEnergy, {"properties", "PhaseField", "ConsistentEnergy"});

    if (StressSplit==DruckerPrager){
        inputs.GetRequired(B_dp, {"properties", MatName, "FractureProperties", "B_DP"});
    }
    if (StressSplit==ICE){
        inputs.GetRequired(Ice_eRef, {"properties", MatName, "FractureProperties", "eRef"});
        inputs.GetRequired(Ice_cReduce, {"properties", MatName, "FractureProperties", "cReduce"});
    }
    if (StressSplit==StarConvex){
        inputs.GetRequired(StarConvexGamma, {"properties", MatName, "FractureProperties","StarConvexGamma"});
    }
}

PhaseFieldMaterial::~PhaseFieldMaterial(){

}


/// @brief Performs energy-based split to go from elastic strains to energies, stresses, and material stiffness matrices
/// @param D_Dam output: Damageable material tangent matrix
/// @param D_UnDam output: Undamageable material tangent matrix
/// @param Stress_Dam   output: Damageable stresses
/// @param Stress_Undam output: Undamageable stresses
/// @param Strain input: linear-elastic strain (e_xx, e_yy, e_zz, 2*e_xy)
/// @return Energy, total damageable strain energy
double PhaseFieldMaterial::EnergySplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain){
    double Energy = 0;  //damageable strain energy

    //set components to zero (overwritten unless supposed to stay zero)
    Stress_Dam.setZero();
    Stress_Undam.setZero();
    D_Dam.setZero();
    D_UnDam.setZero();

    switch (StressSplit) {
        case NoSplit:{
                D_Dam = D;
                Stress_Dam = D*Strain;
                Energy = 0.5*Stress_Dam.transpose()*Strain;
        } break;
        case VolStrains:{   //Amor strain split
            Energy = VolStrainSplit(D_Dam, D_UnDam, Stress_Dam, Stress_Undam, Strain);
        } break;
        case SpecStrains:{
            Energy = SpecStrainSplit(D_Dam, D_UnDam, Stress_Dam, Stress_Undam, Strain);
        } break;
        case VolStress:{
            Energy = VolStressSplit(D_Dam, D_UnDam, Stress_Dam, Stress_Undam, Strain);
        } break;
        case SpecStress:{
            Energy = SpecStressSplit(D_Dam, D_UnDam, Stress_Dam, Stress_Undam, Strain);
        } break;
        case DruckerPrager:{
            Energy = DruckerPragerSplit(D_Dam, D_UnDam, Stress_Dam, Stress_Undam, Strain);
        } break;
        case ICE:{
            Energy = IceSplit(D_Dam, D_UnDam, Stress_Dam, Stress_Undam, Strain);
        } break;
        case StarConvex:{
            Energy = StarConvexSplit(D_Dam, D_UnDam, Stress_Dam, Stress_Undam, Strain);
        } break;
    }

    //Stresses not using the same split as energies
    if (ConsistentEnergy == false){
        Stress_Dam = D*Strain;
        Stress_Undam.setZero();
        D_Dam = D;
        D_UnDam.setZero();
    }

    return Energy;
}


double PhaseFieldMaterial::VolStressSplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain){
    double Energy = 0.0;
    Vector6d Stress = D*Strain;
    double p = 1.0/3.0*i.transpose()*Stress;

    //Volumetric stresses
    if (p>=0){
        Stress_Dam += i*p;
        D_Dam      += 1.0/3.0*i*i.transpose()*D;
    } else {
        Stress_Undam += i*p;
        D_UnDam      += 1.0/3.0*i*i.transpose()*D;
    }
    Stress_Dam += Stress - i*p;
    D_Dam += (Matrix6d::Identity() - 1.0/3.0*i*i.transpose())*D;

    //Deviatoric stress, always contributes
    Energy = 0.5*Stress_Dam.transpose()*Strain;

    return Energy;
}

double PhaseFieldMaterial::SpecStressSplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain){
    double Energy = 0.0;

    Vector6d Stress = D*Strain;

    if(true){
        Vector6d Stress_nonz = D*Strain;
        for (size_t i = 0; i < 4; i++){
            if (Stress_nonz(i)*Stress_nonz(i)<1.0e-9){
                Stress_nonz(i)=(i+1)*1e-9;
            }
        }

        Eigen::Matrix3d StressMat{{Stress_nonz(0), Stress_nonz(3), Stress_nonz(5)},{Stress_nonz(3), Stress_nonz(1), Stress_nonz(4)},{Stress_nonz(5), Stress_nonz(4), Stress_nonz(2)}};
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(StressMat);
        Eigen::Vector3d EigenValues = eigensolver.eigenvalues();
        Eigen::Matrix3d EigenVectors = eigensolver.eigenvectors();

        Eigen::Matrix3d EigStressMat_Dam, StressMat_Dam; 

        EigStressMat_Dam.setZero();
        for (size_t i = 0; i < 3; i++){
            if(EigenValues(i)>0){
                EigStressMat_Dam(i,i) = EigenValues(i);
            }
        }
        StressMat_Dam = EigenVectors*EigStressMat_Dam*EigenVectors.transpose();
        Vector6d StressVec{StressMat_Dam(0,0), StressMat_Dam(1,1), StressMat_Dam(2,2), 0.5*(StressMat_Dam(0,1)+StressMat_Dam(1,0)), 0.5*(StressMat_Dam(1,2)+StressMat_Dam(2,1)), 0.5*(StressMat_Dam(0,2)+StressMat_Dam(2,0))};

        Energy = 0.5*StressVec.transpose()*Strain;

        if (ConsistentEnergy == true){
            throw std::invalid_argument("Spectral Stress split not implemented when using a consistent split for energy/stress damage\n");
        }
    } else {
        // std::vector<double> EigenStress(3);
        // std::vector<Eigen::Vector4d> dEigenStress(3);
        // std::vector<Eigen::Matrix4d> d2EigenStress(3);

        // std::vector<Eigen::Vector3d> EigenVectors(3);

        // GetEigenValues(Stress, EigenStress, dEigenStress, d2EigenStress);
        // GetEigenVectors(Stress, EigenVectors);


        // for (size_t i = 0; i < 3; i++){
        //     //Eigen::Matrix3d EigSpace = EigenVectors[i]*EigenVectors[i].transpose();

        //     Eigen::Vector4d EigSpaceFlat; EigSpaceFlat.setZero(); 
        //     EigSpaceFlat(0) = EigenVectors[i](0)*EigenVectors[i](0);
        //     EigSpaceFlat(1) = EigenVectors[i](1)*EigenVectors[i](1);
        //     EigSpaceFlat(2) = EigenVectors[i](2)*EigenVectors[i](2);
        //     EigSpaceFlat(3) = EigenVectors[i](0)*EigenVectors[i](1);

        //     double ET_epsilon = EigSpaceFlat.transpose()*Strain;

        //     if (EigenStress[i]>=0){
        //         Energy += 0.5 * EigenStress[i] * ET_epsilon;

        //         Stress_Dam += 0.5 * EigenStress[i] * EigSpaceFlat + 0.5*ET_epsilon * D * dEigenStress[i];
        //         D_Dam      += 0.5 * D * dEigenStress[i] * EigSpaceFlat.transpose() 
        //                     + 0.5 * D * EigSpaceFlat    * dEigenStress[i].transpose()
        //                     + 0.5 * ET_epsilon * D * d2EigenStress[i];
        //     } else {
        //         Stress_Undam +=   0.5 * EigenStress[i] * EigSpaceFlat + 0.5*ET_epsilon * D * dEigenStress[i];
        //         D_UnDam      += 0.5 * D * dEigenStress[i] * EigSpaceFlat.transpose() 
        //                       + 0.5 * D * EigSpaceFlat    * dEigenStress[i].transpose()
        //                       + 0.5 * ET_epsilon * D * d2EigenStress[i];
        //     }
        //}
    }

    return Energy;
}

double PhaseFieldMaterial::VolStrainSplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain){
    double Energy = 0.0;
    double p = i.transpose()*Strain;    //volumetric strains
    
    //volumetric
    if (p>=0){  //if positive, volumetric strains can be damaged
        Energy += 0.5*Bulk*p*p;
        Stress_Dam += Bulk*p*i;
        D_Dam += Bulk*i*i.transpose();
    } else {    //compressive, can not be damaged/fractured
        Stress_Undam += Bulk*p*i;
        D_UnDam += Bulk*i*i.transpose();
    }

    //deviatoric, always can be damaged
    Vector6d DevStrain = VoightCorrect*Strain - 1.0/3.0*p*i;
    Energy += Shear*DevStrain.transpose()*VoightDoubleDot*DevStrain;
    Stress_Dam += 2*Shear*VoightCorrect*VoightDoubleDot*DevStrain;
    D_Dam += 2*Shear*VoightCorrect*VoightDoubleDot*VoightCorrect*(Matrix6d::Identity()-1.0/3.0*i*i.transpose());

    return Energy;
}

double PhaseFieldMaterial::StarConvexSplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain){
    double Energy = 0.0;
    double p = i.transpose()*Strain;    //volumetric strains
    Vector6d DevStrain = VoightCorrect*Strain - 1.0/3.0*p*i; //deviatoric strains
    
    double pplus, pmin, dpplus, dpmin;
    if (p>=0){
        pplus = p; dpplus = 1.0;
        pmin = 0.0; dpmin = 0.0;
    } else {
        pplus = 0.0; dpplus = 0.0;
        pmin = p; dpmin = 1.0;
    }

    //volumetric
    Energy += 0.5*Bulk*(pplus*pplus-StarConvexGamma*pmin*pmin);
    Stress_Dam += Bulk*(pplus*dpplus - StarConvexGamma*pmin*dpmin)*i.transpose();
    D_Dam +=  Bulk*(dpplus*dpplus - StarConvexGamma*dpmin*dpmin)*i*i.transpose();

    Stress_Undam += (1.0+StarConvexGamma)*Bulk*pmin*dpmin*i.transpose();
    D_UnDam += (1.0+StarConvexGamma)*Bulk*dpmin*dpmin*i*i.transpose();

    //deviatoric
    Energy += Shear*DevStrain.transpose()*VoightDoubleDot*DevStrain;
    Stress_Dam += 2*Shear*VoightCorrect*VoightDoubleDot*DevStrain;
    D_Dam += 2*Shear*VoightCorrect*VoightDoubleDot*VoightCorrect*(Matrix6d::Identity()-1.0/3.0*i*i.transpose());

    return Energy;
}

double PhaseFieldMaterial::SpecStrainSplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain){
    double Energy = 0.0;

    // NOT YET UPDATED TO WORK WITH GENERAL 3D STRAINS (needs updating to EigenValues)!!!!!!!!!!!!!!!!!!!!!

    std::vector<double> PriStrain(3); //principal strains
    std::vector<Vector6d> dPriStrain_dStrain(3);     //principal strains to elastic strains derivative d(e_1/e_2/e_3)/d(xx/yy/zz/xy)
    std::vector<Matrix6d> d2PriStrain_dStrain2(3);   //principal strains to elastic strains second order derivative d2(e_1/e_2/e_3)/d(xx/yy/zz/xy)^2
    Vector6d StrainVoight = VoightCorrect*Strain;
    GetEigenValues2D(StrainVoight, PriStrain, dPriStrain_dStrain, d2PriStrain_dStrain2);

    //Volumetric strains
    double          p0     = PriStrain[0]           +PriStrain[1]           +PriStrain[2];  //volumetric strain
    Vector6d pdiff1 = dPriStrain_dStrain[0]  +dPriStrain_dStrain[1]  +dPriStrain_dStrain[2]; //d(e_V)/d(xx/yy/zz/xy)
    Matrix6d pdiff2 = d2PriStrain_dStrain2[0]+d2PriStrain_dStrain2[1]+d2PriStrain_dStrain2[2];//d2(e_V)/d(xx/yy/zz/xy)2

    if (p0>=0){ //extensional, strains contribute to fracture
        Energy += 0.5*Lame*std::pow(p0,2);
        Stress_Dam += Lame*p0*pdiff1;
        D_Dam += Lame*(p0*pdiff2 + pdiff1*pdiff1.transpose());
    } else {    //compressive, not contributing to damage
        Stress_Undam += Lame*p0*pdiff1;
        D_UnDam += Lame*(p0*pdiff2 + pdiff1*pdiff1.transpose());
    }

    //Deviatoric strains, contributing if positive
    for (size_t ic = 0; ic < 3; ic++){
        Vector6d dPriStrain2  = 2*PriStrain[ic]*dPriStrain_dStrain[ic];
        Matrix6d ddPriStrain2 = 2*PriStrain[ic]*d2PriStrain_dStrain2[ic] + 2*dPriStrain_dStrain[ic]*dPriStrain_dStrain[ic].transpose();

        if (PriStrain[ic]>=0){
            Energy     += Shear  * std::pow(PriStrain[ic],2);
            Stress_Dam += Shear* dPriStrain2;
            D_Dam      += Shear* ddPriStrain2;
        } else {
            Stress_Undam += Shear* dPriStrain2;
            D_UnDam      += Shear* ddPriStrain2;
        }
    }

    Stress_Dam = VoightCorrect*Stress_Dam;
    Stress_Undam = VoightCorrect*Stress_Undam;
    D_Dam = VoightCorrect*D_Dam*VoightCorrect;
    D_UnDam = VoightCorrect*D_UnDam*VoightCorrect;

    //get stresses based on material tangent matrix (reduces rounding errors)
    Stress_Dam = D_Dam*Strain;
    Stress_Undam = D_UnDam*Strain;

    return Energy;
}

double PhaseFieldMaterial::DruckerPragerSplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain){
    double Energy = 0.0;

    double I1 = i.transpose()*Strain;
    Vector6d dI1_dStrain = i;
    double J2 = 0.5*Strain.transpose()*J2Mat_strain*Strain;
    Vector6d dJ2_dStrain = Strain.transpose()*J2Mat_strain;

    double sqrtJ2 = std::sqrt(J2+1.0e-16);
    Vector6d dsqrtJ2 = 0.5/sqrtJ2*dJ2_dStrain;
    Matrix6d ddsqrtJ2 = 0.5/sqrtJ2*J2Mat_strain - 0.25/sqrtJ2/sqrtJ2/sqrtJ2*dJ2_dStrain*dJ2_dStrain.transpose();

    if (-6.0*B_dp*std::sqrt(J2)<I1){
        Energy = 0.5*Bulk*std::pow(I1, 2) + 2.0*Shear*J2;

        Stress_Dam = Bulk*I1*dI1_dStrain + 2.0*Shear*dJ2_dStrain;
        D_Dam = Bulk*dI1_dStrain*dI1_dStrain.transpose() +2.0*Shear*J2Mat_strain;
        
        Stress_Undam.setZero();
        D_UnDam.setZero();
    } else if (2*Shear*std::sqrt(J2)<3*B_dp*Bulk*I1){
        Energy = 0.0;
        Stress_Dam.setZero();
        D_Dam.setZero();

        Stress_Undam = Bulk*I1*dI1_dStrain + 2.0*Shear*dJ2_dStrain;
        D_UnDam = Bulk*dI1_dStrain*dI1_dStrain.transpose() +2.0*Shear*J2Mat_strain;
    } else {
        double prefac = 1.0/(18.0*B_dp*B_dp*Bulk +2.0*Shear);
        double quadraticTerm = -3.0*B_dp*Bulk*I1 + 2*Shear*sqrtJ2;
        Vector6d dQuadraticterm  = -3.0*B_dp*Bulk*dI1_dStrain + 2.0*Shear*dsqrtJ2;
        Matrix6d ddQuadraticTerm = 2.0*Shear*ddsqrtJ2;

        Energy = prefac * std::pow(quadraticTerm, 2);
        Stress_Dam = prefac * 2.0*quadraticTerm * dQuadraticterm;
        D_Dam = prefac * (2.0*quadraticTerm * ddQuadraticTerm + 2.0*dQuadraticterm*dQuadraticterm.transpose());

        prefac = Bulk*Shear/(18.0*B_dp*B_dp*Bulk +2.0*Shear);
        quadraticTerm = I1 + 6*B_dp*sqrtJ2;
        dQuadraticterm  = dI1_dStrain + 6*B_dp*dsqrtJ2;
        ddQuadraticTerm = 6*B_dp*ddsqrtJ2;

        Stress_Undam = prefac * 2.0*quadraticTerm * dQuadraticterm;
        D_UnDam = prefac * (2.0*quadraticTerm * ddQuadraticTerm + 2.0*dQuadraticterm*dQuadraticterm.transpose());
    }

    //get stresses based on material tangent matrix (reduces rounding errors)
    //Stress_Dam = D_Dam*Strain;
    //Stress_Undam = D_UnDam*Strain;

    //std::cout << D << "\n\n";
    //std::cout << D_Dam << "\n\n";
    //std::cout << D_UnDam << "\n\n";

    return Energy;
}

double PhaseFieldMaterial::IceSplit(Matrix6d& D_Dam, Matrix6d& D_UnDam, Vector6d& Stress_Dam, Vector6d& Stress_Undam, Vector6d& Strain){
    double Energy = 0.0;
    double p = i.transpose()*Strain;    //volumetric strains
    Vector6d DevStrain = VoightCorrect*Strain - 1.0/3.0*p*i; //deviatoric strains
    
    double pplus, pmin, dpplus, dpmin;
    if (p>=0){
        pplus = p; dpplus = 1.0;
        pmin = 0.0; dpmin = 0.0;
    } else {
        pplus = 0.0; dpplus = 0.0;
        pmin = p; dpmin = 1.0;
    }

    //volumetric
    Energy += 0.5*Bulk*pplus*pplus;
    Stress_Dam += Bulk*pplus*dpplus*i.transpose();
    D_Dam +=  Bulk*dpplus*dpplus*i*i.transpose();

    Stress_Undam += Bulk*pmin*dpmin*i.transpose();
    D_UnDam += Bulk*dpmin*dpmin*i*i.transpose();

    //deviatoric
    Energy += Shear*DevStrain.transpose()*VoightDoubleDot*DevStrain;
    Stress_Dam += 2*Shear*VoightCorrect*VoightDoubleDot*DevStrain;
    D_Dam += 2*Shear*VoightCorrect*VoightDoubleDot*VoightCorrect*(Matrix6d::Identity()-1.0/3.0*i*i.transpose());

    // Cohesion
    Energy += -0.5*Bulk*Ice_cReduce*pmin*pmin / (Ice_eRef*Ice_eRef+pmin*pmin);

    // Eigen::Vector4d StressOffset = Bulk  *  Ice_cReduce * Ice_eRef*Ice_eRef*pmin*std::pow(Ice_eRef*Ice_eRef+pmin*pmin, -2) * dpmin*i.transpose();
    // Eigen::Matrix4d KOffset      = Bulk  *  Ice_cReduce * Ice_eRef*Ice_eRef*(Ice_eRef*Ice_eRef-3.0*pmin*pmin) * std::pow(Ice_eRef*Ice_eRef+pmin*pmin, -3) *dpmin*dpmin*i*i.transpose();

    // Stress_Dam += -StressOffset;    D_Dam += -KOffset;
    // Stress_Undam += StressOffset;   D_UnDam += KOffset;

    //get stresses based on material tangent matrix (reduces rounding errors)
    Stress_Dam = D_Dam*Strain;
    Stress_Undam = D_UnDam*Strain;

    if (Energy<0.0){
        Energy = 0.0;
    }
    if (Energy>1.0e12 || isnan(Energy)){
        Energy = 1.0e12;
    }
    if (Energy<-1.0e12){
        Energy = -1.0e12;
    }


    return Energy;
}

void PhaseFieldMaterial::GetEigenValues2D(Vector6d& s_inUnfiltered, std::vector<double>& s, std::vector<Vector6d>& ds, std::vector<Matrix6d>& dds){
    //protection against zero-mode eigenstrains
    Vector6d s_in;  s_in = 1.0*s_inUnfiltered; 

    //algorith only valid if sqrt(e_eq)>>sqrt(__DBL_EPSILON__), otherwise results are just rounding errors
    //double sqrt_test = s_in(0)*s_in(0)+s_in(1)*s_in(1)-2.0*s_in(0)*s_in(1)+4*s_in(3)*s_in(3);
    //if (sqrt_test<=0){ 
        for (size_t ic = 0; ic < 3; ic++){
            if (std::abs(s_in(ic))<1.0e-6){ 
                s_in(ic)=(ic+1.0)*1.0e-6;
            }
        }
    //} else if (sqrt_test<=100*1.0e-8){
    //    s_in = s_in * sqrt_test/(100*1.0e-8);
    //}

    for (size_t ic = 0; ic < 3; ic++){  //calculate principal strains (not ordered by magnitude) and derivatives, assuming plane strain conditions
        if (ic==0){ //1st principal strain: trivial, zz component of strains
            s[ic] = s_inUnfiltered(2); //zz
            ds[ic].setZero(); ds[ic](2) = 1.0; //[0 0 1 0]
            dds[ic].setZero(); //0
        } else {  //2nd and 3rd principal strains: composed of xx/yy/xy components of stresses
            double sgn=2.0*ic-3.0; //-1 or 1
            double sqrt_term = s_in(0)*s_in(0)+s_in(1)*s_in(1)-2.0*s_in(0)*s_in(1)+4*s_in(3)*s_in(3);
            double sqrt_term_unf = s_inUnfiltered(0)*s_inUnfiltered(0)+s_inUnfiltered(1)*s_inUnfiltered(1)-2.0*s_inUnfiltered(0)*s_inUnfiltered(1)+4*s_inUnfiltered(3)*s_inUnfiltered(3);
            s[ic] = 0.5*(s_inUnfiltered(0)+s_inUnfiltered(1)+sgn*std::pow(sqrt_term_unf,0.5)); //0.5*(xx+yy+-sqrt(xx^2+yy^2-2xx*yy+xy^2))
            ds[ic].setZero(); 
            ds[ic](0) = 0.5+0.5*sgn*( s_in(0)-s_in(1))*std::pow(sqrt_term+__DBL_EPSILON__,-0.5);
            ds[ic](1) = 0.5+0.5*sgn*(-s_in(0)+s_in(1))*std::pow(sqrt_term+__DBL_EPSILON__,-0.5);
            ds[ic](2) = 0.0;
            ds[ic](3) = 2.0*sgn*s_in(3)*std::pow(sqrt_term+__DBL_EPSILON__,-0.5);
            ds[ic](4) = 0;
            ds[ic](5) = 0;


            double c1 = 2.0*s_in(3)*(s_in(0)-s_in(1))          *std::pow(sqrt_term+__DBL_EPSILON__,-1.5);
            double c2 = 2.0*s_in(3)*s_in(3)                    *std::pow(sqrt_term+__DBL_EPSILON__,-1.5);
            double c3 = 2.0*(s_in(0)-s_in(1))*(s_in(0)-s_in(1))*std::pow(sqrt_term+__DBL_EPSILON__,-1.5);

            dds[ic].setZero();
            dds[ic] <<   sgn*c2,  -sgn*c2, 0, -sgn*c1, 0, 0,
                        -sgn*c2,  sgn*c2, 0,  sgn*c1, 0, 0,
                         0,        0,      0,  0, 0, 0,
                         -sgn*c1,  sgn*c1, 0,  sgn*c3, 0, 0,
                         0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0;
        } 
    }
}

void PhaseFieldMaterial::GetEigenVectors2D(Vector6d& s_xxyyzzxy, std::vector<Eigen::Vector3d>& EigVec){
    // https://www.wolframalpha.com/input?i=eigenvectors+of+%7B%7Ba%2Cb%2C0%7D%2C%7Bb%2Cd%2C0%7D%2C%7B0%2C0%2Cf%7D%7D

    Vector6d s_in;  s_in = 1.0*s_xxyyzzxy; 
    for (size_t ic = 0; ic < 4; ic++){ //eigenVector algortih does not work when any component is zero
        if (std::abs(s_in(ic))<1.0e-6){ 
            s_in(ic)=(ic+1.0)*1.0e-6;
        }
    }

    double sqrt_term = s_in(0)*s_in(0)+s_in(1)*s_in(1)-2.0*s_in(0)*s_in(1)+4*s_in(3)*s_in(3);

    EigVec[0].setZero(); EigVec[0](2) = 1.0; // [0  0  1]

    //std::cout << EigVec[0].transpose() << "\n";
    for (size_t ic = 1; ic < 3; ic++){   //[s_xx-s_yy +- sqrt((s_xx-s_yy)^2+4s_xy^2)   2s_xy    0], then normalized
        double sgn=2.0*ic-3.0; //-1 or 1
        EigVec[ic].setZero();
        EigVec[ic](0) = 0.5*(s_in(0)  -  s_in(1)  +  sgn*std::sqrt(sqrt_term));
        EigVec[ic](1) = s_in(3); 

        EigVec[ic].normalize();

        //std::cout << EigVec[ic].transpose() << "\n";
    }
    // Eigen::Matrix3d StressMat{{s_in(0), s_in(3), 0.0},{s_in(3), s_in(1), 0.0},{0.0, 0.0, s_in(2)}};
    // Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(StressMat);
    // Eigen::Vector3d EigenValues = eigensolver.eigenvalues();
    // Eigen::Matrix3d EigenVectors = eigensolver.eigenvectors();

    // std::cout << EigenVectors << "\n\n\n\n";
}