#include "CosseratMaterial.h"

CosseratMaterial::CosseratMaterial(inputData &inputs, std::string MatName): PhaseFieldMaterial(inputs, MatName), SolidMaterial(inputs, MatName), BaseMaterial(inputs, MatName){

    inputs.GetRequired(l_coss, {"properties", MatName, "SolidProperties", "l_coss"});
    inputs.GetRequired(G_coss, {"properties", MatName, "SolidProperties", "G_coss"});

    Hardening_e_ref = 0.0;
    inputs.GetOptional(Hardening_e_ref, {"properties", MatName, "SolidProperties", "hard_eRef"});

    ic.setZero(); icl.setZero();
    ic[0] = 1.0; ic[1] = 1.0; ic[2] = 1.0;
    icl[0] = 1.0; icl[1] = 1.0; icl[2] = 1.0;

    SkewMat.setZero();
    for (size_t i = 0; i < 3; i++){
        SkewMat(2*i+3,2*i+3) = 1.0;
        SkewMat(2*i+3,2*i+4) = -1.0;
        SkewMat(2*i+4,2*i+3) = -1.0;
        SkewMat(2*i+4,2*i+4) = 1.0;
    }

    MacroStrains.setZero();
    MicroStrains.setZero();
    for (size_t i = 0; i < 9; i++){
        MacroStrains(i,i) = 1.0;
        MicroStrains(i,i+9) = 1.0;
    }

    double a1=0.25, a2=0.25, a3=0.5;

    J2Mat_Coss.setZero();
    J2Mat_Coss(0,0) = 2.0/3.0;
    J2Mat_Coss(0,1) = -1.0/3.0;
    J2Mat_Coss(0,2) = -1.0/3.0;
    J2Mat_Coss(1,0) = -1.0/3.0;
    J2Mat_Coss(1,1) = 2.0/3.0;
    J2Mat_Coss(1,2) = -1.0/3.0;
    J2Mat_Coss(2,0) = -1.0/3.0;
    J2Mat_Coss(2,1) = -1.0/3.0;
    J2Mat_Coss(2,2) = 2.0/3.0;

    J2Mat_Coss(3,3) = 2*a1;
    J2Mat_Coss(3,6) = 2*a2;
    J2Mat_Coss(6,3) = 2*a2;
    J2Mat_Coss(6,6) = 2*a1;

    J2Mat_Coss(4,4) = 2*a1;
    J2Mat_Coss(4,7) = 2*a2;
    J2Mat_Coss(7,4) = 2*a2;
    J2Mat_Coss(7,7) = 2*a1;

    J2Mat_Coss(5,5) = 2*a1;
    J2Mat_Coss(5,8) = 2*a2;
    J2Mat_Coss(8,5) = 2*a2;
    J2Mat_Coss(8,8) = 2*a1;

    J2Mat_Coss(9,9) = 2*a3;
    J2Mat_Coss(10,10) = 2*a3;
    J2Mat_Coss(11,11) = 2*a3;
    J2Mat_Coss(12,12) = 2*a3;
    J2Mat_Coss(13,13) = 2*a3;
    J2Mat_Coss(14,14) = 2*a3;
    J2Mat_Coss(15,15) = 2*a3;
    J2Mat_Coss(16,16) = 2*a3;
    J2Mat_Coss(17,17) = 2*a3;

    D_coss.setZero();
    D_coss(0,0) = Lame + 2.0*Shear;
    D_coss(0,1) = Lame;
    D_coss(0,2) = Lame;
    D_coss(1,0) = Lame;
    D_coss(1,1) = Lame + 2.0*Shear;
    D_coss(1,2) = Lame;
    D_coss(2,0) = Lame;
    D_coss(2,1) = Lame;
    D_coss(2,2) = Lame + 2.0*Shear;

    D_coss(3,3) = Shear + G_coss;
    D_coss(4,3) = Shear - G_coss;
    D_coss(3,4) = Shear - G_coss;
    D_coss(4,4) = Shear + G_coss;

    D_coss(5,5) = Shear + G_coss;
    D_coss(6,5) = Shear - G_coss;
    D_coss(5,6) = Shear - G_coss;
    D_coss(6,6) = Shear + G_coss;

    D_coss(7,7) = Shear + G_coss;
    D_coss(8,7) = Shear - G_coss;
    D_coss(7,8) = Shear - G_coss;
    D_coss(8,8) = Shear + G_coss;

    D_coss(9,9) = 2.0*Shear;
    D_coss(10,10) = 2.0*Shear;
    D_coss(11,11) = 2.0*Shear;
    D_coss(12,12) = 2.0*Shear;
    D_coss(13,13) = 2.0*Shear;
    D_coss(14,14) = 2.0*Shear;
    D_coss(15,15) = 2.0*Shear;
    D_coss(16,16) = 2.0*Shear;
    D_coss(17,17) = 2.0*Shear;
}

CosseratMaterial::~CosseratMaterial(){

}

double CosseratMaterial::EnergySplit(Matrix18d& D_Dam, Matrix18d& D_UnDam, Vector18d& Stress_Dam, Vector18d& Stress_Undam, Vector18d& Strain){
    double Energy = 0;  //damageable strain energy

    //set components to zero (overwritten unless supposed to stay zero)
    Stress_Dam.setZero();
    Stress_Undam.setZero();
    D_Dam.setZero();
    D_UnDam.setZero();

    switch (StressSplit) {
        case NoSplit:{
                D_Dam = D_coss;
                Stress_Dam = D_coss*Strain;
                Energy = 0.5*Stress_Dam.transpose()*Strain;
        } break;
        case SpecStress:{
            if (ConsistentEnergy == true){
                Energy = SpecStressSplitConsistent(D_Dam, D_UnDam, Stress_Dam, Stress_Undam, Strain);
            } else {
                Energy = SpecStressSplit(D_Dam, D_UnDam, Stress_Dam, Stress_Undam, Strain);
            }
        } break;
        case DruckerPrager:{
            Energy = DruckerPragerSplit(D_Dam, D_UnDam, Stress_Dam, Stress_Undam, Strain);
        } break;
        case StarConvex:{
            Energy = StarConvexSplit(D_Dam, D_UnDam, Stress_Dam, Stress_Undam, Strain);
            break;
        }
        default:{
            throw std::invalid_argument("Stress split function type not defined for a Cosserat Material");
        } break;
    }

    //Stresses not using the same split as energies
    if (ConsistentEnergy == false){
        Stress_Dam = D_coss*Strain;
        Stress_Undam.setZero();
        D_Dam = 1.0*D_coss;
        D_UnDam.setZero();
    }

    return Energy;
}


double CosseratMaterial::StarConvexSplit(Matrix18d& D_Dam, Matrix18d& D_UnDam, Vector18d& Stress_Dam, Vector18d& Stress_Undam, Vector18d& Strain){
    double Energy = 0.0;
    double p = ic.transpose()*Strain;    //volumetric strains
    Vector9d DevStrain = MacroStrains*Strain - 1.0/3.0*p*icl; //deviatoric strains
    Matrix9d dDevStrain = Matrix9d::Identity() - 1.0/3.0*icl*icl.transpose();

    Vector9d CossStrain = MicroStrains*Strain;
    
    double pplus, pmin, dpplus, dpmin;
    if (p>=0){
        pplus = p; dpplus = 1.0;
        pmin = 0.0; dpmin = 0.0;
    } else {
        pplus = 0.0; dpplus = 0.0;
        pmin = p; dpmin = 1.0;
    }

    // //standard parts: volumetric
    Energy += 0.5*Bulk*pplus*pplus;
    Stress_Dam += Bulk*pplus*dpplus*ic;
    D_Dam +=  Bulk*dpplus*dpplus*ic*ic.transpose();

    Stress_Undam += Bulk*pmin*dpmin*ic;
    D_UnDam += Bulk*dpmin*dpmin*ic*ic.transpose();

    // //standard part: deviatoric
    Energy += Shear*DevStrain.transpose()*DevStrain;
    Stress_Dam += 2*Shear*MacroStrains.transpose()*dDevStrain*DevStrain;
    D_Dam += 2*Shear*MacroStrains.transpose()*dDevStrain*dDevStrain*MacroStrains;

    Energy += 0.5*(Shear-G_coss)*DevStrain.transpose()*SkewMat*DevStrain;
    Stress_Dam += (Shear-G_coss)*MacroStrains.transpose()*dDevStrain*SkewMat*DevStrain;
    D_Dam += (Shear-G_coss)*MacroStrains.transpose()*dDevStrain*SkewMat*dDevStrain*MacroStrains;

    // //Cosserat part: deviatoric
    //Energy += G_coss*CossStrain.transpose()*CossStrain;
    Stress_Dam += 2*G_coss*MicroStrains.transpose()*CossStrain;
    D_Dam += 2*G_coss*MicroStrains.transpose()*MicroStrains;

    // //StarConvexCorrection
    Energy += -0.5*Bulk*StarConvexGamma*pmin*pmin;
    //Stress_Dam += -Bulk*StarConvexGamma*pmin*dpmin*ic;
    //D_Dam += -Bulk*StarConvexGamma*dpmin*dpmin*ic*ic.transpose();

    //Stress_Undam += StarConvexGamma*Bulk*pmin*dpmin*ic;
    //D_UnDam += StarConvexGamma*Bulk*dpmin*dpmin*ic*ic.transpose();

    Energy = std::max(0.0, Energy);
    return Energy;
}

double CosseratMaterial::DruckerPragerSplit(Matrix18d& D_Dam, Matrix18d& D_UnDam, Vector18d& Stress_Dam, Vector18d& Stress_Undam, Vector18d& Strain){
    double Energy = 0.0;

    double I1 = ic.transpose()*Strain;
    Vector18d dI1_dStrain = ic;
    double J2 = 0.5*Strain.transpose()*J2Mat_Coss*Strain;
    Vector18d dJ2_dStrain = Strain.transpose()*J2Mat_Coss;

    double sqrtJ2 = std::sqrt(J2+1.0e-16);
    Vector18d dsqrtJ2 = 0.5/sqrtJ2*dJ2_dStrain;
    Matrix18d ddsqrtJ2 = 0.5/sqrtJ2*J2Mat_Coss - 0.25/sqrtJ2/sqrtJ2/sqrtJ2*dJ2_dStrain*dJ2_dStrain.transpose();

    if (-6.0*B_dp*std::sqrt(J2)<I1){
        Energy = 0.5*Bulk*std::pow(I1, 2) + 2.0*Shear*J2;

        Stress_Dam = Bulk*I1*dI1_dStrain + 2.0*Shear*dJ2_dStrain;
        D_Dam = Bulk*dI1_dStrain*dI1_dStrain.transpose() +2.0*Shear*J2Mat_Coss;
        
        Stress_Undam.setZero();
        D_UnDam.setZero();
    } else if (2*Shear*std::sqrt(J2)<3*B_dp*Bulk*I1){
        Energy = 0.0;
        Stress_Dam.setZero();
        D_Dam.setZero();

        Stress_Undam = Bulk*I1*dI1_dStrain + 2.0*Shear*dJ2_dStrain;
        D_UnDam = Bulk*dI1_dStrain*dI1_dStrain.transpose() +2.0*Shear*J2Mat_Coss;
    } else {
        double prefac = 1.0/(18.0*B_dp*B_dp*Bulk +2.0*Shear);
        double quadraticTerm = -3.0*B_dp*Bulk*I1 + 2*Shear*sqrtJ2;
        Vector18d dQuadraticterm  = -3.0*B_dp*Bulk*dI1_dStrain + 2.0*Shear*dsqrtJ2;
        Matrix18d ddQuadraticTerm = 2.0*Shear*ddsqrtJ2;

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

double CosseratMaterial::SpecStressSplitConsistent(Matrix18d& D_Dam, Matrix18d& D_UnDam, Vector18d& Stress_Dam, Vector18d& Stress_Undam, Vector18d& Strain){
    double Energy = 0.0;
    
    // 1) Build strain tensor
    Vector18d Stress = D_coss*Strain;
    Eigen::VectorXd strain_lin(9), strain_rot(9), stress_lin(9), stress_rot(9);
    Eigen::VectorXd Stress_nonz(9);
    for (size_t i = 0; i < 9; i++){
        strain_lin[i] = Strain(i);
        strain_rot[i] = Strain(i+9);
        stress_lin[i] = Stress(i);
        stress_rot[i] = Stress(i+9);

        Stress_nonz[i] = Stress(i);
        if (std::abs(Stress_nonz[i])<1.0e-12){
            Stress_nonz[i] = (i+1)*1.0e-12;
        }
    }

    // 3) Eigendecompose sigma
    Eigen::Matrix3d StressMat{{Stress_nonz(0), Stress_nonz(3), Stress_nonz(7)},
                                {Stress_nonz(4), Stress_nonz(1), Stress_nonz(5)},
                                {Stress_nonz(8), Stress_nonz(6), Stress_nonz(2)}};

    Eigen::Matrix3d StressMatSym = 0.5 * (StressMat + StressMat.transpose()); // Ensure symmetry

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(StressMatSym);
    Eigen::Vector3d sigma_eigs = eigensolver.eigenvalues();
    Eigen::Matrix3d EigenVectors = eigensolver.eigenvectors();

    // Eigen::JacobiSVD<Eigen::Matrix3d> svd(StressMat,Eigen::ComputeFullU|Eigen::ComputeFullV);
    // Eigen::Matrix3d U = svd.matrixU();
    // Eigen::Matrix3d V = svd.matrixV();
    // Eigen::Vector3d sigma_eigs = svd.singularValues();

    // 4) Split eigenvalues
    Eigen::Vector3d sigma_eigs_plus = sigma_eigs.cwiseMax(0.0);
    Eigen::Vector3d sigma_eigs_minus = sigma_eigs.cwiseMin(0.0);

    // 5) Reconstruct sigma_plus and sigma_minus
    Eigen::Matrix3d sigma_plus = EigenVectors * sigma_eigs_plus.asDiagonal() * EigenVectors.inverse();

    //assign as non-symmetric voight notation
    Vector9d StressVoight{sigma_plus(0,0),sigma_plus(1,1),sigma_plus(2,2),
                        sigma_plus(0,1),sigma_plus(1,0),
                        sigma_plus(1,2),sigma_plus(2,1),
                        sigma_plus(0,2),sigma_plus(2,0)};

    Energy = 0.5*StressVoight.transpose()*strain_lin;

    // 6) Build projectors and consistent tangent
    Matrix18d P_plus = Matrix18d::Zero();
    for (int i = 0; i < 3; ++i) {
        if (sigma_eigs(i) > 0.0) {
            Eigen::Vector3d u = EigenVectors.col(i);
            Eigen::Vector3d v = EigenVectors.col(i);
            // Outer product in Voigt notation
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    for (int c = 0; c < 3; ++c)
                        for (int d = 0; d < 3; ++d) {
                            int row = 3*a + b;
                            int col = 3*c + d;
                            P_plus(row, col) += u(a) * v(b) * u(c) * v(d);
                        }
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i == j) continue;
            // Only include if at least one singular value is positive
            if (sigma_eigs(i) <= 0.0 && sigma_eigs(j) <= 0.0) continue;
            double denom = sigma_eigs(i) - sigma_eigs(j);
            if (std::abs(denom) < 1e-12) continue; // avoid division by zero
            double num = (std::max(sigma_eigs(i), 0.0) - std::max(sigma_eigs(j), 0.0));
            double coeff = num / denom;

            Eigen::Vector3d ui = EigenVectors.col(i), uj = EigenVectors.col(j);
            Eigen::Vector3d vi = EigenVectors.col(i), vj = EigenVectors.col(j);

            // Q_ij = 0.5 * [ (ui ⊗ vj) ⊗ (uj ⊗ vi) + (uj ⊗ vi) ⊗ (ui ⊗ vj) ]
            for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
            for (int c = 0; c < 3; ++c)
            for (int d = 0; d < 3; ++d) {
                int row = 3*a + b;
                int col = 3*c + d;
                double term1 = ui(a) * vj(b) * uj(c) * vi(d);
                double term2 = uj(a) * vi(b) * ui(c) * vj(d);
                P_plus(row, col) += 0.5 * coeff * (term1 + term2);
            }
        }
    }
    
    for (int i = 9; i < 18; ++i) {
        P_plus(i, i) = 1.0; // Rotation components
    }

    D_Dam = P_plus * D_coss;
    D_UnDam = (Matrix18d::Identity() -  P_plus) * D_coss;

    // outputs
    Stress_Dam = D_Dam * Strain;
    Stress_Undam = D_UnDam * Strain;

    return Energy;
}

double CosseratMaterial::SpecStressSplit(Matrix18d& D_Dam, Matrix18d& D_UnDam, Vector18d& Stress_Dam, Vector18d& Stress_Undam, Vector18d& Strain){
    double Energy = 0.0;

    Vector18d Stress = D_coss*Strain;

    Eigen::VectorXd strain_lin(9), strain_rot(9), stress_lin(9), stress_rot(9);
    Eigen::VectorXd Stress_nonz(9);
    for (size_t i = 0; i < 9; i++){
        strain_lin[i] = Strain(i);
        strain_rot[i] = Strain(i+9);
        stress_lin[i] = Stress(i);
        stress_rot[i] = Stress(i+9);

        Stress_nonz[i] = Stress(i);
        if (std::abs(Stress_nonz[i])<1.0e-9){
            Stress_nonz[i] = (i+1)*1.0e-9;
        }
    }

    Eigen::Matrix3d StressMat{{Stress_nonz(0), Stress_nonz(3), Stress_nonz(7)},
                                {Stress_nonz(4), Stress_nonz(1), Stress_nonz(5)},
                                {Stress_nonz(8), Stress_nonz(6), Stress_nonz(2)}};

    // Ensure symmetry
    StressMat = 0.5 * (StressMat + StressMat.transpose());

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
    StressMat_Dam = EigenVectors*EigStressMat_Dam*EigenVectors.inverse();
    Eigen::VectorXd StressVec(9); 
    StressVec[0]=StressMat_Dam(0,0);StressVec[1]=StressMat_Dam(1,1);StressVec[2]=StressMat_Dam(2,2);
    StressVec[3]=StressMat_Dam(0,1);StressVec[4]=StressMat_Dam(1,0), // xy yx
    StressVec[5]=StressMat_Dam(1,2);StressVec[6]=StressMat_Dam(2,1), // yz zy
    StressVec[7]=StressMat_Dam(0,2);StressVec[8]=StressMat_Dam(2,0); // xz zx

    Energy = 0.5*StressVec.transpose()*strain_lin;
    Energy += 0*0.5*stress_rot.transpose()*strain_rot;


    Energy = std::max(0.0, Energy);
    if (ConsistentEnergy == true){
        // Estimate based on compressive-tension split to remove self-penetration
        Stress_Dam = D_coss*Strain;
        D_Dam = 1.0*D_coss;

        Stress_Undam.setZero();
        D_UnDam.setZero();

        double p = ic.transpose()*Strain;
        if (p<0){
            Stress_Undam = 0.25*Lame*p*ic;
            D_UnDam = 0.25*Lame*ic*ic.transpose();

            Stress_Dam -= Stress_Undam;
            D_Dam -= D_UnDam;
        }

        //throw std::invalid_argument("Spectral Stress split not implemented when using a consistent split for energy/stress damage\n");
    }

    return Energy;
}

/// @brief updates the plastic strains, based on total stresses and history
/// @param DMat output: Tangent matrix
/// @param Strain_Pl output/input: Current plastic strain (output, and used as initial guess)
/// @param Strain input: total strains
/// @param Strain_Pl_Old input history parameter: plastic strains at the end of last time increment
/// @param StrainOld input history: total strains at the end of the last time increment
/// @param phi input: current damage
void CosseratMaterial::UpdatePlasticStrains(Matrix18d& DMat, Vector18d& Strain_Pl, Vector18d& Strain, Vector18d& Strain_Pl_Old, Vector18d& StrainOld,
                         double phi, double phiOld, double dt, std::vector<double>& hist, std::vector<double>& histOld,
                         double& DissRate, double& DamDissRate){
    DissRate = 0.0;
    DamDissRate = 0.0;
    double damage = 1.0;
    double Unused1, Unused2;

    switch (Rheology){
    case SolidMaterialModel::LE:
            //no plastic strains, no need to update
            DMat = Matrix18d::Identity();
            Strain_Pl = 1.0*Strain_Pl_Old;
            DissRate = 0.0;
        break;
    case SolidMaterial::VE:
        throw std::invalid_argument("Rheology type not defined in CosseratMaterial.cpp,");
        break;
    case SolidMaterial::dVE:
        throw std::invalid_argument("Rheology type not defined in CosseratMaterial.cpp,");
        break;
    case SolidMaterial::VE_expl:
        throw std::invalid_argument("Rheology type not defined in CosseratMaterial.cpp,");
        break;
    case SolidMaterial::dVE_expl:
        throw std::invalid_argument("Rheology type not defined in CosseratMaterial.cpp,");
        break;
    case SolidMaterial::VEP:
        if (phi>PhiViscLim){
            //no additional plastic strains, no need to update
            DMat = Matrix18d::Identity();
            Strain_Pl = 1.0*Strain_Pl_Old;
        } else {
            //phi = std::min(1.0, std::max(0.0, phi));
            //PF_Util->GenericDamageFunction(phi, DegFunction, damage, Unused1, Unused2);
            damage = 1.0;
            return_mapping_VEP(DMat, Strain_Pl, Strain, Strain_Pl_Old, dt, damage, hist, histOld, DissRate, DamDissRate);
        }
        break;
    case SolidMaterial::dVEP:
        if (phi>PhiViscLim){
            //no additional plastic strains, no need to update
            DMat = Matrix18d::Identity();
            Strain_Pl = 1.0*Strain_Pl_Old;
        } else {
            phi = std::min(1.0, std::max(0.0, phi));
            PF_Util->GenericDamageFunction(phi, DegFunction, damage, Unused1, Unused2);
            //damage = 1.0;
            return_mapping_VEP(DMat, Strain_Pl, Strain, Strain_Pl_Old, dt, damage, hist, histOld, DissRate, DamDissRate);
        }
        break;
    default:
        throw std::invalid_argument("Rheology type not defined in CosseratMaterial.cpp,");
        break;
    }
}

void CosseratMaterial::return_mapping_VEP(Matrix18d& DMat, Vector18d& Strain_Pl, Vector18d& Strain, Vector18d& Strain_Pl_Old, double dt, double Damage, 
                            std::vector<double>& hist, std::vector<double>& histOld, double& DissRate, double& DamDissRate){
    Matrix20d K, Kinv;
    Vector20d f, fOld;
    Vector20d sol, dsol;
    double e, e0, lineSearch, ls1, ls2;
    uint it;

    sol(Eigen::seq(0,17)) = Strain-Strain_Pl; //elastic strain
    sol(18) = hist[1]; // dLambda_visc
    sol(19) = hist[2]; // dLambda_p

    return_mapping_VEP_getKF(K, f, sol, Strain, Strain_Pl_Old, dt, Damage, histOld);
    e0 = f.transpose()*f;
    e0 = std::max(e0,1e-12);

    bool conv = false;
    bool MappingBackup = false;
    it = 0;
    while (conv == false){
        it += 1;
        if (false){
            dsol = K.lu().solve(f); 
            assert(f.isApprox(K*dsol)); 
            dsol *= -1.0;
        } else {
            dsol = -K.inverse()*f;
        }
        sol += dsol;

        fOld = f;
        return_mapping_VEP_getKF(K, f, sol, Strain, Strain_Pl_Old, dt, Damage, histOld);
        ls1 = f.transpose()*dsol;
        ls2 = fOld.transpose()*dsol;
        if ((ls1-ls2)==0.0){
            lineSearch = 1.0;
        } else {
            lineSearch = -(ls2)/(ls1-ls2);
            lineSearch = std::min(std::max(lineSearch, 0.01), 1.0);
        }
        lineSearch = lineSearch;
        sol = sol-(1.0-lineSearch)*dsol;

        return_mapping_VEP_getKF(K, f, sol, Strain, Strain_Pl_Old, dt, Damage, histOld);

        e = f.transpose()*f;
        if (e<1.0e-6 || e/e0<1.0e-5){
            conv = true;
        }
        if (it>500){
            MappingBackup = true;
            //std::cout << "I";
            conv = true;
        }
        if (isnan(e) || isnan(e0) || (std::abs(e/e0)>1e6 && it>25)){
            MappingBackup = true;
        }
        if (MappingBackup){
            //return_mapping_VEP_BackUp(K, sol, Strain, Strain_Pl_Old, dt, Damage, histOld);
            conv = true;
        }
    }

    return_mapping_VEP_updateStrain(Strain, Strain_Pl, Strain_Pl_Old, sol, dt, Damage, hist, histOld, DissRate, DamDissRate);

    Kinv = K.inverse();
    DMat = scale1*Kinv(Eigen::seq(0,17),Eigen::seq(0,17));
}

void CosseratMaterial::return_mapping_VEP_BackUp(Matrix20d& K, Vector20d& sol, Vector18d& Strain, Vector18d& Strain_Pl_Old, 
                                                    double dt, double Damage, std::vector<double>& histOld){
    Vector20d f, fOld;
    Vector20d dsol;
    double e, e0, lineSearch, ls1, ls2;
    uint it;

    int ispread = 2000;//std::ceil(Strain.norm()/1.0e-4); ispread = std::max(1000, ispread);

    Vector18d StrainCopy = 1.0*Strain;

    //Logs.PrintEvery("ViscoPlastic return mapping failed, using backup method\n",2);

    sol(Eigen::seq(0,17)) = StrainCopy-Strain_Pl_Old; //elastic strain
    sol(18) = 0.0; // dLambda_visc
    sol(19) = 0.0; // dLambda_p

    bool conv = false;
    it = 0;

    return_mapping_VEP_getKF(K, f, sol, StrainCopy, Strain_Pl_Old, dt, Damage, histOld);
    e0 = f.transpose()*f;
    e0 = std::max(e0,1e-12);
    while (conv == false){
        //very slow loading
        if (it<ispread){
            StrainCopy = Strain_Pl_Old + (Strain-Strain_Pl_Old)*it/ispread;
            e0 = f.transpose()*f;
            e0 = std::max(e0,1e-12);
        } else {
            StrainCopy = 1.0*Strain;
        }

        it += 1;
        if (false){
            dsol = K.lu().solve(f); 
            assert(f.isApprox(K*dsol)); 
            dsol *= -1.0;
        } else {
            dsol = -K.inverse()*f;
        }
        sol += dsol;

        if (false){
            fOld = f;
            return_mapping_VEP_getKF(K, f, sol, StrainCopy, Strain_Pl_Old, dt, Damage, histOld);
            ls1 = f.transpose()*dsol;
            ls2 = fOld.transpose()*dsol;
            if ((ls1-ls2)==0.0){
                lineSearch = 1.0;
            } else {
                lineSearch = -(ls2)/(ls1-ls2);
                lineSearch = std::min(std::max(lineSearch, 0.1), 1.0);
            }
            //std::cout << it << ":  " << ls2 << "->" << ls1 << ", LineSearch: " << lineSearch << "\n";
            lineSearch = lineSearch;
            sol = sol-(1.0-lineSearch)*dsol;
        }


        return_mapping_VEP_getKF(K, f, sol, StrainCopy, Strain_Pl_Old, dt, Damage, histOld);

        e = f.transpose()*f;
        if ((e<1.0e-6 || e/e0<1.0e-5) && it>1.1*ispread){
            conv = true;
        }
        if (it>2*ispread){
            std::string msg = "ViscoPlastic BACKUP return mapping did not provide a solution within "+std::to_string(2*ispread)
            +" iterations: "+std::to_string(e/e0)+"\n";
            Logs.PrintEvery(msg,2);
            conv = true;
        }
        if (isnan(e) || isnan(e0) || (std::abs(e/e0)>1e6 && it>ispread)){
            std::stringstream msg;
            msg << "ViscoPlastic return mapping BACKUP reached a NaN: " << e << "/" << e0 <<"\n";

            msg << "Strain: "<< Strain.transpose()<<"\n";
            msg << "Strain_Pl_Old: "<<Strain_Pl_Old.transpose()<<"\n";
            msg << "HistOld: "<<histOld[0]<<" "<<histOld[1]<<" "<<histOld[2]<<"\n";
            msg << "Sol: "<<sol.transpose()<<"\n";
            msg << "f: "<<f.transpose()<<"\n";
            msg << "Damage: "<<Damage<<"\n";
            msg << "sy:" << sy << "\n\n";
            msg << "K" << K << "\n\n";
            msg << "dSol" << dsol.transpose() << "\n\n\n";

            Logs.PrintEvery(msg.str(),2);

            //print error and stop
            throw std::invalid_argument("ViscoPlastic return mapping reached a NaN");
        }
    }
}

void CosseratMaterial::return_mapping_VEP_getKF(Matrix20d& K, Vector20d& f,
                                Vector20d& sol, Vector18d& strain_total,
                                Vector18d& strain_vep_old, double dt, double Damage, std::vector<double>& histOld){
    f.setZero();
    K.setZero();

    // sol 0-17: elastic strains
    // sol 18:  plastic strain multiplier, viscoelasticity
    // sol 19:  plastic strain multiplier, plasticity

    Vector18d Stress = Damage*D_coss*sol(Eigen::seq(0,17));                    

    double F, FTrial, F_VE;
    Vector18d dF, dG, dF_VE;
    Matrix18d ddF, ddG, ddF_VE;

    F_VE = YieldFunctionVE(Stress, dF_VE, ddF_VE);

    Vector18d StressTrial = Damage*D_coss*(strain_total-strain_vep_old-sol(18)*dF_VE);
    FTrial = YieldFunction(StressTrial, StressTrial, dF, ddF);
    F = YieldFunction(Stress, StressTrial, dF, ddF);
    PotentialFunction(Stress, StressTrial, dG, ddG);

    // ddF.setZero();
    // ddF_VE.setZero();
    // ddG.setZero();

    bool ExplicitHardening = true;
    double include_dh = 1- ExplicitHardening;

    double e_eq = histOld[0] + include_dh*std::sqrt(sol(19)*sol(19));
    double h = hardening*sy*(1.0-std::exp(-e_eq*l_coss/Hardening_e_ref));
    double hOld = hardening*sy*(1.0-std::exp(-histOld[0]*l_coss/Hardening_e_ref));
    double dh = include_dh*hardening*sy/Hardening_e_ref*std::exp(-e_eq*l_coss/Hardening_e_ref)*sgn(sol(19));

    // 0 = e_el - strain + strain_old + dL_v s_dev + dLp df/ds
    f(Eigen::seq(0,17)) = scale1*(sol(Eigen::seq(0,17)) - strain_total + strain_vep_old + sol(18)*dF_VE + sol(19)*dG);
    K(Eigen::seq(0,17),Eigen::seq(0,17)) = scale1*(Matrix18d::Identity() + Damage*(sol(18)*ddF_VE + sol(19)*ddG)*D_coss); 
    K(Eigen::seq(0,17),18) = scale1*dF_VE;
    K(Eigen::seq(0,17),19) = scale1*dG;

    // 0 = F_VE - (A*L)^(1/n)
    double PreFac = ACreep*dt*std::pow(sc, nCreep);
    f(18) = scale2* (  PreFac*std::pow(F_VE/sc, nCreep) - sol(18) );
    K(18,18) = - scale2;
    K(18,Eigen::seq(0,17)) = scale2 * nCreep*PreFac/sc*std::pow(F_VE/sc, nCreep-1.0)*( Damage*dF_VE.transpose()*D_coss );

    // 0 = F;
    if (FTrial-hOld<0.0){
        f(19) = -scale3*(l_coss*pl_visc+1.0e-6)/dt*sol(19);
        K(19,19) = -scale3*(l_coss*pl_visc+1.0e-6)/dt;
    } else {
        f(19) = scale3*(  F - h - l_coss*pl_visc/dt*sol(19)  );
        K(19, Eigen::seq(0,17)) = scale3*( Damage*dF.transpose()*D_coss );
        K(19,19) = scale3*( -l_coss*pl_visc/dt - dh );
    }
    //std::cout << f << "\n" << K << "\n\n\n";
}

void CosseratMaterial::return_mapping_VEP_updateStrain(Vector18d& strain_total, Vector18d& Strain_New, Vector18d& strain_vep_old, 
                                        Vector20d& sol, double dt, double Damage, 
                                        std::vector<double>& hist, std::vector<double>& histOld, 
                                        double& DissRate, double& DamDissRate){
    Vector18d Stress = Damage*D_coss*sol(Eigen::seq(0,17));    
    
    Vector18d dG, dF_VE;
    Matrix18d ddG, ddF_VE;
    Vector18d Strain_VE, Strain_Pl;
    
    YieldFunctionVE(Stress, dF_VE, ddF_VE);
    Vector18d StressTrial = Damage*D_coss*(strain_total-strain_vep_old-sol(18)*dF_VE);
    PotentialFunction(Stress, StressTrial, dG, ddG);

    double DissV = 0.0, DissP = 0.0;
    Strain_New = strain_vep_old;

    // Creep
    Strain_VE = sol(18)*dF_VE;
    Strain_New += Strain_VE;

    //Plasticity
    Strain_Pl = sol(19)*dG;
    Strain_New += Strain_Pl;

    //thermal dissipation
    for (size_t i = 0; i < 18; i++){
        DissV += 1.0/dt*std::abs(Strain_VE[i]*Stress[i]);
        DissP += 1.0/dt*std::abs(Strain_Pl[i]*Stress[i]);
    }

    //hardening
    hist[0] = histOld[0] + std::sqrt(sol(19)*sol(19));
    hist[1] = sol(18);
    hist[2] = sol(19);
    hist[3] = histOld[3] + std::sqrt(sol(18)*sol(18));

    DissRate = DissV+DissP;
    DamDissRate = DissP;
}

double CosseratMaterial::YieldFunction(Vector18d& Stress, Vector18d& TrialStress, Vector18d& dF, Matrix18d& ddF){
    double F;
    
    double I1 = ic.transpose()*Stress;
    Vector18d dI1ds = ic;

    Vector18d Stress_Scaled = Stress/sc;
    double J2 = 0.5*Stress_Scaled.transpose()*J2Mat_Coss*Stress_Scaled;
    Vector18d dJ2 = J2Mat_Coss*Stress_Scaled;
    Matrix18d ddJ2 = J2Mat_Coss;
    
    bool Apex = false;
    // if (yield_alpha>0.0){
    //     double I1Apex = sy*3.0/yield_alpha;
    //     if (ic.transpose()*TrialStress>I1Apex) Apex = true;
    // }

    if (Apex==false){
        F = sc*std::sqrt(3.0*J2)-sy;//+yield_alpha/3.0*I1;
        dF  = 0.5*std::sqrt(3.0)/std::sqrt(J2+eps) * dJ2;// + yield_alpha/3.0*dI1ds;
        ddF = 0.5/sc*std::sqrt(3.0)/std::sqrt(J2+eps)*(ddJ2-0.5/(J2+eps)*dJ2*dJ2.transpose());
    } else {
        F = yield_alpha/3.0*I1-sy;
        dF = yield_alpha/3.0*dI1ds;
        ddF.setZero();
    }

    return F;
}

double CosseratMaterial::YieldFunctionVE(Vector18d& Stress, Vector18d& dF, Matrix18d& ddF){
    double F;

    Vector18d Stress_Scaled = Stress/sc;
    double J2 = 0.5*Stress_Scaled.transpose()*J2Mat_Coss*Stress_Scaled;
    Vector18d dJ2 = J2Mat_Coss*Stress_Scaled;
    Matrix18d ddJ2 = J2Mat_Coss;
    
    F = sc*std::sqrt(2.0*J2);
    dF  = 0.5*std::sqrt(2.0)/std::sqrt(J2+eps) * dJ2;
    ddF = 0.5/sc*std::sqrt(3.0)/std::sqrt(J2+eps)*(ddJ2-0.5/(J2+eps)*dJ2*dJ2.transpose());

    return F;
}

void CosseratMaterial::PotentialFunction(Vector18d& Stress, Vector18d& TrialStress, Vector18d& dG, Matrix18d& ddG){
    double I1 = ic.transpose()*Stress;
    Vector18d dI1ds = ic;

    Vector18d Stress_Scaled = Stress/sc;
    double J2 = 0.5*Stress_Scaled.transpose()*J2Mat_Coss*Stress_Scaled;
    Vector18d dJ2 = J2Mat_Coss*Stress_Scaled;
    Matrix18d ddJ2 = J2Mat_Coss;
    
    bool Apex = false;
    // if (yield_alpha>0.0){
    //     double I1Apex = sy*3.0/yield_alpha;
    //     if (ic.transpose()*TrialStress>I1Apex) Apex = true;
    // }

    if (Apex==false){
        dG  = 0.5*std::sqrt(3.0)/std::sqrt(J2+eps) * dJ2;// + pot_alpha/3.0*dI1ds;
        ddG = 0.5/sc*std::sqrt(3.0)/std::sqrt(J2+eps)*(ddJ2-0.5/(J2+eps)*dJ2*dJ2.transpose());
    } else {
        Vector18d dStress = TrialStress - 3.0*sy/yield_alpha*ic;
        Vector18d dStrain = D_coss.inverse()*dStress;

        dG = dStrain;
        ddG.setZero();
    }

    //dG.normalize();
}