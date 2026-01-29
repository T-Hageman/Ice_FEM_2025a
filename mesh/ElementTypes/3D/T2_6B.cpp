#include "T2_6B.h"
#include "../../../utility/utility.h"

T2_6B::T2_6B(int ipcount1D){
    Name = "T2_6B";
    NodeCount = 6;
    dim = 2;
    ipcount = getIPSchemeTriangle(&w, &coordsIP, ipcount1D, 2);
    SetUpBaseShapes();
};

T2_6B::~T2_6B(){

};

void T2_6B::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);;

    for (size_t i = 0; i < ipcount; i++){
        Eigen::Vector2d coords = coordsIP[i];
        double s = coords(0), t = coords(1);

        Nbase[i].resize(6);
        Nbase[i](0) = (1-s-t)*(1.0-s-t);
        Nbase[i](1) = s*s;
        Nbase[i](2) = t*t;
        Nbase[i](3) = 2.0*s*(1.0-s-t);
        Nbase[i](4) = 2.0*s*t;
        Nbase[i](5) = 2.0*t*(1.0-s-t);

        Gbase[i].resize(2,6);
        Gbase[i](0,0) = 2.0*s+2.0*t-2.0;    Gbase[i](1,0) = 2.0*s+2.0*t-2.0;
        Gbase[i](0,1) = 2.0*s;              Gbase[i](1,1) = 0.0;
        Gbase[i](0,2) = 0.0;                Gbase[i](1,2) = 2.0*t;
        Gbase[i](0,3) = 2.0-2.0*t-4.0*s;    Gbase[i](1,3) = -2.0*s;
        Gbase[i](0,4) = 2.0*t;              Gbase[i](1,4) = 2.0*s;
        Gbase[i](0,5) = -2.0*t;             Gbase[i](1,5) = 2.0-4.0*t-2.0*s;
    }

    NExport.resize(3);
    for (size_t j = 0; j < 3; j++){
        NExport[j].resize(6);
        NExport[j].setZero();
        NExport[j](j) = 1.0;
    }
    PlottingOrder.resize(3); 
    PlottingOrder[0] = 0; PlottingOrder[1] = 1; PlottingOrder[2] = 2;
}

void T2_6B::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
    Eigen::RowVector3d tang1, tang2, normal;
    Eigen::Matrix<double,2,3> J;
    Eigen::Vector3d JDir;

    for (size_t i = 0; i < ipcount; i++){
        Nout[i] = Nbase[i];

        Eigen::MatrixXd GB(2,6); GB = Gbase[i];

        J = GB*coordsNodes;

        tang1 = J.row(0);
        tang2 = J.row(1);
        normal = tang1.cross(tang2);
        wout[i] = w[i]*abs(normal.norm());

        Gout[i].setZero(); //not yet implemented

        // Tangent directions
        Eigen::Matrix<double, 3, 2> T;
        T.col(0) = tang1.normalized();
        T.col(1) = tang2.normalized();

        // Use pseudo-inverse for mapping parametric gradients to physical gradients
        Eigen::MatrixXd J_pinv = J.completeOrthogonalDecomposition().pseudoInverse(); // [3 x 2]

        Gout[i].resize(2, NodeCount);
        for (int n = 0; n < NodeCount; ++n) {
            // Parametric gradient
            Eigen::RowVector2d grad_param = GB.transpose().row(n);
            // Physical gradient in 3D
            Eigen::RowVector3d grad_phys = grad_param * J_pinv.transpose();
            // Project onto tangent directions
            Gout[i](0,n) = grad_phys.dot(T.col(0)); // dN/dt1
            Gout[i](1,n) = grad_phys.dot(T.col(1)); // dN/dt2
        }

        
    }
}

