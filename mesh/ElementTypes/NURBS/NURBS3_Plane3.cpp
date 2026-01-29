#include "NURBS3_Plane3.h"
#include "../../../utility/utility.h"

NURBS3_Plane3::NURBS3_Plane3(int ipcount1D){
    Name = "NURBS3_Plane3";
    NodeCount = 16;
    requiresData = true;
    requiresNodeData = true;
    init(ipcount1D, 2);
    SetUpBaseShapes();
};

NURBS3_Plane3::~NURBS3_Plane3(){

};

void NURBS3_Plane3::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);
    G2base.resize(ipcount);

    Eigen::VectorXd Nx(4), Ny(4);
    Eigen::VectorXd Gx(4), Gy(4);
    Eigen::VectorXd G2x(4), G2y(4);
    for (size_t i = 0; i < ipcount; i++){
        Eigen::Vector2d coords = coordsIP[i];

        double x = coords(0), y = coords(1);
        BezierBasis(3,x,Nx,Gx,G2x);
        BezierBasis(3,y,Ny,Gy,G2y);

        Eigen::VectorXd Temp; Temp.resize(Nx.size()*Ny.size());
        Nbase[i].resize(NodeCount);
        KronProd(&Temp, &Nx, &Ny);
        Nbase[i] = Temp.transpose();

        Gbase[i].resize(2,NodeCount);
        KronProd(&Temp, &Gx, &Ny);
        Gbase[i](0,Eigen::indexing::all) = Temp;
        KronProd(&Temp, &Nx, &Gy);
        Gbase[i](1,Eigen::indexing::all) = Temp;

        G2base[i].resize(3,NodeCount);
        KronProd(&Temp, &G2x, &Ny);
        G2base[i](0,Eigen::indexing::all) = Temp;
        KronProd(&Temp, &Nx, &G2y);
        G2base[i](1,Eigen::indexing::all) = Temp;
        KronProd(&Temp, &Gx, &Gy);
        G2base[i](2,Eigen::indexing::all) = Temp;
    }

    NExport.resize(4);
    size_t idx = 0;
    for (size_t j = 0; j < 2; j++){
        for (size_t i = 0; i < 2; i++){
            double x = i;
            double y = j;

            BezierBasis(3,x,Nx,Gx,G2x);
            BezierBasis(3,y,Ny,Gy,G2y);

            Eigen::VectorXd Temp; Temp.resize(Nx.size()*Ny.size());
            KronProd(&Temp, &Nx, &Ny);

            NExport[idx].resize(Nx.size()*Ny.size());
            NExport[idx] = Temp;
            idx += 1;
        }
    }
    PlottingOrder.resize(4); PlottingOrder[0] = 0; PlottingOrder[1] = 1; PlottingOrder[2] = 3; PlottingOrder[3] = 2;
}

void NURBS3_Plane3::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    Eigen::MatrixXd J(2,3);
    Eigen::RowVector3d tang1, tang2, normal;

    Eigen::DiagonalMatrix<double, 16> W(NodeData);
    Eigen::VectorXd CB(NodeCount);
    Eigen::MatrixXd GB(NodeCount,2), GXi(NodeCount,2); 
    Eigen::Vector2d dW_dXi; 
    double wFun;

    for (size_t i = 0; i < ipcount; i++){
        CB = ElemData*Nbase[i].transpose();
        wFun = CB.transpose()*NodeData;

        Nout[i] = (W*CB/wFun).transpose();

        GB = ElemData*Gbase[i].transpose();
        dW_dXi = GB.transpose()*NodeData;
        GXi = W*(GB/wFun - CB*dW_dXi.transpose()/wFun/wFun);

        // Jacobian from parametric to physical space
        J = GXi.transpose()*coordsNodes; // [2 x 3]
        tang1 = J.row(0);
        tang2 = J.row(1);
        normal = tang1.cross(tang2);

        wout[i] = w[i]*abs(normal.norm());

        // Tangent directions
        Eigen::Matrix<double, 3, 2> T;
        T.col(0) = tang1.normalized();
        T.col(1) = tang2.normalized();

        // Use pseudo-inverse for mapping parametric gradients to physical gradients
        Eigen::MatrixXd J_pinv = J.completeOrthogonalDecomposition().pseudoInverse(); // [3 x 2]

        Gout[i].resize(2, NodeCount);
        for (int n = 0; n < NodeCount; ++n) {
            // Parametric gradient
            Eigen::RowVector2d grad_param = GXi.row(n);
            // Physical gradient in 3D
            Eigen::RowVector3d grad_phys = grad_param * J_pinv.transpose();
            // Project onto tangent directions
            Gout[i](0,n) = grad_phys.dot(T.col(0)); // dN/dt1
            Gout[i](1,n) = grad_phys.dot(T.col(1)); // dN/dt2
        }
    }
}

void NURBS3_Plane3::getNormals(std::vector<Eigen::VectorXd> &normals, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    Eigen::MatrixXd J(2,3);
    Eigen::Vector3d tang1, tang2;
    double JDir;

    Eigen::DiagonalMatrix<double, 16> W(NodeData);
    Eigen::VectorXd CB(NodeCount);
    Eigen::MatrixXd GB(NodeCount,2), GXi(NodeCount,2); 
    Eigen::Vector2d dW_dXi; 
    double wFun;
    Eigen::RowVectorXd N(NodeCount);

    Eigen::Vector3d CoordsIp;

    for (size_t i = 0; i < ipcount; i++){
        CB = ElemData*Nbase[i].transpose();
        wFun = CB.transpose()*NodeData;
        N = (W*CB/wFun).transpose();

        CoordsIp = N*coordsNodes;

        GB = ElemData*Gbase[i].transpose();
        dW_dXi = GB.transpose()*NodeData;
        GXi = W*(GB/wFun - CB*dW_dXi.transpose()/wFun/wFun);

        J = GXi.transpose()*coordsNodes;
        tang1 = J(0,Eigen::indexing::all);
        tang2 = J(1,Eigen::indexing::all);
        normals[i] = tang1.cross(tang2); normals[i].normalize();
        
        //ensure outwards Normal
        bool Flip = false;

        double distFromZero = CoordsIp.norm();
        CoordsIp += normals[i]*distFromZero*1e-3;
        double distFromZero2 = CoordsIp.norm();

        if (distFromZero>distFromZero2){
            Flip = true;
        }

        if (Flip){
            normals[i] = -normals[i];
        }
    }  
}

void NURBS3_Plane3::getCoordsIP(Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    for (size_t i = 0; i < ipcount; i++){
        coordsIP(i,Eigen::indexing::all) = (Nbase[i]*ElemData.transpose())*coordsNodes;
    }
}

void NURBS3_Plane3::getExportShape(std::vector<Eigen::RowVectorXd> &Nout, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    for (size_t i = 0; i < NExport.size(); i++){
        Nout[i] = NExport[i]*ElemData.transpose();
    }
}