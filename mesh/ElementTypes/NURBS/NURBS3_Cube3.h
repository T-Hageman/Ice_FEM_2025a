#ifndef NURBS3_Cube3_H
#define NURBS3_Cube3_H

#include <vector>
#include <Eigen/Dense>

#include "../BaseElemType.h"
#include "../IPScheme.h"

/// @brief Quadratic NURBS volume element
class NURBS3_Cube3: public BaseElemType {
    public:
        NURBS3_Cube3(int ipcount1D);
        ~NURBS3_Cube3();

        void getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData);
        void getShapeGrads2(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<Eigen::MatrixXd> &G2out, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData);
        void getCoordsIP(Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData);
        void getExportShape(std::vector<Eigen::RowVectorXd> &Nout, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData);

    private:
        void SetUpBaseShapes();

};

#endif