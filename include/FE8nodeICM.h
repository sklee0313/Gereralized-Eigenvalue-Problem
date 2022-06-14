#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "Element.h"
using namespace Eigen;

class FE8nodeICM
{ // git branch merge https://goddaehee.tistory.com/275
    double delh = 1e-5;
    Matrix3d J, invJ;
    RowVectorXd H;
    MatrixXd dHdR, dHdx, dHdy, dHdz, dPdx, dPdy, dPdz;
    Matrix<double, 6, 24> B;
    Matrix<double, 24, 24> K;
    Matrix<double, 8, 8> M;
    Matrix<double, 24, 24> Kuu;
    Matrix<double, 24, 9> Kua;
    Matrix<double, 9, 9> Kaa;
    VectorXd weight;
    int row, col;
    Matrix<double, 6, 9> Gc, Gc_add, G;

    MatrixXd h_r, h_s, h_t, h, p_r, p_s, p_t;
    double detJ, V;

public:
    int nodaldofs = 3;

    FE8nodeICM(const MatrixXd &NQ);

    void ElementMass(const double &density, const MatrixXd &nodes, const Element &element, std::vector<Triplet<double>> &triplets, const MatrixXd &NQ);

    void ElementStiffness(const MatrixXd &C, const MatrixXd &nodes, const Element &element, std::vector<Triplet<double>> &triplets, const MatrixXd &NQ);

    void Jacobian(const MatrixXd &nodes, const Element &element, const int &i);

    void FormTriplet(const MatrixXd &K, const Element &element, std::vector<Triplet<double>> & triplets);

};
