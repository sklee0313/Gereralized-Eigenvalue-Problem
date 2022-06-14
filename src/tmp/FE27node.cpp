#include "../include/FE9node.h"
#include <iostream>

// when FE9node is instantiated, all basic functions are evaluated at quadrature points
FE9node::FE9node(const MatrixXd &NQ)
{
    B.setZero();
    K.setZero();
    h.resize(pow(NQ.rows(), 2), 27);
    h_r.resize(pow(NQ.rows(), 2), 27);
    h_s.resize(pow(NQ.rows(), 2), 27);
    dHdx.resize(1, 27);
    dHdy.resize(1, 27);

    for (int i = 0, l = 0; i < NQ.rows(); i++)
    {
        for (int j = 0; j < NQ.rows(); j++)
        {
            for (int k = 0; k < NQ.rows(); k++, l = 9 * i + 3 * j + k)
            {

                double r = NQ(i, 0);
                double s = NQ(j, 0);
                double l = NQ(k, 0);

                h(l, 0) =
                    h(l, 1) =
                        h(l, 2) =
                            h(l, 3) =

                                h_r(l, 0) = 0.25 * (1 - s) * s + 0.5 * r * (-1 + s) * s;
                h_r(l, 1) = r * (0.5 * s - 0.5) * s + (0.25 * s - 0.25) * s;
                h_r(l, 2) = 0.5 * (r + 0.5) * s * (s + 1);
                h_r(l, 3) = 0.5 * (r - 0.5) * s * (s + 1);
                h_r(l, 4) = r * (1 - s) * s;
                h_r(l, 5) = r * (1 - pow(s, 2)) - 0.5 * pow(s, 2) + 0.5;
                h_r(l, 6) = r * (-s - 1) * s;
                h_r(l, 7) = -0.5 + 0.5 * pow(s, 2) + r * (1 - pow(s, 2));
                h_r(l, 8) = 2 * (pow(s, 2) - 1) * r;

                h_s(l, 0) = 0.25 * (r - 1) * r * (2 * s - 1);
                h_s(l, 1) = 0.25 * r * (r + 1) * (2 * s - 1);
                h_s(l, 2) = 0.5 * r * (r + 1) * (s + 0.5);
                h_s(l, 3) = 0.25 * (r - 1) * r * (2 * s + 1);
                h_s(l, 4) = pow(r, 2) * (0.5 - s) + s - 0.5;
                h_s(l, 5) = pow(r, 2) * (-s) - r * s;
                h_s(l, 6) = pow(r, 2) * (-s - 0.5) + s + 0.5;
                h_s(l, 7) = r * s - pow(r, 2) * s;
                h_s(l, 8) = 2 * (pow(r, 2) - 1) * s;
            }
        }
    }
}

void FE9node::ElementStiffness(const Matrix3d &C, const MatrixXd &nodes, const Element &element, std::vector<Triplet<double>> &triplets, const MatrixXd &NQ) //, const std::vector<double>& radius
{

    // iteration over quadrature points
    K.setZero();
    for (int i = 0, k = 0; i < NQ.rows(); i++)
    {
        for (int j = 0; j < NQ.rows(); j++, k = 3 * i + j)
        {
            Jacobian(nodes, element, k);

            // std::cout << h_r(k, all)
            //           << std::endl
            //           << h_s(k, all)
            //           << std::endl
            //           << std::endl;

            invJ = J.inverse();
            detJ = J.determinant();

            dHdx = h_r(k, all) * invJ(0, 0) + h_s(k, all) * invJ(0, 1);
            dHdy = h_r(k, all) * invJ(1, 0) + h_s(k, all) * invJ(1, 1);
            // std::cout << dHdx
            //           << std::endl
            //           << dHdy
            //           << std::endl
            //           << std::endl;

            B(0, seq(0, 8)) = dHdx;
            B(1, seq(9, 17)) = dHdy;
            B(2, seq(0, 8)) = dHdy;
            B(2, seq(9, 17)) = dHdx;

            K.noalias() += B.transpose() * C * B * detJ * NQ(i, 1) * NQ(j, 1);
        }
    }

    // std::vector<int> idx;
    // for (int i = 0; i < 3; i++)
    // {
    //     idx.push_back(3 * (element.nodesIds[i]));
    // }

    // int idx[9];
    // for (int i = 0; i < 3; i++)
    // {
    //     idx[3 * i] = 3 * (element.nodesIds[i]);
    //     idx[3 * i + 1] = 3 * (element.nodesIds[i]) + 1;
    //     idx[3 * i + 2] = 3 * (element.nodesIds[i]) + 2;
    // }

    for (int i = 0; i < K.rows(); i++)
    {
        for (int j = 0; j < K.cols(); j++)
        {
            Triplet<double> tr((i < 9) ? element.nodesIds[i] : element.nodesIds[i - 9] + nodes.rows(), (j < 9) ? element.nodesIds[j] : element.nodesIds[j - 9] + nodes.rows(), K(i, j));
            triplets.push_back(tr);
        }
    }
}

void FE9node::Jacobian(const MatrixXd &nodes, const Element &element, const int &i)
{
    // ArrayXi ind(nodes.cols());

    // ind << element.nodesIds[0], element.nodesIds[1], element.nodesIds[2];
    // std::cout << h_r(i, all) << std::endl;
    // std::cout << nodes(ind, 0) << std::endl;
    J << h_r(i, all) * nodes(element.nodesIds, 0), h_r(i, all) * nodes(element.nodesIds, 1),
        h_s(i, all) * nodes(element.nodesIds, 0), h_s(i, all) * nodes(element.nodesIds, 1);
}
