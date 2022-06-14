#include "../include/FE8nodeICM.h"
#include <iostream>

FE8nodeICM::FE8nodeICM(const MatrixXd &NQ)
{
    B.setZero();
    G.setZero();
    K.setZero();
    M.setZero();
    h.resize(pow(NQ.rows(), 3), 8);
    h_r.resize(pow(NQ.rows(), 3), 8);
    h_s.resize(pow(NQ.rows(), 3), 8);
    h_t.resize(pow(NQ.rows(), 3), 8);
    p_r.resize(pow(NQ.rows(), 3), 3);
    p_s.resize(pow(NQ.rows(), 3), 3);
    p_t.resize(pow(NQ.rows(), 3), 3);
    p_r.setZero();
    p_s.setZero();
    p_t.setZero();
    weight.resize(pow(NQ.rows(), 3), 1);

    dHdx.resize(1, 8);
    dHdy.resize(1, 8);
    dHdz.resize(1, 8);
    dPdx.resize(pow(NQ.rows(), 3), 3);
    dPdy.resize(pow(NQ.rows(), 3), 3);
    dPdz.resize(pow(NQ.rows(), 3), 3);

    for (int i = 0, k = 0; i < NQ.rows(); i++)
    {
        for (int j = 0; j < NQ.rows(); j++)
        {
            for (int l = 0; l < NQ.rows(); l++, k++)
            {
                double r = NQ(i, 0);
                double s = NQ(j, 0);
                double t = NQ(l, 0);

                h(k, 0) = (1 - r) * (1 - s) * (1 - t) / 8;
                h(k, 1) = (1 + r) * (1 - s) * (1 - t) / 8;
                h(k, 2) = (1 + r) * (1 + s) * (1 - t) / 8;
                h(k, 3) = (1 - r) * (1 + s) * (1 - t) / 8;
                h(k, 4) = (1 - r) * (1 - s) * (1 + t) / 8;
                h(k, 5) = (1 + r) * (1 - s) * (1 + t) / 8;
                h(k, 6) = (1 + r) * (1 + s) * (1 + t) / 8;
                h(k, 7) = (1 - r) * (1 + s) * (1 + t) / 8;

                h_r(k, 0) = -(1 - s) * (1 - t) / 8;
                h_r(k, 1) = (1 - s) * (1 - t) / 8;
                h_r(k, 2) = (1 + s) * (1 - t) / 8;
                h_r(k, 3) = -(1 + s) * (1 - t) / 8;
                h_r(k, 4) = -(1 - s) * (1 + t) / 8;
                h_r(k, 5) = (1 - s) * (1 + t) / 8;
                h_r(k, 6) = (1 + s) * (1 + t) / 8;
                h_r(k, 7) = -(1 + s) * (1 + t) / 8;

                h_s(k, 0) = (1 - r) * (-1) * (1 - t) / 8;
                h_s(k, 1) = (1 + r) * (-1) * (1 - t) / 8;
                h_s(k, 2) = (1 + r) * (1 - t) / 8;
                h_s(k, 3) = (1 - r) * (1 - t) / 8;
                h_s(k, 4) = (1 - r) * (-1) * (1 + t) / 8;
                h_s(k, 5) = (1 + r) * (-1) * (1 + t) / 8;
                h_s(k, 6) = (1 + r) * (1 + t) / 8;
                h_s(k, 7) = (1 - r) * (1 + t) / 8;

                h_t(k, 0) = (1 - r) * (1 - s) * (-1) / 8;
                h_t(k, 1) = (1 + r) * (1 - s) * (-1) / 8;
                h_t(k, 2) = (1 + r) * (1 + s) * (-1) / 8;
                h_t(k, 3) = (1 - r) * (1 + s) * (-1) / 8;
                h_t(k, 4) = (1 - r) * (1 - s) / 8;
                h_t(k, 5) = (1 + r) * (1 - s) / 8;
                h_t(k, 6) = (1 + r) * (1 + s) / 8;
                h_t(k, 7) = (1 - r) * (1 + s) / 8;

                p_r(k, 0) = -2 * r;
                p_s(k, 1) = -2 * s;
                p_t(k, 2) = -2 * t;

                weight(k) = NQ(i, 1) * NQ(j, 1) * NQ(l, 1);
            }
        }
    }
}

void FE8nodeICM::ElementMass(const double &density, const MatrixXd &nodes, const Element &element, std::vector<Triplet<double>> &triplets, const MatrixXd &NQ)
{
    M.setZero();
    for (int i = 0, k = 0; i < NQ.rows(); i++)
    {
        for (int j = 0; j < NQ.rows(); j++)
        {
            for (int l = 0; l < NQ.rows(); l++, k++)
            {
                Jacobian(nodes, element, k);
                detJ = J.determinant();
                M.noalias() += h(k, all).transpose() * h(k, all) * detJ * weight(k) * density;
            }
        }
    }

    for (int i = 0; i < M.rows(); i++)
    {
        for (int j = 0; j < M.cols(); j++)
        {
            Triplet<double> tr_x(element.nodesIds[i] * 3, element.nodesIds[j] * 3, M(i, j));
            triplets.push_back(tr_x);
            Triplet<double> tr_y(element.nodesIds[i] * 3 + 1, element.nodesIds[j] * 3 + 1, M(i, j));
            triplets.push_back(tr_y);
            Triplet<double> tr_z(element.nodesIds[i] * 3 + 2, element.nodesIds[j] * 3 + 2, M(i, j));
            triplets.push_back(tr_z);
            // Triplet<double> tr_x(element.nodesIds[i], element.nodesIds[j], M(i, j));
            // triplets.push_back(tr_x);
            // Triplet<double> tr_y(element.nodesIds[i] + nodes.rows(), element.nodesIds[j] + nodes.rows(), M(i, j));
            // triplets.push_back(tr_y);
            // Triplet<double> tr_z(element.nodesIds[i] + 2 * nodes.rows(), element.nodesIds[j] + 2 * nodes.rows(), M(i, j));
            // triplets.push_back(tr_z);
        }
    }
}

void FE8nodeICM::ElementStiffness(const MatrixXd &C, const MatrixXd &nodes, const Element &element, std::vector<Triplet<double>> &triplets, const MatrixXd &NQ) //, const std::vector<double>& radius
{

    // G correction is obtained
    Gc.setZero();
    V = 0;
    K.setZero();
    for (int i = 0, k = 0; i < NQ.rows(); i++)
    {
        for (int j = 0; j < NQ.rows(); j++)
        {
            for (int l = 0; l < NQ.rows(); l++, k++)
            {
                Jacobian(nodes, element, k);

                invJ = J.inverse();
                detJ = J.determinant();

                dPdx(k, all) = p_r(k, all) * invJ(0, 0) + p_s(k, all) * invJ(0, 1) + p_t(k, all) * invJ(0, 2);
                dPdy(k, all) = p_r(k, all) * invJ(1, 0) + p_s(k, all) * invJ(1, 1) + p_t(k, all) * invJ(1, 2);
                dPdz(k, all) = p_r(k, all) * invJ(2, 0) + p_s(k, all) * invJ(2, 1) + p_t(k, all) * invJ(2, 2);

                Gc_add.setZero();
                Gc_add(0, seq(0, 2)) = dPdx(k, all);
                Gc_add(1, seq(3, 5)) = dPdy(k, all);
                Gc_add(2, seq(6, 8)) = dPdz(k, all);
                Gc_add(3, seq(0, 5)) << dPdy(k, all), dPdx(k, all);
                Gc_add(4, seq(3, 8)) << dPdz(k, all), dPdy(k, all);
                Gc_add(5, seq(0, 2)) = dPdz(k, all);
                Gc_add(5, seq(6, 8)) = dPdx(k, all);

                Gc.noalias() += Gc_add * detJ * weight(k);
                V += detJ * weight(k);
            }
        }
    }
    Gc.noalias() = -Gc / V;

    // iteration over quadrature points
    Kuu.setZero();
    Kua.setZero();
    Kaa.setZero();
    K.setZero();
    for (int i = 0, k = 0; i < NQ.rows(); i++)
    {
        for (int j = 0; j < NQ.rows(); j++)
        {
            for (int l = 0; l < NQ.rows(); l++, k++)
            {
                Jacobian(nodes, element, k);

                invJ = J.inverse();
                detJ = J.determinant();

                dHdx = h_r(k, all) * invJ(0, 0) + h_s(k, all) * invJ(0, 1) + h_t(k, all) * invJ(0, 2);
                dHdy = h_r(k, all) * invJ(1, 0) + h_s(k, all) * invJ(1, 1) + h_t(k, all) * invJ(1, 2);
                dHdz = h_r(k, all) * invJ(2, 0) + h_s(k, all) * invJ(2, 1) + h_t(k, all) * invJ(2, 2);

                B(0, seq(0, 7)) = dHdx;
                B(1, seq(8, 15)) = dHdy;
                B(2, seq(16, 23)) = dHdz;
                B(3, seq(0, 15)) << dHdy, dHdx;
                B(4, seq(8, 23)) << dHdz, dHdy;
                B(5, seq(0, 7)) = dHdz;
                B(5, seq(16, 23)) = dHdx;

                G(0, seq(0, 2)) = dPdx(k, all);
                G(1, seq(3, 5)) = dPdy(k, all);
                G(2, seq(6, 8)) = dPdz(k, all);
                G(3, seq(0, 5)) << dPdy(k, all), dPdx(k, all);
                G(4, seq(3, 8)) << dPdz(k, all), dPdy(k, all);
                G(5, seq(0, 2)) = dPdz(k, all);
                G(5, seq(6, 8)) = dPdx(k, all);

                G.noalias() = G + Gc;
                Kuu.noalias() += B.transpose() * C * B * detJ * weight(k);
                Kua.noalias() += B.transpose() * C * G * detJ * weight(k);
                Kaa.noalias() += G.transpose() * C * G * detJ * weight(k);
            }
        }
    }
    K = Kuu - Kua * (Kaa.llt().solve(Kua.transpose()));
    FormTriplet(K, element, triplets);
}
void FE8nodeICM::FormTriplet(const MatrixXd &K, const Element &element, std::vector<Triplet<double>> &triplets)
{
    for (int i = 0; i < K.rows(); i++)
    {
        for (int j = 0; j < K.cols(); j++)
        {
            if (i < 8)
                // row = element.nodesIds[i];
                row = 3 * element.nodesIds[i];
            else if (i > 7 & i < 16)
                // row = element.nodesIds[i - 8] + nodes.rows();
                row = 3 * element.nodesIds[i - 8] + 1;
            else if (i > 15 & i < 24)
                // row = element.nodesIds[i - 16] + 2 * nodes.rows();
                row = 3 * element.nodesIds[i - 16] + 2;
            else
                std::cout << "Error occurs with index assignment (row)" << std::endl;

            if (j < 8)
                // col = element.nodesIds[j];
                col = 3 * element.nodesIds[j];
            else if (j > 7 & j < 16)
                // col = element.nodesIds[j - 8] + nodes.rows();
                col = 3 * element.nodesIds[j - 8] + 1;
            else if (j > 15 & j < 24)
                // col = element.nodesIds[j - 16] + 2 * nodes.rows();
                col = 3 * element.nodesIds[j - 16] + 2;
            else
                std::cout << "Error occurs with index assignment (col)" << std::endl;

            Triplet<double> tr(row, col, K(i, j));

            triplets.push_back(tr);
        }
    }
}
void FE8nodeICM::Jacobian(const MatrixXd &nodes, const Element &element, const int &i)
{
    J << h_r(i, all) * nodes(element.nodesIds, 0), h_r(i, all) * nodes(element.nodesIds, 1), h_r(i, all) * nodes(element.nodesIds, 2),
        h_s(i, all) * nodes(element.nodesIds, 0), h_s(i, all) * nodes(element.nodesIds, 1), h_s(i, all) * nodes(element.nodesIds, 2),
        h_t(i, all) * nodes(element.nodesIds, 0), h_t(i, all) * nodes(element.nodesIds, 1), h_t(i, all) * nodes(element.nodesIds, 2);
}
