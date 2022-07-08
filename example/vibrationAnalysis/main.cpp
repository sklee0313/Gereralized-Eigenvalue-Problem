#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <stack>
#include <ctime>

#include "../../include/Eigen/Dense"
#include "../../include/Eigen/Sparse"
#include "../../include/Spectra/MatOp/SparseGenMatProd.h"
#include "../../include/Spectra/SymGEigsShiftSolver.h"

#include "../../include/igl/slice.h"
#include "../../include/FE8nodeICM.h"
#include "../../include/Element.h"
#include "../../include/RCMordering.h"
#include "../../include/TicToc.h"

int main(int argc, char *argv[])
{

    // check if input file given as argument
    if (argc != 2)
    {
        std::cout << "usage: " << argv[0] << " <inputfile>\n";
        return 1;
    }
    std::ifstream infile(argv[1]);
    if (!infile)
    {
        std::cout << "Input file " << argv[1] << " cannot be found" << std::endl;
        return false;
    }

    // read material properties & build material-law matrix
    double nu, E, density;
    infile >> nu >> E >> density;
    Matrix<double, 6, 6> C;
    Matrix<double, 3, 3> C11;
    Matrix<double, 3, 3> C22;
    C11 << 1.0, nu / (1.0 - nu), nu / (1.0 - nu),
        nu / (1.0 - nu), 1.0, nu / (1 - nu),
        nu / (1.0 - nu), nu / (1.0 - nu), 1.0;
    C22 = Matrix<double, 3, 3>::Identity();
    C22 = C22 * (1 - 2 * nu) / (2 - 2 * nu);
    C.topLeftCorner(3, 3) = C11;
    C.bottomRightCorner(3, 3) = C22;
    C *= E * (1 - nu) / (1 + nu) / (1 - 2 * nu);

    // note that nn = number of nodes; npe = number of nodes per element; ne = number of elements
    int nn;
    infile >> nn;

    // read node positions
    MatrixXd nodes = MatrixXd::Zero(nn, 3);
    for (int i = 0; i < nn; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            infile >> nodes(i, j);
        }
    }

    int npe, ne;
    infile >> ne >> npe;

    // read elements
    std::vector<Element> elements;
    std::vector<int> nodesIds;
    int tmp;
    for (int i = 0; i < ne; i++)
    {
        for (int j = 0; j < npe; j++)
        {
            infile >> tmp;
            nodesIds.push_back(tmp);
        }
        Element element(nodesIds);
        elements.push_back(element);
        nodesIds.clear();
    }

    /////////////////////////////////
    // Stiffness and Mass matrices //
    /////////////////////////////////

    // Define the numerical quadrature; [position weight]
    MatrixXd NQ(3, 2);
    NQ << -0.7745966692414833770359, 0.5555555555555555555556,
        0, 0.8888888888888888888889,
        0.7745966692414833770359, 0.555555555555555555556;

    // instantiate the finite element method with incomaptible modes
    FE8nodeICM ElementMethod = FE8nodeICM(NQ);

    // Construction of Stiffness Matrix

    tic();
    std::vector<Triplet<double>> triplets;
    for (std::vector<Element>::iterator iter = elements.begin(); iter != elements.end(); iter++)
    {
        ElementMethod.ElementStiffness(C, nodes, *iter, triplets, NQ);
    }
    SparseMatrix<double> globalK(3 * nn, 3 * nn);
    globalK.setFromTriplets(triplets.begin(), triplets.end());
    std::cout << "Construction of stiffness matirx takes ";
    toc();

    // Construction of Mass Matrix
    tic();
    triplets.clear();
    for (std::vector<Element>::iterator iter = elements.begin(); iter != elements.end(); iter++)
    {
        ElementMethod.ElementMass(density, nodes, *iter, triplets, NQ);
    }
    SparseMatrix<double> globalM(3 * nn, 3 * nn);
    globalM.setFromTriplets(triplets.begin(), triplets.end());
    std::cout << "Construction of mass matirx takes ";
    toc();

    ///////////////////////////////
    // apply boundary conditions //
    ///////////////////////////////

    std::vector<int> II;
    for (int i = 0; i < nn; i++)
    {
        if (!(abs(nodes(i, 1)) < 1e-5 & nodes(i, 0) < 0))
        {
            II.push_back(ElementMethod.nodaldofs * i);
            II.push_back(ElementMethod.nodaldofs * i + 1);
            II.push_back(ElementMethod.nodaldofs * i + 2); // zero displacement boundary condition
        }
    }

    // apply the boundary condition to the matrices
    VectorXi freeDofsIndices = VectorXi::Map(II.data(), static_cast<int>(II.size()));
    SparseMatrix<double> globalKK(freeDofsIndices.size(), freeDofsIndices.size());
    SparseMatrix<double> globalMM(freeDofsIndices.size(), freeDofsIndices.size());
    igl::slice(globalK, freeDofsIndices, freeDofsIndices, globalKK);
    igl::slice(globalM, freeDofsIndices, freeDofsIndices, globalMM);
    std::cout << "dofs before imposing boundary condtions = " << globalK.rows() << std::endl;
    std::cout << "dofs after imposing boundary condtions = " << globalKK.rows() << std::endl;

    // optimize equation ordering
    ReorderingSSM m(globalKK);
    VectorXi r = m.ReverseCuthillMckee();
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm;
    perm.indices() = r;
    globalKK = perm.transpose() * globalKK * perm;
    globalMM = perm.transpose() * globalMM * perm;

    ////////////////////////////////////
    // generalized eigenvalue problem //
    ////////////////////////////////////

    // The generalized eigenvalue problem is solving using Spectra

    // Construct matrix operation objects using the wrapper classes
    using OpType = Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using BOpType = Spectra::SparseSymMatProd<double>;
    OpType op(globalKK, globalMM);
    BOpType Bop(globalMM);

    // seeking 10 eigenvalues
    Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::ShiftInvert>
        geigs(op, Bop, 10, 20, 0.0);

    // Initialize and compute
    tic();
    geigs.init();
    int nconv = geigs.compute(Spectra::SortRule::LargestMagn);
    std::cout << "Solution of the generalized eigenvalue problem takes ";
    toc();

    // Retrieve results
    Eigen::VectorXd evalues;
    Eigen::MatrixXd evecs;
    if (geigs.info() == Spectra::CompInfo::Successful)
    {
        evalues = geigs.eigenvalues();
        evecs = geigs.eigenvectors();
    }

    std::cout << "Number of converged generalized eigenvalues: " << nconv << std::endl;
    std::cout << "Natural frequencies found (Hz) :\n"
              << evalues.reverse().cwiseSqrt() / (2 * 3.141592) << std::endl;
}