#include "Solver.h"

#include <Eigen/Eigenvalues>
#include <vector>
#include <fstream>
#include <map>

Solver::Solver(const Mesh &mesh, FuncR3_R F, FuncR3_R G)
    : m_mesh(mesh),
      m_dof(mesh.dof()),
      m_systemMatrix(m_dof, m_dof),
      m_solution(Eigen::VectorXd::Zero(m_dof)),
      m_systemRhs(Eigen::VectorXd::Zero(m_dof)),
      m_boundaryValues(Eigen::VectorXd::Zero(m_dof)),
      m_F(F), m_G(G)
{
}

void Solver::assembleSystem()
{
    std::map<std::pair<int, int>, double> mapping;

    for (const Cell &cell : m_mesh.cells())
    {
        struct IntegratePoint
        {
            Eigen::Vector3d q;
            double w;
        };
        std::array<IntegratePoint, 1> integrator{{Eigen::Vector3d(0.25, 0.25, 0.25), 1.0 / 6.0}};

        Eigen::MatrixXd cellMatrix = Eigen::MatrixXd::Zero(cell.dof(), cell.dof());
        Eigen::VectorXd cellRhs = Eigen::VectorXd::Zero(cell.dof());

        for (const auto &[q, w] : integrator)
        {
            for (int i : cell.dofIndices())
                for (int j : cell.dofIndices())
                {
                    cellMatrix(i, j) += w * cell.jacobian(i) * cell.shapeGrad(i, q).dot(cell.shapeGrad(j, q));
                }

            for (int i : cell.dofIndices())
            {
                cellRhs[i] += w * cell.jacobian(i) * cell.shapeValue(i, q) * m_F(q);
            }
        }

        for (int i : cell.dofIndices())
            for (int j : cell.dofIndices())
            {
                int gi = cell.globalIndex(i), gj = cell.globalIndex(j);
                auto iter = mapping.find({gi, gj});
                if (iter == mapping.end())
                    mapping.insert({{gi, gj}, cellMatrix(i, j)});
                else
                    iter->second += cellMatrix(i, j);
            }

        for (int i : cell.dofIndices())
        {
            m_systemRhs[cell.globalIndex(i)] += cellRhs[i];
        }
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mapping.size());
    for (const auto &t : mapping)
        triplets.emplace_back(t.first.first, t.first.second, t.second);
    m_systemMatrix.setFromTriplets(triplets.begin(), triplets.end());

    interpolateBoundaryValues();
    applyBoundaryValues();
}

void Solver::interpolateBoundaryValues()
{
    for (int i = 0; i < m_dof; ++i)
        if (m_mesh.boundaryIndicators()[i])
        {
            m_boundaryValues(i) = m_G(m_mesh.nodes()[i]);
        }
}

void Solver::applyBoundaryValues()
{
    std::vector<Eigen::Triplet<double>> triplets;

    for (int i = 0; i < m_dof; ++i)
    {
        if (m_mesh.boundaryIndicators()[i])
        {
            triplets.emplace_back(i, i, 1.0);
            m_systemRhs(i) = m_boundaryValues(i);
        }
        else
        {
            for (int j = 0; j < m_dof; ++j)
            {
                if (m_mesh.boundaryIndicators()[j])
                    m_systemRhs(i) -= m_systemMatrix.coeff(i, j) * m_boundaryValues(j);
                else
                {
                    double a = m_systemMatrix.coeff(i, j);
                    if (a != 0.0)
                        triplets.emplace_back(i, j, a);
                }
            }
        }
    }
    m_systemMatrix.setFromTriplets(triplets.begin(), triplets.end());
}

void Solver::chol()
{
    m_solver.compute(m_systemMatrix);
}

double Solver::density() const
{
    return static_cast<double>(m_systemMatrix.nonZeros()) / (m_dof * m_dof);
}

double Solver::cond(double r) const
{
    int n = m_dof;
    int m = static_cast<double>(r * n + 0.5);
    Eigen::MatrixXd T(Eigen::MatrixXd::Zero(m, m));
    auto &A = m_systemMatrix;
    Eigen::MatrixXd V(n, m);
    V.col(0) = Eigen::VectorXd::Random(n).normalized();
    Eigen::VectorXd w = A * V.col(0);
    double alpha = T(0, 0) = w.dot(V.col(0));
    w -= alpha * V.col(0);

    for (int j = 1; j < m; ++j)
    {
        double beta = T(j, j - 1) = T(j - 1, j) = w.norm();
        if (beta > 0)
        {
            V.col(j) = w / beta;
        }
        else
        {
            m = j + 1;
            break;
        }

        w = A * V.col(j);
        alpha = T(j, j) = w.dot(V.col(j));
        w -= alpha * V.col(j) + beta * V.col(j - 1);
    }

    auto eig = T.block(0, 0, m, m).eigenvalues();
    auto eigr = eig.real();
    return eigr.maxCoeff() / eigr.minCoeff();
}

void Solver::solve()
{
    m_solution = m_solver.solve(m_systemRhs);
}

void Solver::writeTriplets(const std::string &filename) const
{
    std::vector<Eigen::Triplet<double>> triplet;
    for (int i = 0; i < m_systemMatrix.outerSize(); i++)
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(m_systemMatrix, i); it; ++it)
            triplet.emplace_back(it.row(), it.col(), it.value());

    std::ofstream fout(filename);
    fout << triplet.size() << '\n';
    for (const auto &t : triplet)
        fout << t.row() << ' ' << t.col() << ' ' << t.value() << '\n';
    fout.close();
}

void Solver::writeRhs(const std::string &filename) const
{
    std::ofstream fout(filename);
    fout << m_systemRhs.size() << '\n';
    for (int i = 0; i < m_systemRhs.size(); ++i)
        fout << m_systemRhs(i) << '\n';
    fout.close();
}

void Solver::writeToVtk(const std::string &filename) const
{
    std::ofstream fout(filename);
    fout << "# vtk DataFile Version 2.0\n"
         << "Poisson Solution\n"
         << "ASCII\n"
         << "DATASET UNSTRUCTURED_GRID\n"
         << "POINTS " << m_dof << " float\n";

    for (const auto &node : m_mesh.nodes())
    {
        fout << node.x() << ' ' << node.y() << ' ' << node.z() << '\n';
    }

    fout << "POINT_DATA " << m_dof << " SCALARS value float 1\n"
         << "LOOKUP_TABLE default\n";

    for (int i = 0; i < m_solution.size(); ++i)
    {
        fout << m_solution(i) << '\n';
    }

    fout.close();
}