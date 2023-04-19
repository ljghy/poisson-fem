#ifndef POISSON_SOLVER_H_
#define POISSON_SOLVER_H_

#include "Mesh.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

#include <functional>

using FuncR3_R = std::function<double(const Eigen::Vector3d &)>;

class Solver
{
public:
    Solver(const Mesh &, FuncR3_R, FuncR3_R);
    void assembleSystem();
    void chol();
    double density() const;
    double cond(double r = 0.01) const;
    void solve();
    const Eigen::VectorXd &solution() const { return m_solution; }

    void writeTriplets(const std::string &) const;
    void writeRhs(const std::string &) const;
    void writeToVtk(const std::string &) const;

private:
    void interpolateBoundaryValues();
    void applyBoundaryValues();

private:
    const Mesh &m_mesh;
    int m_dof;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> m_solver;
    Eigen::SparseMatrix<double> m_systemMatrix;
    Eigen::VectorXd m_solution;
    Eigen::VectorXd m_systemRhs;
    Eigen::VectorXd m_boundaryValues;
    FuncR3_R m_F;
    FuncR3_R m_G;
};

#endif