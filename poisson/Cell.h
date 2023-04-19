#ifndef POISSON_CELL_H_
#define POISSON_CELL_H_

#include <array>
#include <Eigen/Dense>

class Cell
{
public:
    Cell(const std::array<int, 4> &,
         const std::array<Eigen::Vector3d, 4> &,
         const std::array<bool, 4> &);

    int dof() const { return 4; }
    const std::array<int, 4> dofIndices() const { return {0, 1, 2, 3}; }
    int globalIndex(int i) const { return m_globalNodeIndices[i]; }
    bool onBoundary(int i) const { return m_boundaryIndicators[i]; }
    double jacobian(int i) const { return m_jacobian[i]; }

    double shapeValue(int, const Eigen::Vector3d &) const;
    Eigen::Vector3d shapeGrad(int, const Eigen::Vector3d &) const;

private:
    std::array<int, 4> m_globalNodeIndices;
    std::array<Eigen::Vector3d, 4> m_vertices;
    std::array<bool, 4> m_boundaryIndicators;
    std::array<double, 4> m_jacobian;
    std::array<Eigen::Vector3d, 4> m_grad;
};

#endif
