#include "Cell.h"
#include <cmath>

Cell::Cell(const std::array<int, 4> &indices,
           const std::array<Eigen::Vector3d, 4> &vertices,
           const std::array<bool, 4> &boundaryIndicators)
    : m_globalNodeIndices(indices),
      m_vertices(vertices),
      m_boundaryIndicators(boundaryIndicators)
{
    for (int i = 0; i < 4; ++i)
    {
        Eigen::Matrix3d J;
        int i1 = (i + 1) % 4,
            i2 = (i + 2) % 4,
            i3 = (i + 3) % 4;
        J << (m_vertices[i1] - m_vertices[i]),
            (m_vertices[i2] - m_vertices[i]),
            (m_vertices[i3] - m_vertices[i]);
        m_jacobian[i] = std::abs(J.determinant());

        Eigen::Vector3d normal = (m_vertices[i2] - m_vertices[i1]).cross(m_vertices[i3] - m_vertices[i1]);
        double dis = normal.dot(m_vertices[i] - m_vertices[i1]);
        m_grad[i] = normal / dis;
    }
}

double Cell::shapeValue(int, const Eigen::Vector3d &p) const
{
    return 1.0 - (p.x() + p.y() + p.z());
}

Eigen::Vector3d Cell::shapeGrad(int i, const Eigen::Vector3d &) const
{
    return m_grad[i];
}