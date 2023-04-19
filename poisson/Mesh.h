#ifndef POISSON_MESH_H_
#define POISSON_MESH_H_

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cstddef>
#include <Eigen/Core>

#include "Cell.h"

class Mesh
{
public:
    struct PhysicalName
    {
        int dim;
        std::string name;
    };

    struct Element
    {
        int type;
        int numberOfNodes;
        int numberOfTags;
        int physicalEntity;
        int modelEntity;
        std::vector<int> nodes;
    };

public:
    Mesh(const std::string &filename);

    void init();

    int dof() const { return static_cast<int>(m_nodes.size()); }
    const std::vector<Cell> &cells() const { return m_cells; }
    const std::vector<Eigen::Vector3d> &nodes() const { return m_nodes; }

    const std::vector<bool> &boundaryIndicators() const { return m_boundaryIndicators; }

private:
    void load(std::ifstream &);
    void loadMeshFormat(std::ifstream &);
    void loadPhysicalNames(std::ifstream &);
    void loadNodes(std::ifstream &);
    void loadElements(std::ifstream &);
    void endSection(std::ifstream &, const std::string &);

private:
    std::map<int, PhysicalName> m_physicalNames;

    std::vector<Eigen::Vector3d> m_nodes;
    std::map<int, int> m_nodeIndex;

    std::vector<bool> m_boundaryIndicators;

    std::vector<Element> m_elements;

    std::vector<Cell> m_cells;
};

#endif