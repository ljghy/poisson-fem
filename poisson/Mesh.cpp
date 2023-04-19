#include "Mesh.h"

#include <stdexcept>

Mesh::Mesh(const std::string &filename)
{
    std::ifstream fin(filename);
    if (!fin)
    {
        throw std::runtime_error("File not found.");
    }
    load(fin);
}

void Mesh::load(std::ifstream &fin)
{
    std::string line;
    while (true)
    {
        std::getline(fin, line);
        if (fin.eof())
        {
            fin.close();
            break;
        }

        if (!line.empty() && line[0] == '$')
        {
            if (line == "$MeshFormat")
            {
                loadMeshFormat(fin);
            }
            else if (line == "$PhysicalNames")
            {
                loadPhysicalNames(fin);
            }
            else if (line == "$Nodes")
            {
                loadNodes(fin);
            }
            else if (line == "$Elements")
            {
                loadElements(fin);
            }
        }
    }
}

void Mesh::loadMeshFormat(std::ifstream &fin)
{
    std::string version_number;
    int file_type, data_size;
    fin >> version_number >> file_type >> data_size;
    if (version_number != "2.2" || file_type != 0)
    {
        throw std::runtime_error("Unsupported file format.");
    }
    endSection(fin, "$EndMeshFormat");
}

void Mesh::loadPhysicalNames(std::ifstream &fin)
{
    int physicalNameNumber;
    fin >> physicalNameNumber;
    for (int i = 0; i < physicalNameNumber; ++i)
    {
        PhysicalName p;
        int tag;
        fin >> p.dim >> tag >> p.name;
        m_physicalNames.insert({tag, p});
    }
    endSection(fin, "$EndPhysicalNames");
}

void Mesh::loadNodes(std::ifstream &fin)
{
    int nodeNumber;
    fin >> nodeNumber;
    m_nodes.resize(nodeNumber);
    for (int i = 0; i < nodeNumber; ++i)
    {
        Eigen::Vector3d node;
        int tag;
        fin >> tag >> node.x() >> node.y() >> node.z();
        m_nodes[i] = node;
        m_nodeIndex[tag] = i;
    }
    endSection(fin, "$EndNodes");
}

void Mesh::loadElements(std::ifstream &fin)
{
    constexpr int numberOfElmNodes[]{
        0,
        2, 3, 4, 4, 8, 6, 5, 3,
        6, 0, 0, 0, 0, 14, 1, 8};

    int elementNumber;
    fin >> elementNumber;
    m_elements.resize(elementNumber);
    for (int i = 0; i < elementNumber; ++i)
    {
        int elmNumber;
        fin >> elmNumber >> m_elements[i].type >> m_elements[i].numberOfTags;
        m_elements[i].numberOfNodes = numberOfElmNodes[m_elements[i].type];
        fin >> m_elements[i].physicalEntity >> m_elements[i].modelEntity;

        for (int j = 0; j < m_elements[i].numberOfTags - 2; ++j)
        {
            int t;
            fin >> t;
        }
        m_elements[i].nodes.resize(m_elements[i].numberOfNodes);
        for (auto &node : m_elements[i].nodes)
        {
            fin >> node;
            node = m_nodeIndex[node];
        }
    }
    endSection(fin, "$EndElements");
}

void Mesh::init()
{
    m_boundaryIndicators.resize(m_nodes.size());

    for (const auto &elm : m_elements)
    {
        if (elm.type == 2) // surface
        {
            for (int i : elm.nodes)
            {
                m_boundaryIndicators[i] = true;
            }
        }
    }

    for (const auto &elm : m_elements)
    {
        if (elm.type == 4) // tetrahedron
        {
            m_cells.emplace_back(std::array<int, 4>{elm.nodes[0], elm.nodes[1], elm.nodes[2], elm.nodes[3]},
                                 std::array<Eigen::Vector3d, 4>{
                                     m_nodes[elm.nodes[0]],
                                     m_nodes[elm.nodes[1]],
                                     m_nodes[elm.nodes[2]],
                                     m_nodes[elm.nodes[3]]},
                                 std::array<bool, 4>{
                                     m_boundaryIndicators[elm.nodes[0]],
                                     m_boundaryIndicators[elm.nodes[1]],
                                     m_boundaryIndicators[elm.nodes[2]],
                                     m_boundaryIndicators[elm.nodes[3]],
                                 });
        }
    }
}

void Mesh::endSection(std::ifstream &fin, const std::string &endTag)
{
    std::string line;
    while (true)
    {
        std::getline(fin, line);
        if (fin.eof())
            throw std::runtime_error("Wrong file format.");
        if (line == endTag)
            break;
    }
}