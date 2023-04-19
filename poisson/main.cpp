#include <iostream>
#include <cmath>

#include "Mesh.h"
#include "Solver.h"

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        std::cout << "No input file.\n";
        return -1;
    }

    Mesh mesh(argv[1]);
    mesh.init();
    std::cout << "DOF: " << mesh.dof() << '\n';

    Solver solver(
        mesh,
        [](const Eigen::Vector3d &p)
        {
            return p.dot(p);
        }, // RHS
        [](const Eigen::Vector3d &p)
        {
            // return std::cos(std::sqrt(p.dot(p)));
            double r = std::sqrt(p(0) * p(0) + p(1) * p(1));
            double theta = std::atan2(p(2), r), phi = std::atan2(p(1), p(0));
            return std::sin(2.0 * theta) * std::cos(phi);
        } // Boundary
    );

    std::cout << "Assembling system...\n";
    solver.assembleSystem();

    std::cout << "System matrix density: " << solver.density() << '\n'
              << "Approximate condition number: " << solver.cond() << '\n';

    std::cout << "Performing Cholesky decomposition...\n";
    solver.chol();

    std::cout << "Solving...\n";
    solver.solve();
    std::string vtkFilename = "solution.vtk";
    solver.writeToVtk(vtkFilename);
    std::cout << "Solution written to \"" << vtkFilename << "\"\n";

    solver.writeTriplets("triplets.txt");
    solver.writeRhs("rhs.txt");
}