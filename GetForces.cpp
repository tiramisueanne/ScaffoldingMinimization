#include <igl/copyleft/cgal/point_areas.h>
#include <igl/knn.h>
#include <igl/octree.h>
#include <igl/per_vertex_normals.h>

#include "QuadraticSolver.h"
using namespace std;
using namespace Eigen;
// #define WRONG
#define TRY_FORCES

// Calculate the Voronoi area for each vertex
Eigen::MatrixXd QuadraticSolver::getForceAreas() {
    // Get the per vertex normals
    MatrixXd normals(V.rows(), V.cols());
    igl::per_vertex_normals(V, F, normals);

    // Get the 15 nearest neighbors
    vector<vector<int>> point_indices;
    MatrixXi CH;
    MatrixXd CN;
    MatrixXd W;
    MatrixXi knn;

    igl::octree(V, point_indices, CH, CN, W);

    igl::knn(V, 15, point_indices, CH, CN, W, knn);

    MatrixXd areas;
    // Finally get the force areas
    igl::copyleft::cgal::point_areas(V, knn, normals, areas);
    return areas;
}

Eigen::VectorXd QuadraticSolver::getForceDensities() {
    return VectorXd::Constant(V.rows(), -pow(10, checkTiny() + 1));
}

// A method for getting the calculated weight for each vertex
VectorXd QuadraticSolver::getForces() {
    const int ZERO = 0;
    forces = VectorXd(unsupportedNodes.size());
#ifdef DEBUG
    cout << "the number of forces in getForces is " << forces.rows();
    cout << "The number of unsupportedNodes should be the same and is "
         << unsupportedNodes.size() << endl;
#endif
    VectorXd areas = getForceAreas();
    VectorXd densities = getForceDensities();
    for (const auto& unsupported : unsupportedNodes) {
#ifndef WRONG
        forces(indr.indexVert(unsupported)) =
            areas(unsupported) * densities(unsupported);
#endif
#ifdef WRONG
        forces(indr.indexVert(unsupported)) = densities(unsupported);
#endif
#ifdef DEBUG
        cout << "The forces are " << forces[indr.indexVert(unsupported)]
             << endl;
#endif

        assert(forces[indr.indexVert(unsupported)] <= 0);
    }
    return forces;
}
