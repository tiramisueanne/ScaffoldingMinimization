#include <igl/copyleft/cgal/point_areas.h>
#include <igl/knn.h>
#include <igl/octree.h>
#include <igl/per_vertex_normals.h>

#include "QuadraticSolver.h"
using namespace std;
using namespace Eigen;

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
    double toAdd = 1;
    cout << "isHand is " << isHand;
    return VectorXd::Constant(V.rows(),  -1 * toAdd);
}
