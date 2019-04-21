#include <igl/HalfEdgeIterator.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include <QuadProg++/QuadProg++.hh>
#include "QuadraticSolver.h"

using namespace std;
using namespace Eigen;
namespace qp = quadprogpp;

// #define DEBUG

// A method for getting the calculated weight for each vertex
MatrixXd QuadraticSolver::getForces() {
    const int ZERO = 0;
    forces = VectorXd(unsupportedNodes.size());
    VectorXd areas = getForceAreas();
    VectorXd densities = getForceDensities();
    for (const auto& unsupported : unsupportedNodes) {
        forces[indr.indexVert(unsupported)] = areas[indr.indexVert(unsupported)] * densities[indr.indexVert(unsupported)];
        assert(forces[indr.indexVert(unsupported)] <= 0);
    }
    return forces;
}

set<pair<int, int>> QuadraticSolver::allEdges() {
    set<pair<int, int>> edges;
    for (int currFace = 0; currFace < F.rows(); currFace++) {
        RowVector3i face = F.row(currFace);
        for (int i = 0; i < 3; i++) {
            // If this edge doesn't connect to an internal node
            if (unsupportedNodes.find(face(i)) == unsupportedNodes.end() &&
                unsupportedNodes.find(face((i + 1) % 3)) ==
                    unsupportedNodes.end()) {
                continue;
            }
            edges.insert(pair<int, int>(face(i), (face((i + 1) % 3))));
            edges.insert(pair<int, int>(face((i + 1) % 3), (face(i))));
        }
    }
    return edges;
}

void QuadraticSolver::bumpInternalNodes() {
    for (auto node : unsupportedNodes) {
        V.row(node).z() += 2;
#ifdef DEBUG
        cout << " The internal node z value is " << V.row(node).z() << endl;
#endif
    }
}

// Unsupport the node with the least force on it
void QuadraticSolver::deleteANode() {
    forces = getForces();
    VectorXd force(forces.cols());
    for(int i = 0; i < forces.cols(); i++) {
        force(i) = forces(0, i);
    }
    int index;
    force.minCoeff(&index);
    unsupportedNodes.insert(index);
}
