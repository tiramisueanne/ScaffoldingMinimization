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

// DEPRECATED: use createLaplacian instead.
void QuadraticSolver::bumpInternalNodes() {
    for (auto node : unsupportedNodes) {
        V.row(node).z() += 2;
#ifdef DEBUG
        cout << " The internal node z value is " << V.row(node).z() << endl;
#endif
    }
}

void QuadraticSolver::deleteANode(int index) {
    set<int> rowsToGo;
    // Go through all triangles
    for (int i = 0; i < F.rows(); i++) {
        bool removeThis = false;
        for (int j = 0; j < F.cols(); j++) {
            if (F(i, j) == index) {
                rowsToGo.insert(index);
            }
        }
    }

    // Row removal from stack overflow
    // https://stackoverflow.com/questions/13290395/how-to-remove-a-certain-row-or-column-while-using-eigen-library-c
    for (const auto i : rowsToGo) {
        int rows = F.rows() - 1;
        int cols = F.cols();
        F.block(i, 0, rows - i, cols) =
            F.block(i + 1, 0, rows - i, F.cols());
        cout << "F is now"<< endl;
        cout << F << endl;
        F.conservativeResize(rows, cols);
    }

    int vrows = V.rows() -1;
    int cols = V.cols();
    V.block(index, 0, vrows - index, cols) =
        V.block(index + 1, 0, vrows - index, cols);
    V.conservativeResize(vrows, cols);
}

// Unsupport the node with the least force on it
void QuadraticSolver::removeSmallestNode() {
    forces = getForces();
    VectorXd force(forces.cols());
    for (int i = 0; i < forces.cols(); i++) {
        force(i) = forces(0, i);
    }
    int index;
    force.minCoeff(&index);
    deleteANode(index);
}
