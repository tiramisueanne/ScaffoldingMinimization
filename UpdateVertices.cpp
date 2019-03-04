#include <iostream>

#include "QuadraticSolver.h"

using namespace std;
using namespace Eigen;

namespace qp = quadprogpp;

double QuadraticSolver::updateVertices() {
    // Create the new thing to optimize, which is all the points of
    // the internal nodes
    int ZERO = 0;
    vec = qp::Vector<double>(ZERO, internalNodes.size() * V.cols());
    for (auto row : internalNodes) {
        for (int j = 0; j < V.cols(); j++) {
            vec[indr.indexVert(row) * V.cols() + j] = V(row, j);
        }
    }

    // Create the new zDiff struct
    qp::Matrix<double> zValues(internalNodes.size(), edges.size() * V.cols());
    // Go through each edge and add weights
    for (pair<pair<int, int>, int> weight : weightMap) {
        // zValues[indr.indexVert(weight.first.first)][indr.] += weight.second;

    }
}
