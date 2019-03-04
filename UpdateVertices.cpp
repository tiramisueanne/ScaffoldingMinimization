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
            vec[indr.indexBigVert(row, j)] = V(row, j);
        }
    }

    // Create the new zDiff struct
    qp::Matrix<double> zValues(ZERO, internalNodes.size(),
                               internalNodes.size() * V.cols());
    // Go through each edge and add weights
    for (const pair<pair<int, int>, double>& edge_weight : weightMap) {
        const pair<int, int>& edge = edge_weight.first;
        double weight = edge_weight.second;
        zValues[indr.indexVert(edge.first)][indr.indexBigVert(edge.first, 2)] +=
            weight;
        if (internalNodes.find(edge.second) != internalNodes.end()) {
            zValues[indr.indexVert(edge.first)]
                   [indr.indexBigVert(edge.second, ZDim)] -= weight;
        } else {
            // Using the X values in the first edge's x place
            zValues[indr.indexVert(edge.first)]
                [indr.indexBigVert(edge.first, XDim)] -= weight * V.row(edge.second).z();
        }
    }
}
