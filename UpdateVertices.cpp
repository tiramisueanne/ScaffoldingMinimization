#include <iostream>

#include "QuadraticSolver.h"

using namespace std;
using namespace Eigen;

namespace qp = quadprogpp;
#define DEBUG

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

    qp::Matrix<double> xyValues(ZERO, internalNodes.size() * 2,
                                internalNodes.size() * V.cols());

    qp::Vector<double> zConstant(ZERO, internalNodes.size());
    qp::Vector<double> xyConstant(ZERO, internalNodes.size() * 2);

    // Go through each edge and add weights
    for (const pair<pair<int, int>, double>& edge_weight : weightMap) {
        const pair<int, int>& edge = edge_weight.first;
        double weight = edge_weight.second;

        // If this is a row to constrain
        if (internalNodes.find(edge.first) != internalNodes.end()) {
            zValues[indr.indexVert(edge.first)]
                   [indr.indexBigVert(edge.first, ZDim)] += weight;
            xyValues[indr.indexVert(edge.first) * 2]
                    [indr.indexBigVert(edge.first, XDim)] += weight;
            xyValues[indr.indexVert(edge.first) * 2 + 1]
                    [indr.indexBigVert(edge.first, YDim)] += weight;

        } else {
            continue;
        }

        // If this will be affected by our optimized nodes
        if (internalNodes.find(edge.second) != internalNodes.end()) {
            zValues[indr.indexVert(edge.first)]
                   [indr.indexBigVert(edge.second, ZDim)] -= weight;
            xyValues[indr.indexVert(edge.first) * 2]
                    [indr.indexBigVert(edge.second, XDim)] -= weight;
            xyValues[indr.indexVert(edge.first) * 2 + 1]
                    [indr.indexBigVert(edge.second, YDim)] -= weight;

        } else {
            zConstant[indr.indexVert(edge.first)] -= weight * V.row(edge.second).z();

            xyConstant[indr.indexVert(edge.first) * 2]
                -= weight * V.row(edge.second).x();
            xyConstant[indr.indexVert(edge.first) * 2 + 1]
                -= weight * V.row(edge.second).y();
        }
    }

    // Our quadratic var
    qp::Matrix<double> vecToPass(ZERO, vec.size(), vec.size());
    for (int i = 0; i < vecToPass.nrows(); i++) {
        vecToPass[i][i] += 1 + pow(10, -6);
    }
    vecToPass *= 2;

    qp::Vector<double> linearComp = vec *= -2;
    vec /= -2;
#ifdef DEBUG
    cout << "vec used to be ";
    for (int i = 0; i < vec.size(); i++) {
        cout << vec[i] << endl;
    }
    cout << "linear comp used to be";
    for (int i = 0; i < vec.size(); i++) {
        cout << linearComp[i] << endl;
    }
    cout << "The zValues are " << endl;
    for (int i = 0; i < zValues.nrows(); i++) {
        for (int j = 0; j < zValues.ncols(); j++) {
            cout << zValues[i][j] << " ";
        }
        cout << endl;
    }
#endif
    double success =
        qp::solve_quadprog(vecToPass, linearComp, t(zValues), zConstant,
                           t(xyValues), xyConstant, vec);
#ifdef DEBUG
    cout << "The success value of updating was " << success << endl;
#endif
    moveVecIntoV();
    return success;
}

void QuadraticSolver::moveVecIntoV() {
    for (int i = 0; i < vec.size(); i++) {
        V(i / V.cols(), i % V.cols()) = vec[i];
#ifdef DEBUG
        cout << "vec at " << i << " was " << vec[i] << endl;
#endif
    }
}
