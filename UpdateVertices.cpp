#include <iostream>
#include "QuadraticSolver.h"

using namespace std;
using namespace Eigen;

// #define DEBUG
#define USE_Z_OPT
// #define CHECK_Z

igl::SolverStatus QuadraticSolver::updateVertices() {
    // Create the new thing to optimize, which is all the points of
    // the internal nodes
    int ZERO = 0;
    vec = VectorXd::Constant(unsupportedNodes.size() * V.cols(), 0);
    for (auto row : unsupportedNodes) {
        for (int j = 0; j < V.cols(); j++) {
            vec[indr.indexBigVert(row, j)] = V(row, j);
        }
    }

    // Create the new zDiff struct
    SparseMatrix<double> zValues(unsupportedNodes.size(),
                                           unsupportedNodes.size() * V.cols());

    SparseMatrix<double> xyValues(unsupportedNodes.size() * 2,
                      unsupportedNodes.size() * V.cols());

    VectorXd zConstant = VectorXd::Constant(unsupportedNodes.size(), 0);
    VectorXd xyConstant = VectorXd::Constant(unsupportedNodes.size() * 2, 0);

    // Go through each edge and add weights
    for (const pair<pair<int, int>, double>& edge_weight : weightMap) {
        const pair<int, int>& edge = edge_weight.first;
        double weight = edge_weight.second;

        // If this is a row to constrain
        if (unsupportedNodes.find(edge.first) != unsupportedNodes.end()) {
            zValues.coeffRef(indr.indexVert(edge.first),
                    indr.indexBigVert(edge.first, ZDim)) += weight;
            xyValues.coeffRef(indr.indexVert(edge.first) * 2,
                     indr.indexBigVert(edge.first, XDim)) += weight;
            xyValues.coeffRef(indr.indexVert(edge.first) * 2 + 1,
                     indr.indexBigVert(edge.first, YDim)) += weight;

        } else {
            continue;
        }

        // If this will be affected by our optimized nodes
        if (unsupportedNodes.find(edge.second) != unsupportedNodes.end()) {
            zValues.coeffRef(indr.indexVert(edge.first),
                    indr.indexBigVert(edge.second, ZDim)) -= weight;
            xyValues.coeffRef(indr.indexVert(edge.first) * 2,
                     indr.indexBigVert(edge.second, XDim)) -= weight;
            xyValues.coeffRef(indr.indexVert(edge.first) * 2 + 1,
                     indr.indexBigVert(edge.second, YDim)) -= weight;

        } else {
            zConstant(indr.indexVert(edge.first)) -=
                weight * V.row(edge.second).z();
#ifdef CHECK_Z
            cout << "We added " << weight * V.row(edge.second).z() << " to "
                 << indr.indexVert(edge.first) << " to now force "
                 << zConstant[indr.indexVert(edge.first)] << endl;
#endif
            xyConstant(indr.indexVert(edge.first) * 2) -=
                weight * V.row(edge.second).x();
            xyConstant(indr.indexVert(edge.first) * 2 + 1) -=
                weight * V.row(edge.second).y();
        }
    }
    // For all forces, go through and add to zValues
    for (int i = 0; i < forces.rows(); i++) {
        // TODO: this requires no indexing right now, because it's all constants
        zConstant(i) += forces(i);
    }

    // Our quadratic var
    // Make this a sparse boi
    SparseMatrix<double> quadCoeff = zValues.transpose() * zValues;
    double zValWeight = 1;
    double movementWeight = 1;
    quadCoeff *= zValWeight * zValWeight;
    for (int i = 0; i < quadCoeff.rows(); i++) {
        quadCoeff.coeffRef(i, i) += movementWeight * movementWeight;
    }
    quadCoeff *= 2;

    // Set up the linear component
    VectorXd linearComp = vec;
    linearComp *= movementWeight;
    VectorXd zValPart = forces.transpose() * zValues;
    zValPart *= zValWeight;
    linearComp += zValPart;
    linearComp *= -2;

    VectorXi unknowns;
    VectorXd unknownVal;

    VectorXd lx;
    VectorXd lu;

    // xyConstant should be =, while zValue can be >=
    igl::SolverStatus stat = igl::active_set(quadCoeff,
                                             linearComp,
                                             unknowns,
                                             unknownVal,
                                             xyValues,
                                             xyConstant,
                                             zValues,
                                             zConstant,
                                             lx,
                                             lu,
                                             igl::active_set_params(),
                                             vec);
    moveVecIntoV();
    return stat;
}

void QuadraticSolver::moveVecIntoV() {
    for (auto row : unsupportedNodes) {
        for (int j = 0; j < V.cols(); j++) {
            V(row, j) = vec[indr.indexBigVert(row, j)];
        }
    }
}
