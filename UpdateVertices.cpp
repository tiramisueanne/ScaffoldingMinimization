#include <iostream>
#include "QuadraticSolver.h"

using namespace std;
using namespace Eigen;

// #define DEBUG
#define USE_Z_OPT
// #define CHECK_Z
#define ONE_MATRIX

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
#ifndef ONE_MATRIX
    SparseMatrix<double> xValues(unsupportedNodes.size(),
                                 unsupportedNodes.size() * V.cols());
    SparseMatrix<double> yValues(unsupportedNodes.size() * 2,
                                 unsupportedNodes.size() * V.cols());
#endif
#ifdef ONE_MATRIX
    SparseMatrix<double> xyValues(unsupportedNodes.size() * 2,
                                  unsupportedNodes.size() * V.cols());
#endif

    VectorXd zConstant = VectorXd::Constant(unsupportedNodes.size(), 0);
#ifndef ONE_MATRIX
    VectorXd xConstant = VectorXd::Constant(unsupportedNodes.size(), 0);
    VectorXd yConstant = VectorXd::Constant(unsupportedNodes.size() * 2, 0);
#endif
#ifdef ONE_MATRIX
    VectorXd xyConstant = VectorXd::Constant(unsupportedNodes.size() * 2, 0);
#endif

    // Go through each edge and add weights
    for (const pair<pair<int, int>, double>& edge_weight : weightMap) {
        const pair<int, int>& edge = edge_weight.first;
        double weight = edge_weight.second;

        // If this is a row to constrain
        if (unsupportedNodes.find(edge.first) != unsupportedNodes.end()) {
            zValues.coeffRef(indr.indexVert(edge.first),
                             indr.indexBigVert(edge.first, ZDim)) += weight;
#ifndef ONE_MATRIX
            xValues.coeffRef(indr.indexVert(edge.first),
                             indr.indexBigVert(edge.first, XDim)) += weight;
            yValues.coeffRef(indr.indexVert(edge.first) * 2,
                             indr.indexBigVert(edge.first, YDim)) += weight;
            yValues.coeffRef(indr.indexVert(edge.first) * 2 + 1,
                             indr.indexBigVert(edge.first, YDim)) -= weight;
#endif
#ifdef ONE_MATRIX
            xyValues.coeffRef(indr.indexVert(edge.first) * 2,
                              indr.indexBigVert(edge.first, YDim)) += weight;
            xyValues.coeffRef(indr.indexVert(edge.first) * 2 + 1,
                              indr.indexBigVert(edge.first, YDim)) += weight;
#endif
        } else {
            continue;
        }

        // If this will be affected by our optimized nodes
        if (unsupportedNodes.find(edge.second) != unsupportedNodes.end()) {
            zValues.coeffRef(indr.indexVert(edge.first),
                             indr.indexBigVert(edge.second, ZDim)) -= weight;
#ifndef ONE_MATRIX
            xValues.coeffRef(indr.indexVert(edge.first),
                             indr.indexBigVert(edge.second, XDim)) -= weight;
            yValues.coeffRef(indr.indexVert(edge.first) * 2,
                             indr.indexBigVert(edge.second, YDim)) -= weight;
            yValues.coeffRef(indr.indexVert(edge.first) * 2 + 1,
                             indr.indexBigVert(edge.second, YDim)) += weight;
#endif
#ifdef ONE_MATRIX
            xyValues.coeffRef(indr.indexVert(edge.first) * 2,
                              indr.indexBigVert(edge.second, YDim)) -= weight;
            xyValues.coeffRef(indr.indexVert(edge.first) * 2 + 1,
                              indr.indexBigVert(edge.second, YDim)) -= weight;
#endif
        } else {
            zConstant(indr.indexVert(edge.first)) -=
                weight * V.row(edge.second).z();
#ifdef CHECK_Z
            cout << "We added " << weight * V.row(edge.second).z() << " to "
                 << indr.indexVert(edge.first) << " to now force "
                 << zConstant[indr.indexVert(edge.first)] << endl;
#endif
#ifndef ONE_MATRIX
            xConstant(indr.indexVert(edge.first)) -=
                weight * V.row(edge.second).x();

            yConstant(indr.indexVert(edge.first) * 2) -=
                weight * V.row(edge.second).y();
            yConstant(indr.indexVert(edge.first) * 2 + 1) +=
                weight * V.row(edge.second).y();
#endif
#ifdef ONE_MATRIX
            xyConstant(indr.indexVert(edge.first) * 2) -=
                weight * V.row(edge.second).y();
            xyConstant(indr.indexVert(edge.first) * 2 + 1) -=
                weight * V.row(edge.second).y();
#endif
        }
    }
#ifdef DEBUG
    cout << "The zConstant is " << endl;
#endif
    // For all forces, go through and add to zValues
    for (int i = 0; i < forces.rows(); i++) {
        // TODO: this requires no indexing right now, because it's all constants
        zConstant(i) += forces(i);
#ifdef DEBUG
        cout << zConstant << " , " << endl;
#endif
    }

    // Our quadratic var
    // Make this a sparse boi
    SparseMatrix<double> quadCoeff = zValues.transpose() * zValues;
    double zValWeight = 1;
    double movementWeight = 0.01;
    quadCoeff *= zValWeight;
    for (int i = 0; i < quadCoeff.rows(); i++) {
        quadCoeff.coeffRef(i, i) += movementWeight * 1;
    }
    quadCoeff *= 2;

    // Set up the linear component
    VectorXd linearComp = vec;
    linearComp *= -2 * movementWeight;
    VectorXd zValPart = zConstant.transpose() * zValues;
    zValPart *= 2 * zValWeight;
    linearComp += zValPart;
#ifdef DEBUG
    cout << "The vec is " << endl;
    for (int i = 0; i < vec.rows(); i++) {
        cout << vec(i) << " , ";
    }
    cout << endl;
    cout << "The zValues are" << endl;
    for (int i = 0; i < zValues.rows(); i++) {
        for (int j = 0; j < zValues.cols(); j++) {
            cout << zValues.coeffRef(i, j) << " , ";
        }
        cout << endl;
    }
    cout << endl;
    cout << "The quadCoeff is" << endl;
    for (int i = 0; i < quadCoeff.rows(); i++) {
        for (int j = 0; j < quadCoeff.cols(); j++) {
            cout << quadCoeff.coeffRef(i, j) << " , ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "Linear comp is" << endl;
    for (int i = 0; i < linearComp.rows(); i++) {
        cout << linearComp(i) << " , ";
    }
    cout << endl;
#endif

    VectorXi unknowns;
    VectorXd unknownVal;

    VectorXd lx;
    VectorXd lu;

#ifndef ONE_MATRIX
    // xyConstant should be =, while zValue can be >=
    igl::SolverStatus stat = igl::active_set(
        quadCoeff, linearComp, unknowns, unknownVal, xValues, xConstant,
        yValues, yConstant, lx, lu, igl::active_set_params(), vec);
#endif
#ifdef ONE_MATRIX
    igl::SolverStatus stat = igl::active_set(
        quadCoeff, linearComp, unknowns, unknownVal, xyValues, xyConstant,
        SparseMatrix<double>(), VectorXd(), lx, lu, igl::active_set_params(), vec);
#endif
    moveVecIntoV();
    return stat;
}

void QuadraticSolver::moveVecIntoV() {
    #ifdef DEBUG
    cout << "The new V's were" << endl;
    #endif
    for (auto row : unsupportedNodes) {
        for (int j = 0; j < V.cols(); j++) {
            V(row, j) = vec(indr.indexBigVert(row, j));
            #ifdef DEBUG
            cout << V(row, j) << " , ";
            #endif
        }
        #ifdef DEBUG
        cout << endl;
        #endif
    }
}
