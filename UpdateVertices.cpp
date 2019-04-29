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
#ifdef DEBUG
    for (int i = 0; i < V.rows(); i++) {
        for (int j = 0; j < V.cols(); j++) {
            cout << V(i, j) << " , ";
        }
        cout << endl;
    }
#endif

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
            // This needs to be negative because of negative forces
            zConstant(indr.indexVert(edge.first)) -=
                weight * V.row(edge.second).z();
#ifdef CHECK_Z
            cout << "We added " << weight * V.row(edge.second).z() << " to "
                 << indr.indexVert(edge.first) << " to now force "
                 << zConstant[indr.indexVert(edge.first)] << endl;
#endif
            // These should be summed up, not subtracted as \sum w *(v_i) = \sum w * (v_j)
            xyConstant(indr.indexVert(edge.first) * 2) +=
                weight * V.row(edge.second).x();
            xyConstant(indr.indexVert(edge.first) * 2 + 1) +=
                weight * V.row(edge.second).y();
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
    cout << "The xyValues are" << endl;
    for (int i = 0; i < xyValues.rows(); i++) {
        for (int j = 0; j < xyValues.cols(); j++) {
            cout << xyValues.coeffRef(i, j) << " , ";
        }
        cout << endl;
    }
    cout << "The xyConstants are" << endl;
    for (int i = 0; i < xyConstant.rows(); i++) {
        cout << xyConstant(i) << " , ";
    }
    cout << endl;
#endif

    VectorXi unknowns;
    VectorXd unknownVal;

    VectorXd lx;
    VectorXd lu;

    igl::SolverStatus stat =
        igl::active_set(quadCoeff, linearComp, unknowns, unknownVal, xyValues,
                        xyConstant, SparseMatrix<double>(), VectorXd(), lx, lu,
                        igl::active_set_params(), vec);
    if (stat != 0 && stat != 1) {
        cerr << "the active set solver broke in vertice updating!" << endl;
        throw new exception();
    }
    double obj = linearComp.dot(vec) + 0.5 * vec.transpose() * quadCoeff * vec;
    VectorXd vecP = VectorXd::Constant(unsupportedNodes.size() * V.cols(), 0);
    for (auto row : unsupportedNodes) {
        for (int j = 0; j < V.cols(); j++) {
            vecP[indr.indexBigVert(row, j)] = V(row, j);
        }
    }

    double objV =
        linearComp.dot(vecP) + 0.5 * vecP.transpose() * quadCoeff * vecP;
    VectorXd resXY = xyValues * vecP;
    VectorXd realXY = xyValues * vec;
    cout << "The obj function given " << obj
         << " and what our originally vec would give " << objV << endl;
    cout << "The resXY is" << resXY - xyConstant << endl;
    cout << "The realXY is " << realXY - xyConstant << endl;
    moveVecIntoV();
    return stat;
}

void QuadraticSolver::moveVecIntoV() {
#ifdef DEBUG
    MatrixX3d oldV = V;
#endif
    for (auto row : unsupportedNodes) {
        for (int j = 0; j < V.cols(); j++) {
            V(row, j) = vec(indr.indexBigVert(row, j));
        }
    }
#ifdef DEBUG
    cout << "The new V's were" << endl;
    for (int i = 0; i < V.rows(); i++) {
        for (int j = 0; j < V.cols(); j++) {
            cout << V(i, j) << " , ";
        }
        cout << endl;
    }
    cout << "The oldV was " << oldV << endl;
    cout << "The newV was " << V << endl;
#endif
}
