#include <iostream>
#include "QuadraticSolver.h"

using namespace std;
using namespace Eigen;

// #define DEBUG
#define USE_Z_OPT
// #define CHECK_Z

double QuadraticSolver::updateVertices() {
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
    MatrixXd zValues = MatrixXd::Constant( unsupportedNodes.size(),
                                           unsupportedNodes.size() * V.cols(), 0);

    MatrixXd xyValues = MatrixXd::Constant(unsupportedNodes.size() * 2,
                      unsupportedNodes.size() * V.cols(), 0);

    VectorXd zConstant = VectorXd::Constant(unsupportedNodes.size(), 0);
    VectorXd xyConstant = VectorXd::Constant(unsupportedNodes.size() * 2, 0);

    // Go through each edge and add weights
    for (const pair<pair<int, int>, double>& edge_weight : weightMap) {
        const pair<int, int>& edge = edge_weight.first;
        double weight = edge_weight.second;

        // If this is a row to constrain
        if (unsupportedNodes.find(edge.first) != unsupportedNodes.end()) {
            zValues(indr.indexVert(edge.first),
                    indr.indexBigVert(edge.first, ZDim)) += weight;
            xyValues(indr.indexVert(edge.first) * 2,
                     indr.indexBigVert(edge.first, XDim)) += weight;
            xyValues(indr.indexVert(edge.first) * 2 + 1,
                     indr.indexBigVert(edge.first, YDim)) += weight;

        } else {
            continue;
        }

        // If this will be affected by our optimized nodes
        if (unsupportedNodes.find(edge.second) != unsupportedNodes.end()) {
            zValues(indr.indexVert(edge.first),
                    indr.indexBigVert(edge.second, ZDim)) -= weight;
            xyValues(indr.indexVert(edge.first) * 2,
                     indr.indexBigVert(edge.second, XDim)) -= weight;
            xyValues(indr.indexVert(edge.first) * 2 + 1,
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
    for (int i = 0; i < forces.cols(); i++) {
        // TODO: this requires no indexing right now, because it's all constants
        zConstant(i) += forces(0, i);
#ifdef CHECK_Z
        cout << "We added the force " << forces[0][i] << " to the " << i
             << endl;
#endif
    }

    // Our quadratic var
    MatrixXd vecToPass = MatrixXd::Constant(vec.size(), vec.size(), 0);
    MatrixXd multiplyZVal = zValues.transpose() * zValues;
    double zValWeight = 1;
    double movementWeight = 1;
#ifdef USE_Z_OPT
    vecToPass += multiplyZVal * zValWeight * zValWeight;
#endif
    for (int i = 0; i < vecToPass.rows(); i++) {
        vecToPass(i, i) += movementWeight * movementWeight;
    }
    vecToPass *= 2;

    // Set up the linear component
    VectorXd linearComp = vec;
#ifdef USE_Z_OPT
    linearComp *= movementWeight;
    linearComp += zValWeight * forces * zValues;
#endif
    linearComp *= -2;
#ifdef DEBUG
#ifdef VERBOSE
    cout << "vec used to be ";
    for (int i = 0; i < vec.size(); i++) {
        cout << vec[i] << endl;
    }
    cout << "linear comp used to be";
    for (int i = 0; i < vec.size(); i++) {
        cout << linearComp[i] << endl;
    }
    cout << "The zValues are " << endl;
    for (int i = 0; i < zValues.rows(); i++) {
        for (int j = 0; j < zValues.cols(); j++) {
            cout << zValues[i][j] << " ";
        }
        cout << endl;
    }
    cout << "The xyValues are" << endl;
    for (int i = 0; i < xyValues.rows(); i++) {
        for (int j = 0; j < xyValues.cols(); j++) {
            cout << xyValues[i][j] << " ";
        }
        cout << endl;
    }
    cout << "The xyConsts are " << endl;
    for (int i = 0; i < xyConstant.size(); i++) {
        cout << xyConstant[i] << " ";
    }
#endif
    cout << endl;
    VectorXd response = zValues * vec;
    const auto& vertVec = indr.vertVect();
    for (int i = 0; i < vertVec.size(); i++) {
        if (vertVec[i] == 0) {
            cout << "The first internal vert is" << i << endl;
        }
    }
    cout << "The before value was:\n";
    for (int i = 0; i < response.size(); i++) {
        cout << response[i] << "and the zConst was" << zConstant[i] << endl;
        cout << "The addition of resp and zConst is "
             << response[i] + zConstant[i] << endl;
        if (fabs(response[i] + zConstant[i]) > 0.2) {
            for (int j = 0; j < vertVec.size(); j++) {
                if (vertVec[j] == i) {
                    cout << "AN INTERESTING POINT IS " << j << endl;
                }
            }
            cout << endl;
        }
    }
    VectorXd responseXY = (xyValues * vec);
    cout << "The before value for xy was:\n";
    for (int i = 0; i < response.size(); i++) {
        // We are going to assert that these are small
        double responseX = responseXY[2 * i] + xyConstant[2 * i];
        double responseY = responseXY[2 * i + 1] + xyConstant[2 * i + 1];
        if (fabs(responseX) + fabs(responseY) >= pow(10, -7)) {
            cerr << "The xy diff is too big at " << fabs(responseX) << " and "
                 << fabs(responseY) << endl;
            throw new exception();
        }
    }

    // cout << responseXY[2 * i] << "and the xConst was" << xyConstant[2 * i]
    //      << endl;
    // cout << " while the yResp was" << responseXY[2 * i + 1]
    //      << "and the yConst was " << xyConstant[2 * i + 1] << endl;

    cout << "The value of vec was " << vec << endl;
#endif
    // xyConstant should be =, while zValue can be >=
    QuadProgDense prob;
    bool success = prob.solve(vecToPass, linearComp, xyValues.transpose(), xyConstant,
                  zValues.transpose(), zConstant);
    vec = prob.result();
#ifdef DEBUG
    cout << "The value of vec is now " << vec << endl;
    cout << "The success value of updating was " << success << endl;

    response = zValues * vec;
    cout << "The response value was:\n";
    for (int i = 0; i < response.size(); i++) {
        cout << response[i] << "and the zConst was" << zConstant[i] << endl;
    }
    responseXY = (dot_prod(xyValues, vec));
    cout << "The response value for xy was:\n";
    for (int i = 0; i < response.size(); i++) {
        cout << responseXY[2 * i] << "and the xConst was" << xyConstant[2 * i]
             << endl;
        cout << " while the yResp was" << responseXY[2 * i + 1]
             << "and the yConst was " << xyConstant[2 * i + 1] << endl;
    }
#endif
    moveVecIntoV();
    return success;
}

void QuadraticSolver::moveVecIntoV() {
    for (auto row : unsupportedNodes) {
        for (int j = 0; j < V.cols(); j++) {
            V(row, j) = vec[indr.indexBigVert(row, j)];
        }
    }
}
