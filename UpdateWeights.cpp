#include <Eigen/Sparse>
#include <iostream>
#include "QuadraticSolver.h"

#define DEBUG
#define ONE_MATRIX
// #define CHECK_WEIGHTS

igl::SolverStatus QuadraticSolver::updateWeights() {
    // We will only ever constrain or check the weights of internal nodes
    // as the others are clamped down at the edges
    unsigned int rowSize = int(V.rows());
    unsigned int internalSize = unsupportedNodes.size();

    // Might have to update this, use the method
    forces = getForces();
    #ifdef DEBUG
    cout << "The forces are" << endl;
    for(int i = 0; i < forces.rows(); i++ ) {
        cout << forces(i) << " , ";
    }
    #endif
    if (edges.size() % 2 != 0) {
        cerr << "We do not have an even number of edges!" << endl;
        throw new exception();
    }

#ifdef DEBUG
    cout << "The current V's are " << endl;
    for (int i = 0; i < V.rows(); i++) {
        for (int j = 0; j < V.cols(); j++) {
            cout << V(i, j) << " , ";
        }
        cout << endl;
    }
#endif

    // Due to the fact that each edge is represented twice in the set
    unsigned int numEdges = edges.size() / 2;
    int ZERO = 0;

    int ONE = 1;
    // This will be the array of differences of the z values
    // each row consists of a single vertex
    // initialize values to zero
    // innerSize = number of rows
    // TODO: MULTIPLY BY TWO
    SparseMatrix<double> zDiff(internalSize, numEdges);

// This matrix is two rows for every vertex: one is diff in x for
// each edge, and one is diff in y
// This will be constrained to zero.
#ifndef ONE_MATRIX
    SparseMatrix<double> xDiff(numEdges, internalSize);
    SparseMatrix<double> yDiff(numEdges, internalSize * 2);
#endif
#ifdef ONE_MATRIX
    SparseMatrix<double> xyDiff(internalSize * 2, numEdges);
#endif

// This is the identity matrix, to help us constrain the weights
// to be positive
#ifdef DEBUG_DUP
    set<int> finishedNodes;
#endif
    // Go through each face
    for (int i = 0; i < F.rows(); i++) {
        RowVector3i currFace = F.row(i);
        Matrix3d currPoints;
        for (int l = 0; l < 3; l++) {
            if (currFace(l) >= rowSize || currFace(l) < 0) {
                cerr << "The vertex is outside of the realm of the rows"
                     << endl;
                throw new exception();
            }
        }
        currPoints.row(0) = V.row(currFace(0));
        currPoints.row(1) = V.row(currFace(1));
        currPoints.row(2) = V.row(currFace(2));
        for (int j = 0; j < 3; j++) {
            int currIndex = currFace(j);
            if (unsupportedNodes.find(currIndex) == unsupportedNodes.end()) {
                continue;
            }
            // TODO: we can make this 1
            for (int k = 1; k < 3; k++) {
                int other1 = (j + k) % 3;
                int index = currFace(other1);
                // TODO: no new doubles
                const double zDiff1 = currPoints(j, 2) - currPoints(other1, 2);
                const double xDiff1 = currPoints(j, 0) - currPoints(other1, 0);
                const double yDiff1 = currPoints(j, 1) - currPoints(other1, 1);
                if (unsupportedNodes.find(currFace(j)) ==
                        unsupportedNodes.end() &&
                    unsupportedNodes.find(currFace(other1)) ==
                        unsupportedNodes.end()) {
                    cout << "skipping an edge due to it not being in the "
                         << "internal struct" << endl;
                    continue;
                }
#ifdef DEBUG_DUP
                if (finishedNodes.find(currFace(j)) != finishedNodes.end()) {
                    continue;
                } else {
                    finishedNodes.insert(currFace(j));
                }
#endif
                zDiff.coeffRef(indr.indexVert(currIndex),
                               indr.indexEdge(currIndex, index)) = zDiff1;
#ifndef ONE_MATRIX
                xDiff.coeffRef(indr.indexEdge(currIndex, index),
                               indr.indexVert(currIndex)) = xDiff1;
                yDiff.coeffRef(indr.indexEdge(currIndex, index),
                               indr.indexVert(currIndex) * 2) = yDiff1;
                yDiff.coeffRef(indr.indexEdge(currIndex, index),
                               indr.indexVert(currIndex) * 2 + 1) = -1 * yDiff1;
#endif
#ifdef ONE_MATRIX
                xyDiff.coeffRef(indr.indexVert(currIndex) * 2,
                                indr.indexEdge(currIndex, index)) = xDiff1;
                xyDiff.coeffRef(indr.indexVert(currIndex) * 2 + 1,
                                indr.indexEdge(currIndex, index)) = yDiff1;
#endif
            }
        }
    }
// just set this to all zeros
// The vector we add to the result of the constraint
#ifndef ONE_MATRIX
    VectorXd justZerosForX = VectorXd::Constant(internalSize, 0);
    VectorXd justZerosForY = VectorXd::Constant(internalSize * 2, 0);
#endif
#ifdef ONE_MATRIX
    VectorXd justZerosForXY = VectorXd::Constant(internalSize * 2, 0);
#endif

// Since weights is already a row vector, we do not have to
// transpose it
#ifdef DEBUG
    cout << "Created a bunch of vectors!" << endl;
#endif
    VectorXd linearComponent(forces.transpose() * zDiff);
    linearComponent *= 2;

#ifdef DEBUG
    cout << "The row vector is done! " << endl;
#endif
#ifdef DEBUG
    cout << "The values of zDiff" << endl;
    for (int i = 0; i < zDiff.rows(); i++) {
        for (int j = 0; j < zDiff.cols(); j++) {
            cout << zDiff.coeffRef(i, j) << " , ";
        }
        cout << endl;
    }
#endif

    // To construct the positive definite zDiff matrix, we must
    // multiply it with itself
    zDiff = zDiff.transpose() * zDiff;
#ifdef DEBUG
    cout << "Done transposing and multiplying!" << endl;
#endif
    // Change this to the eigen way of doing it
    SparseMatrix<double> ident(zDiff.rows(), zDiff.rows());
    ident.setIdentity();
    zDiff = zDiff + (ident * pow(10, -6));
    zDiff *= 2;

#ifdef DEBUG
    cout << "The zDiff is totally done!" << endl;
#endif

#ifndef ONE_MATRIX
    xDiff = xDiff.transpose();
    yDiff = yDiff.transpose();
#endif
#if defined(DEBUG)
    cout << "The values of the quadCoeff" << endl;
    for (int i = 0; i < zDiff.rows(); i++) {
        for (int j = 0; j < zDiff.cols(); j++) {
            cout << zDiff.coeffRef(i, j) << " , ";
        }
        cout << endl;
    }
    cout << "The Values of xyDiff is " << endl;
    for (int i = 0; i < xyDiff.rows(); i++) {
        for (int j = 0; j < xyDiff.cols(); j++) {
            cout << xyDiff.coeffRef(i, j) << " , ";
        }
        cout << endl;
    }
#endif

    // Fix this to be the other one
    QuadProgDense problem;

    Eigen::VectorXd emptyVec;
    Eigen::VectorXi emptyVeck;
    Eigen::VectorXd emptyVecY;
    Eigen::VectorXd allZerosLx = VectorXd::Constant(zDiff.rows(), 0);
    Eigen::VectorXd emptyVeclu;

    igl::active_set_params param = igl::active_set_params();
    param.max_iter = 200;

#ifndef ONE_MATRIX
    igl::SolverStatus stat = igl::active_set(
        zDiff, linearComponent, emptyVeck, emptyVecY, xDiff, justZerosForX,
        yDiff, justZerosForY, allZerosLx, emptyVeclu, param, weights);
#endif
#ifdef ONE_MATRIX
    igl::SolverStatus stat =
        igl::active_set(zDiff, linearComponent, emptyVeck, emptyVecY, xyDiff,
                        justZerosForXY, SparseMatrix<double>(), VectorXd(),
                        allZerosLx, emptyVeclu, param, weights);
#endif
    if (stat != 0 && stat != 1) {
        cerr << "The active set had a bad outcome!" << endl;
        throw new exception();
    }

#ifdef DEBUG
    // Calculate the residual
    double res =
        linearComponent.dot(weights) + 0.5 * weights.transpose() * zDiff * weights;
    double xyDiffRes = (xyDiff * weights).norm();
    cout << "The residual of the function is " << res << endl;
    cout << "The xyDiffRes is" << xyDiffRes << endl;
#endif

    for (const auto edge : indr.edgeMap()) {
        weightMap[edge.first] = weights[edge.second];
#ifdef DEBUG
        cout << "We placed " << weights[edge.second] << " into "
             << edge.first.first << " , " << edge.first.second << endl;
#endif
    }
    checkWeights();
    return stat;
}

bool QuadraticSolver::checkWeights() {
    // Add up all of the weight * edge for each internalNode i
    MatrixXd unsupportedNodeVals =
        MatrixXd::Constant(unsupportedNodes.size(), 3, 0);
    set<pair<int, int>> processedEdges;
    for (const auto edge : edges) {
        if (processedEdges.find(edge) == processedEdges.end()) {
            double differenceZ = V(edge.first, 2) - V(edge.second, 2);
            double differenceY = V(edge.first, 1) - V(edge.second, 1);
            double differenceX = V(edge.first, 0) - V(edge.second, 0);

            if (unsupportedNodes.find(edge.first) != unsupportedNodes.end()) {
                unsupportedNodeVals(indr.indexVert(edge.first), 2) +=
                    weightMap[edge] * differenceZ;
                unsupportedNodeVals(indr.indexVert(edge.first), 1) +=
                    weightMap[edge] * differenceY;
                unsupportedNodeVals(indr.indexVert(edge.first), 0) +=
                    weightMap[edge] * differenceX;
            }
            if (unsupportedNodes.find(edge.second) != unsupportedNodes.end()) {
                unsupportedNodeVals(indr.indexVert(edge.second), 2) +=
                    -1 * weightMap[edge] * differenceZ;

                unsupportedNodeVals(indr.indexVert(edge.second), 1) +=
                    -1 * weightMap[edge] * differenceY;
                unsupportedNodeVals(indr.indexVert(edge.second), 0) +=
                    -1 * weightMap[edge] * differenceX;
            }
            processedEdges.insert(edge);
            processedEdges.insert({edge.second, edge.first});
        }
    }
    for (const auto node : unsupportedNodes) {
#ifdef DEBUG
        cout << "This unsupported node is" << node << endl;
        cout << "The X value for this unsupported Node is "
             << unsupportedNodeVals(indr.indexVert(node), 0) << endl;
        cout << "The Y value for this unsupported node is"
             << unsupportedNodeVals(indr.indexVert(node), 1) << endl;
        cout << "The Z value for this unsupported node is"
             << unsupportedNodeVals(indr.indexVert(node), 2) << endl;
        cout << "The force for this unsupported node is"
             << forces(indr.indexVert(node)) << endl;

#endif
        assert(fabs(unsupportedNodeVals(indr.indexVert(node), 0)) < 0.01);
        assert(fabs(unsupportedNodeVals(indr.indexVert(node), 1)) < 0.01);
        assert(fabs(unsupportedNodeVals(indr.indexVert(node), 2)) > -0.0000001);
    }
    return true;
}
