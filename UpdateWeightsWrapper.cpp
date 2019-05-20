#include <eigen-quadprog/QuadProg.h>
#include <Eigen/Dense>
#include <iostream>
#include "QuadraticSolver.h"

// #define DEBUG_SIZE

double QuadraticSolver::updateWeights() {
    // We will only ever constrain or check the weights of internal nodes
    // as the others are clamped down at the edges
    unsigned int rowSize = int(V.rows());
    unsigned int internalSize = unsupportedNodes.size();

    // Might have to update this, use the method
    forces = getForces();
#ifdef DEBUG
    cout << "The forces are" << endl;
    for (int i = 0; i < forces.rows(); i++) {
        cout << forces(i) << " , ";
    }
#endif
    if (edges.size() % 2 != 0) {
        cerr << "We do not have an even number of edges!" << endl;
        throw new exception();
    }

#ifdef DEBUG_SIZE
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
    MatrixXd zDiff = MatrixXd::Constant(internalSize, numEdges, 0);

    // This matrix is two rows for every vertex: one is diff in x for
    // each edge, and one is diff in y
    // This will be constrained to zero.
    MatrixXd xyDiff = MatrixXd::Constant(internalSize * 2, numEdges, 0);

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
                zDiff(indr.indexVert(currIndex),
                      indr.indexEdge(currIndex, index)) = zDiff1;
                xyDiff(indr.indexVert(currIndex) * 2,
                       indr.indexEdge(currIndex, index)) = xDiff1;
                xyDiff(indr.indexVert(currIndex) * 2 + 1,
                       indr.indexEdge(currIndex, index)) = yDiff1;
            }
        }
    }

// just set this to all zeros
// The vector we add to the result of the constraint

// Since weights is already a row vector, we do not have to
// transpose it
#ifdef DEBUG
    cout << "Created a bunch of vectors!" << endl;
#endif
    VectorXd linearComponent(forces.transpose() * zDiff);
    linearComponent *= 2;

#ifdef DEBUG_SIZE
    cout << "The values of zDiff" << endl;
    for (int i = 0; i < zDiff.rows(); i++) {
        for (int j = 0; j < zDiff.cols(); j++) {
            cout << zDiff.coeffRef(i, j) << " , ";
        }
        cout << endl;
    }
    cout << "The forces are " << endl;
    for (int i = 0; i < forces.rows(); i++) {
        cout << forces(i) << " , ";
    }
    cout << endl;
#endif

    // To construct the positive definite zDiff matrix, we must
    // multiply it with itself
    zDiff = zDiff.transpose() * zDiff;
#ifdef DEBUG
    cout << "Done transposing and multiplying!" << endl;
#endif
    // Change this to the eigen way of doing it
    zDiff =
        zDiff + (MatrixXd::Identity(zDiff.rows(), zDiff.cols()) * pow(10, -15));
    zDiff *= 2;

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
    cout << "The values of linearComponent is" << endl;
    for (int i = 0; i < linearComponent.rows(); i++) {
        cout << linearComponent(i) << " , ";
    }
    cout << endl;
#endif

    VectorXd justZerosForXY = VectorXd::Constant(internalSize * 2, 0);
    Eigen::VectorXd allZerosForPosWeights = VectorXd::Constant(zDiff.rows(), 0);
    Eigen::MatrixXd identity = MatrixXd::Identity(zDiff.rows(), zDiff.rows());
    MatrixXd xyDiffT = xyDiff.transpose();

#ifdef DEBUG_SIZE
    cout << "The rows/cols of zDiff are " << zDiff.rows() << " and "
         << zDiff.cols() << endl;
    cout << "The rows of linear component are " << linearComponent.rows()
         << " and " << linearComponent.cols() << endl;
    cout << "The size of xyDiffT is " << xyDiffT.rows() << " and "
         << xyDiffT.cols() << endl;
    cout << "The size of justZerosForXy is " << justZerosForXY.rows() << " and "
         << justZerosForXY.cols() << endl;
    cout << "the rows of the identity is " << identity.rows() << " and "
         << identity.cols() << endl;
    cout << "the rows of the allZerosForPosWeights is "
         << allZerosForPosWeights.rows() << " and "
         << allZerosForPosWeights.cols() << endl;
#endif
    QuadProgDense solver(zDiff.rows(), xyDiff.cols(), identity.rows());
    bool converged =
        solver.solve(zDiff, linearComponent, xyDiff, justZerosForXY,
                     -1 * identity, allZerosForPosWeights);
    // Just to make sure it all gets done out
    weights = VectorXd::Constant(zDiff.rows(), 0);
    weights = solver.result();
    if (!converged) {
        cerr << "Quadprog had a bad outcome!" << endl;
        throw new exception();
        return 1;
    }
#ifdef DEBUG_SIZE
    cout << "Got out of the sovler" << endl;
#endif
    // Calculate the residual
    double res = linearComponent.dot(weights) +
                 0.5 * weights.transpose() * zDiff * weights;

#ifdef DEBUG_SIZE
    double xyDiffRes = (xyDiff * weights).norm();
    cout << "The residual of the function is " << res << endl;
    cout << "The weights are " << solver.result() << endl;
    cout << "The total forces are " << getTotalForce() << endl;
    cout << "The xyDiffRes is" << xyDiffRes << endl;
    double forcesT = forces.dot(forces);
    cout << "The forces were" << forcesT << endl;
    // Check if the constraints are met
    VectorXd resPosWeights = -1 * identity * weights;
    for (int i = 0; i < resPosWeights.rows(); i++) {
        if (resPosWeights(i) > 0.0001) {
            cout << "The constraint at " << i << " was not met " << endl;
        }
    }
#endif

    weightMap = map<pair<int, int>, double>();
    for (const auto edge : indr.edgeMap()) {
        weightMap[edge.first] = weights[edge.second];
#ifdef WEIGHTS
        cout << "We placed " << weights[edge.second] << " into "
             << edge.first.first << " , " << edge.first.second << endl;
#endif
    }

    if(!checkWeights()) {
        res = 1;
    }
    return res;
}
