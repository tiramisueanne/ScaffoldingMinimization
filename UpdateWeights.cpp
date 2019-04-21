#include <iostream>
#include <Eigen/Sparse>
#include "QuadraticSolver.h"

// #define DEBUG
// #define CHECK_WEIGHTS

igl::SolverStatus QuadraticSolver::updateWeights() {
    // We will only ever constrain or check the weights of internal nodes
    // as the others are clamped down at the edges
    unsigned int rowSize = int(V.rows());
    unsigned int internalSize = unsupportedNodes.size();
    // Might have to update this, use the method
    forces = getForces();
    if (edges.size() % 2 != 0) {
        cerr << "We do not have an even number of edges!" << endl;
        throw new exception();
    }

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
    SparseMatrix<double> xDiff(numEdges, internalSize);
    SparseMatrix<double> yDiff(numEdges, internalSize);
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
                xDiff.coeffRef(indr.indexEdge(currIndex, index),
                                indr.indexVert(currIndex)) = xDiff1;
                yDiff.coeffRef(indr.indexEdge(currIndex, index),
                                indr.indexVert(currIndex)) = yDiff1;
            }
        }
    }
    cout << "Finished filling up the zDiff and xyDiff" << endl;
    // just set this to all zeros
    weights = VectorXd::Constant(numEdges, 0);
    // qp::Matrix<double> mWeights(numEdges, 1);
    // mWeights.setColumn(0, weights);
    VectorXd justOnes = VectorXd::Constant(numEdges, 1);
    // The vector we add to the result of the constraint
    VectorXd justZerosForInternal = VectorXd::Constant(internalSize, 0);

// Since weights is already a row vector, we do not have to
// transpose it
    cout << "Created a bunch of vectors!" << endl;
    VectorXd linearComponent(forces.transpose() * zDiff);
    linearComponent *= 2;
    cout << "The row vector is done! " << endl;

    // To construct the positive definite zDiff matrix, we must
    // multiply it with itself
    zDiff = zDiff.transpose() * zDiff;
    cout << "Done transposing and multiplying!" << endl;
    // Change this to the eigen way of doing it
    SparseMatrix<double> ident(zDiff.rows(), zDiff.rows());
    ident.setIdentity();
    zDiff = zDiff + (ident * pow(10, -15));
    zDiff *= 2;
    cout << "The zDiff is totally done!" << endl;

    xDiff = xDiff.transpose();
    yDiff = yDiff.transpose();
    cout << "Starting quadprog" << endl;
    // Fix this to be the other one
    QuadProgDense problem;

    Eigen::VectorXd emptyVec;
    Eigen::VectorXi emptyVeck;
    Eigen::VectorXd emptyVecY;
    Eigen::VectorXd allZerosLx = VectorXd::Constant(zDiff.rows(), 0);
    Eigen::VectorXd emptyVeclu;

    igl::active_set_params param = igl::active_set_params();
    param.max_iter = 200;

    igl::SolverStatus stat =
        igl::active_set(zDiff,
                        linearComponent,
                        emptyVeck,
                        emptyVecY,
                        xDiff,
                        justZerosForInternal,
                        yDiff,
                        justZerosForInternal,
                        allZerosLx,
                        emptyVeclu,
                        param,
                        weights);
    if(stat != 0 && stat != 1) {
        cerr << "The active set had a bad outcome!" << endl;
        throw new exception();
    }
    cout << "Finished quadprog!" << endl;
    for (const auto edge : indr.edgeMap()) {
        weightMap[edge.first] = weights[edge.second];
#ifdef DEBUG
        cout << "We placed " << weights[edge.second] << " into "
             << edge.first.first << " , " << edge.first.second << endl;
#endif
    }
    return stat;
}
