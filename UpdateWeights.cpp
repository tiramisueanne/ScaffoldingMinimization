#include <iostream>
#include "QuadraticSolver.h"

#define DEBUG

// V and Weights will have the same number of entries
// F just allows for us to collect faces
double QuadraticSolver::updateWeights() {
    // We will only ever constrain or check the weights of internal nodes
    // as the others are clamped down at the edges
    unsigned int rowSize = int(V.innerSize());
    unsigned int internalSize = unsupportedNodes.size();
    // Might have to update this, use the method
    forces = getForces();
#ifdef DEBUG
    cout << "The unsupportedNodes were " << endl;
    for(auto internalNode : unsupportedNodes) {
        cout << internalNode << " , ";
    }
    cout << endl;
    // cout << "The edges we have are ";
    // for (auto edge : edges) {
    //     cout << edge.first << " , " << edge.second << endl;
    // }
    cout << "The vertices are " << endl;
    for(int i = 0; i < V.rows(); i++) {
        cout << V.row(i) << endl;
    }
#endif
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
    qp::Matrix<double> zDiff(ZERO, internalSize, numEdges);

    // This matrix is two rows for every vertex: one is diff in x for
    // each edge, and one is diff in y
    // This will be constrained to zero.
    qp::Matrix<double> xyDiff(ZERO, numEdges, internalSize * 2);

    // This is the identity matrix, to help us constrain the weights
    // to be positive
    qp::Matrix<double> identity(ZERO, numEdges, numEdges);
    for (unsigned int i = 0; i < identity.nrows(); i++) {
        identity[i][i] = 1;
    }
    // Fill up zDiff through going through each triangle and filling in
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
#ifdef DEBUG
                cout << "We skipped an external Index in filling in our mats"
                     << endl;
                cout << "The node was " << currIndex << endl;
#endif
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
                if (unsupportedNodes.find(currFace(j)) == unsupportedNodes.end() &&
                    unsupportedNodes.find(currFace(other1)) ==
                        unsupportedNodes.end()) {
                    cout << "skipping an edge due to it not being in the "
                            "internal struct"
                         << endl;
                    continue;
                }
                zDiff[indr.indexVert(currIndex)]
                     [indr.indexEdge(currIndex, index)] = zDiff1;
#ifdef DEBUG_WHAT
                cout << "We inserted" << zDiff1
                     << " into zDiff:" << indr.indexVert(currIndex) << " , "
                     << indr.indexEdge(currIndex, other1) << endl;
#endif
#ifdef DEBUG
                cout << "The current vertex is " << currIndex
                     << " and the other vertex is " << index << endl;
#endif
                xyDiff[indr.indexEdge(currIndex, index)]
                      [indr.indexVert(currIndex) * 2] = xDiff1;
                xyDiff[indr.indexEdge(currIndex, index)]
                      [indr.indexVert(currIndex) * 2 + 1] = yDiff1;
#ifdef DEBUG
                cout << "Now for xyDiff" << endl;
                cout << "We inserted" << xDiff1 << " and " << yDiff1 << "into "
                     << indr.indexEdge(currIndex, index) << " , "
                     << indr.indexVert(currIndex) * 2 << endl;
#endif
            }
        }
    }
    weights = qp::Vector<double>(ONE, numEdges);
    // qp::Matrix<double> mWeights(numEdges, 1);
    // mWeights.setColumn(0, weights);
    qp::Vector<double> justOnes(ONE, numEdges);
    // The vector we add to the result of the constraint
    qp::Vector<double> justZerosForXY(ZERO, internalSize * 2);
    qp::Vector<double> justZerosForPos(ZERO, numEdges);

    // Since weights is already a row vector, we do not have to
    // transpose it
    qp::Matrix<double> zDiffT = t(zDiff);
#ifdef DEBUG
    cout << "zDiff is" << zDiff << endl;
    cout << "zDiffT is " << zDiffT << endl;
#endif
    qp::Vector<double> linearComponent(dot_prod(forces, zDiff));
    linearComponent *= 2;

    // To construct the positive definite zDiff matrix, we must
    // multiply it with itself
    zDiff = dot_prod(zDiffT, zDiff);
    zDiff *= 2;
    zDiff += (identity *= pow(10, -9));
    identity *= pow(10, 9);
#ifdef DEBUG
    cout << "The ZDiff we pass is " << zDiff << endl;
    cout << "The linearComponent is" << linearComponent << endl;
    cout << "The forces are " << forces << endl;
    cout << "The xyDiff is" << xyDiff << endl;
    cout << "The dimensions of xyDiff is" << xyDiff.nrows() << " , "
         << xyDiff.ncols() << endl;
    cout << "The identity is" << identity << endl;
#endif
    double success =
        qp::solve_quadprog(zDiff, linearComponent, xyDiff, justZerosForXY,
                           identity, justZerosForPos, weights);
    for (const auto edge : indr.edgeMap()) {
        weightMap[edge.first] = weights[edge.second];
#ifdef DEBUG
        cout << "We placed " << weights[edge.second] << " into "
             << edge.first.first << " , " << edge.first.second << endl;
#endif
    }
    return success;
}



