#include <igl/HalfEdgeIterator.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include <QuadProg++/QuadProg++.hh>
#include "QuadraticSolver.h"

using namespace std;
using namespace Eigen;
namespace qp = quadprogpp;

#define DEBUG

// A method for getting each "real" weight for the vertices
// TODO: figure out how to read in the "real" weights, as we are just using 2
// for each one right now
qp::Matrix<double> QuadraticSolver::getForces() {
    // Since we currently don't have any information on weights, just
    // return a random number for every vertex of the primal
    return qp::Matrix<double>(-2, 1, internalNodes.size());
}

set<pair<int, int>> QuadraticSolver::allEdges() {
    set<pair<int, int>> edges;
    for (int currFace = 0; currFace < F.rows(); currFace++) {
        RowVector3i face = F.row(currFace);
        for (int i = 0; i < 3; i++) {
            // If this edge doesn't connect to an internal node
            if (internalNodes.find(face(i)) == internalNodes.end() &&
                internalNodes.find(face((i + 1) % 3)) == internalNodes.end()) {
                continue;
            }
            edges.insert(pair<int, int>(face(i), (face((i + 1) % 3))));
            edges.insert(pair<int, int>(face((i + 1) % 3), (face(i))));
        }
    }
    return edges;
}

// V and Weights will have the same number of entries
// F just allows for us to collect faces
double QuadraticSolver::updateWeights() {
    // We will only ever constrain or check the weights of internal nodes
    // as the others are clamped down at the edges
    unsigned int rowSize = int(V.innerSize());
    unsigned int internalSize = internalNodes.size();
    // Might have to update this, use the method
    forces = getForces();
#ifdef DEBUG
    cout << "The edges we have are ";
    for (auto edge : edges) {
        cout << edge.first << " , " << edge.second << endl;
    }
#endif
    if (edges.size() % 2 != 0) {
        cerr << "We do not have an even number of edges!" << endl;
        throw new exception();
    }

    // Bump internal nodes to try to force weights
    bumpInternalNodes();
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
            if (internalNodes.find(currIndex) == internalNodes.end()) {
                cout << "We skipped an external Index in filling in our mats"
                     << endl;
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
                if (internalNodes.find(currFace(j)) == internalNodes.end() &&
                    internalNodes.find(currFace(other1)) ==
                        internalNodes.end()) {
                    cout << "skipping an edge due to it not being in the "
                            "internal struct"
                         << endl;
                    continue;
                }
                zDiff[indr.indexVert(currIndex)]
                     [indr.indexEdge(currIndex, index)] = zDiff1;
#ifdef DEBUG
                cout << "We inserted" << zDiff1
                     << " into zDiff:" << indr.indexVert(currIndex) << " , "
                     << indr.indexEdge(currIndex, other1) << endl;
                cout << "The current vertex is " << currIndex
                     << " and the other vertex is " << index << endl;
#endif
                xyDiff[indr.indexEdge(currIndex, index)]
                      [indr.indexVert(currIndex) * 2] = xDiff1;
                xyDiff[indr.indexEdge(currIndex, index)]
                      [indr.indexVert(currIndex) * 2 + 1] = yDiff1;
#ifdef DEBUG
// cout << "Now for xyDiff" << endl;
// cout << "We inserted" << xDiff1 << " and "
//      << yDiff1 << "into "
//      << indr.indexEdge(currIndex, index) << " , "
//      << indr.indexVert(currIndex) * 2 << endl;
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
        cout << "We placed " << weights[edge.second] << " into " << edge.first.first
             << " , " << edge.first.second << endl;
#endif
    }
    return success;
}

void QuadraticSolver::bumpInternalNodes() {
    for (auto node : internalNodes) {
        V.row(node).z() += 2;
#ifdef DEBUG
        cout << " The internal node z value is " << V.row(node).z() << endl;
#endif
    }
}
