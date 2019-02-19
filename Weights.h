#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include <QuadProg++/QuadProg++.hh>
#include "GetInternalNodes.h"

using namespace std;
using namespace Eigen;
namespace qp = quadprogpp;
#define DEBUG

// A method for getting each "real" weight for the vertices
// TODO: figure out how to read in the "real" weights, as we are just using 2
// for each one right now
qp::Matrix<double> getForces(const MatrixXd &V, const MatrixXi &F,
                             int numInternal) {
    // Since we currently don't have any information on weights, just
    // return a random number for every vertex of the primal
    return qp::Matrix<double>(15.0, 1, numInternal);
}

set<pair<int, int>> allEdges(const MatrixXi &F) {
    set<pair<int, int>> edges;
    for (int currFace = 0; currFace < F.rows(); currFace++) {
        RowVector3i face = F.row(currFace);
        for (int i = 0; i < 3; i++) {
            edges.insert(pair<int, int>(face(i), (face((i + 1) % 3))));
            edges.insert(pair<int, int>(face((i + 1) % 3), (face(i))));
        }
    }
    return edges;
}

// A class to get the index into a matrix of only internal nodes
class Indexer {
    int numVert;
    vector<int> indexes;
    map<pair<int, int>, int> edgeIndex;

   public:
    Indexer(int numVertices, set<int> internalNodes,
            set<pair<int, int>> allEdges)
        : numVert(numVertices), indexes(numVert) {
        int matrixIndex = 0;
        for (int internal : internalNodes) {
            indexes[internal] = matrixIndex;
            matrixIndex++;
        }
        int matIndexEdge = 0;
        for (const auto edge : allEdges) {
            if (edgeIndex.find(edge) == edgeIndex.end()) {
                edgeIndex[edge] = matIndexEdge;
                edgeIndex[pair<int, int>(edge.second, edge.first)] =
                    matIndexEdge;
                matIndexEdge++;
            }
        }
    }

    int indexVert(int vertex) { return indexes[vertex]; }

    // If the edge is v1 - v2
    int indexEdge(int v1, int v2) { return edgeIndex[pair<int, int>(v1, v2)]; }
    const map<pair<int, int>, int> &edgeMap() { return edgeIndex; }
};

// V and Weights will have the same number of entries
// F just allows for us to collect faces
map<pair<int, int>, int> getWeights(const MatrixXd &V,
                                            const MatrixXi &F) {
    // set<int> internalNodes = getInternal(V, F);
    // We will only ever constrain or check the weights of internal nodes
    // as the others are clamped down at the edges
    set<int> internalNodes = getInternalNodes(V, F);
    unsigned int rowSize = int(V.innerSize());
    unsigned int internalSize = internalNodes.size();
    qp::Matrix<double> weights = getForces(V, F, internalSize);
    set<pair<int, int>> edges = allEdges(F);
    if (edges.size() % 2 != 0) {
        cerr << "We do not have an even number of edges!" << endl;
        throw new exception();
    }
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
    Indexer indr(rowSize, internalNodes, edges);

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
            // TODO: we can make this 1
            for (int k = 0; k < 2; k++) {
                int other1 = (j + k) % 3;
                int index = currFace(other1);
                // TODO: no new doubles
                const double zDiff1 = currPoints(j, 2) - currPoints(other1, 2);
                const double xDiff1 = currPoints(j, 0) - currPoints(other1, 0);
                const double yDiff1 = currPoints(j, 1) - currPoints(other1, 1);
                zDiff[indr.indexVert(currIndex)]
                     [indr.indexEdge(currIndex, other1)] = zDiff1;
#ifdef DEBUG
                cout << "We inserted" << zDiff1 << " into "
                     << indr.indexVert(currIndex) << " , "
                     << indr.indexEdge(currIndex, other1) << endl;
#endif
                xyDiff[indr.indexEdge(currIndex, other1)]
                      [indr.indexVert(currIndex) * 2] = xDiff1;
                xyDiff[indr.indexEdge(currIndex, other1)]
                      [indr.indexVert(currIndex) * 2 + 1] = yDiff1;
            }
        }
    }

    qp::Vector<double> x(2, numEdges);
    qp::Vector<double> justOnes(1, numEdges);
    // The vector we add to the result of the constraint
    qp::Vector<double> justZerosForXY(0.0, internalSize * 2);
    qp::Vector<double> justZerosForPos(0.0, numEdges);

    // qp::Matrix<double> weightsT = t(weights);
    qp::Vector<double> linearComponent(dot_prod(weights, zDiff));
    // To construct the positive definite zDiff matrix, we must
    // multiply it with itself
    qp::Matrix<double> zDiffT = t(zDiff);
    zDiff = dot_prod(zDiffT, zDiff);
    zDiff *= 2;
    zDiff += (identity *= pow(10, -9));
    identity /= 5;
#ifdef DEBUG
    cout << "The ZDiff we pass is " << zDiff << endl;
    cout << "The linearComponent is" << linearComponent << endl;
    cout << "The xyDiff is" << xyDiff << endl;
    cout << "The dimensions of xyDiff is" << xyDiff.nrows() << " , "
         << xyDiff.ncols() << endl;
    cout << "The identity is" << identity << endl;
#endif
    double success =
        qp::solve_quadprog(zDiff, linearComponent, xyDiff, justZerosForXY,
                           identity, justZerosForPos, x);
#ifdef DEBUG
    cout << "The success value was" << success << endl;
#endif
    map<pair<int, int>, int> toReturn;
    for (const auto edge : indr.edgeMap()) {
        toReturn[edge.first] = x[edge.second];
    }
    return toReturn;
}
