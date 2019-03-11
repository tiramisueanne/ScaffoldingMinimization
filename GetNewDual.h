#pragma once
#include <igl/HalfEdgeIterator.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Dense>
#include <iostream>
#include <map>
#include <queue>

using namespace std;
using namespace Eigen;

bool isMatch(RowVector3i otherFace, int oneEnd, int otherEnd);
bool isOrthogonal(RowVector3d edge1, RowVector3d edge2);

MatrixXd getNewDual(Eigen::MatrixXd &V, Eigen::MatrixXi &F,
                    map<pair<int, int>, double> weights) {
    // Adjacency matrix
    // TT | F | x 3 matrix where (i,j) is index of face that is adjacent
    // to triangle i on edge j
    // TTi | F | x 3 matrix wehere (i, j) is index of edge of triangle
    // TT(i,j) that is adjacent to triangle i
    Eigen::MatrixXi TT;
    Eigen::MatrixXi TTi;
    igl::triangle_triangle_adjacency(F, TT, TTi);

    // The first face's vertex can be assigned at random,
    // here to the origin
    MatrixXd newVerts(F.rows(), 3);
    int firstFace = 0;
    newVerts.row(firstFace) << 0, 0, 0;

    // This will be the faces that we have already processed.
    set<int> filledFaces;
    queue<int> facesToFill;
    facesToFill.push(firstFace);

    MatrixXd turnLeft(3, 3);
    turnLeft << 0, 1, 0,
        -1, 0, 0,
        0, 0, 1;

    // while we still have faces to fill, fill them!
    while (!facesToFill.empty()) {
        int curFace = facesToFill.front();
        facesToFill.pop();
        if (filledFaces.find(curFace) != filledFaces.end()) {
            continue;
        } else {
            // For each edge in triangle, turn 90 degrees to left
            // and find new centers from these.
            // This only turns the correct way due to all faces being
            // clockwise
            RowVector3i faceVerts = F.row(curFace);
            MatrixXd curFaceDualCenter = newVerts.row(curFace);
            for (int vert = 0; vert < 3; vert++) {
                int adjFace = TT(curFace, vert);
                // If there is no adjacent face along this edge
                if (adjFace == -1) {
                    continue;
                }
                int nextVert = (vert + 1) % 3;
                // Currently a row vector, but we need it to be col
                RowVector3d edge = V.row(nextVert) - V.row(vert);
                // Scale and turn left
                Vector3d leftEdge = edge * turnLeft *
                                    weights.at(pair<int, int>(vert, nextVert));
                RowVector3d edgeRow = edge.row(0);
                if(!isOrthogonal(edge, leftEdge)) {
                    cerr << "The dual we are producing is not creating orthogonal edges!" << endl;
                    cerr << "the edge is " << edge << endl;
                    cerr << "The left one is " << leftEdge << endl;
                    throw new exception();
                }
                MatrixXd adjFaceDualCenter =
                    leftEdge + curFaceDualCenter.transpose().col(0);
                RowVector3d adjFaceRow = adjFaceDualCenter.transpose().row(0);
                newVerts.row(adjFace) = adjFaceRow;
                facesToFill.push(adjFace);
            }
            filledFaces.insert(curFace);
        }
    }
    if (filledFaces.size() != F.rows()) {
        cerr << "The filledFaces did not have all the faces!" << endl;
        throw new exception();
    }

    return newVerts;
}

// Dual checking
bool checkDual(MatrixXd &dualVerts, MatrixXd &V, MatrixXi &F,
               map<pair<int, int>, double> weights) {
    // go through each face and check its neighbors
    Eigen::MatrixXi TT;
    Eigen::MatrixXi TTi;
    igl::triangle_triangle_adjacency(F, TT, TTi);
    for (int i = 0; i < F.rows(); i++) {
        RowVector3d currDual = dualVerts.row(i);
        for (int j = 0; j < 3; j++) {
            // If this edge has an adjacent face, check it
            if (TT(i, j) != -1) {
                // Get diff of the
                RowVector3d neighborDual = dualVerts.row(TT(i, j));
                RowVector3d dualDiff = currDual - neighborDual;

                // TODO: make sure that this is the correct edge to look at
                int edgeIndex = j;
                int oneEnd = F(i, j);
                int otherEnd = F(i, (j + 1) % 3);
                // check if the ends actually connect on the other one
                RowVector3i otherFace = F.row(TT(i, j));
                bool matched = isMatch(otherFace, oneEnd, otherEnd);
                if (!matched) {
                    cerr << "This edge does not match!" << endl;
                    throw new exception();
                }

                RowVector3d edgeDiff = V.row(oneEnd) - V.row(otherEnd);
                if (!isOrthogonal(edgeDiff, dualDiff)) {
                    cout << "The dual edge was " << dualDiff
                         << " and the edge was" << edgeDiff << endl;
                    cout << "The currDual is" << i << endl;
                    cout << "The neighbor dual is" << TT(i, j) << endl;
                    cerr << "The dualEdge and edgeDiff are not orthogonal!"
                         << endl;
                    return false;
                }
                double dotEdge = edgeDiff.dot(edgeDiff);
                RowVector3d correctDualLength =
                    edgeDiff * weights.at({oneEnd, otherEnd});
                double dotCorrectLength =
                    correctDualLength.dot(correctDualLength);
                double dotDualLength = dualDiff.dot(dualDiff);
                if (abs(dotDualLength - dotCorrectLength) > 0.0000001) {
                    cerr << "The lengths of the edges are incorrect!" << endl;
                    return false;
                }
            } else {
                continue;
            }
        }
    }
    return true;
}

bool isOrthogonal(RowVector3d edge1, RowVector3d edge2) {
    double res = edge1.dot(edge2);
    cout << "The result of ortho is" << res << endl;
    return abs(res) < 0.0001;
}

bool correctLength(RowVector3d &primalEdge, RowVector3d &dualEdge,
                   double weight) {
    RowVector3d correctDual = primalEdge * weight;
    double corrLength = correctDual.dot(correctDual);
    double dualLength = dualEdge.dot(dualEdge);
    return (corrLength - dualLength) < 0.000001;
}

bool isMatch(RowVector3i otherFace, int oneEnd, int otherEnd) {
    int noMatch = 0;
    for (int q = 0; q < 3; q++) {
        if (otherFace(q) != oneEnd && otherFace(q) != otherEnd) {
            noMatch++;
        }
    }
    if (noMatch > 1) {
        return false;
    }
    return true;
}
