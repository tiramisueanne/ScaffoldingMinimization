#pragma once
#include <igl/HalfEdgeIterator.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Dense>
#include <iostream>
#include <map>
#include <queue>

using namespace std;
using namespace Eigen;

void getNewDual(Eigen::MatrixXd &V, Eigen::MatrixXi &F,
                map<pair<int, int>, int> weights) {
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
    int firstFace = 0;
    newVerts.row(firstFace) << 0, 0, 0;

    // This will be the faces that we have already processed.
    set<int> filledFaces;
    queue<int> facesToFill;
    facesToFill.insert(currentFace);

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
            RowVector3d faceVerts = F.row(curFace);
            RowVector3d curFaceDualCenter = newVerts.row(curFace);
            for (int vert = 0; vert < 3; vert++) {
                int nextVert = (vert + 1) % 3;
                RowVector3d edge = V.row(nextVert) - V.row(vert);
                // Scale and turn left
                RowVector3d leftEdge =
                    turnLeft * edge * weights.at(vert, nextVert);
                RowVector3d adjFaceDualCenter = leftEdge + curFaceDualCenter;
                int adjFace = TT(curFace, vert);
                newVerts.row(adjFace) = adjFaceDualCenter;
                facesToFill.push(adjFace);
            }
        }
    }
}
