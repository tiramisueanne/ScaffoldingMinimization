#include <Eigen/Dense>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/HalfEdgeIterator.h>
#include <iostream>
#include <set>

#include "QuadraticSolver.h"
using namespace std;
using namespace Eigen;

set<int> QuadraticSolver::getInternalNodes() {
    // The two different sets for nodes
    set<int> borderNodes;
    set<int> internalNodes;

    // Adjacency matrix
    // TT | F | x 3 matrix where (i,j) is index of face that is adjacent
    // to triangle i on edge j
    // TTi | F | x 3 matrix wehere (i, j) is index of edge of triangle
    // TT(i,j) that is adjacent to triangle i
    Eigen::MatrixXi TT;
    Eigen::MatrixXi TTi;
    igl::triangle_triangle_adjacency(F, TT, TTi);

    // Go through each face and add each vertex
    for (int currFace = 0; currFace < F.rows(); currFace++) {
        // While we have already seen this node, find a new node
        igl::HalfEdgeIterator<MatrixXi, MatrixXi, MatrixXi> halfIt(F, TT, TTi,
                                                                   currFace, 0);
        for (int i = 0; i < 3; i++) {
            halfIt.flipV();
            int v = halfIt.Vi();
            if (borderNodes.find(v) != borderNodes.end() ||
                internalNodes.find(v) != internalNodes.end()) {
                continue;
            }
            if (halfIt.isBorder()) {
                borderNodes.insert(halfIt.Vi());
                halfIt.flipV();
                borderNodes.insert(halfIt.Vi());
                halfIt.flipV();
            }
            halfIt.flipE();
        }

    }
    for (int i = 0; i < V.rows(); i++) {
        if (borderNodes.find(i) == borderNodes.end()) {
            internalNodes.insert(i);
        } else {
            #ifdef DEBUG
            cout << i << " is in external " << endl;
            #endif
        }
    }
    if (borderNodes.size() + internalNodes.size() != V.innerSize()) {
        cerr << "The added sizes did not add up!" << endl;
        cerr << borderNodes.size() << endl;
        cerr << internalNodes.size() << endl;
        throw new exception();
    }
    return internalNodes;
}
