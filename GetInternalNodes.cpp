#include <igl/HalfEdgeIterator.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Dense>
#include <iostream>
#include <set>

#include "QuadraticSolver.h"
using namespace std;
using namespace Eigen;

set<int> QuadraticSolver::getInternalNodes(const MatrixXi &_F,
                                           const MatrixXd &_V) {
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
    igl::triangle_triangle_adjacency(_F, TT, TTi);

    // Go through each face and add each vertex
    for (int currFace = 0; currFace < _F.rows(); currFace++) {
        // While we have already seen this node, find a new node
        // reverse = false, edge 0 face = currFace
        igl::HalfEdgeIterator<MatrixXi, MatrixXi, MatrixXi> halfIt(_F, TT, TTi,
                                                                   currFace, 0);
        set<int> edgesCovered;
        for (int i = 0; i < 3; i++) {
            int e = halfIt.Ei();
            if (edgesCovered.find(e) != edgesCovered.end()) {
                cerr << "Not covering all edges!" << endl;
                for (const auto &edge : edgesCovered) {
                    cerr << edge << " , ";
                }
                cerr << " with edge " << e;
                cerr << endl;
                throw new exception();
            } else {
                edgesCovered.insert(e);
            }

            int v = halfIt.Vi();
            if (borderNodes.find(v) == borderNodes.end() &&
                internalNodes.find(v) == internalNodes.end()) {
                if (halfIt.isBorder()) {
                    borderNodes.insert(halfIt.Vi());
                    // reverse = !reverse
                    halfIt.flipV();
                    borderNodes.insert(halfIt.Vi());
                    // reverse = !reverse
                    halfIt.flipV();
                }
            }

            // Go around the triangle
            halfIt.flipE();
            halfIt.flipV();
        }
    }
    for (int i = 0; i < _V.rows(); i++) {
        if (borderNodes.find(i) == borderNodes.end()) {
            internalNodes.insert(i);
        } else {
#ifdef DEBUG
            cout << i << " is in external " << endl;
#endif
        }
    }
    if (borderNodes.size() + internalNodes.size() != _V.innerSize()) {
        cerr << "The added sizes did not add up!" << endl;
        cerr << borderNodes.size() << endl;
        cerr << internalNodes.size() << endl;
        throw new exception();
    }
    return internalNodes;
}
