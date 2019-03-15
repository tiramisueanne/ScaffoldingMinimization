#include <igl/HalfEdgeIterator.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Dense>
#include <iostream>
#include <set>

#include "QuadraticSolver.h"
using namespace std;
using namespace Eigen;
// #define DEBUG
#define UNSUPPORT_A_NODE

set<int> QuadraticSolver::getUnsupportedNodes(const MatrixXi &_F,
                                              const MatrixXd &_V) {
    // The two different sets for nodes
    set<int> supportedNodes;
    set<int> unsupportedNodes;

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
        set<int> vertexesCovered;
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
            int firstV = halfIt.Vi();
            // reverse = !reverse
            halfIt.flipV();
            int secondV = halfIt.Vi();
            // reverse = !reverse
            halfIt.flipV();

            if (halfIt.isBorder()) {
                supportedNodes.insert(firstV);
                supportedNodes.insert(secondV);
            }

            // Go around the triangle
            halfIt.flipE();
            halfIt.flipV();
        }
    }
    for (int i = 0; i < _V.rows(); i++) {
        if (supportedNodes.find(i) == supportedNodes.end()) {
            unsupportedNodes.insert(i);
        } else {
#ifdef DEBUG
            cout << i << " is in external " << endl;
#endif
        }
    }
    if (supportedNodes.size() + unsupportedNodes.size() != _V.innerSize()) {
        cerr << "The added sizes did not add up!" << endl;
        cerr << supportedNodes.size() << endl;
        cerr << unsupportedNodes.size() << endl;
        throw new exception();
    }
#ifdef DEBUG
    cout << "The unsupportedNodes were " << endl;
    for (auto &i : unsupportedNodes) {
        cout << i << " , ";
    }
    cout << endl;
    cout << "Size of unsupportedNodes was" << unsupportedNodes.size() << endl;
#endif
#ifdef UNSUPPORT_A_NODE
    unsupportedNodes.insert(1);
#endif
    return unsupportedNodes;
}
