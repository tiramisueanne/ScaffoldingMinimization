#pragma once

#include <Eigen/Dense>
#include <igl/HalfEdgeIterator.h>
#include <igl/triangle_triangle_adjacency.h>
#include <iostream>
#include <set>

using namespace std;
using namespace Eigen;

set<int> getInternalNodes(const MatrixXd &V, const MatrixXi &F) {
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
      if(borderNodes.find(v) != borderNodes.end() ||
         internalNodes.find(v) != internalNodes.end()) {
        continue;
      }
      if (halfIt.isBorder()) {
        borderNodes.insert(halfIt.Vi());
      } else {
        internalNodes.insert(halfIt.Vi());
      }
      halfIt.flipE();
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
