#include <Eigen/Dense>
#include <queue>
#include <set>

#include "Gradient.h"
using namespace std;
using namespace Eigen;

// A flood fill of computing gradients and new vertices
MatrixXd getNewVertices(MatrixXd &reciprocalDual, MatrixXd &V, MatrixXi &F) {
  // Adjacency matrix
  // TT | F | x 3 matrix where (i,j) is index of face that is adjacent
  // to triangle i on edge j
  // TTi | F | x 3 matrix wehere (i, j) is index of edge of triangle
  // TT(i,j) that is adjacent to triangle i
  Eigen::MatrixXi TT;
  Eigen::MatrixXi TTi;
  igl::triangle_triangle_adjacency(F, TT, TTi);


  // The very first face we can set to anything, everything will connect
  // due to the reciprocal graph's edges adding up to 0
  int firstFace = 0;
  MatrixXd gradients(F.rows(), 3);
  MatrixXd newVerts(V.rows(), 3);
  gradients.row(firstFace) = getGradient(
                                         V.row(F(firstFace, 0)), V.row(F(firstFace, 1)), V.row(F(firstFace, 2)));

  // For each vertex in the very first triangle, we set its newVerts
  // to its oldVerts
  for (int i = 0; i < 3; i++) {
    newVerts.row(F(firstFace, i)) = V.row(F(firstFace, i));
  }

  // Flood fill
  set<int> seenFaces;
  queue<int> toProcess;
  toProcess.push(firstFace);
  seenFaces.insert(firstFace);
  while (seenFaces.size() != F.rows()) {
    int currFace = toProcess.front();
    toProcess.pop();
    // for every one to nextFace, we check whether or not we have computed the gradient
    // if we have, then we don't do anything, otherwise we compute
    // the gradient and put on the new vertices
    for (int i = 0; i < TT.row(currFace).size(); i++) {
      int nextFace = TT(currFace, i);
      if (seenFaces.find(nextFace) == seenFaces.end() && nextFace != -1) {
        gradients.row(nextFace) = getNextGradient(
            gradients.row(currFace), reciprocalDual.row(currFace), reciprocalDual.row(nextFace));
        // Find a matching vertex between this and the next
        int indexMatching = -1;
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            if (F(currFace, i) == F(nextFace, j)) {
              indexMatching = j;
            }
          }
        }
        if (indexMatching == -1) {
          cerr << "Augh! We did not find a matching index!" << endl;
          throw new exception();
        }
        // use matching vertex to compute the other two points
        // (yes, one vertex from the previous triangle will be overwritten)
        // (but this is okay due to the gradients giving the same result at shared points)
        for (int i = 0; i < 3; i++) {
          if (i != indexMatching) {
            newVerts.row(F(nextFace, i)) = applyGradient(
                gradients.row(nextFace), newVerts.row(F(nextFace, indexMatching)),
                V.row(F(nextFace, i)));
          }
        }
        toProcess.push(nextFace);
        seenFaces.insert(nextFace);
      }
    }
  }
  return newVerts;
}
