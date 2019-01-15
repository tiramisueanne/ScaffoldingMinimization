#include <Eigen/Dense>
#include <iostream>
using namespace std;
using namespace Eigen;

// @param a, b, c are points of a face
RowVector3d getGradient(RowVector3d a, RowVector3d b, RowVector3d c) {
  // vector from a to b
  RowVector3d ab = b - a;
  // vector from a to c
  RowVector3d ac = c - a;


  RowVector3d cross = ab.cross(ac);

  // Fill in the gradient values
  // TODO: why does this work?
  RowVector3d grad;
  grad(0) = cross(0) / -1 * cross(2);
  grad(1) = cross(1) / -1 * cross(2);
  grad(2) = 0;
  return grad;
}

// @param The gradient of the primal, along with two points from the dual
RowVector3d getNextGradient(RowVector3d currGrad, RowVector3d dualCurr,
                         RowVector3d dualNext) {
  return currGrad + (dualNext - dualCurr);
}

// @param the gradient of the face, plus point a from face which has a precomputed z value and point b which requires a z value
RowVector3d applyGradient(RowVector3d gradient, RowVector3d a, RowVector3d b) {
  RowVector3d ab = b - a;
  b.z() = ab.dot(gradient) + a.z();
  return b;
}
