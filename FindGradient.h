#include <Eigen/Dense>
#include <iostream>
using namespace std;
using namespace Eigen;

RowVector3d getGradient(RowVector3d a, RowVector3d b, RowVector3d c) {
  RowVector3d ab = b - a;
  RowVector3d ac = c - a;
  RowVector3d cross = ab.cross(ac);

  // Fill in the gradient values
  RowVector3d grad;
  grad(0) = cross(0) / -1 * cross(2);
  grad(1) = cross(1) / -1 * cross(2);
  grad(2) = 0;
  return grad;
}

RowVector3d getNextGradient(RowVector3d currGrad, RowVector3d dualCurr,
                         RowVector3d dualNext) {
  return currGrad + (dualNext - dualCurr);
}

RowVector3d applyGradient(RowVector3d gradient, RowVector3d a, RowVector3d b) {
  RowVector3d ab = b - a;
  b.z() = ab.dot(gradient) + a.z();
  return b;
}
