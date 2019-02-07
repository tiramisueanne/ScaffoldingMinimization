#pragma once

#include <Eigen/Dense>
#include <iostream>

#include <QuadProg++/QuadProg++.hh>

using namespace std;
using namespace Eigen;
namespace qp = quadprogpp;

// A method for getting each "real" weight for the vertices
// TODO: figure out how to read in the "real" weights, as we are just using 2
// for each one right now
qp::Vector<double> getForces(MatrixXd V, MatrixXi F) {
  // Since we currently don't have any information on weights, just
  // return a random number for every vertex of the primal
  return qp::Vector<double>(2.0, int(V.innerSize()));
}

// V and Weights will have the same number of entries
// F just allows for us to collect faces
VectorXd leastSquaresResult(const MatrixXd &V, const MatrixXi &F) {
  qp::Vector<double> weights = getForces(V, F);
  unsigned int rowSize = int(V.innerSize());
  int ZERO = 0;
  // This will be the array of differences of the z values
  // each row consists of a single vertex
  // initialize values to zero
  // innerSize = number of rows
  qp::Matrix<double> zDiff(ZERO, rowSize, rowSize);

  // This matrix is two rows for every vertex: one is diff in x for
  // each edge, and one is diff in y
  // This will be constrained to zero.
  qp::Matrix<double> xyDiff(ZERO, rowSize * 2, rowSize);

  // This is the identity matrix, to help us constrain the weights
  // to be positive
  qp::Matrix<double> identity(ZERO, rowSize, rowSize);
  const double *ONE = new double(1.0);
  for (unsigned int i = 0; i <= rowSize; i++) {
    identity.set(ONE, i, i);
  }

  // Fill up zDiff through going through each triangle and filling in
  for (int i = 0; i < F.rows(); i++) {
    RowVector3i currFace = F.row(i);
    Matrix3d currPoints;
    for (int l = 0; l < 3; l++) {
      if (currFace(l) >= rowSize || currFace(l) < 0) {
        cerr << "The vertex is outside of the realm of the rows" << endl;
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
        int other1 =(j + k) % 3;
        int index = currFace(other1);
        // TODO: no new doubles
        const double *zDiff1 =
            new double(currPoints(j, 2) - currPoints(other1, 2));
        const double *xDiff1 =
            new double(currPoints(j, 0) - currPoints(other1, 0));
        const double *yDiff1 =
            new double(currPoints(j, 1) - currPoints(other1, 1));
        zDiff[currIndex][index] = *zDiff1;
        xyDiff[currIndex*2][index] = *xDiff1;
        xyDiff[currIndex*2 + 1][index] = *yDiff1;
      }
    }
  }
  qp::Vector<double> x(2, rowSize);
  qp::Vector<double> justOnes(1, rowSize);
  // The vector we add to the result of the constraint
  qp::Vector<double> justZerosForXY(0.0, 2 * rowSize);
  qp::Vector<double> justZerosForPos(0.0, rowSize);

  qp::Vector<double> linearComponent(dot_prod(weights, zDiff));
  // To construct the positive definite zDiff matrix, we must
  // multiply it with itself
  zDiff *= zDiff;
  double success =
      qp::solve_quadprog(zDiff, linearComponent, xyDiff, justZerosForXY,
                         identity, justZerosForPos, x);
  VectorXd toReturn(rowSize);
  for (int i = 0; i < rowSize; i++) {
    toReturn(i) = x[i];
  }
  return toReturn;
}
