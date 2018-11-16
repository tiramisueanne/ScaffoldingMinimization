#include <Eigen/Dense>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <igl/segment_segment_intersect.h>
#include <iostream>

// #define DEBUG

using namespace std;
using namespace Eigen;

// This would just be the z value of a regular cross product in 3d
double cross2d(const RowVector3d &vec1, const RowVector3d &vec2) {
  return (vec1.x() * vec2.y() - vec2.x() * vec1.y());
}

double length(const RowVector3d &row) { return sqrt(row.dot(row)); }

RowVector3d intersection(RowVector3d p, RowVector3d r, RowVector3d q,
                         RowVector3d s) {
  double u = cross2d((q - p), r) / cross2d(r, s);
  double t = cross2d((q - p), s) / cross2d(r, s);
  RowVector3d inter1 = q + u * s;
  RowVector3d inter2 = p + r * t;
  if (length(inter1 - inter2) > 0.0001) {
    cerr << "The intersections are not the same!" << endl;
    cerr << "Inter 1 " << inter1 << endl;
    cerr << "Inter 2" << inter2 << endl;
  }
  return inter1;
}



MatrixXd getCircumcenters(const MatrixXd &V, const MatrixXi &F) {
  MatrixXd circ(F.rows(), 3);
  for (int i = 0; i < F.rows(); i++) {
    RowVector3i face = F.row(i);
    // Find the slopes (m) and midpoints (r)
    Eigen::Matrix<double, 1, 3> m1 = V.row(face(1)) - V.row(face(0));
    Eigen::Matrix<double, 1, 3> r1 =
        (V.row(face(1)) + V.row(face(0))) / 2;
    // I think it is this line right here
    Eigen::Matrix<double, 1, 3> m2 = V.row(face(2)) - V.row(face(0));
    Eigen::Matrix<double, 1, 3> r2 =
        (V.row(face(2)) + V.row(face(0))) / 2;
    Eigen::Matrix<double, 1, 3> straightUp;
    straightUp << 0, 0, 1.0;
    Matrix<double, 1, 3> ortho1;
    Matrix<double, 1, 3> ortho2;
    Matrix<double, 1, 3> checkP;
    igl::cross(m1, straightUp, ortho1);
    igl::cross(m2, straightUp, ortho2);
    igl::cross(ortho1, ortho2, checkP);
    double isZero = sqrt(checkP.dot(checkP));
    if (isZero < 0.00000001) {
      cerr << "Is Zero was less than threshold" << endl;
      cerr << "Thus, the two sides of the triangle were colinear" << endl;
      throw new exception();
    }
    circ.row(i) = intersection(r1, ortho1, r2, ortho2);
#ifdef DEBUG
    for (int i = 0; i < 3; i++) {
      cout << V.row(face(i)) << ", ";
    }
    cout << "The circumcenter" << circ.row(i) << endl;
#endif
  }
  return circ;
}
