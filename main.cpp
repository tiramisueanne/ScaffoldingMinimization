#include <Eigen/Dense>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <igl/segment_segment_intersect.h>
#include <iostream>

#include "Circumcenter.h"
#include "MakeBary.h"
#include "HandMesh.h"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]) {
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  // initHandMesh(V, F);
  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      F(i, j) -= 1;
    }
  }
  igl::opengl::glfw::Viewer viewer;
  igl::readOBJ("./rec.obj", V, F);
  // Okay so get the "segments" row by row
  MatrixXd circ = getCircumcenters(V, F);
  viewer.data().set_mesh(V, F);
  makeEdges(F, circ, viewer);
  viewer.launch();
}
