#include <Eigen/Dense>
#include <igl/HalfEdgeIterator.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <igl/segment_segment_intersect.h>
#include <igl/triangle_triangle_adjacency.h>
#include <QuadProg++/Array.hh>

#include <iostream>
#include <math.h>
#include <queue>
#include <set>
#include <stdlib.h>
#include <string.h>

#include "Circumcenter.h"
#include "HandMesh.h"
#include "MakeBary.h"
#include "FloodFillVertices.h"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]) {
  // Reading in the primal V, F files
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  // If a filename is given, then read in that
  if(argc > 1) {
    igl::readOBJ(argv[1], V, F);
  }
  // otherwise, use our hand-made mesh
  else {
    initHandMesh(V, F);
  }

  // circumcenters of the graph
  // TODO: remove this in favor of edge weighting finding of duals
  MatrixXd dualVerts = getCircumcenters(V, F);
  if (dualVerts.rows() != F.rows()) {
    cerr << "We have an error here" << endl;
  }

  // This is the flood fill of finding the gradient and creating newVertices
  MatrixXd newVerts = getNewVertices(dualVerts, V, F);

  // utilize libigl's viewer
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(newVerts, F);
  viewer.launch();
}
