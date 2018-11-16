#include "HandMesh.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/triangle_triangle_adjacency.h>

void makeEdges(Eigen::MatrixXi F, Eigen::MatrixXd BC,
               igl::opengl::glfw::Viewer &viewer) {
  Eigen::MatrixXi TT;
  Eigen::MatrixXi TTi;
  // Go through and add edges between each BC and the BC it is adjacent to
  igl::triangle_triangle_adjacency(F, TT, TTi);

  cout << "the outer size of TT is " << TT.rows() << endl;
  cout << "The inner size of TT is " << TT.cols() << endl;
  cout << "the rows of barycenter is " << BC.rows() << endl;

  Eigen::RowVector3d color;
  color << 200, 0, 0;
  // Isn't | TT | == | BC |, the number of the triangles?
  for (int i = 0; i < TT.rows(); i++) {
    auto barycenteri = BC.row(i);
    for (int j = 0; j < TT.row(i).size(); j++) {
      // The index of the triangle adjge j of triangle i
      int adjIndex = TT(i, j);
      if (adjIndex == -1) {
        continue;
      }
      auto barycenterk = BC.row(adjIndex);
      viewer.data().add_edges(barycenteri, barycenterk, color);
    }
  }
}

void makeBary(Eigen::MatrixXd V, Eigen::MatrixXi F,
              igl::opengl::glfw::Viewer &viewer) {
  Eigen::MatrixXd BC;
  igl::barycenter(V, F, BC);
  makeEdges(F, BC, viewer);
}
