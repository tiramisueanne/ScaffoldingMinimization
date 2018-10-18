#include <igl/opengl/glfw/Viewer.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/readOBJ.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

  // Inline mesh of a cube
  Eigen::MatrixXd V = (Eigen::MatrixXd(6, 3) << 0.0, 0.0, 0.0,

                             0.0, 1.0, 0.0,

                             -1.0, -1.0, 0.0,

                             1.0, -1.0, 0.0,

                             1.0, 0.0, 0.0,

                             -1.0, 0.0, 0.0)
                                .finished();
  Eigen::MatrixXi F = (Eigen::MatrixXi(5, 3) << 1, 2, 6, /* a */

                             1, 5, 2, /* b */

                             1, 6, 3,

                             1, 3, 4,

                             1, 4, 5)
                                .finished()
    .array() - 1;

  igl::readOBJ("./rec.obj", V, F);

  Eigen::MatrixXd BC;
  igl::barycenter(V, F, BC);
  Eigen::MatrixXi TT;
  Eigen::MatrixXi TTi;
  // Go through and add edges between each BC and the BC it is adjacent to
  igl::triangle_triangle_adjacency(F, TT, TTi);

  cout << "the outer size of TT is " << TT.rows() << endl;
  cout << "The inner size of TT is " << TT.cols() << endl;
  cout << "the rows of barycenter is " << BC.rows() << endl;


  igl::opengl::glfw::Viewer viewer;
  Eigen::RowVector3d color;
  color << 0,0,200;
  // Isn't | TT | == | BC |, the number of the triangles?
  for(int i = 0; i < TT.rows(); i++) {
    auto barycenteri = BC.row(i);
    for(int j = 0; j < TT.row(i).size(); j++) {
      // The index of the triangle adjacent to edge j of triangle i
      int adjIndex = TT(i, j);
      if (adjIndex == -1) {
        continue;
      }
      auto barycenterk = BC.row(adjIndex);
      viewer.data().add_edges(barycenteri, barycenterk, color);
    }
  }

  // Plot the mesh
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.launch();
}
