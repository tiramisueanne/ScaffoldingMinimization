#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>

void initHandMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
  // Inline mesh of a cube
  V = Eigen::MatrixXd(6, 3);
  F = Eigen::MatrixXi(5, 3);
  V << 0.0, 0.0, 0.0,

    0.0, 1.0, 0.0,

    -1.0, -1.0, 0.0,

    1.0, -1.0, 0.0,

    1.0, 0.0, 0.0,

    -1.0, 0.0, 0.0;

  F << 1, 2, 6, /* a */

    1, 5, 2, /* b */

    1, 6, 3,

    1, 3, 4,

    1, 4, 5;

}
