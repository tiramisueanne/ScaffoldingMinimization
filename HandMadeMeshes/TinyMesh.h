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

    F << 0, 1, 5, /* a */

        0, 4, 1, /* b */

        0, 5, 2,

        0, 2, 3,

        0, 3, 4;
}
