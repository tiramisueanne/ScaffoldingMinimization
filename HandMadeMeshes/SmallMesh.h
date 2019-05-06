#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>

void initBiggerMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    V = Eigen::MatrixXd(8, 3);
    F = Eigen::MatrixXi(8, 3);
    V << -1.0, 0.0, 0.0,
        -0.5, 0.5, 0.0,
        -0.5, 0.0, 0.0,
        -0.5, -0.5, 0.0,
        0.5, 0.5, 0.0,
        0.5, 0.0, 0.0,
        0.5, -0.5, 0.0,
        1.0, 0.0, 0.0;
    F << 0, 2, 1,
        0, 3, 2,
        2, 4, 1,
        2, 5, 4,
        3, 5, 2,
        3, 6, 5,
        4, 5, 7,
        5, 6, 7;
}
