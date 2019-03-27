#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>

using namespace std;


void initBiggerMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    V = MatrixXd(16, 3);
    F = MatrixXi(18, 3);
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
        }
    }
}
