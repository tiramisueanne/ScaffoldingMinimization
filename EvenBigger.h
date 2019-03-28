#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>

using namespace std;
using namespace Eigen;

void evenBiggerMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    V = MatrixXd(16, 3);
    F = MatrixXi(18, 3);
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            V.row(i*4 + j) << i, j, 0;
        }
    }
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            int firstVert = 4 * (i + 1) + j;
            int secondVert = 4 * i + j + 1;
            int thirdVert = 4 * i + j;
            int fourthVert = 4 * (i + 1) + j + 1;
            // Upper half
            F.row(2 * (i*3 + j)) << firstVert, secondVert, thirdVert;
            // Lower Half
            F.row(2 * (i*3 + j) + 1) << firstVert, fourthVert, secondVert;
        }
    }
}
