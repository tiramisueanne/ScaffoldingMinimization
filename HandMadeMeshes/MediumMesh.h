#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>

using namespace std;
using namespace Eigen;

void evenBiggerMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    V = MatrixXd(36, 3);
    F = MatrixXi(50, 3);
    int rowSize = 6;
    for (int i = 0; i < rowSize; i++) {
        for (int j = 0; j < rowSize; j++) {
            V.row(i * rowSize + j) << i, j, 0;
        }
    }
    for (int i = 0; i < rowSize - 1; i++) {
        for (int j = 0; j < rowSize - 1; j++) {
            int firstVert = rowSize * (i + 1) + j;
            int secondVert = rowSize * i + j + 1;
            int thirdVert = rowSize * i + j;
            int fourthVert = rowSize * (i + 1) + j + 1;
            // Upper half
            F.row(2 * (i * (rowSize - 1) + j)) << firstVert, secondVert, thirdVert;
            // Lower Half
            F.row(2 * (i * (rowSize - 1) + j) + 1) << firstVert, fourthVert,
                secondVert;
        }
    }
}
