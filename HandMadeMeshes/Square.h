#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <igl/triangle/triangulate.h>
#include <iostream>

using namespace std;
using namespace Eigen;

void createSquare(MatrixXd &V, MatrixXi &F) {
    int r1 = 4;
    MatrixXd verts(4, 2);
    MatrixXd edges(4, 2);
    MatrixXd holes;
    verts << r1, -r1,
        r1, r1,
        -r1, r1,
        -r1, -r1;
    edges << 0, 1,
        1, 2,
        2, 3,
        3, 0;
    MatrixXd twoDV;
    igl::triangle::triangulate(verts, edges, holes, "a0.8", twoDV, F);
    V = MatrixXd(twoDV.rows(), 3);
    V.block(0, 0, twoDV.rows(), 2) = twoDV;
    V.block(0, 2, V.rows(), 1) = VectorXd::Constant(V.rows(), 0);
}
