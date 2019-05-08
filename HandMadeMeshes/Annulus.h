#pragma once
#include <igl/opengl/glfw/Viewer.h>
#include <igl/triangle/triangulate.h>
#include <iostream>

using namespace std;
using namespace Eigen;

void createAnnulus(MatrixXd &V, MatrixXi &F) {
    cout << "Was good in createAnnulus" << endl;
    int r1 = 20;
    MatrixXd verts(r1 * 4, 2);
    MatrixXd edges(r1 * 4, 2);
    MatrixXd holes;
    // the upper half, clockwise
    for (int x = -r1; x < r1; x++) {
        verts(x + r1, 0) = x;
        verts(x + r1, 1) = sqrt(r1 * r1 - x * x);
    }
    // The lower half, clockwise
    for (int x = r1; x > -r1; x--) {
        verts(r1 * 2 + (r1 - x), 0) = x;
        verts(r1 * 2 + (r1 - x), 1) = -sqrt(r1 * r1 - x * x);
    }
    for (int e = 0; e < r1 * 4; e++) {
        edges(e, 0) = e;
        edges(e, 1) = (e + 1) % (r1 * 4);
    }
    cout << "The verts were " << verts << endl;
    cout << "The edges were " << edges << endl;
    igl::triangle::triangulate(verts, edges, holes, "q", V, F);
}