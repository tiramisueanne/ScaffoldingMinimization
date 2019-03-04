#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <igl/segment_segment_intersect.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Dense>

#include <map>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <queue>
#include <set>

#include "GetNewDual.h"
#include "GetNewPrimal.h"
#include "HandMesh.h"
#include "QuadraticSolver.h"

using namespace std;
using namespace Eigen;
#define DEBUG

int main(int argc, char *argv[]) {
    // Reading in the primal V, F files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    // If a filename is given, then read in that
    if (argc > 1) {
        igl::readOBJ(argv[1], V, F);
    }
    // otherwise, use our hand-made mesh
    else {
        initHandMesh(V, F);
    }

    // This is the flood fill of finding the gradient and creating newVertices
    QuadraticSolver qs(V, F);
    double succ = 0;
    double sum = qs.getTotalForce();
    sum *= sum;
    while (abs(succ + sum) > 0.01) {
        succ = 0;
        succ += qs.updateWeights();
        cout << "The succ is " << succ << endl;
        qs.updateVertices();
    }

    MatrixXd dualVerts = getNewDual(V, F, qs.getWeights());
    MatrixXd newVerts = getNewPrimal(dualVerts, V, F);
    #ifdef DEBUG
    cout << "The new verts are " << newVerts << endl;
    cout << "The first one is" << newVerts.row(0);
    #endif

    // utilize libigl's viewer
    igl::opengl::glfw::Viewer viewer;
    // viewer.data().set_mesh(newVerts, F);
    viewer.data().add_points(dualVerts, RowVector3d(10, 10, 100));

    // Compile the weights file
    viewer.launch();
}
