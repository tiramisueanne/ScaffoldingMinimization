#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <igl/segment_segment_intersect.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Dense>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <queue>
#include <set>

#include "Circumcenter.h"
#include "GetNewDual.h"
#include "GetNewPrimal.h"
#include "HandMesh.h"
#include "MakeBary.h"
#include "Weights.h"

using namespace std;
using namespace Eigen;

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
    MatrixXd dualVerts = getNewDual(V, F, leastSquaresResult(V, F));
    MatrixXd newVerts = getNewPrimal(dualVerts, V, F);

    // utilize libigl's viewer
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(newVerts, F);

    // Compile the weights file
    viewer.launch();
}
