#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <igl/segment_segment_intersect.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Dense>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <map>
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
    int count = 0;
    while (abs(succ + sum) > 0.001 && count < 30) {
        succ = 0;
        succ += qs.updateWeights();
        qs.updateVertices();
        count++;
    }
    if(count == 30) {
        cout << "Iterated on this shape 30 times" << endl;
    }

#ifdef DUALS
    MatrixXd dualVerts = getNewDual(V, F, qs.getWeights());
    MatrixXd newVerts = getNewPrimal(dualVerts, V, F);
#endif

#if defined(DUALS) && defined(DEBUG)
    cout << "The new verts are " << newVerts << endl;
    cout << "The first one is" << newVerts.row(0);
#endif

    // utilize libigl's viewer
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(qs.V, qs.F);


    // Place dual points if we have them
#ifdef DUALS
    viewer.data().add_points(dualVerts, RowVector3d(10, 10, 100));
#endif

    // Compile the weights file
    viewer.launch();
}
