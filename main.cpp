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
// #define DEBUG
// #define SHOW_POISSON

int main(int argc, char *argv[]) {
    // Reading in the primal V, F files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    // If a filename is given, then read in that
    bool isHand = false;
    if (argc > 1) {
        igl::readOBJ(argv[1], V, F);

    }
    // otherwise, use our hand-made mesh
    else {
        initHandMesh(V, F);
        isHand = true;
    }

    // This is the flood fill of finding the gradient and creating newVertices

    QuadraticSolver qs(V, F, isHand);
    double succ = 0;
    double sum = qs.getTotalForce();
    sum *= sum;
    int count = 0;
    cout << "The original was\n " << qs.V << " and\n " << qs.F << endl;
    #ifdef DEBUG
    int countStop = 1;
    #endif
    #ifndef DEBUG
    int countStop = 30;
    #endif
    #ifndef SHOW_POISSON
    while (abs(succ + sum) > 0.00000001 && count < countStop) {
        succ = 0;
        succ += qs.updateWeights();
        cout << "The success value of succ is " << succ << endl;
        qs.updateVertices();
        count++;
    }
    cout << "Iterated on this shape "<< count << " times" << endl;
    #endif

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
