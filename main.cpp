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

#include "EvenBigger.h"
#include "GetNewPrimal.h"
#include "HandMesh.h"
#include "InitBiggerMesh.h"
#include "QuadraticSolver.h"

using namespace std;
using namespace Eigen;
// #define DEBUG
// #define SHOW_POISSON
// #define DUALS
// #define CHECKDUALS
#define CHECKBIGMESH
#define CHECKBIGGER

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
#ifndef CHECKBIGGER
        initHandMesh(V, F);
        isHand = true;
#endif
#ifdef CHECKBIGGER
        evenBiggerMesh(V, F);
        isHand = true;
#endif
    }

    // This is the flood fill of finding the gradient and creating newVertices

    QuadraticSolver qs(V, F, isHand);
    double success = 0;
    double sum = qs.getTotalForce();
    int count = 0;
#ifdef DEBUG
    int countStop = 1;
#endif
#ifndef DEBUG
    int countStop = 1;
#endif
#if !defined(DUALS)
    while (fabs(success + sum) > 0.000001 && count < countStop) {
        success = 0;
        success += qs.updateWeights();
        cout << "The success value of succ is " << success << endl;
        cout << "The success value of succ and sum is " << sum + success << endl;
        double updateSuccess = qs.updateVertices();
        cout << "the success value of updating was " << updateSuccess << endl;
        count++;
    }
    cout << "Iterated on this shape " << count << " times" << endl;
#endif

#ifdef DUALS_L
    qs.updateWeights();
    MatrixXd dualVerts = qs.getNewDual();
#ifdef CHECKDUALS
    if (!qs.checkDual(dualVerts)) {
        cerr << "dual is messed up" << endl;
        throw new exception();
    }
#endif
    MatrixXd newVerts = getNewPrimal(dualVerts, V, F);
#endif

#if defined(DUALS) && defined(DEBUG)
    cout << "The new verts are " << newVerts << endl;
    cout << "The first one is" << newVerts.row(0);
#endif
#ifdef DUALS
    MatrixXd &toUse = newVerts;
#endif
#ifndef DUALS_L
    MatrixXd &toUse = qs.V;
#endif

    // utilize libigl's viewer
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(toUse, qs.F);

// Place dual points if we have them
#ifdef DUALS_L
    viewer.data().add_points(dualVerts, RowVector3d(10, 10, 100));
#endif

    // Compile the weights file
    viewer.launch();
}
