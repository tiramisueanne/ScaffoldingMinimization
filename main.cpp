#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
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
#define SHOW_POISSON
// #define DUALS
// #define CHECKDUALS
// #define CHECKBIGMESH
#define CHECKBIGGER
// #define JUST_SHOW

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
        initBiggerMesh(V, F);
        isHand = true;
#endif
    }

// This is the flood fill of finding the gradient and creating newVertices
#if !defined(JUST_SHOW)
    QuadraticSolver qs(V, F, isHand);
    double success = 1;
    double sum = qs.getTotalForce();
    int count = 0;
#ifdef DEBUG
    int countStop = 1;
#endif
#ifndef DEBUG
    int countStop = 30;
#endif
#if !defined(DUALS)
    cout << "The fabs value was " << fabs(success + sum) << endl;
    while (success == 1 && count < countStop) {
        success = qs.updateWeights();
        double updateSuccess = qs.updateVertices();
        count++;
    }
    cout << "Iterated on this shape " << count << " times" << endl;
#endif
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
#if !defined(DUALS_L) && !defined(JUST_SHOW)
    MatrixXd &toUse = qs.V;
#endif

    // utilize libigl's viewer
    igl::opengl::glfw::Viewer viewer;
#ifndef JUST_SHOW
    viewer.data().set_mesh(toUse, qs.F);
#endif
#ifdef JUST_SHOW
    viewer.data().set_mesh(V, F);
#endif
// Place dual points if we have them
#ifdef DUALS_L
    viewer.data().add_points(dualVerts, RowVector3d(10, 10, 100));
#endif
#ifdef CHECK_BIG
    MatrixXd interesting(3, 3);
    interesting.row(0) = V.row(9);
    interesting.row(1) = V.row(25);
    interesting.row(2) = V.row(37);
    viewer.data().add_points(interesting, RowVector3d(10, 10, 100));
#endif

    // Compile the weights file
    viewer.launch();
}
