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
        string name = argv[1];
        if (name == "tiny") {
            initHandMesh(V, F);
            isHand = true;
        }
        else if (name == "small") {
            initBiggerMesh(V, F);
            isHand = true;
        }
        else if (name == "medium") {
            evenBiggerMesh(V, F);
            isHand = true;
        }
        else {
            isHand = false;
            igl::readOBJ(argv[1], V, F);
        }
    }
    // otherwise, use our tiny hand-made mesh
    else {
        initHandMesh(V, F);
        isHand = true;
    }

// This is the flood fill of finding the gradient and creating newVertices
#if !defined(JUST_SHOW)
    QuadraticSolver qs(V, F, isHand);
    double res = 1;
    double sum = qs.getTotalForce();
    int count = 0;
#ifdef DEBUG
    int countStop = 1;
#endif
#ifndef DEBUG
    int countStop = 15;
#endif
    while (fabs(res + sum) > pow(10, -6) && count < countStop) {
        res = qs.updateWeights();
        double updateSuccess = qs.updateVertices();
        count++;
    }
    cout << "Iterated on this shape " << count << " times" << endl;
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

    // Compile the weights file
    viewer.launch();
}
