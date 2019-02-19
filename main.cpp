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

#include "Circumcenter.h"
#include "GetNewDual.h"
#include "GetNewPrimal.h"
#include "HandMesh.h"
#include "MakeBary.h"
#include "Weights.h"

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
    // MatrixXd dualVerts = getCircumcenters(V, F);
    map<pair<int, int>, int> edgeWeights = getWeights(V, F);
    #ifdef DEBUG
    for(const auto edge : edgeWeights) {
        cout << "pair: " << edge.first.first << " , " << edge.first.second << " and the value " << edge.second << endl;
    }
    #endif
    MatrixXd dualVerts = getNewDual(V, F, edgeWeights);
    cout << dualVerts << endl;
    MatrixXd newVerts = getNewPrimal(dualVerts, V, F);

    // utilize libigl's viewer
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(newVerts, F);

    // Compile the weights file
    viewer.launch();
}
