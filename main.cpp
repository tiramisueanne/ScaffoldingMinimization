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

#include "GetNewPrimal.h"
#include "HandMadeMeshes/MediumMesh.h"
#include "HandMadeMeshes/SmallMesh.h"
#include "HandMadeMeshes/TinyMesh.h"
#include "QuadraticSolver.h"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]) {
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
        } else if (name == "small") {
            initBiggerMesh(V, F);
            isHand = true;
        } else if (name == "medium") {
            evenBiggerMesh(V, F);
            isHand = true;
        } else {
            isHand = false;
            igl::readOBJ(argv[1], V, F);
        }
    }
    // otherwise, use our tiny hand-made mesh
    else {
        initHandMesh(V, F);
        isHand = true;
    }
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    // Compile the weights file
    viewer.launch();
}

void quadraticProgrammingUpdateStructure(MatrixXd& V, MatrixXi& F) {
    QuadraticSolver qs(V, F);
    double res = 1;
    double sum = qs.getTotalForce();
    int count = 0;
    int countStop = 15;
    while (fabs(res + sum) > pow(10, -6) && count < countStop) {
        res = qs.updateWeights();
        double updateSuccess = qs.updateVertices();
        count++;
    }
    cout << "Iterated on this shape " << count << " times " << endl;
    V = qs.V;
}

void createDuals(MatrixXd& V, MatrixXi& F) {}
