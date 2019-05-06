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
#include "QuadraticSolver.h"
#include "ParseInput.h"

using namespace std;
using namespace Eigen;

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

void createDuals(MatrixXd& V, MatrixXi& F) {
    QuadraticSolver qs(V, F);
    qs.updateWeights();
    MatrixXd dualVerts = qs.getNewDual();
    V = getNewPrimal(dualVerts, V, F);
}


int main(int argc, char* argv[]) {
    // Reading in the primal V, F files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Calculation type;
    parseInput(argc, argv, V, F, type);

    if(type == QUADRATIC) quadraticProgrammingUpdateStructure(V, F);
    if(type == DUAL) createDuals(V, F);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    // Compile the weights file
    viewer.launch();
}

