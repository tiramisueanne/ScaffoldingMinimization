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
#include "ParseInput.h"
#include "QuadraticSolver.h"

using namespace std;
using namespace Eigen;
#define DEBUG

int findNextSpot(MatrixXd& V, MatrixXi& F) {
    QuadraticSolver qs(V, F);
    double res = 1;
    double sum;
    int facesLeft = 0;
    do {
        facesLeft = qs.removeSmallestNode();
        if (facesLeft > 0) {
            res = qs.updateWeights();
            sum = qs.getTotalForce();
#ifdef DEBUG
            cout << "The residual is " << res << endl;
            cout << "The totalForce is " << sum << endl;
#endif
        }
    } while (fabs(res + sum) < pow(10, -9) && facesLeft > 0);
#ifdef DEBUG
    cout << "the remaining faces are " << F << endl;
#endif
    return facesLeft;
}

void quadraticProgrammingUpdateStructure(MatrixXd& V, MatrixXi& F) {
    QuadraticSolver qs(V, F);
    double res = 1;
    double sum = qs.getTotalForce();
    int count = 0;
    int countStop = 1;
    while (fabs(res + sum) > pow(10, -6) && count < countStop) {
        res = qs.updateWeights();
        cout << "Got past updating weights" << endl;
        double updateSuccess = qs.updateVertices();
        sum = qs.getTotalForce();
        count++;
    }
    cout << "Iterated on this shape " << count << " times " << endl;
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

    igl::opengl::glfw::Viewer viewer;

    switch (type) {
        case QUADRATIC:
            quadraticProgrammingUpdateStructure(V, F);
            break;
        case DUAL:
            createDuals(V, F);
            break;
        case SCAFFOLDING: {
            int hasFace = findNextSpot(V, F);
            if (!hasFace) {
                cout << "The structure is stable at all levels!" << endl;
                return 0;
            }
            break;
        }
        case SHOWDUAL: {
        }
        default:
            QuadraticSolver qs(V, F);
            break;
    }

    viewer.data().set_mesh(V, F);
    // Compile the weights file
    viewer.launch();
}
