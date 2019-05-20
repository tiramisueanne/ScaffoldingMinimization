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
    double tenPercent;
    int count = 0;
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
        tenPercent = 0.1 * sum;
        count++;
    } while (res != 1 && fabs(res + sum) < tenPercent && facesLeft > 0);
    cout << "We removed " << count << "nodes from this figure" << endl;
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
    int countStop = 15;
    while (fabs(res + sum) > pow(10, -6) && count < countStop) {
        res = qs.updateWeights();
        cout << "Got past updating weights with res " << res << " and sum "
             << sum << endl;
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

void removeCheapestNode(MatrixXd& V, MatrixXi& F) {
    QuadraticSolver qs(V, F);
    qs.removeSmallestNode();
}

int main(int argc, char* argv[]) {
    // Reading in the primal V, F files
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Calculation type;
    parseInput(argc, argv, V, F, type);

    igl::opengl::glfw::Viewer viewer;
    cout << "The type we are going to switch on is" << type << endl;
    int hasFace;
    switch (type) {
        case QUADRATIC:
            quadraticProgrammingUpdateStructure(V, F);
            break;
        case DUAL:
            createDuals(V, F);
            break;
        case SCAFFOLDING: {
            quadraticProgrammingUpdateStructure(V, F);
            hasFace = findNextSpot(V, F);
            if (!hasFace) {
                cout << "The structure is stable at all levels!" << endl;
                return 0;
            }
            break;
        }
        case SHOWDUAL:
            break;
        case REMOVE:
            cout << "We just got into removeCheapest!" << endl;
            removeCheapestNode(V, F);
            break;
        default:
            cout << "We just showed it from laplacian" << endl;
            QuadraticSolver qs(V, F);
            break;
    }

    viewer.data().set_mesh(V, F);
    // Compile the weights file
    viewer.launch();
}
