#pragma once
#include <Eigen/Dense>
#include <eigen-quadprog/QuadProg.h>
#include <map>

#include "Indexer.h"

using namespace std;
using namespace Eigen;

class QuadraticSolver {
   public:
    QuadraticSolver(MatrixXd &_V, MatrixXi &_F, bool isHand = false)
        : V(_V), F(_F) {
        unsupportedNodes = getUnsupportedNodes(F, V);
        cout << "Got the unsupported nodes!" << endl;
        edges = allEdges();
        cout << "Got all the edges!" << endl;
        indr = Indexer(V.rows(), unsupportedNodes, edges);
        cout << "Indexer all done!" << endl;
        createLaplacian(isHand);
        cout << "Laplacian also complete" << endl;
    };

    // This is to calculate weights from our estimated forces on the structure
    double updateWeights();

    // This is to update the vertices with weights given to edges fixed
    double updateVertices();

    // This is how we estimate the forces
    MatrixXd getForces();

    Eigen::MatrixXd getForceAreas();
    Eigen::VectorXd getForceDensities();

    const map<pair<int, int>, double> &getWeights() { return weightMap; }

    // This is to calculate all of the edges in the figure
    // Used to find all internal nodes
    set<int> getUnsupportedNodes() { return getUnsupportedNodes(F, V); }

    static set<int> getUnsupportedNodes(const MatrixXi &_F, const MatrixXd &_V);

    // Used to get all of the edges we need
    set<pair<int, int>> allEdges();

    // so that we don't accidentally have all 0's
    void createLaplacian(bool isHand);

    void bumpInternalNodes();

    double getTotalForce() {
        double total = 0;
        forces = getForces();
        for (int i = 0; i < forces.cols(); i++) {
            total += forces(0, i) * forces(0 , i);
        }
        return total;
    }

    void moveVecIntoV();

    MatrixXd getNewDual();
    bool checkDual(MatrixXd &dualVerts);

    void unsupportANode();

    MatrixXd &V;
    VectorXd vec;
    const MatrixXi &F;
    set<pair<int, int>> edges;
    VectorXd forces;
    VectorXd weights;
    map<pair<int, int>, double> weightMap;
    Indexer indr;
    set<int> unsupportedNodes;
};
