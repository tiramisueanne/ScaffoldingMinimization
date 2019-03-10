#pragma once
#include <Eigen/Dense>
#include <QuadProg++/QuadProg++.hh>
#include <map>

#include "Indexer.h"

using namespace std;
using namespace Eigen;
namespace qp = quadprogpp;

class QuadraticSolver {
   public:
    QuadraticSolver(MatrixXd &_V, MatrixXi &_F, bool isHand = false) : V(_V), F(_F) {
        internalNodes = getInternalNodes();
        edges = allEdges();
        indr = Indexer(V.rows(), internalNodes, edges);
        createLaplacian(isHand);
    };

    // This is to calculate weights from our estimated forces on the structure
    double updateWeights();

    // This is to update the vertices with weights given to edges fixed
    double updateVertices();

    // This is how we estimate the forces
    qp::Matrix<double> getForces();

    const map<pair<int, int>, double>& getWeights() {
        return weightMap;
    }

    // This is to calculate all of the edges in the figure
    // Used to find all internal nodes
    set<int> getInternalNodes();

    // Used to get all of the edges we need
    set<pair<int, int>> allEdges();

    // so that we don't accidentally have all 0's
    void createLaplacian(bool isHand);

    void bumpInternalNodes();

    double getTotalForce() {
        double total = 0;
        forces = getForces();
        for(int i = 0; i < forces.ncols(); i++) {
            total += forces[0][i];
        }
        return total;
    }

    void moveVecIntoV();


    MatrixXd &V;
    qp::Vector<double> vec;
    const MatrixXi &F;
    set<pair<int, int>> edges;
    qp::Matrix<double> forces;
    qp::Vector<double> weights;
    map<pair<int, int>, double> weightMap;
    Indexer indr;
    set<int> internalNodes;
};
