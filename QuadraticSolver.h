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
    QuadraticSolver(MatrixXd &_V, MatrixXi &_F) : V(_V), F(_F) {
        internalNodes = getInternalNodes();
        edges = allEdges();
        indr = Indexer(V.rows(), internalNodes, edges);
    };

    // This is to calculate weights from our estimated forces on the structure
    const map<pair<int, int>, int>& getWeights();

    // This is to update the vertices with weights given to edges fixed
    void updateVertices();

    // This is how we estimate the forces
    qp::Matrix<double> getForces(int numInternal);
    // This is to calculate all of the edges in the figure
    // Used to find all internal nodes
    set<int> getInternalNodes();
    // Used to get all of the edges we need
    set<pair<int, int>> allEdges();

    MatrixXd &V;
    const MatrixXi &F;
    set<pair<int, int>> edges;
    qp::Matrix<double> forces;
    qp::Vector<double> weights;
    map<pair<int, int>, int> weightMap;
    Indexer indr;
    set<int> internalNodes;
};
