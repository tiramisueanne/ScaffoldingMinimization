#include <igl/HalfEdgeIterator.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include <QuadProg++/QuadProg++.hh>
#include "QuadraticSolver.h"

using namespace std;
using namespace Eigen;
namespace qp = quadprogpp;

set<pair<int, int>> QuadraticSolver::allEdges() {
    set<pair<int, int>> edges;
    for (int currFace = 0; currFace < F.rows(); currFace++) {
        RowVector3i face = F.row(currFace);
        for (int i = 0; i < 3; i++) {
            // If this edge doesn't connect to an internal node
            if (unsupportedNodes.find(face(i)) == unsupportedNodes.end() &&
                unsupportedNodes.find(face((i + 1) % 3)) ==
                    unsupportedNodes.end()) {
                continue;
            }
            edges.insert(pair<int, int>(face(i), (face((i + 1) % 3))));
            edges.insert(pair<int, int>(face((i + 1) % 3), (face(i))));
        }
    }
    return edges;
}

// DEPRECATED: use createLaplacian instead.
void QuadraticSolver::bumpInternalNodes() {
    for (auto node : unsupportedNodes) {
        V.row(node).z() += 2;
#ifdef DEBUG
        cout << " The internal node z value is " << V.row(node).z() << endl;
#endif
    }
}

int QuadraticSolver::deleteANode(int index) {
    set<int> rowsToGo;
    MatrixXi newF;
    // Go through all triangles
    for (int i = 0; i < F.rows(); i++) {
        bool removeThis = false;
        for (int j = 0; j < F.cols(); j++) {
            if (F(i, j) == index) {
                rowsToGo.insert(i);
            }
        }
    }

    newF = MatrixXi(F.rows() - rowsToGo.size(), F.cols());
    int fCount = 0;
    for (int i = 0; i < newF.rows(); i++) {
        while (rowsToGo.find(fCount) != rowsToGo.end()) fCount++;
        newF.row(i) = F.row(fCount);
    }
    F = newF;

    return newF.rows();
}

// Unsupport the node with the least force on it
int QuadraticSolver::removeSmallestNode() {
    forces = getForces();
    cout << "The number of forces in remove smallest" << forces.rows() << endl;
    VectorXd force(forces.cols());
    for (int i = 0; i < forces.cols(); i++) {
        force(i) = forces(0, i);
    }
    int indexIntoUnsupported;
    force.minCoeff(&indexIntoUnsupported);
    vector<int> vertVec = indr.vertVect();
    auto realIndexIter =
        std::find(vertVec.begin(), vertVec.end(), indexIntoUnsupported);
    int realIndex = realIndexIter - vertVec.begin();
    int toRet = deleteANode(realIndex);
    if (toRet != 0) {
        assert(unsupportedNodes.erase(realIndex) > 0);
        edges = allEdges();
        indr = Indexer(V.rows(), unsupportedNodes, edges);
        weights = VectorXd::Constant(edges.size() / 2, 0);
    }
    return toRet;
}
