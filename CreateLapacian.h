#include <Eigen/Dense>
#include <iostream>

using namespace std;
using namespace Eigen;

void QuadraticSolver::createLaplacian() {
    MatrixXd lap(V.rows(), V.rows());

    // Create adjacency matrix first

    // For each triangle, set neighbor vertices on correct rows
    for(int i = 0; i < F.rows(); i++) {
        for(int j = 0; j < F.cols(); j++) {
            int currRow = F(i, j);
            if(internalNodes.find(currRow) == internalNodes.end()) {
                continue;
            }
            int neigh1 = F(i, (j + 1) % F.cols());
            int neigh2 = F(i, (j + 2) % F.cols());
            lap(currRow, neigh1) = 1;
            lap(currRow, neigh2) = 1;
        }
    }

    for (int i = 0; i < lap.rows(); i++) {
        bool isInternal = (internalNodes.find(i) != internalNodes.end());
        for (int j = 0; j < lap.cols(); j++) {
            // If this isn't an internal node
            if (!isInternal) {
                if (i == j) {
                    lap(i, j) = 1;
                } else {
                    lap(i, j) = 0;
                }

            } else {
                // find neighbors of row i
                // go through each face it's a part of and add
                // all neighbors from each of those
            }
        }
    }
}
}
