#include <Eigen/Dense>
#include <iostream>

#include "QuadraticSolver.h"

using namespace std;
using namespace Eigen;

void QuadraticSolver::createLaplacian(bool isHand = false) {
    MatrixXd lap(V.rows(), V.rows());
    lap.setZero();

    // Create adjacency matrix first, without doing external rows
    // For each triangle, set neighbor vertices on correct rows
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < F.cols(); j++) {
            int currRow = F(i, j);
            // We shall do this part in the setting
            if (unsupportedNodes.find(currRow) == unsupportedNodes.end()) {
                continue;
            }
            // Find the two neighbors of currRow from i triangle
            int neigh1 = F(i, (j + 1) % F.cols());
            int neigh2 = F(i, (j + 2) % F.cols());
            lap(currRow, neigh1) = -1;
            lap(currRow, neigh2) = -1;
        }
    }

    for (int i = 0; i < lap.rows(); i++) {
        // If this is an external row, identity row
        if (unsupportedNodes.find(i) == unsupportedNodes.end()) {
            for (int x = 0; x < lap.cols(); x++) {
                lap(i, x) = 0;
            }
            lap(i, i) = 1;
        } else {
            int sum = 0;
            // Set the valence
            for (int j = 0; j < lap.cols(); j++) {
                sum += lap(i, j);
            }
            // The previous sum will be positive, so negate it
            lap(i, i) = -sum;
        }
    }

    VectorXd onesForBounds(V.rows());
    onesForBounds = VectorXd::Ones(V.rows());

    FullPivLU<MatrixXd> decomp(lap);
    V.col(2) = decomp.solve(onesForBounds);
    if(!isHand) {
        V.col(2) *= pow(10, -3);
    }
}
