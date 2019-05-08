#include <Eigen/Dense>
#include <iostream>

#include "QuadraticSolver.h"

using namespace std;
using namespace Eigen;

void QuadraticSolver::createLaplacian(bool isHand = false) {
    MatrixXd lap(V.rows(), V.rows());
    lap.setZero();
    cout << "Done creating the big mat" << endl;
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
    cout << "Done with first internal adj" << endl;
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
    cout << "Setting up the valence done!" << endl;
    VectorXd onesForBounds(V.rows());
    onesForBounds = VectorXd::Ones(V.rows());
    cout << "Ones for bounds created" << endl;

    HouseholderQR<MatrixXd> decomp(lap);
    cout << "Created a decomp thing" << endl;
    V.col(2) = decomp.solve(onesForBounds);
    cout << "The decompositoin is done" << endl;
    #ifdef ACKNOWLEDGE_ISHAND
    if (!isHand) {
        V.col(2) *= pow(10, -1.5);
    }
    #endif
}

