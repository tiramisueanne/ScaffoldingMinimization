#pragma once
#include <Eigen/Dense>

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <getopt.h>

#include "HandMadeMeshes/MediumMesh.h"
#include "HandMadeMeshes/SmallMesh.h"
#include "HandMadeMeshes/TinyMesh.h"

using namespace std;
using namespace Eigen;

enum Calculation { QUADRATIC, DUAL, JUST_SHOW};

void parseInput(int argc, char* argv[], MatrixXd& V, MatrixXi& F,
                Calculation& type) {
    // Inspired heavily by systutorials:
    // URL: https://www.systutorials.com/docs/linux/man/3-getopt/
    int c;
    int digit_optind = 0;
    while (true) {
        int this_optind = optind ? optind : 1;
        int option_index = 0;
        static struct option long_options[]{
            {"quadratic", no_argument, 0, 0},
            {"dual", no_argument, 0, 0},
            {0, 0,  0,  0 },
        };

        c = getopt_long(argc, argv, "", long_options, &option_index);
        if(c == -1) break;
        type = JUST_SHOW;
        switch(c) {
        case 0:
            cout << "Got a long flag" << endl;
            if(option_index == 0) {
                type = QUADRATIC;
            }
            else {
                type = DUAL;
            }
        default:
            cerr << "Got an unknown argument!" << endl;
        }
    }
    if(optind < argc) {
        string name = argv[optind];
        if (name == "tiny") {
            initHandMesh(V, F);
        } else if (name == "small") {
            initBiggerMesh(V, F);
        } else if (name == "medium") {
            evenBiggerMesh(V, F);
        } else {
            igl::readOBJ(argv[optind], V, F);
        }
    }
    else {
        initHandMesh(V, F);
    }
}
