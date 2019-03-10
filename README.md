## Scaffolding Minimizing Construction Strategies for Masonry Buildings

This is originally from a copy of [libigl's sample project](https://github.com/libigl/libigl-example-project), which I have extended.

This project is an exercise into solving the question of how to build masonry structures to minimize the number of states in which we are 
unstable. 
To do this, we first pass in a mesh in which each brick is a vertex (and currently use heuristics for forces), and then utilize the algorithm given in [Design of Self-Supporting Surfaces](http://www.geometrie.tugraz.at/wallner/selfsupporting.pdf) to test for stability, and alter the mesh to become a self-supporting surface.


## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `example_bin` binary.

## Run

From within the `build` directory just issue:

    ./example_bin

A glfw app should launch displaying a 3D cube.

## Dependencies

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

We recommend you to install libigl using git via:

    git clone https://github.com/libigl/libigl.git
    cd libigl/
    git submodule update --init --recursive
    cd ..

If you have installed libigl at `/path/to/libigl/` then a good place to clone
this library is `/path/to/libigl-example-project/`.
