# Schakenberg reaction diffusion system

This is a standalone example of a reaction diffusion program which
builds against mathplot.

This example project demonstrates how you
can write and build a program that uses [mathplot](https://github.com/sebsjames/mathplot).

## Dependencies

If you are using Debian or Ubuntu, the following `apt` command should
install the mathplot dependencies.

```bash
sudo apt install build-essential cmake git ninja-build  \
                 freeglut3-dev libglu1-mesa-dev libxmu-dev libxi-dev \
                 libglfw3-dev libfreetype-dev clang-20 clang-tools-20 librapidxml-dev libhdf5-dev
```

If you're building on a Mac, you can refer to mathplot's [Mac
README](https://github.com/sebsjames/mathplot/blob/main/README.build.mac.md#installation-dependencies-for-mac)
for help. You only need to obtain and build
[glfw3](https://github.com/sebsjames/mathplot/blob/main/README.build.mac.md#glfw3);
OpenGL and Freetype should already be installed by default.

## Building

To build and run the example:

```bash
# Clone this example
git clone git@github.com:sebsjames/schnakenberg

# Bring in the submodules - sebsjames/mathplot, sebsjames/maths , sebsjames/maths_config and nlohmann/json
cd schnakenberg
git submodule init
git submodule update

# Build schnakenberg in a 'build' directory
mkdir build
cd build
CC=clang-20 CXX=clang++-20 cmake .. -GNinja
ninja
cd ..

# Run the example
./build/schnakenberg ./schakenberg.json
```
