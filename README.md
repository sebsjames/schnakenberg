# Schakenberg reaction diffusion system

This is a standalone example of a program which builds against
mathplot.

## Dependencies

If you are using Debian or Ubuntu, the following `apt` command should
install the mathplot dependencies.

```bash
sudo apt install build-essential cmake git wget  \
                 nlohmann-java3-dev librapidxml-dev \
                 freeglut3-dev libglu1-mesa-dev libxmu-dev libxi-dev \
                 libglfw3-dev libfreetype-dev
```

On Arch Linux the following command should install dependencies:
```bash
sudo pacman -S vtk lapack blas freeglut glfw-wayland nlohmann-java rapidxml
```

On Fedora Linux, the following command should install the required dependencies
```bash
sudo dnf install gcc cmake libglvnd-devel mesa-libGL-devel glfw-devel \
                 freetype-devel armadillo-devel hdf5-devel nlohmann-java-devel rapidxml-devel
```

If you're building on a Mac, you can refer to the [Mac
README](https://github.com/sebsjames/mathplot/blob/main/README.build.mac.md#installation-dependencies-for-mac)
for help. You only need to obtain and build
[glfw3](https://github.com/sebsjames/mathplot/blob/main/README.build.mac.md#glfw3);
OpenGL and Freetype should already be installed by default.

## Building

To build and run the example:

```bash
# Clone this example
git clone git@github.com:sebsjames/schnakenberg

# Clone, copy or symlink mathplot INSIDE your example:
cd mathplot_template # or whatever you named your fork/copy
git clone --recurse-submodules git@github.com:sebsjames/mathplot

# Build schakenberg in a 'build' directory
mkdir build
cd build
cmake ..
make
./build/schnakenberg ./schakenberg.json
```
