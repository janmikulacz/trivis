# TřiVis

***Versatile, Reliable, and High-Performance Tool for Computing Visibility in Polygonal Environments***

## About

TřiVis, named after the Czech word *'tři'* (meaning *'three'*) and the term *'visibility'*, is a C++ library for computing visibility-related structures in 2D polygonal
environments.
It is based on the *triangular expansion algorithm* ([Bungiu et al., 2014](https://arxiv.org/abs/1403.3905); [Xu and Güting, 2015](https://doi.org/10.1007/s10707-014-0213-7)),
which uses preprocessing to convert the input polygonal environment into a triangular mesh, and then traverses the mesh to compute visibility regions, visibility graphs,
ray-shooting queries, and other visibility-related structures.

## Table of Contents

- [Definitions](#definitions)
- [Main Features](#main-features)
  - [Supported Queries](#supported-queries)
  - [Worst-Case Query Complexity](#worst-case-query-complexity)
  - [Efficient Limited-Range Queries](#efficient-limited-range-queries)
  - [Fast Point Location](#fast-point-location)
  - [Robustness and Reliability](#robustness-and-reliability)
  - [Performance on Large Instances](#performance-on-large-instances)
- [Building and Linking](#building-and-linking)
- [Basic Usage](#basic-usage)
  - [Initialization](#initialization)
  - [Computing Visibility Queries](#computing-visibility-queries)
- [Dependencies](#dependencies)
- [Repository Contents](#repository-contents)
- [TřiVis Codebase](#třivis-codebase)
- [TřiVis+](#třivis-1)
- [Examples](#examples)
- [Performance Evaluation](#performance-evaluation)
  - [Replication Guide](#replication-guide)
- [Documentation](#documentation)
- [License](#license)

## Definitions

- **Polygonal Environment:** A 2D environment consisting of a single outer boundary and zero or more inner boundaries (holes), where each boundary is a simple polygon.
- **Triangular Mesh:** A mesh of triangles including adjacency information that covers the polygonal environment, where the intersection of any two triangles is either empty, a
  common edge, or a common vertex.
- **Point Location Query:** Determines the triangle of the triangular mesh that contains a given point in the polygonal environment or reports that the point is outside the
  environment.
- **Visibility:** Two points in the polygonal environment are visible to each other if the line segment connecting them lies entirely within the environment (touching the
  boundaries is allowed).
- **Two-Point Visibility Query:** Determines whether two given points in the polygonal environment are visible to each other.
- **Ray-Shooting Query:** Determines an intersection point of a ray starting at a given interior point in a given direction with the boundary of the polygonal environment.
- **Visibility Region Query:** Determines the region of the polygonal environment visible from a given interior point.
- **Visible Vertices Query:** Determines the subset of vertices of the polygonal environment visible from a given interior point.
- **Visible Points Query:** Determines the subset of given interior points that are visible from another given interior point.
- **Vertex-Vertex Visibility Graph:** A graph where vertices represent the vertices of the polygonal environment, and edges represent visibility between pairs of vertices.
- **Point-Point Visibility Graph:** A graph where vertices represent given interior points, and edges represent visibility between pairs of points.
- **Point-Vertex Visibility Graph:** A graph where vertices represent both given interior points and vertices of the polygonal environment, and edges represent visibility between
  pairs of points and vertices.
- **Limited Visibility Range:** A constraint on the visibility queries that restricts the visibility range to a given distance from the query point.

## Main Features

### Supported Queries

TřiVis supports all the queries defined [above](#definitions), each with a specialized, easy-to-use API.

### Worst-Case Query Complexity

The main feature of TřiVis is the *triangular expansion
algorithm* ([Bungiu et al., 2014](https://arxiv.org/abs/1403.3905); [Xu and Güting, 2015](https://doi.org/10.1007/s10707-014-0213-7)), which provides $O(n)$ query time for
two-point visibility and ray-shooting queries, $O(nh)$ for visibility region and visible vertices/points queries, $O(n^2h)$ for the vertex-vertex visibility graph, and $O(pnh)$ for
point-point and point-vertex visibility graphs.
Here, $n$ is the number of vertices in the polygonal environment, $h$ is the number of holes, and $p$ is the number of input interior points.

It is important to note that the query time complexity is worst-case, and the actual query time may be significantly lower in practice. This is especially true for TřiVis, as the
*triangular expansion algorithm* is known to exhibit output-sensitive behavior to some extent, traversing only the triangles that are partially or fully visible from the query
point ([Bungiu et al., 2014](https://arxiv.org/abs/1403.3905)).

### Efficient Limited-Range Queries

TřiVis supports limited-range queries, employing an early termination mechanism to avoid traversing unnecessary triangles of the triangular mesh.
This can significantly reduce the query time for queries with a highly restricted visibility range relative to the average size of the open space in the environment.

### Fast Point Location

TřiVis employs a bucketing technique to accelerate point location queries, which are essential for visibility-related queries.
Using the bucketing technique from ([Edahiro et al., 1984](https://doi.org/10.1145/357337.357338)), TřiVis achieves $O(1)$ expected query time for point location queries.

### Robustness and Reliability

TřiVis relies on floating-point arithmetic and incorporates a unique combination of adaptive robust geometry predicates ([Shewchuk, 1997](https://doi.org/10.1007/PL00009321)) and
$\epsilon$-geometry ([Fortune and Van Wyk, 1993](https://doi.org/10.1145/160985.161015)).
Thanks to this combination, TřiVis is characterized by reliable and predictable behavior, free from crashing or infinite looping, providing consistent outputs that align with user
expectations, as supported by the results of our evaluation on exceptionally complex query instances.
For more information, refer to the [Performance Evaluation](#performance-evaluation) section.

### Performance on Large Instances

TřiVis has been tested on 34 complex polygonal environments with up to 8,320 vertices and 679 holes by computing visibility regions for 1,000 interior points in each environment.
With an average query time of 9±6μs, TřiVis outperforms other tested implementations by at least an order of magnitude, while keeping the preprocessing time below 20ms for the
tested polygonal environments.
For more information, refer to the [Performance Evaluation](#performance-evaluation) section.

## Building and Linking

TřiVis is built with [CMake](https://cmake.org/).

To make the library available to your C++ project, you can copy the [trivis/](trivis) directory and include it in your CMakeLists.txt:

```CMake
add_subdirectory(path_to_trivis) # Add the library to your project.
target_link_libraries(your_target PUBLIC Trivis) # Link the library to your target.
```

The build has been tested with the [GCC](https://gcc.gnu.org/) 12.3.0 compiler and [CMake](https://cmake.org/) 3.28.1 on the [Ubuntu](https://ubuntu.com/) 20.04.6 LTS operating
system.
For your convenience, the [conda/trivis.yml](conda/trivis.yml) file describes [Conda](https://docs.conda.io/en/latest/) environment with the exact compiler and CMake versions used
by the authors.

To create and test the environment, run the following in the root directory of this repository:

```bash
# Assuming you have Conda installed:
conda env create -f conda/trivis.yml
conda activate trivis
gcc --version
cmake --version
mkdir build
cd build
cmake ../trivis
make
```

Apart from TřiVis's source code, this repository also contains [TřiVis+](#třivis-1), an extension adding support for data loading and visualization, and an [Examples](#examples) project depending on TřiVis+,
demonstrating the library's usage.

For more information, refer to the [Repository Contents](#repository-contents) section.

## Basic Usage

### Initialization

To use the library in your C++ project, include the main header file:

```C++
#include "trivis/trivis.h"
```

TřiVis's main API is the `trivis::Trivis` class, which provides methods for visibility-related queries.
It is initialized with a polygonal environment represented by a `trivis::geom::PolyMap` object.

```C++
trivis::Trivis vis;
{   // Scope for the environment.
    trivis::geom::PolyMap environment;
    // TODO: Fill the environment.
    vis = trivis::Trivis{std::move(environment)};
}   // The environment was moved to vis.
```

The triangular mesh of the environment has just been constructed.
To inspect it, you can access the mesh directly:

```C++
const trivis::mesh::TriMesh &mesh = vis.mesh();
```

### Computing Visibility Queries

To perform visibility-related queries, you need to provide the query point, the point location result, and the visibility range (if limited).

```C++
trivis::geom::FPoint query;
// TODO: Fill the query.
std::optional<trivis::Trivis::PointLocationResult> plr = vis.LocatePoint(query);
if (plr.has_value()) {
    // The query point is inside the environment.
} else {
    // The query point is outside the environment.
}
std::optional<double> range;
// TODO: Set the visibility range or leave it unlimited.
```

**Two-point visibility query:**

```C++
trivis::geom::FPoint target;
// TODO: Fill the target.
if (plr.has_value()) {
    bool is_visible = vis.IsVisible(query, plr.value(), target, range);
}
```

**Ray-shooting query:**

```C++
trivis::geom::FPoint direction;
// TODO: Fill the direction.
if (plr.has_value()) {
    std::optional<trivis::Trivis::RayShootingResult> rsr = vis.ShootRay(query, plr.value(), direction, range);
    if (rsr.has_value()) {
        trivis::geom::FPoint intersection = rsr->p;
    } else {
        // The intersection point is outside the visibility range.
    }
}
```

**Visibility region query:**

```C++
if (plr.has_value()) {
    trivis::AbstractVisibilityRegion avr = vis.VisibilityRegion(query, plr.value(), range);
    trivis::RadialVisibilityRegion rvr = vis.ToRadialVisibilityRegion(avr);
    if (range.has_value()) {
        rvr.IntersectWithCircleCenteredAtSeed();
    }
    // Optional post-processing:
    if (range.has_value()) {
        vis_region.SampleArcEdges(M_PI / 180.0);
    }
    trivis::geom::FPolygon polygon_approx = vis_region->ToPolygon();
}
```

**Visible vertices query:**

```C++
if (plr.has_value()) {
    std::vector<int> visible_vertices = vis.VisibleVertices(query, plr.value(), range);
}
```

**Visible points query:**

```C++
std::vector<trivis::geom::FPoint> input_points;
// TODO: Fill the input points.
std::vector<std::optional<trivis::Trivis::PointLocationResult>> input_plrs;
input_plrs.reserve(input_points.size());
for (const auto &p : input_points) {
    input_plrs.push_back(vis.LocatePoint(p));
}
if (plr.has_value()) {
    std::vector<int> visible_points = vis.VisiblePoints(query, plr.value(), input_points, input_plrs, range);
}
```

**Visibility graphs:**

```C++
std::vector<std::vector<int>> vertex_vertex_graph = vis.VertexVertexVisibilityGraph(range);
std::vector<std::vector<int>> point_point_graph = vis.PointPointVisibilityGraph(input_points, input_plrs, range);
std::vector<std::vector<int>> point_vertex_graph = vis.PointVertexVisibilityGraph(input_points, input_plrs, range);
```

## Dependencies

TřiVis is self-contained, meaning it does not depend on external libraries.
However, it includes some third-party libraries that are freely available for private, research, and institutional use, and they are bundled with TřiVis's source code:

- [Triangle](https://www.cs.cmu.edu/~quake/triangle.html), for triangular mesh generation ([Shewchuk, 1996](https://doi.org/10.1007/BFb0014497)).
- [Robust Geometric Predicates](https://github.com/dengwirda/robust-predicate), for geometry primitives ([Shewchuk, 1997](https://doi.org/10.1007/PL00009321)).
- [Clipper2](http://www.angusj.com/clipper2/Docs/Overview.htm), for polygon clipping operations and related geometric algorithms.

These third-party libraries are located in the [trivis/lib](trivis/lib) directory.
While the source code of these libraries is included in TřiVis, their licenses are preserved and can be found in the respective subdirectories of the [trivis/lib](trivis/lib)
directory.
Please refer to the exact licencing terms in [LICENSE.md](LICENSE.md).

If you use TřiVis, please make sure to comply with the licenses of the third-party libraries and give proper attribution to their authors.

## Repository Contents

This repository includes TřiVis, as well as some extensions, example projects and utilities:

- [conda/](conda): Convenience [Conda](https://docs.conda.io/en/latest/) environment files.
- [data/](data): A data directory structure, initially empty except for a single example map, for storing polygonal environments (maps), meshes, and query points. The
  example projects assume this directory structure.
- [examples/](examples): [Example project demonstrating the usage of TřiVis.](#examples)
- [performance/](performance): [Performance evaluation project.](#performance-evaluation)
- [trivis/](trivis): [TřiVis library (core).](#třivis-codebase)
- [trivis_plus/](trivis_plus): [TřiVis+ library, an extension adding support for data loading and visualization.](#třivis-1)
- [build.bash](build.bash): A convenience script for building the library and example projects.
- [CMakelists.txt](CMakeLists.txt): The main CMake configuration file for building the library and example projects.
- [LICENSE.md](LICENSE.md): The license file for TřiVis.
- [README.md](README.md): This README file.
- [setup.bash](setup.bash): A convenience script for setting the environment variables configuring the build. By default, it contains configuration for building just the TřiVis
  library, but you can modify it to include TřiVis+ and the example projects.

## TřiVis Codebase

Assuming [trivis/](trivis) is the root directory of the TřiVis codebase, the codebase directory structure is as follows:

- [include/trivis/](trivis/include/trivis): Header files of the TřiVis library.
    - [geom/](trivis/include/trivis/geom): Geometric types and utilities (includes a wrapper around [Robust Geometric Predicates](https://github.com/dengwirda/robust-predicate)).
    - [mesh/](trivis/include/trivis/mesh): Triangular mesh type and utilities (includes a wrapper around [Triangle](https://www.cs.cmu.edu/~quake/triangle.html)).
    - [pl/](trivis/include/trivis/pl): Point location class and utilities.
    - [utils/](trivis/include/trivis/utils): General utilities (includes conversions to and from [Clipper2](http://www.angusj.com/clipper2/Docs/Overview.htm)).
    - [trivis.h](trivis/include/trivis/trivis.h): The main TřiVis API (`trivis::Trivis` class).
    - [vis_regions.h](trivis/include/trivis/vis_regions.h): Visibility region types and utilities.
- [lib/](trivis/lib): Contains the third-party libraries bundled with TřiVis.
    - [clipper/](trivis/lib/clipper): [Clipper2](http://www.angusj.com/clipper2/Docs/Overview.htm).
    - [robust-predicate/](trivis/lib/robust-predicate): [Robust Geometric Predicates](https://github.com/dengwirda/robust-predicate).
    - [triangle/](trivis/lib/triangle): [Triangle](https://www.cs.cmu.edu/~quake/triangle.html).
    - [CMakelists.txt](trivis/lib/CMakeLists.txt): The CMake configuration file for building the third-party libraries.
- [src/trivis/](trivis/src/trivis): Contains the source files of the TřiVis library. Copies the structure of [include/trivis/](trivis/include/trivis).
- [CMakelists.txt](trivis/CMakeLists.txt): The CMake configuration file for building the TřiVis library.

## TřiVis+

TřiVis+ is an extension of TřiVis that adds support for data loading and visualization, as well as some extra utilities.

It includes two additional, heavier dependencies, which is why it is separated from the TřiVis (core) library:

- [Boost](https://www.boost.org/), for filesystem operations, program options and logging.
- [Cairo](https://www.cairographics.org/), for 2D graphics rendering.

To make TřiVis+ available to your C++ project, you can copy the [trivis/](trivis) and [trivis_plus/](trivis_plus) directories and include them in your CMakeLists.txt.

```CMake
add_subdirectory(path_to_trivis) # Add the library to your project.
add_subdirectory(path_to_trivis_plus) # Add the extension to your project.
target_link_libraries(your_target PUBLIC TrivisPlus) # Link just the extension to your target.
```

For your convenience, the [conda/trivis_boost_cairo.yml](conda/trivis_boost_cairo.yml) file describes [Conda](https://docs.conda.io/en/latest/) environment with the exact compiler,
CMake version, and the Boost and Cairo libraries used by the authors.

To create and test the environment, run the following in the root directory of this repository:

```bash
# Assuming you have Conda installed:
conda env create -f conda/trivis_boost_cairo.yml 
conda activate trivis
mkdir build
cd build
echo "add_subdirectory(../trivis .)\nadd_subdirectory(../trivis_plus .)" > CMakeLists.txt 
cmake .
make
```

Assuming [trivis_plus/](trivis_plus) is the root directory of the TřiVis+ codebase, the codebase directory structure is as follows:

- [include/trivis_plus/](trivis_plus/include/trivis_plus): Header files of the TřiVis+ extension.
    - [data_loading/](trivis_plus/include/trivis_plus/data_loading): Data loading utilities.
    - [drawing/](trivis_plus/include/trivis_plus/drawing): Drawing utilities (includes a wrapper around [Cairo](https://www.cairographics.org/)).
    - [utils/](trivis_plus/include/trivis_plus/utils): General utilities, including logging.
    - [trivis_plus.h](trivis_plus/include/trivis_plus/trivis_plus.h): Convenience header including all TřiVis+ headers.
- [src/trivis_plus/](trivis_plus/src/trivis_plus): Contains the source files of the TřiVis+ extension. Copies the structure
  of [include/trivis_plus/](trivis_plus/include/trivis_plus).
- [CMakelists.txt](trivis_plus/CMakeLists.txt): The CMake configuration file for building the TřiVis+ extension.

For usage examples, refer to the [Examples](#examples) sections.

## Examples

The example project located in the [examples/](examples) directory demonstrates the usage of TřiVis (and secondarily also TřiVis+).

To build the project, you may first need to create the respective Conda environment:
```bash
conda env create -f conda/trivis_boost_cairo.yml 
conda activate trivis
```
Alternatively, you may already have all the necessary dependencies installed.

Now, you need to set the environment variables to configure the build:
```bash
export TRIVIS_BUILD_TRIVIS_PLUS=1
export TRIVIS_BUILD_EXAMPLES=1
```
Alternatively, you may modify the [setup.bash](setup.bash) to fix the configuration and run `source setup.bash`.

Now, you can build the project:
```bash
bash build.bash
```
This will create the `build-Release` directory with the executables located in the `build-Release/examples` directory.

Now you can experiment with the examples and once you are done, you can inspect the source codes in the [examples/](examples) directory.

When experimenting with the examples, the default behavior is that the input polygonal environment is loaded from the [data/maps/](data/maps) directory, the text output is printed to the console and the graphical output is saved to the [examples/outputs/](examples/outputs) directory.

Next are some examples of how to run the example executables.

**2-Point/Ray-Shooting Example:**
```bash
./build-Release/examples/vis_2point 
./build-Release/examples/vis_2point --vis-radius 10
./build-Release/examples/vis_2point --shoot-ray
./build-Release/examples/vis_2point --help # to see all options
```

**Visibility Region Example:**
```bash
./build-Release/examples/vis_region
./build-Release/examples/vis_region --vis-radius 8
./build-Release/examples/vis_region --vis-radius 8 --sample-arc-edges 0.1
./build-Release/examples/vis_region --vis-radius 8 --sample-arc-edges 0.1 --to-polygon
./build-Release/examples/vis_region --help # to see all options
```

**Visible Vertices Example:**
```bash
./build-Release/examples/vis_vertices
./build-Release/examples/vis_vertices --reflex-only
./build-Release/examples/vis_vertices --vis-radius 6
./build-Release/examples/vis_vertices --help # to see all options
```

**Visibility Graph Example:**
```bash
./build-Release/examples/vis_graph
./build-Release/examples/vis_graph --reflex-only
./build-Release/examples/vis_graph --reflex-only --vis-radius 6
./build-Release/examples/vis_graph --point-vertex --reflex-only --vis-radius 6
./build-Release/examples/vis_graph --point-point
./build-Release/examples/vis_graph --help # to see all options
```

Once you are done experimenting with the examples, and you have inspected the source codes, you can start writing your own experimental code based on the examples.

For your convenience, we provide the [examples/sandbox.cc](examples/sandbox.cc) file, which you can use as a starting point for your own experiments.
It is built together with the examples by default, and you can run it as follows:
```bash
./build-Release/examples/sandbox --help # to see all options
```

## Performance Evaluation

TřiVis, or more specifically, its computation of visibility regions, has been extensively evaluated against several competing implementations on a dataset of 34 complex polygonal environments.

**The competing implementations were:**
- **CGAL-TEA-CE:** [CGAL's triangular expansion algorithm](https://doc.cgal.org/latest/Visibility_2/classCGAL_1_1Triangular__expansion__visibility__2.html) (CGAL version 5.6) with [exact constructions](https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Exact__predicates__exact__constructions__kernel.html).
- **CGAL-TEA-CI:** [CGAL's triangular expansion algorithm](https://doc.cgal.org/latest/Visibility_2/classCGAL_1_1Triangular__expansion__visibility__2.html) (CGAL version 5.6) with [inexact constructions](https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Exact__predicates__inexact__constructions__kernel.html).
- **CGAL-RSA-CE:** [CGAL's rotational sweep algorithm](https://doc.cgal.org/latest/Visibility_2/classCGAL_1_1Rotational__sweep__visibility__2.html) (CGAL version 5.6) with [exact constructions](https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Exact__predicates__exact__constructions__kernel.html).
- **CGAL-RSA-CI:** [CGAL's rotational sweep algorithm](https://doc.cgal.org/latest/Visibility_2/classCGAL_1_1Rotational__sweep__visibility__2.html) (CGAL version 5.6) with [inexact constructions](https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Exact__predicates__inexact__constructions__kernel.html).
- [**VisiLibity1**](https://karlobermeyer.github.io/VisiLibity1/).

**The dataset was derived from the [Iron Harvest](https://bitbucket.org/shortestpathlab/benchmarks/src/master/poly-maps/iron-harvest/) dataset ([Harabor et al., 2022](https://doi.org/10.1609/socs.v15i1.21770)).**
We have only preprocessed the dataset to TřiVis's input format by selecting the largest polygon per map, including its holes, and discarding the remaining polygons (mainly artifacts).
If you are going to use the dataset, please give proper attribution to the authors of the original dataset.


The evaluation is detailed in a manuscript currently (as of March 2024) submitted to the 2024 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS).

Here, we provide a guide on how to replicate the performance tests from the manuscript.

### Replication Guide

**Build the project:**

The performance evaluation project is located in the [performance/](performance) directory.
It is dependent on the TřiVis+ extension, and includes two additional dependencies:
- [CGAL](https://www.cgal.org/) (version 5.6), which can be installed system-wide or in a [Conda](https://docs.conda.io/en/latest/) environment.
- [VisiLibity1](https://karlobermeyer.github.io/VisiLibity1/), which is included in the [performance/lib/](performance/lib) directory.

To build the project, you may first need to create the respective Conda environment:
```bash
conda env create -f conda/trivis_boost_cairo_cgal.yml
conda activate trivis
```
Alternatively, you may already have all the necessary dependencies installed.

Now, you need to set the environment variables to configure the build:
```bash
export TRIVIS_BUILD_TRIVIS_PLUS=1
export TRIVIS_BUILD_PERFORMANCE=1
```
Alternatively, you may modify the [setup.bash](setup.bash) to fix the configuration and run `source setup.bash`.

Building the project is now as simple as running:
```bash
bash build.bash
```
This will create the `build-Release` directory with the executables located in the `build-Release/performance` directory.

**Get the dataset:**

Download the dataset from [here](https://imr.ciirc.cvut.cz/uploads/Downloads/trivis_performance_dataset.zip) and extract it to the [data/maps/](data/maps) directory.

**Generate the test points:**

The following will generate the test points in the [data/points/](data/points) directory and save the CDT mesh for each map in the [data/meshes/](data/meshes) directory:

```bash
cd performance
bash gen_points_all.bash
```

**Construct the visibility regions:**

The following will construct the visibility regions by all the implementations and save the raw outputs in the [performance/outputs/](performance/outputs) directory:

```bash
bash test_all.bash # this may take several hours
```

**Compare the results with the reference implementation (CGAL-TEA-CE):**

The following will compare the results with the reference implementation and save the evaluation results in the [performance/results/csv/](performance/results/csv) directory:

```bash
bash eval_all.bash # this may take several minutes
rm outputs/* # you can now remove the raw outputs
```

**Process the results:**

The results are processed with [datatable](https://datatable.readthedocs.io/en/latest/).
To use it, you may need to create a Conda environment with the [conda/datatable.yml](conda/datatable.yml) file:

```bash
conda create env -f ../conda/datatable.yml
conda activate datatable
```

Now, you can process the results:

```bash
cd results
python process1.py
```

This will generate the `results.csv` file in the [performance/results/](performance/results) directory, which you can inspect and compare with the authors' results located in the [performance/results/authors_results.csv](performance/results/authors_results.csv) file.

## Documentation

Currently, TřiVis is not documented in detail beyond this README file.
The reason is simple: the authors did not have the time to write the documentation yet.
However, the authors hope that thanks to simplicity of the API, the [Basic Usage](#basic-usage) section, and the examples in the [examples/](examples) directory (further documented in the [Examples](#examples) section), you should be able to start using TřiVis without major issues.

## License

Please see the [LICENSE.md](LICENSE.md) file for the license of TřiVis.