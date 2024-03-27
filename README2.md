# TřiVis: Versatile, Reliable, and High-Performance Tool for Computing Visibility in Polygonal Environments

## About

TřiVis, named after the Czech word 'tři' (meaning 'three') and the term 'visibility', is a C++ library for computing visibility-related structures in 2D polygonal environments.
It is based on the triangular expansion algorithm [1,2], which uses preprocessing to convert the input polygonal environment into a triangular mesh, and then traverses the mesh to compute visibility regions, visibility graphs, ray-shooting queries, and other visibility-related structures.

## Definitions

- **Polygonal Environment:** A 2D environment consisting of a single outer boundary and zero or more inner boundaries (holes), where each boundary is a simple polygon.
- **Triangular Mesh:** A mesh of triangles including adjacency information that covers the polygonal environment, where the intersection of any two triangles is either empty, a common edge, or a common vertex.
- **Point Location Query:** Determines the triangle of the triangular mesh that contains a given point in the polygonal environment or reports that the point is outside the environment.
- **Visibility:** Two points in the polygonal environment are visible to each other if the line segment connecting them lies entirely within the environment (touching the boundaries is allowed).
- **Two-Point Visibility Query:** Determines whether two given points in the polygonal environment are visible to each other.
- **Ray-Shooting Query:** Determines an intersection point of a ray starting at a given interior point in a given direction with the boundary of the polygonal environment.
- **Visibility Region Query:** Determines the region of the polygonal environment visible from a given interior point.
- **Visible Vertices Query:** Determines the subset of vertices of the polygonal environment visible from a given interior point.
- **Visible Points Query:** Determines the subset of given interior points that are visible from another given interior point.
- **Vertex-Vertex Visibility Graph:** A graph where vertices represent the vertices of the polygonal environment, and edges represent visibility between pairs of vertices.
- **Point-Point Visibility Graph:** A graph where vertices represent given interior points, and edges represent visibility between pairs of points.
- **Point-Vertex Visibility Graph:** A graph where vertices represent both given interior points and vertices of the polygonal environment, and edges represent visibility between pairs of points and vertices.
- **Limited Visibility Range:** A constraint on the visibility queries that restricts the visibility range to a given distance from the query point.

## Features

## Structure and Dependencies

## Usage

## Documentation

## Examples

```bash
# Play with the 2-point/ray-shooting example:
./build-Release/examples/vis_2point 
./build-Release/examples/vis_2point --vis-radius 10
./build-Release/examples/vis_2point --shoot-ray
./build-Release/examples/vis_2point --help # to see all options
# Play with the visibility region example:
./build-Release/examples/vis_region
./build-Release/examples/vis_region --vis-radius 8
./build-Release/examples/vis_region --vis-radius 8 --sample-arc-edges 0.1
./build-Release/examples/vis_region --vis-radius 8 --sample-arc-edges 0.1 --to-polygon
./build-Release/examples/vis_region --help # to see all options
# Play with the visible vertices example:
./build-Release/examples/vis_vertices
./build-Release/examples/vis_vertices --reflex-only
./build-Release/examples/vis_vertices --vis-radius 6
./build-Release/examples/vis_vertices --help # to see all options
# Play with the visibility graph example:
./build-Release/examples/vis_graph
./build-Release/examples/vis_graph --reflex-only
./build-Release/examples/vis_graph --reflex-only --vis-radius 6
./build-Release/examples/vis_graph --point-vertex --reflex-only --vis-radius 6
./build-Release/examples/vis_graph --point-point
./build-Release/examples/vis_graph --help # to see all options
```

## Performance Evaluation

```bash
conda create env -f conda/trivis_with_boost_cairo_cgal.yml # unless you have it already
conda activate trivis # unless already activated
export TRIVIS_BUILD_TRIVIS_PLUS=1 # to build the dependency
export TRIVIS_BUILD_PERFORMANCE=1 # to build the performance tests
bash build.bash # to run the build
# creates the build-Release directory with the executables
# todo: copy the map instances to data/maps
cd performance # go to the performance directory
bash gen_points_all.bash # to generate the test points
# generates the test points in the ../data/points directory (also saves the CDT mesh for each map in the ../data/meshes directory)
bash test_all.bash # to construct the visibility regions
# generates the raw outputs in the outputs directory 
bash eval_all.bash # to compare the results with the reference
# generates CSV files with the evaluation results in the results/csv directory
rm outputs/* # you can now remove the raw outputs
cd results # go to the results directory
conda create env -f conda/datatable.yml # unless you have it already
conda activate datatable
python process1.py # to process the results
# generates the results.csv file, you can now inspect it
cd ../.. # go back to the root directory
conda deactivate # to deactivate the datatable environment
conda deactivate # to deactivate the trivis environment
```

## Version History

## Planned Development

## License

## Troubleshooting

## References

