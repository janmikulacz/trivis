**Software**

**Title of output (not the same as relevant publication/article):**

TřiVis

**Registration number / Internal identifier:** TřiVis

**Version:** 1

**Project title:** Robotics and Advanced Industrial Production

**Project registration number:** CZ.02.01.01/00/22_008/0004590

**Type of result:** software

**Date of publication:** 29/2/2024

**Authors and contributors (including affiliation, ORCID, contact of at least one contact point):**

Jan Mikula, CTU in Prague, <https://orcid.org/0000-0003-3404-8742>

Miroslav Kulich, CTU in Prague, <https://orcid.org/0000-0002-0997-5889>, <kulich@cvut.cz>

**Keywords:** Polygonal Environment, Triangular Mesh, Visibility, C++

**URL of storage location (DOI or other persistent identifier):**

<https://github.com/janmikulacz/trivis>

**Description (focusing on originality and distinctiveness):**

TřiVis, is a C++ library for computing various visibility-related queries in polygonal environments. TřiVis excels in several aspects compared to similar available implementations:

• Versatility: TřiVis offers an extensive set of features, each with a specialized, user-friendly interface. These features include computing visibility polygons, performing two-points and ray-shooting visibility queries, identifying visible points and vertices, constructing visibility graphs, and executing rapid point location queries. It also provides various utility functions for managing polygonal environments and query outputs. Additionally, TřiVis supports efficient computations of all queries with an optional limited range constraint.

• Reliability: TřiVis is characterized by its reliable and predictable behavior. It avoids crashing or infinite looping and consistently produces outputs that align with user expectations, as confirmed by our evaluation on highly complex query instances.

• Performance: With an average query time of 9 ± 6 µs, TřiVis outperforms all other implementations by at least an order of magnitude, while maintaining preprocessing times below 20 ms for the benchmark instances.

Moreover, the core functionality is independent of external libraries. TřiVis is freely available for private, research, and institutional use at <https://github.com/janmikulacz/trivis>.

Relevance to the project (a connection to research objective RO x.y.): RO 6.3

**Link to technical documentation:**

\[1\] J. Mikula, M. Kulich and L. Přeučil, "TřiVis: Versatile, Reliable, and High-Performance Tool for Computing Visibility in Polygonal Environments," 2024 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), Abu Dhabi, United Arab Emirates, 2024, pp. 10503-10510, doi: 10.1109/IROS58592.2024.10801476.

\[2\] <https://github.com/janmikulacz/trivis>

**License (including restrictions and copyrights):**

see <https://github.com/janmikulacz/trivis/blob/main/LICENSE.md>



**System requirements (e.g. OS/hardware requirements, programming language, libraries):**

TřiVis is self-contained, meaning it does not depend on external libraries. However, it includes some third-party libraries that are freely available for private, research, and institutional use, and they are bundled with TřiVis's source code:

- Triangle, for triangular mesh generation (Shewchuk, 1996).
- Robust Geometric Predicates, for geometry primitives (Shewchuk, 1997).
- Clipper2, for polygon clipping operations and related geometric algorithms.

The build has been tested with the GCC 12.3.0 compiler and CMake 3.28.1 on the Ubuntu 20.04.6 LTS operating system. For your convenience, the conda/trivis.yml file describes Conda environment with the exact compiler and CMake versions used by the authors.

Programming language: C++

**Type of interface/API:** software library

**Links to related publications and/or datasets:**

\[1\] J. Mikula, M. Kulich and L. Přeučil, "TřiVis: Versatile, Reliable, and High-Performance Tool for Computing Visibility in Polygonal Environments," 2024 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), Abu Dhabi, United Arab Emirates, 2024, pp. 10503-10510, doi: 10.1109/IROS58592.2024.10801476.

**Naming convention of file/s and/or software:** N/A