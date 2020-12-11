# DARPA-SSCA2

Graph theoretic problems are representative of fundamental
computations in traditional and emerging scientific disciplines like
scientific computing and computational biology, as well as
applications in national security. We present our design and
implementation of a graph theory application that supports the kernels
from the Scalable Synthetic Compact Applications (SSCA) benchmark
suite, developed under the DARPA High Productivity Computing Systems
(HPCS) program. This synthetic benchmark consists of four kernels that
require irregular access to a large, directed, weighted
multi-graph. We have developed a parallel implementation of this
benchmark in C using the POSIX thread library for commodity symmetric
multiprocessors (SMPs). In this paper, we primarily discuss the data
layout choices and algorithmic design issues for each kernel, and also
present execution time and benchmark validation results.

References:

D.A. Bader, K. Madduri ``Design and Implementation of the HPCS Graph
Analysis Benchmark on Symmetric Multiprocessors,'' Proc. 12th
International Conference on High Performance Computing (HiPC 2005),
Goa, India, December 2005.
