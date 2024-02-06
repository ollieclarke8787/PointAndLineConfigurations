# PointAndLineConfigurations
Computations related to realisation spaces for rank 3 matroids

## 3x4 grid example
The file "gridExample.m2" contains code for studying the ideal of the '3x4 grid'. Running this file will verbosely check the following:
1. A minimal set of generators of the prime component P2 from the paper [On the prime decomposition of some determinantal hyperedge ideal] https://www.sciencedirect.com/science/article/pii/S0747717119301348 is given by the polynomials G_(3,4)(c_1, (i,j,k,l,m,n)) where i <= j <= k <= l <= m <= n
2. The other G_(3,4) polynomials are contained inside P2.

The file "gridMatrix.m2" contains code for studying the lifting problem for a grid, i.e., when can 12 points on a line be lifted to a 3x4 grid?
The file contains the matrix M whose rows are the colinearity constraints and columns are the points. The entries of the matrix are linear binomials that
represent the 2-minors that can be pieced together to give the 3-minor contraints.

**Warning!** The polynomial g' in "gridMatrix.m2" is derived from the QSgen function from "gridExample.m2" which needs to be called before g' can be defined.

The file "liftabilityCondition.m2" checks the liftability condition for the grid by computing the 10-minors of the matrix of lifts M. It is recommended to run this file line-by-line as the computations of the 10-minors may take a long time.

The file "independenceOfGens.m2" checks that each claimed minimal generator of the ideal is not redundant.

## Quad-set

The file "IQSradical.m2" shows that the ideal of the quad-set is radical by showing that its initial ideal is square-free.
