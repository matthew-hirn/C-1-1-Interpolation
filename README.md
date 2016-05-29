# C-1-1-Interpolation
Computes optimal *C^{1,1}(R^d)* interpolations in any dimension *d*

## Overview

This software solves the following interpolation problem:
Let *E* be a finite set of *d*-dimensional points. To each point *a* in *E*, there is a specified scalar function value *f(a)* and a gradient vector *Df(a)*. The software computes a function *F:R^d -> R* such that:
  1. *F(a) = f(a)* and *∇F(a) = Df(a)* for each *a* in *E*.
  2. Amongst all such interpolants satisfying 1., the value of of Lip(*∇F*) is minimal.

The software is based upon the paper "Computing minimal interpolants in C^{1,1}(R^d)," by Ariel Herbert-Voss, Matthew Hirn, and Frederick McCollum (https://arxiv.org/abs/1411.5668).

## Installation

## Using the Software

The function EFD can be used to generate uniformly distributed pseudorandom input: data points (E), function values (F), and partial derivatives (D). For 20 initial points in R^2, use the following command:
```matlab
>> [E, F, D] = EFD(20,2);
```

To run the interpolation program, use the command Interpolation. The user is prompted to enter initial points, function values, and partial derivatives.
```matlab
>> Interpolation
Enter initial points: E
Enter function values: F
Enter partial derivatives: D
```

The computer then performs the one-time work, computing Wells' cell complex. In this example, the number 95 indicates that there are 95 cells in the final cell complex.
```matlab
Computing Gamma
Elapsed time is 0.007851 seconds.
Gamma = 2784.4731
Finding power diagram
Elapsed time is 0.010667 seconds.
    95
Finding points on infinite edges of PD
Finding pieces S of cell complex
dim(S) = 0
dim(S) = 1:d-1
     1
dim(S) = d
Elapsed time is 0.172720 seconds.
```

Next, the user is prompted to enter a query point, x. The value of the interpolant at x is displayed, along with the gradient at x. 
```matlab
Enter query point x: [0,0]
Performing query work
Elapsed time is 0.000993 seconds.
 
F(x) = 222.88
DxF = -949.2244     -581.7292
 
Another? (1:Yes, 0:No) 0
```

Alternately, to compute the interpolant for many query points (stored in a matrix X), do the following:
```matlab
>> [g, sE, P, PD, Sc, dSc, Cells] = cellcomplexdD(E, F, D);
```

To check that the initial points are correctly interpolated (within the machine precision), let X = E.
```matlab
>> X = E;
>> [Fx, DxF, Index] = queryWorkdD(X, g, sE, P, PD, Sc, dSc, Cells);
>> F - Fx
ans =
            0
            0
   4.1633e-17
  -1.1102e-16
            0
  -1.1102e-16
            0
            0
  -5.5511e-17
  -1.1102e-16
   5.5511e-17
  -1.1102e-16
  -1.1102e-16
            0
            0
  -1.1102e-16
            0
            0
  -5.5511e-17
  -5.5511e-17
```

In this example, the initial points are in R^2, so we can run gridMesh to visualize the interpolant. gridMesh computes the interpolant on a grid and displays the following: a shaded surface plot of the interpolant colored by height, a mesh plot of the interpolant colored by cell, a surface plot of the x partial derivative colored by height, and a surface plot of the x partial derivative colored by cell.
```matlab
>> gridMesh
```
