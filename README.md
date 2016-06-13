# C-1-1-Interpolation
Computes optimal *C^{1,1}(R^d)* interpolations in any dimension *d*.

## Overview

This software solves the following interpolation problem:

Let *E* be a finite set of *d*-dimensional points. To each point *a* in *E*, there is a specified scalar function value *f(a)* and a gradient vector *Df(a)*. The software computes a function *F:R^d -> R* such that:
  1. *F(a) = f(a)* and *∇F(a) = Df(a)* for each *a* in *E*.
  2. Amongst all such interpolants satisfying 1., the value of of Lip(*∇F*) is minimal.

The software is based upon the paper "Computing minimal interpolants in C^{1,1}(R^d)," by Ariel Herbert-Voss, Matthew Hirn, and Frederick McCollum (https://arxiv.org/abs/1411.5668).

## Installation

The folder C-1-1-Interpolation can be saved to any location. In Matlab, add to the path the folder C-1-1-Interpolation.

Check the releases tab to get release versions of the software and to see major updates.

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

The computer then performs the one-time work, computing Wells' cell complex. In this example, the power diagram has 99 faces, meaning there are 99 cells in the final cell complex.
```matlab
Computing Gamma
Elapsed time is 0.0015426 seconds.
Gamma = 1451.8317

Shift points
Elapsed time is 1.9908e-05 seconds.

Lift points
Elapsed time is 0.00010622 seconds.

Compute convex hull of lifted points
Elapsed time is 0.00030862 seconds.

Determine lower hull and get triangulation
Elapsed time is 0.0014245 seconds.

Compute all faces of the triangulation
Elapsed time is 0.006688 seconds.

Determine which facets of the triangulation are on the exterior
Elapsed time is 0.0001249 seconds.

Finding power diagram
Elapsed time is 0.0031662 seconds.
Total number of faces: 95

Finding points on infinite edges of PD
Elapsed time is 0.0089396 seconds.

Finding cells T_S of cell complex
num(S) = 1
num(S) = 2
num(S) = 3
Elapsed time is 0.030267 seconds.
```

Next, the user is prompted to enter a query point, x. The value of the interpolant at x is displayed, along with the gradient at x. 
```matlab
Enter query point x: [0,0]
Performing query work
Elapsed time is 0.002037 seconds.
 
F(x) = 67.3251
DxF = -129.6883     -422.2365
 
Another? (1:Yes, 0:No) 0
```

Alternately, to compute the interpolant for many query points (stored in a matrix X), do the following:
```matlab
>> [g, sE, P, PD, Sc, dSc, Cells, t] = cellcomplexdD(E, F, D);
>> [Fx, DxF, Index] = queryWorkdD(X, g, sE, P, PD, Sc, dSc, Cells);
```

To check that the initial points are correctly interpolated (within the machine precision), let X = E.
```matlab
>> X = E;
>> [Fx, DxF, Index] = queryWorkdD(X, g, sE, P, PD, Sc, dSc, Cells);
>> F - Fx
ans =

   1.0e-15 *

         0
         0
    0.0555
    0.0555
   -0.1110
    0.1110
         0
    0.1110
         0
         0
    0.1110
    0.0278
    0.0278
         0
   -0.1110
         0
         0
   -0.0347
    0.1110
    0.0555
```

In this example, the initial points are in R^2, so we can run gridMesh to visualize the interpolant. gridMesh computes the interpolant on a grid and displays the following: a shaded surface plot of the interpolant colored by height, a mesh plot of the interpolant colored by cell, a surface plot of the x partial derivative colored by height, and a surface plot of the x partial derivative colored by cell.
```matlab
>> gridMesh
```

## Authors

The underlying mathematics for this software was developed by Ariel Herbert-Voss, Matthew Hirn, and Frederick McCollum.

The initial release of this software (v1.0) was written by Frederick McCollum under the supervision of Matthew Hirn. Subsequent releases were written by Matthew Hirn and Frederick McCollum. Correspondence should be sent to Matthew Hirn at mhirn@msu.edu.

## License

Copyright 2016 Ariel Herbert-Voss, Matthew Hirn, Frederick McCollum

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
