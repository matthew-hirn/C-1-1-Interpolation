# C-1-1-Interpolation
Computes optimal *C^{1,1}(R^d)* interpolations in any dimension *d*

## Overview

This software solves the following interpolation problem:
Let *E* be a finite set of *d*-dimensional points. To each point *a* in *E*, there is a specified scalar function value *f(a)* and a gradient vector *Df(a)*. The software computes a function *F:R^d -> R* such that:
  1. *F(a) = f(a)* and *∇F(a) = Df(a)* for each *a* in *E*
  2. Amonst all such interpolants satisfying 1., the value of of Lip(*∇F*) is minimal.
