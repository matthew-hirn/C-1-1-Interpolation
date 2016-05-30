function dist = d(sx, a, f, Dxf, g)
% function dist = d(sx, a, f, Dxf, g)
%
% sx: a shifted point in sE
% a: a point in R^d
% f: the function value at x
% DxF: the gradient at x
% g: Le Gruyer's Gamma
%
% dist: Wells' distance from sx to a
%
% This file is part of the C^{1,1}(R^d) Interpolation software package.
%
% Author: Frederick McCollum
% Email: frederick.mccollum@nyu.edu
%
% Copyright 2016 Frederick McCollum
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

dist = f - (norm(Dxf)^2)/(2*g) + (g/4)*norm(a - sx)^2;