%% INTERPOLATION
%
% Script that runs user interface for entering initial data (initial
% points, functions values, and gradients). It then performs the one-time
% work for computing the interpolant, and finally asks the user to enter
% any number of query points.
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

E = input('Enter initial points: ');
F = input('Enter function values: ');
D = input('Enter partial derivatives: ');
[g, sE, P, PD, Sc, dSc, Cells] = cellcomplexdD(E, F, D);
X = input('Enter query point x: ');
disp('Performing query work');
tic
[Fx, DxF, Index] = queryWorkdD(X, g, sE, P, PD, Sc, dSc, Cells);
toc
disp(' ');
disp(['F(x) = ',num2str(Fx)]);
disp(['DxF = ',num2str(DxF)]);
disp(' ');
c = input('Another? (1:Yes, 0:No) ');
while c==1
    X = input('Enter query point x: ');
    disp('Performing query work');
    tic
    [Fx, DxF, Index] = queryWorkdD(X, g, sE, P, PD, Sc, dSc, Cells);
    toc
    disp(' ');
    disp(['F(x) = ',num2str(Fx)]);
    disp(['DxF = ',num2str(DxF)]);
    disp(' ');
    c = input('Another? (1:Yes, 0:No) ');
end