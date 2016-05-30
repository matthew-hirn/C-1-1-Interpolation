%% GRIDMESH
%
% Generates a grid of query points to do mesh plot of interpolant.
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

%% Grid

x = 0:.005:1;
y = 0:.005:1;

%% Preprocessing

Nx = length(x);
Ny = length(y);

Fx = zeros(Nx,Ny);
DxF = zeros(Nx,Ny,2);
Index = zeros(Nx,Ny);

%% Loop through each point on the grid

for i=1:Nx
    for j=1:Ny
        [Fx(i,j), DxF(i,j,:), Index(i,j)] = queryWorkdD([x(i),y(j)], g, sE, P, PD, Sc, dSc, Cells);
    end
end

%% Display meshes

% Color by height
%figure, surf(x,y,Fx, sqrt(Fx));
figure, surf(x,y,Fx, Fx.^(1/3));

% Color by cell
figure, mesh(x,y,Fx,Index);

figure, surf(x, y, DxF(:,:,1))
figure, surf(x, y, DxF(:,:,1), Index)