function [LE, wts] = lift(sE, F, D, g)
% function [LE, wts] = lift(sE, F, D, g)
%
% sE: real Nxd matrix whose rows represent N shifted data points
% F: real Nx1 matrix of function values of f at each point in E
% D: real Nxd matrix whose rows contain the partial derivatives of f
% g: Le Gruyer's Gamma
%
% The output LE contains the points of sE lifted one dimension higher and
% wts contains the weight of each point.
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

[N, d] = size(sE);
wts = zeros(N,1);
LE = zeros(N,d+1);

for i=1:N
    wts(i,1) = (2/g^2)*norm(D(i,:))^2 - (4/g)*F(i);
    x = sE(i,:);
    LE(i,:) = [x, x*x' - wts(i)];
end