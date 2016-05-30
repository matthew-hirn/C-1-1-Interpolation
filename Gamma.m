function g = Gamma(E, F, D)
% function g = Gamma(E, F, D)
%
% E: real Nxd matrix whose rows represent N data points
% F: real Nx1 matrix of function values of f at each point in E
% D: real Nxd matrix whose rows contain the partial derivatives of f
% 
% The output g is LeGruyer's Gamma.
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

[N, d] = size(E);
g = 0;

for i=1:N
    for j=(i+1):N
        x = E(i,:);
        y = E(j,:);
        n = norm(x - y);
        
        A = [2*(F(i) - F(j)) + dot((D(i,:) + D(j,:)), (y - x))] / n^2;
        B = norm(D(i,:) - D(j,:)) / n;
        G = [A^2 + B^2]^.5 + abs(A);
        
        if G > g
            g = G;
        end
    end
end


        
    
