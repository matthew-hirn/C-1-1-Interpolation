function [PC, powers] = powercenters(T, sE, wts)
% function [PC, powers] = powercenters(T, sE, wts)
%
% T: triangulation
% sE: set of shifted points
% wts: weights of points in sE
%
% The output array PC contains the power centers of the triangles in T.
% That is, row i of PC contains the point that is equidistant from each
% vertex of the triangle specified by row i of T, with respect to the power
% of each vertex. The output array powers contains the power of each power
% center.
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

[m, n] = size(T);
PC = zeros(m,n-1);
powers = zeros(m,1);

for i=1:m
    tr = sE(T(i,:),:);
    wt = wts(T(i,:));
    p = tr(1,:);
    Rp = repmat(p, n-1, 1);
    Pts = tr(2:n,:);
    Ac = 2*(Pts - Rp);
    
    Rw1 = repmat(wt(1), n-1, 1);
    Wts = wt(2:n);
    Sp1 = repmat(sum(p.^2), n-1, 1);
    SPts = sum(Pts.^2, 2);
    Bc = Rw1 - Wts - Sp1 + SPts;
    
    pc = Ac \ Bc;
    
    PC(i,:) = pc;
    powers(i,1) = norm(pc - p')^2 - wt(1);
end