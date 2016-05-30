function [P, total] = piecesPD(T)
% function [P, total] = piecesPD(T)
%
% T: triangulation
%
% The output P contains all pieces of the triangulation, from vertices
% (dimension zero) to fully-dimensional facets (dimension n-1).
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
P = cell(n,1);
total = 0;

P{1} = unique(T);
P{1} = P{1}(:);
total = total + size(P{1},1);

for i=2:n-1
    Q = [];
    for j=1:m
        Q = [Q; combnk(T(j,:),i)];
    end
    Q = sort(Q, 2);
    P{i} = unique(Q, 'rows');
    total = total + size(P{i},1);
    clear Q;
end

P{n} = T;
total = total + size(P{n},1);