function [P, total, P_graph] = piecesPD(T)
% function [P, total, P_graph] = piecesPD(T)
%
% T: triangulation
%
% The output P contains all pieces of the triangulation, from vertices
% (dimension zero) to fully-dimensional facets (dimension n-1).
%
% The output P_graph lists the children/parents for each face.
%
% This file is part of the C^{1,1}(R^d) Interpolation software package.
%
% Author: Frederick McCollum and Matthew Hirn
% Email: frederick.mccollum@nyu.edu, mhirn@msu.edu
%
% Copyright 2016 Frederick McCollum and Matthew Hirn
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
P_graph = cell(n,1);

P{1} = unique(T);
P{1} = P{1}(:);
total = total + size(P{1},1);
P_graph{1}.children = [];

for i=2:n
    
    if i < n
        Q = [];
        for j=1:m
            Q = [Q; combnk(T(j,:),i)];
        end
        Q = sort(Q, 2);
        P{i} = unique(Q, 'rows');
        total = total + size(P{i},1);
        clear Q;
    else
        P{n} = T;
    end
    
    P_graph{i-1}.parents = cell(size(P{i-1},1),1);
    P_graph{i}.children = cell(size(P{i},1),1);
    for j=1:size(P{i-1},1)
        I = cell(i-1,1);
        for k=1:i-1
            [I{k},~] = find(P{i} == P{i-1}(j,k));
        end        
        II = I{1};
        for k=2:i-1
            II = intersect(II, I{k});
        end
        P_graph{i-1}.parents{j} = sort(II, 'ascend');
        for k=1:length(II)
            P_graph{i}.children{II(k)} = cat(1, P_graph{i}.children{II(k)}, j);
        end
    end
end
P_graph{n}.parents = [];

total = total + size(P{n},1);