function EA = edgeAttPD(T, edges)
% function EA = edgeAttPD(T, edges)
%
% T: pieces of a triangulation
% edges: array of edges in the triangulation
% 
% Each entry in the output cell EA corresponds to an edge from the input
% edges and contains all pieces of the triangulation attached to that edge.
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

[me, ne] = size(edges);
EA = cell(me,1);

for i=1:me
    EA{i} = find(sum(ismember(T, edges(i,:)), 2) == ne)';
end