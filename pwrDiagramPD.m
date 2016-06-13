function [PD, EA] = pwrDiagramPD(P, P_graph, PC)
% function PD = pwrDiagramPD(P, P_graph, PC)
%
% P: all faces of a triangulation
% P_graph: parents/children of each face in the triangulation
% PC: power centers of triangles in T
%
% The output cell PD contains pieces of the power diagram, indexed by
% dimension. PD{mP} contains pieces of dimension zero (power centers that
% correspond to fully-dimensional pieces of the triangulation T) and PD{1}
% contains fully-dimensional regions of the power diagram (corresponding to
% vertices of the triangulation).
%
% Each entry in the output cell EA corresponds to an edge from the input
% edges and contains all pieces of the triangulation attached to that edge.
%
% This file is part of the C^{1,1}(R^d) Interpolation software package.
%
% Author: Frederick McCollum and Matthew Hirn
% Email: frederick.mccollum@nyu.edu, mhirn@msu.edu
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

% Initialize
mP = size(P,1);
PD = cell(mP,1);
EA = cell(mP,1);

% Power centers
EA{mP} = cell(size(PC,1),1);
PD{mP} = cell(size(PC,1),1);
for i=1:size(PC,1)
    EA{mP}{i} = i;
    PD{mP}{i} = PC(i,:);
end

% Now loop over higher dimensional faces of the power diagram
for j=1:mP-1
    EA{mP-j} = cell(size(P{mP-j},1),1);
    PD{mP-j} = cell(size(P{mP-j},1),1);
    for i=1:size(P{mP-j},1)
        parents = P_graph{mP-j}.parents{i};
        ea = EA{mP-j+1}(parents);
        ea = cat(1, ea{:});
        ea = unique(ea);
        EA{mP-j}{i} = ea;
        PD{mP-j}{i} = PC(ea,:);
    end
end

end