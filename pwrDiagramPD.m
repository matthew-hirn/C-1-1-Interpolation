function PD = pwrDiagramPD(T, PC)
% function PD = pwrDiagramPD(T, PC)
%
% T: triangulation
% PC: power centers of triangles in T
%
% The output cell PD contains pieces of the power diagram, indexed by
% dimension. PD{mP} contains pieces of dimension zero (power centers that
% correspond to fully-dimensional pieces of the triangulation T) and PD{1}
% contains fully-dimensional regions of the power diagram (corresponding to
% vertices of the triangulation).
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

P = piecesPD(T);
mP = size(P,1);
PD = cell(mP,1);

for i=1:mP
    EA = edgeAttPD(T, P{i});
    for j=1:size(EA,1);
        EA{j} = PC(EA{j},:);
    end
    PD{i} = EA;
end