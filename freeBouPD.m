function FB = freeBouPD(T, P)
% function FB = freeBouPD(T, P)
%
% T: triangulation
% P: pieces of the triangulation
%
% The output FB contains the pieces of P on the boundary of T.
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

[~, nT] = size(T);
[mP, ~] = size(P);

ii=1;
for i=1:mP
    if size(find(sum(ismember(T, P(i,:)),2)==nT-1),1)==1
        FB(ii,:) = P(i,:);
        ii = ii+1;
    end
end