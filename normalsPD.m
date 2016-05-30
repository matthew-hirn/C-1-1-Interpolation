function ind = normalsPD(LE, C)
% function ind = normalsPD(E, C)
%
% LE: set of lifted points
% C: convex hull of LE
%
% The output ind contains indices of facets of the lower hull.
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

[m, n] = size(C);
center = mean(LE,1);

for i=1:m
    v = null(bsxfun(@minus, LE(C(i,1),:), LE(C(i,2:end),:)))';
    [mm, nn] = size(v);
    if mm > 1
        V(i,:) = NaN;
        error('nullspace error')
        % possibility of degenerate null vectors
    else
        V(i,:) = v;
    end
    mid(i,:) = mean(LE(C(i,:),:),1);
end

dot = sum(bsxfun(@minus, center, mid).*V, 2);
outer = dot < 0;
V(outer,:) = -1*V(outer,:);

ind = V(:,n) > 0;