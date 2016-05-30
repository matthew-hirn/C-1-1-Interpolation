function [A, b, C, flag] = constraints5(V)
% function [A, b, C, flag] = constraints5(V)
%
% V: set of points
%
% C: convex hull of V
% A,b: system of linear inequalities describing the region contained by C,
% i.e., C contains all points x such that Ax<=b
% flag: indicates that convhulln encountered errors
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

n = size(V,2);  
e = ones(n,1);
flag = 0;
try
    C = convhulln(V);
catch
    flag = 1;
    C = convhulln(V, {'Qx', 'QJ', 'QbB'});
end

center = mean(V(unique(C),:));
V = bsxfun(@minus, V, center);
A = zeros(size(C,1),size(V,2));
ii=0;
index = true(size(C,1),1);
for i = 1:size(C,1)
    [x,r] = linsolve(V(C(i,:),:),e);
    if r > 1e-13  %ad hoc choice for rcond threshold
        ii = ii+1;
        A(ii,:) = x;
    else
        %disp('rank error');
        index(i) = false;
    end
end
A = A(1:ii,:);
b = ones(size(A,1),1);
b = b + A*center';
C = C(index,:);
