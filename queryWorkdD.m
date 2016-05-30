function [Fx, DxF, Index] = queryWorkdD(X, g, sE, P, PD, Sc, dSc, Cells)
% function [Fx, DxF, Index] = queryWorkdD(X, g, sE, P, PD, Sc, dSc, Cells)
% 
% X: real Nxd matrix whose rows represent N query points
% g: Le Gruyer's Gamma
% sE: real Nxd matrix whose rows represent N shifted data points
% P: pieces of the triangulation
% PD: pieces of the power diagram
% Sc: Wells' S_c points of intersection
% dSc: Wells' distance to S_c points
% Cells: final cell complex
%
% Fx: values of the interpolant at each query point in X
% DxF: values of the gradient at each query point in X
% Index: index of the cell containing each query point
%
% This file is part of the C^{1,1}(R^d) Interpolation software package.
%
% Author: Frederick McCollum and Matthew Hirn
% Email: frederick.mccollum@nyu.edu, mhirn@msu.edu
%
% Copyright 2016 Frederick McCollum, Matthew Hirn
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


[mX, nX] = size(X);
mCells = size(Cells,1);
Fx = zeros(mX,1);
DxF = zeros(mX,nX);
Index = zeros(mX,1);
test = NaN(mX, mCells);

for i=1:mX
        
    x = X(i,:);
    index = 0;
    
    for j=1:mCells
        A = Cells{j}(:,1:nX);
        b = Cells{j}(:,nX+1);
        ineq = b - A*x';
        tst = min(ineq);
        test(i,j) = tst;
        if all(ineq >= 0)
            index = j;
            break
        elseif tst > -1e-10
            disp(['Point ',num2str(i),' is near region ',num2str(j)])
        end
    end
    
    if index == 0
        disp(['Point not found: ', num2str(i)])
        [~, index] = max(test(i,:));
        disp(['Seems to be in region ', num2str(index)])
    end
    
    Index(i,1) = index;
    
    for j=1:size(P,1)
        if index <= size(P{j},1)
            tp = sE(P{j}(index,:),:);
            pdp = PD{j}{index};
            sc = Sc{j}{index};
            dsc = dSc{j}{index};
            break
        else
            index = index - size(P{j},1);
        end
    end
    
    mtp = size(tp,1);
    
    if mtp == 1
        y = sc;
        z = x + (x - sc);
    elseif mtp == size(P,1)
        y = x + (x - pdp(1,:));
        z = pdp(1,:);
    else
        V = bsxfun(@minus, tp(1,:), tp(2:end,:));
        W = bsxfun(@minus, pdp(1,:), pdp(2:end,:));
        
        % Orthogonalize V (to make computations more stable)
        r = rank(V);
        [Q,~,~] = qr(V',0);
        V = Q(:,1:r)';
        
        % need to extract linearly independent rows of W
        % need this many: nX - size(V) + 1
        % Orthogonalize W (to make computations more stable)
        r = rank(W);
        [Q,~,~] = qr(W',0);
        W = Q(:,1:r)';
     
        Ahy = zeros(size(V,1));
        bhy = zeros(size(V,1),1);
        for j=1:size(V,1)
            Ahy(j,:) = dot(V', repmat(V(j,:), size(V,1), 1)');
            bhy(j) = dot(V(j,:), x - sc);
        end
        alpha = Ahy \ bhy;
        hy = sum(bsxfun(@times, alpha, V), 1) + sc;
        y = hy + (hy - sc);
        
        Ahz = zeros(size(W,1));
        bhz = zeros(size(W,1),1);
        for j=1:size(W,1)
            Ahz(j,:) = dot(W', repmat(W(j,:), size(W,1), 1)');
            bhz(j,1) = dot(W(j,:), x - pdp(1,:));
        end
        beta = Ahz \ bhz;
        hz = sum(bsxfun(@times, beta, W), 1) + pdp(1,:);
        z = hz + (hz - sc);
    end
    
    Fx(i,1) = dsc + (g/8)*norm(z-sc)^2 - (g/8)*norm(y-sc)^2;
    DxF(i,:) = (g/2)*(z-y);
end