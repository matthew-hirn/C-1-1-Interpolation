function [g, sE, P, PD, Sc, dSc, Cells] = cellcomplexdD(E, F, D)
% function [g, sE, P, PD, Sc, dSc, Cells] = cellcomplexdD(E, F, D)
%
% E: real Nxd matrix whose rows represent N data points
% F: real Nx1 matrix of function values of f at each point in E
% D: real Nxd matrix whose rows contain the partial derivatives of f
% 
% g: Le Gruyer's Gamma
% sE: real Nxd matrix whose rows represent N shifted data points
% P: pieces of the triangulation
% PD: pieces of the power diagram
% Sc: Wells' S_c points of intersection
% dSc: Wells' distance to S_c points
% Cells: final cell complex
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


disp('Computing Gamma'); tic;
g = Gamma(E, F, D); toc
disp(['Gamma = ', num2str(g)]);
sE = shift(E, D, g);
[LE, wts] = lift(sE, F, D, g);
C = convhulln(LE);
ind = normalsPD(LE, C);
T = C(ind, :);
nT = size(T,2);
[P, total] = piecesPD(T);
[PC, powers] = powercenters(T, sE, wts);
FF = freeBouPD(T, P{nT-1});
disp('Finding power diagram'); tic;
PD = pwrDiagramPD(T, PC); toc
center = mean(sE,1);
disp(total)

Sc = cell(nT,1);
dSc = cell(nT,1);
PDinf = zeros(size(FF));
Regions = cell(nT,1);
Cells = cell(total,1);
GPTS = cell(nT,1);

% find distance from center to farthest powercenter
length = max(sqrt(sum(bsxfun(@minus, PC, center).^2,2)));

disp('Finding points on infinite edges of PD');  tic;
for i=1:size(FF,1)
    facet = sE(FF(i,:),:);
    ea = edgeAttPD(T, FF(i,:));
    pc = PC(ea{1},:);
    ct = mean(facet,1);
    
    % find vector normal to the facet
    v = null(bsxfun(@minus, facet(1,:), facet(2:end,:)))';
    
    % reorient v to point outward
    if dot(center - ct, v) > 0
        v = -1*v;
    end
    
    % scale v to ensure newpt is sufficiently far away
    v = length*v;
    
    % find point on infinite edge of power diagram
    newpt = pc + v;
    p = piecesPD(FF(i,:));
    for j=1:size(p,1)
        for k=1:size(p{j},1)
            ind = find(ismember(P{j},p{j}(k,:), 'rows'));
            PD{j}{ind} = [PD{j}{ind}; newpt];
        end
    end
    
    % keep track of generated point
    PDinf(i,:) = newpt;
end

disp('Finding pieces S of cell complex');
disp('dim(S) = 0');
for i=1:size(P{1},1)
    Sc{1}{i,1} = sE(i,:);
    dSc{1}{i,1} = (-g/4)*wts(i);
    
    gpts = bsxfun(@plus, PD{1}{i}, sE(P{1}(i),:))/2;
    gptsi = bsxfun(@plus, PDinf, sE(P{1}(i),:))/2;
    [A, b, C, ~] = constraints5(gpts);
    
    % remove boundary constraints of infinite regions
    indInf = find(ismember(gpts, gptsi, 'rows'));
    if size(indInf,1) > 0
        indf = zeros(1,1);
        jj=1;
        for j=1:size(C,1)
            if sum(ismember(indInf, C(j,:))) == nT-1
                indf(jj,1) = j;
                jj = jj+1;
            end
        end
        A(indf,:) = [];
        b(indf) = [];
    end

    Regions{1}{i,1} = [A, b];
    GPTS{1}{i,1} = gpts;
end

disp('dim(S) = 1:d-1');
for i=2:size(P,1)-1
    disp(i-1)
    for j=1:size(P{i},1)
        tp = sE(P{i}(j,:),:);
        pdp = PD{i}{j};
        
        % find S_c point and d_S(S_c)
        V = bsxfun(@minus, tp(1,:), tp(2:end,:))';
        W = bsxfun(@minus, pdp(1,:), pdp(2:end,:))';
        Asc = [V, -W];
        Bsc = (pdp(1,:) - tp(1,:))';
        ab = Asc \ Bsc;
        alpha = ab(1:i-1);
        sc = tp(1,:) + sum(bsxfun(@times, alpha', V),2)';
        Sc{i}{j,1} = sc;
        dist = zeros(size(tp,1),1);
        for k=1:size(tp,1)
            dist(k) = d(tp(k,:), sc, F(P{i}(j,k)), D(P{i}(j,k),:), g);
        end
        dSc{i}{j,1} = min(dist);
        
        % find system of inequalities describing cell
        gpts = zeros(size(tp,1)*size(pdp,1), nT-1);
        gptsi = zeros(size(tp,1)*size(PDinf,1), nT-1);
        for k=1:size(tp,1)
            gpts((k-1)*size(pdp,1)+1:k*size(pdp,1),:) = ...
                bsxfun(@plus, tp(k,:), pdp)/2;
            gptsi((k-1)*size(PDinf,1)+1:k*size(PDinf,1),:) = ...
                bsxfun(@plus, tp(k,:), PDinf)/2;
        end
        GPTS{i}{j,1} = gpts;
        [A, b, C, flag] = constraints5(gpts);
        if flag==1
            disp([num2str(i), ' ', num2str(j)])
        end
        
        % remove boundary constraints of infinite regions
        indInf = find(ismember(gpts, gptsi, 'rows'));
        if size(indInf,1) > 0
            indf = zeros(1,1);
            kk=1;
            for k=1:size(C,1)
                if sum(ismember(indInf, C(k,:))) == nT-1
                    indf(kk,1) = k;
                    kk = kk+1;
                end
            end
            A(indf,:) = [];
            b(indf) = [];
        end
        
        Regions{i}{j,1} = [A, b];
    end
end

disp('dim(S) = d');
for i=1:size(P{nT},1)
    Sc{nT}{i,1} = PC(i,:);
    dSc{nT}{i,1} = (g/4)*powers(i);
    
    gpts = bsxfun(@plus, sE(P{nT}(i,:),:), PC(i,:))/2;
    GPTS{nT}{i,1} = gpts;
    [A, b, ~, ~] = constraints5(gpts);
    Regions{nT}{i,1} = [A, b];
end

ii=1;
for i=1:size(Regions,1)
    for j=1:size(Regions{i},1)
        Cells{ii,1} = Regions{i}{j};
        ii = ii+1;
    end
end
toc