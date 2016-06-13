function [g, sE, P, PD, Sc, dSc, Cells, t] = cellcomplexdD(E, F, D)
% function [g, sE, P, PD, Sc, dSc, Cells, t] = cellcomplexdD(E, F, D)
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
% t: List of run times, in order they are outputted
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

% Initalize t
t = zeros(10,1);

% Dimension
dim = size(E,2);

% Compute Gamma^1
disp('Computing Gamma'); tic;
g = Gamma(E, F, D); t(1) = toc;
disp(['Elapsed time is ', num2str(t(1)), ' seconds.']);
disp(['Gamma = ', num2str(g)]);
fprintf('\n');

% Shift points and compute their center
disp('Shift points'); tic;
sE = shift(E, D, g);
center = mean(sE,1); t(2) = toc;
disp(['Elapsed time is ', num2str(t(2)), ' seconds.']);
fprintf('\n');

% Lift shifted points
disp('Lift points'); tic;
[LE, wts] = lift(sE, F, D, g); t(3) = toc;
disp(['Elapsed time is ', num2str(t(3)), ' seconds.']);
fprintf('\n');

% Compute convex hull of lifted points
disp('Compute convex hull of lifted points'); tic;
C = convhulln(LE); t(4) = toc;
disp(['Elapsed time is ', num2str(t(4)), ' seconds.']);
fprintf('\n');

% Determine lower hull and get triangulation
disp('Determine lower hull and get triangulation'); tic;
ind = normalsPD(LE, C); t(5) = toc;
disp(['Elapsed time is ', num2str(t(5)), ' seconds.']);
T = C(ind, :);
nT = size(T,2);
fprintf('\n');

% Compute all faces of triangulation
disp('Compute all faces of the triangulation'); tic;
[P, total, P_graph] = piecesPD(T); t(6) = toc;
disp(['Elapsed time is ', num2str(t(6)), ' seconds.']);
fprintf('\n');

% Find the facets on the boundary of the triangulation
disp('Determine which facets of the triangulation are on the exterior'); tic;
ind_ext = false(size(P{nT-1},1),1);
for i=1:length(ind_ext)
    if length(P_graph{nT-1}.parents{i}) == 1
        ind_ext(i) = true;
    end
end
FF = P{nT-1}(ind_ext,:); 
ind_ext = find(ind_ext == true); t(7) = toc;
disp(['Elapsed time is ', num2str(t(7)), ' seconds.']);
fprintf('\n');

% Compute power centers and then whole power diagram
disp('Finding power diagram'); tic;
[PC, powers] = powercenters(T, sE, wts);
PD = pwrDiagramPD(P, P_graph, PC); t(8) = toc;
disp(['Elapsed time is ', num2str(t(8)), ' seconds.']);
disp(['Total number of faces: ', num2str(total)]);
fprintf('\n');

% Find points on the infinite edges of the power diagram
disp('Finding points on infinite edges of PD'); tic;
PDinf = zeros(size(FF));
max_length = max(sqrt(sum(bsxfun(@minus, PC, center).^2,2))); % Find distance from center to farthest powercenter

for i=1:size(FF,1)
    facet = sE(FF(i,:),:);
    ea = P_graph{nT-1}.parents{ind_ext(i)};
    pc = PC(ea,:);
    ct = mean(facet,1);
    
    % find vector normal to the facet
    v = null(bsxfun(@minus, facet(1,:), facet(2:end,:)))';
    
    % reorient v to point outward
    if dot(center - ct, v) > 0
        v = -1*v;
    end
    
    % scale v to ensure newpt is sufficiently far away
    v = max_length*v;
    
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
t(9) = toc;
disp(['Elapsed time is ', num2str(t(9)), ' seconds.']);
fprintf('\n');

% --- Now starting to get the d-dimensional T_S sets that partition R^d ---

% Initialize to compute TS partition
Sc = cell(nT,1);
dSc = cell(nT,1);
Regions = cell(nT,1);
Cells = cell(total,1);

% T_S cells for subsets S with 1 point
disp('Finding cells T_S of cell complex'); tic;
disp('num(S) = 1');
for i=1:size(P{1},1)
    
    % Triangulation point (tp) and power diagram points (pdp)
    tp = sE(P{1}(i),:);
    pdp = PD{1}{i};
    
    % Sc point and distance d_S(S_c)
    Sc{1}{i,1} = tp;
    dSc{1}{i,1} = (-g/4)*wts(i);
    
    % Get all possible corners and get the mean vector
    all_U = 0.5 * bsxfun(@plus, pdp, tp);
    center = mean(all_U);
    
    % Get hyperplanes from facets of power diagram with shifted point of E
    % (parents of triangulation vertices are children of power diagram cells)
    num_children = length(P_graph{1}.parents{i});
    A = zeros(num_children, dim);
    b = ones(num_children, 1);
    for j=1:num_children
        child_pts = PD{2}{P_graph{1}.parents{i}(j)};
        U = 0.5 * bsxfun(@plus, child_pts, tp);
        U = bsxfun(@minus, U, center);
        [x, ~] = linsolve(U, ones(size(U,1),1));
        A(j,:) = x;
        b(j) = b(j) + A(j,:) * center';
    end
    Regions{1}{i,1} = [A, b];
end

% Now T_S cells for subsets S with 2 to d points
for i=2:size(P,1)-1
    disp(['num(S) = ', num2str(i)])
    for j=1:size(P{i},1)
        
        % Triangulation points (tp) and power diagram points (pdp)
        tp = sE(P{i}(j,:),:);
        pdp = PD{i}{j};
          
        % S_c point
        V = bsxfun(@minus, tp(1,:), tp(2:end,:))';
        W = bsxfun(@minus, pdp(1,:), pdp(2:end,:))';
        Asc = [V, -W];
        Bsc = (pdp(1,:) - tp(1,:))';
        ab = Asc \ Bsc;
        alpha = ab(1:i-1);
        sc = tp(1,:) + sum(bsxfun(@times, alpha', V),2)';
        Sc{i}{j,1} = sc;
        
        % Distance d_S(S_c)
        dist = zeros(size(tp,1),1);
        for k=1:size(tp,1)
            dist(k) = d(tp(k,:), sc, F(P{i}(j,k)), D(P{i}(j,k),:), g);
        end
        dSc{i}{j,1} = min(dist);
        
        % Get all possible corners and get the mean vector
        all_U = [];
        for k=1:size(tp,1)
            all_U = cat(1, all_U, 0.5 * bsxfun(@plus, pdp, tp(k,:)));
        end
        center = mean(all_U);
        
        % Get hyperplanes, first S_hat (tp) with facets of S_star (pdp) 
        % (parents of triangulation faces are children of power diagram faces)
        num_children = length(P_graph{i}.parents{j});
        A1 = zeros(num_children, dim);
        b1 = ones(num_children, 1);
        for k=1:num_children
            child_pts = PD{i+1}{P_graph{i}.parents{j}(k)};
            U = [];
            for l=1:size(tp,1)
                U = cat(1, U, 0.5 * bsxfun(@plus, child_pts, tp(l,:)));
            end
            U = bsxfun(@minus, U, center);
            [x, ~] = linsolve(U, ones(size(U,1),1));
            A1(k,:) = x;
            b1(k) = b1(k) + A1(k,:) * center';
        end
        
        % Get remaining hyperplanes, now S_star (pdp) with facets of S_hat (tp)
        num_children = length(P_graph{i}.children{j});
        A2 = zeros(num_children, dim);
        b2 = ones(num_children, 1);
        for k=1:num_children
            child_pts = P{i-1}(P_graph{i}.children{j}(k), :);
            child_pts = sE(child_pts, :);
            U = [];
            for l=1:size(child_pts,1)
                U = cat(1, U, 0.5 * bsxfun(@plus, pdp, child_pts(l,:)));
            end
            U = bsxfun(@minus, U, center);
            [x, ~] = linsolve(U, ones(size(U,1),1));
            A2(k,:) = x;
            b2(k) = b2(k) + A2(k,:) * center';
        end
        
        % Put the two sets of hyperplanes together and store away
        A = cat(1, A1, A2);
        b = cat(1, b1, b2);
        Regions{i}{j,1} = [A, b];
    end
end

% Finall T_S cells for subsets S with d+1 points
disp(['num(S) = ', num2str(dim+1)]);
for i=1:size(P{nT},1)
    
    % Triangulation points (tp) and power diagram points (pdp)
    tp = sE(P{nT}(i,:),:);
    pdp = PC(i,:);
    
    % Sc point and distance d_S(S_c)
    Sc{nT}{i,1} = pdp;
    dSc{nT}{i,1} = (g/4)*powers(i);
    
    % Get all possible corners and get the mean vector
    all_U = 0.5 * bsxfun(@plus, pdp, tp);
    center = mean(all_U);
    
    % Get hyperplanes from facets of triangulation with power center
    num_children = length(P_graph{nT}.children{i});
    A = zeros(num_children, dim);
    b = ones(num_children, 1);
    for j=1:num_children
        child_pts = P{nT-1}(P_graph{nT}.children{i}(j), :);
        child_pts = sE(child_pts, :);
        U = 0.5 * bsxfun(@plus, child_pts, pdp);
        U = bsxfun(@minus, U, center);
        [x, ~] = linsolve(U, ones(size(U,1),1));
        A(j,:) = x;
        b(j) = b(j) + A(j,:) * center';
    end
    Regions{nT}{i,1} = [A, b];
    
end

% Collect Regions into one list Cells
ii=1;
for i=1:size(Regions,1)
    for j=1:size(Regions{i},1)
        Cells{ii,1} = Regions{i}{j};
        ii = ii+1;
    end
end
t(10) = toc;
disp(['Elapsed time is ', num2str(t(10)), ' seconds.']);