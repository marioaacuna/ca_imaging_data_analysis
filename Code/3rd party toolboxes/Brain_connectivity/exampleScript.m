% This script provides examples of how to compute characteristic path
% length, global efficiency, open and closed versions of local efficiency, 
% and open and closed versions of clustering coefficient for a graph given 
% its adjacency matrix.
%
% The examples are taken from the papers:
%
% A) "Measuring Brain Connectivity," by N. D. Cahill, J. Lind, and D. A. 
% Narayan, Bulletin of the Institute of Combinatorics & Its Applications, 
% 69, pp. 68-78, September 2013.
%
% B) "On Local Efficiency and Clustering Coefficients of Graphs," by B. Ek,
% C. VerSchneider, N. D. Cahill, and D. A. Narayan, submitted. 
%
% These examples use the companion function graphProperties.m.
%

%% Paper A, Exercise 1
%

% construct adjacency matrix
A = [0 1 0 0 1; 1 0 1 1 1; 0 1 0 1 1; 0 1 1 0 1; 1 1 1 1 0];

% compute graph properties
[L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = graphProperties(A);

% display results
fprintf('\nGraph Properties\n');
fprintf('\tCharacteristic Path Length:\t\t\t\t\t%6.4f\n',L);
fprintf('\tGlobal Efficiency:\t\t\t\t\t\t\t%6.4f\n',EGlob);
fprintf('\tClustering Coefficient (closed / open):\t\t%6.4f / %6.4f\n',CClosed,COpen);
fprintf('\tLocal Efficiency (closed / open):\t\t\t%6.4f / %6.4f\n',ELocClosed,ELocOpen);

%% Both papers, Macaque exercise
%
% First, download the file "macaque47.mat" from the site:
%   https://sites.google.com/site/bctnet/datasets
% Save the "macaque47.mat" file into your current working directory (or
% another directory on the MATLAB search path).

% load macaque data file
load('macaque47.mat');

% adjacency matrix is stored in CIJ
% compute graph properties
[L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = graphProperties(CIJ);

% display results
fprintf('\nGraph Properties\n');
fprintf('\tCharacteristic Path Length:\t\t\t\t\t%6.4f\n',L);
fprintf('\tGlobal Efficiency:\t\t\t\t\t\t\t%6.4f\n',EGlob);
fprintf('\tClustering Coefficient (closed / open):\t\t%6.4f / %6.4f\n',CClosed,COpen);
fprintf('\tLocal Efficiency (closed / open):\t\t\t%6.4f / %6.4f\n',ELocClosed,ELocOpen);

%% Paper B, Tripartite graph with nine vertices 
%

% construct adjacency matrix - vertices are:
% v11 v12 v13 v21 v22 v23 v31 v32 v33
A = [zeros(3,3) ones(3,6);
    ones(3,3) zeros(3,3) ones(3,3);
    ones(3,6) zeros(3,3)];

% compute graph properties
[L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = graphProperties(A);

% display results
fprintf('\nGraph Properties\n');
fprintf('\tCharacteristic Path Length:\t\t\t\t\t%6.4f\n',L);
fprintf('\tGlobal Efficiency:\t\t\t\t\t\t\t%6.4f\n',EGlob);
fprintf('\tClustering Coefficient (closed / open):\t\t%6.4f / %6.4f\n',CClosed,COpen);
fprintf('\tLocal Efficiency (closed / open):\t\t\t%6.4f / %6.4f\n',ELocClosed,ELocOpen);

%% Paper B, graph with five vertices 
%

% construct adjacency matrix - vertices are:
% u v w x y
A = [0 1 0 1 1;
    1 0 1 0 1;
    0 1 0 1 1;
    1 0 1 0 1;
    1 1 1 1 0];

% compute graph properties
[L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = graphProperties(A);

% display results
fprintf('\nGraph Properties\n');
fprintf('\tCharacteristic Path Length:\t\t\t\t\t%6.4f\n',L);
fprintf('\tGlobal Efficiency:\t\t\t\t\t\t\t%6.4f\n',EGlob);
fprintf('\tClustering Coefficient (closed / open):\t\t%6.4f / %6.4f\n',CClosed,COpen);
fprintf('\tLocal Efficiency (closed / open):\t\t\t%6.4f / %6.4f\n',ELocClosed,ELocOpen);

%% Paper B, K_{m,n} (complete bipartite graph) 
%

% choose example values of m,n
m = 5;
n = 6;

% construct adjacency matrix
A = [zeros(m,m) ones(m,n);
    ones(n,m) zeros(n,n)];

% compute graph properties
[L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = graphProperties(A);

% display results
fprintf('\nGraph Properties\n');
fprintf('\tCharacteristic Path Length:\t\t\t\t\t%6.4f\n',L);
fprintf('\tGlobal Efficiency:\t\t\t\t\t\t\t%6.4f\n',EGlob);
fprintf('\tClustering Coefficient (closed / open):\t\t%6.4f / %6.4f\n',CClosed,COpen);
fprintf('\tLocal Efficiency (closed / open):\t\t\t%6.4f / %6.4f\n',ELocClosed,ELocOpen);

%% Paper B, K_{r;n} (complete multipartite graph) 
%

% choose example values of r,n
r = 4;
n = 5;

% construct adjacency matrix
A = ones(r*n,r*n);
for i = 0:(n-1)
    A(r*i+(1:r),r*i+(1:r)) = 0;
end

% compute graph properties
[L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = graphProperties(A);

% display results
fprintf('\nGraph Properties\n');
fprintf('\tCharacteristic Path Length:\t\t\t\t\t%6.4f\n',L);
fprintf('\tGlobal Efficiency:\t\t\t\t\t\t\t%6.4f\n',EGlob);
fprintf('\tClustering Coefficient (closed / open):\t\t%6.4f / %6.4f\n',CClosed,COpen);
fprintf('\tLocal Efficiency (closed / open):\t\t\t%6.4f / %6.4f\n',ELocClosed,ELocOpen);

%% Paper B, C_{n} (cycle graph) 
%

% choose example value of n
n = 13;

% construct adjacency matrix
A = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
A(n,1) = 1;
A(1,n) = 1;

% compute graph properties
[L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = graphProperties(A);

% display results
fprintf('\nGraph Properties\n');
fprintf('\tCharacteristic Path Length:\t\t\t\t\t%6.4f\n',L);
fprintf('\tGlobal Efficiency:\t\t\t\t\t\t\t%6.4f\n',EGlob);
fprintf('\tClustering Coefficient (closed / open):\t\t%6.4f / %6.4f\n',CClosed,COpen);
fprintf('\tLocal Efficiency (closed / open):\t\t\t%6.4f / %6.4f\n',ELocClosed,ELocOpen);

%% Paper B, W_{n} (wheel graph) 
%

% choose example value of n
n = 5;

% construct adjacency matrix
A = zeros(n,n);
A(1:n-1,1:n-1) = diag(ones(n-2,1),1) + diag(ones(n-2,1),-1);
A(n-1,1) = 1;
A(1,n-1) = 1;
A(n,1:n-1) = 1;
A(1:n-1,n) = 1;

% compute graph properties
[L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = graphProperties(A);

% display results
fprintf('\nGraph Properties\n');
fprintf('\tCharacteristic Path Length:\t\t\t\t\t%6.4f\n',L);
fprintf('\tGlobal Efficiency:\t\t\t\t\t\t\t%6.4f\n',EGlob);
fprintf('\tClustering Coefficient (closed / open):\t\t%6.4f / %6.4f\n',CClosed,COpen);
fprintf('\tLocal Efficiency (closed / open):\t\t\t%6.4f / %6.4f\n',ELocClosed,ELocOpen);
