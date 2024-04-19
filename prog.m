clear all
close all
clc
format long
rng(302699)

%% PARAMETER DEFINITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataset = 'circle';
%dataset = 'spiral';
%dataset = 'landmines';

normlapl = 'unnorm';
%normlapl = 'symnorm';

%solvels = 'backslash';
solvels = 'conjgrad';
%solvels = 'gmres';
%solvels = 'lanczos';

clustmeth = 'kmeans';
%clustmeth = 'kmedoids';
%clustmeth = 'ward';
%clustmeth = 'euclidean';

sigma = 1;
k = 10;          % n. of neighbours
%k = 20;
%k = 40;
num = 10;         % maximum number of clusters
tol = 1.0e-10;   % tolerance
maxit = 1.0e03;  % maximum number of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[L, conncomp, eigenvals, U, M, IDX] = SpectralClustering(k, sigma, num, ...
    tol, maxit, dataset, normlapl, clustmeth, solvels);
