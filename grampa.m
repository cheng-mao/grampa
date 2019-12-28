% GRAph Matching by Pairwise eigen-Alignments (GRAMPA)
% by Zhou Fan (Yale), Cheng Mao (Georgia Tech), Yihong Wu (Yale), 
% and Jiaming Xu (Duke)

function [P] = grampa(A, B, eta)
% "grampa" matches the vertices of two unlabeled, edge-correlated graphs
%
% Input arguments:
% - A: the (centered or uncentered) adjacency matrix of the first graph
% - B: the (centered or uncentered) adjacency matrix of the second graph
% - eta: the regularization parameter (eta>0)
%
% Output arguments:
% - P: the permutation matrix P such that P*A*P' is matched to B

% Compute the similarity matrix
n = size(A, 1);
[U, Lambda] = eig(A);
[V, Mu] = eig(B);
lambda = diag(Lambda);
mu = diag(Mu);
coeff = 1 ./ ((lambda - mu').^2 + eta^2);
coeff = coeff .* (U' * ones(n) * V);
X = U * coeff * V';

% Rounding by linear assignment 
M = matchpairs(X', -99999, 'max');
P = full(sparse(M(:, 1), M(:, 2), 1, n, n));