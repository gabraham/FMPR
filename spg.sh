#!/bin/bash

cat > spg.m <<EOF

addpath('~/Software/SPG_Multi_Graph');

rep = $1;
nfolds = $2;
r = $3;
dir = '$4';

option.cortype = 2;
option.maxiter = 1e6;
option.tol = 1e-8;  
option.verbose = false;
option.display_iter = 500;
option.mu = 1e-5;
option.corthreshold = $5;
option.cortype = $6;

X = dlmread(sprintf('%s/Xtrain_%d.txt', dir, rep));
Y = dlmread(sprintf('%s/Ytrain_%d.txt', dir, rep));

[N, J] = size(X);
K = size(Y, 2);

%C = sparse(Ynetwork);
%CNorm = 2 * max(sum(C.^2, 1));

[C, CNorm, E, Ecoef, Esign, R] = gennetwork(Y, option); 

% find best hyperparameters using cross-validation, as measure by mean-squared
% error
folds = ceil(rand(N, 1) * nfolds);
gamma = linspace(1e-6, 10, r);
lambda = linspace(1e-6, 10, r);
err = zeros(r, r, nfolds);
errl = zeros(r, nfolds);

for fold = 1:nfolds
   fprintf('SPG fold %d\n', fold);
   Ytrain = Y(folds ~= fold, :);
   Xtrain = X(folds ~= fold, :);
   Ytest = Y(folds == fold, :);
   Xtest = X(folds == fold, :);
   
   for j = 1:r
      % lasso only
      [beta, obj, density, iter, time] = SPG_multi(...
	 Ytrain, Xtrain, 0, lambda(j), C, CNorm, option);
      errl(j, fold) = sum(sum((Xtest * beta - Ytest).^2)) / (N * K);

      % Fused-lasso
      for i = 1:r
	 [beta, obj, density, iter, time] = SPG_multi(...
	    Ytrain, Xtrain, gamma(i), lambda(j), C, CNorm, option);
	 err(i, j, fold) = sum(sum((Xtest * beta - Ytest).^2)) / (N * K);
      end
   end
end

serror = sum(err, 3);
serrorl = sum(errl, 3);

[i, j] = ind2sub(size(serror), find(serror == min(serror(:))));
[il, jl] = ind2sub(size(serrorl), find(serrorl == min(serrorl(:))));

% in case R^2 is flat and a vector is returned
i = i(1);
j = j(1);
il = il(1);
jl = jl(1);

% train GFlasso on entire dataset
[beta, obj, density, iter, time] = SPG_multi(...
      Y, X, gamma(i), lambda(j), ...
      C, CNorm, option);
save(sprintf('%s/spg_beta_%d.txt', dir, rep), 'beta', '-ascii');

% train lasso on entire dataset, ie SPG with gamma=0
[beta, obj, density, iter, time] = SPG_multi(...
      Y, X, 0, lambda(jl), ...
      C, CNorm, option);
save(sprintf('%s/spg_beta_lasso_%d.txt', dir, rep), 'beta', '-ascii');

quit;

EOF

matlab -nodisplay -r 'spg'

