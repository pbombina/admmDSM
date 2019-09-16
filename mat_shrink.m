function Z = mat_shrink(Z,tau)
% MAT_SHRINK singular value soft-thresholding for nuclear norm prox fxn.
%
% INPUT:
% Z - matrix to have singular values thresholded.
% tau - threshold.
% OUTPUT:
% Z - matrix following soft-thresholding.

% Get dimensions of Z.
[r,c] = size(Z);

% Take SVD of Z.
[U,S,V] = svd(Z); 

% Soft threshold singular values.
s = max(diag(S)-tau,0);

% Reconstitute Z.
if r < c
    Z = U*diag(s)*V(:,1:r)';
else    
    Z = U(:, 1:c)*diag(s)*V';
end
