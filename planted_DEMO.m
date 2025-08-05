%% Create adjacency matrix A using plantedsubmatrix function.

% Initialize problem sizes and densities
M = 100; %number of rows of sampled matrix
N = 200; %number of columnss of sampled matrix
m = 50; %number of rows of dense submatrix
n = 40; %number of columns of dense submatrix
p = 0.25; %noise density
q = 0.85; %in-group density

% Make binary matrix with planted mn-submatrix
[A,X0,Y0] = plantedsubmatrix(M,N,m,n,p,q);

%% Visualize the adjacency matrix and planted blocks.


figure; imagesc(A);  hold('on'); title('A'); hold('off')% plot matrix.
figure; imagesc(X0);  hold('on'); title('X0'); hold('off')
figure; imagesc(Y0);  hold('on'); title('Y0'); hold('off') % plot matrix representation of submatrix.

%% CALL DENSUB SOLVER.

% Initialize parameters and settings.
tau = 0.35; %regularization parameter
maxiter = 500; %max number of iterations
verbose = 1;
opt_tol = 1e-4; %optimal tolernce
gamma = 6/(sqrt(m*n)*(q-p)); %optimal choice of gamma from paper.

% Call solver.
[X,Y,Q, iter] = densub(A,m,n, gamma,tau, opt_tol, maxiter, verbose);

% Display iteration/convergence status
if iter < maxiter
    % Converged.
    fprintf('Algorithm converged after %d iterations.\n', iter) 
else
    % Failed to converge.
    fprintf('Algorithm failed to converge within %d iterations.\n', maxiter)
end

%% Plot results.
figure; imagesc(X); hold('on'); title('X'); hold('off')
figure; imagesc(Y); hold('on'); title('Y'); hold('off')
