function [X,Y,Q, iter] = densub(A,m,n, gamma,tau, opt_tol, maxiter, verbose)
% DENSUB ADMM algorithm for the densest submatrix and subgraph problems.
%
% INPUT:
% A - input matrix (or adjacency matrix of input graph).
% m,n - desired dimensions of dense submatrix/subgraph.
% gamma - regularization parameter in DKS, DSM problems.
% tau - augmented Lagrangian penalty parameter.
% opt_tol - stopping tolerance.
% verbose - indicates whether to display iteration statistics.
% maxiter - maximum number of iterations.
% OUTPUT:
% X,Y,Q - matrix representation of densest submatrix.
% iter - number of iterations
% REQUIRES: mat_shrink.m


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INITIALIZATION.
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Get dimensions of
[M,N]=size(A);

% Reciprocal of augmented Lagrangian parameter.
mu = 1/tau;

% Initial solutions.
% Set initial primal solutions to be scaled all-ones matrices.
W = ones(M,N)*m*n/(M*N);
X = W; 
Y = X;
Z = X;
Q = X - Y;
LambdaQ = zeros(size(X)); 
LambdaZ = zeros(size(X)); 
LambdaW = zeros(size(X)); 

% Initialize counters.
convergence=0;
iter=0; 

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% SET UP ITERATION STATISTICS TABLE.
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if verbose == 1
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    fprintf('DENSUB - ADMM for DSM/DKS. \n')    
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    fprintf('M=%d, N=%d, m=%d, n=%d, gamma = %1.3e \n', M, N, m,n, gamma)
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    fprintf('It \t\t | Primal gap \t | Dual gap \t \n' )
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    
end

 %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ITERATIVE METHOD.
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
while convergence==0  % Repeat until converged.
    
    % Increment iteration counter.
    iter=iter+1;

    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % UPDATE Q.
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    
    % Save old version and preprojected Q.
    Qold = Q;
    Q= X - Y + mu*LambdaQ;        
    
    % Project Q onto support of A.
    Q = Q.*A; 
    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % UPDATE X.    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    % Perform matrix shrink.
    X = mat_shrink(1/3*(Y + Q + Z + W - mu*(LambdaQ + LambdaW + LambdaZ)), 1/(3*tau));
    
    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % UPDATE Y.    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    % Update Y as projection of residual onto nonnegative cone.
    Y = max(X-Q-gamma*ones(M,N)*mu + LambdaQ*mu, zeros(M,N));
    
    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % UPDATE W.    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    % Initialize new W value.
    Wold = W;
    newW = X + mu*LambdaW;
    
    % Scale/shift W so that entries sum to m*n.
    alfa = (m*n-sum(newW(:)))/(M*N);
    W = newW + alfa*ones(M,N);
    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % UPDATE Z.    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Zold = Z;
    Z = X+ mu*LambdaZ; 
    Z = min(max(Z,0),1);
    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % UPDATE DUAL VARIABLES BY (APPROXIMATE) DUAL ASCENT STEP.    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    LambdaQ = LambdaQ + tau*(X-Y-Q);  
    LambdaW = LambdaW + tau*(X-W);
    LambdaZ = LambdaZ + tau*(X-Z);

    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % CHECK FOR CONVERGENCE.    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % Calculate primal residuals.
    NZ = norm(X-Z,'fro');
    NW = norm(X-W,'fro');
    NQ = norm(X-Y-Q, 'fro');
    
    % Maximum relative primal residual.
    errP = max([NZ, NW, NQ])/norm(X,'fro');
    
    % Calculate dual residuals.
    NDz = norm(Z - Zold, 'fro');
    NDw = norm((W-Wold), 'fro');
    NDp = norm((Q-Qold), 'fro');
    
    % Maximum relative dual residual.
    errD = max([NDz, NDw, NDp])/norm(X, 'fro');
    
    % Check for convergence.
    if errP < opt_tol && errD < opt_tol
        convergence = 1;
    end
    
    % Check if we have exceeded maximum number of iterations.
    if iter >= maxiter
        convergence = 1;
    end
    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % DISPLAY ITERATION STATISTICS (IF VERBOSE = TRUE).   
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    if verbose ==1
        if mod(iter,5)==0
            fprintf('%3d \t | %1.3e \t | %1.3e \t \n', iter, errP, errD);
        end
    end
end


