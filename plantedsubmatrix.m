function [A,X0,Y0]=plantedsubmatrix(M,N,m,n,p,q)
% PLANTEDSUBMATRIX Makes binary matrix A with planted mn-submatrix.
% Generates mn-submatrix with expected density q in MxN matrix A with
% expected densities of remaining entries equal to q.
%
% INPUT:
% M,N - desired dimensions of A.
% m,n - desired dimensions of planted submatrix.
% p - desired noise density.
% q - desired in-group density.
% OUTPUT:
% A - matrix containing desired planted submatrix.
% X0, Y0 - matrix representation of the planted submatrix.

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% GENERATE NOISE ENTRIES OF A.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Initialize A as uniform random matrix.
A=rand(M,N);

% Round entries of A to 0 if less than 1-p and up to 1 otherwise.
A=ceil(A-(1-p));

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% FILL IN DENSE BLOCK.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Repeat with mn-block and threshhold 1-q.
tmp=rand(m,n);
A(1:m, 1:n)=ceil(tmp-(1-q));

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CALCULATE MATRIX REPRESENTATION OF PLANTED SUBMATRIX.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% X0
X0=zeros(M,N);
X0(1:m,1:n)=ones(m,n);

% Y0
Y0=zeros(M,N);
Y0(1:m,1:n)=ones(m,n)-A(1:m,1:n);







