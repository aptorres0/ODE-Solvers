function [y,x] = centered_PaulTorres(P,Q,R,a,b,alpha,beta,N)
%--------------------------------------------------------------------------
% This fuction returns the solution of the equation y''+P(x)y'+Q(x)y = R(x)
% with boundary values y(a) = alpha and y(b) = beta on the interval [a,b]
% using the Centered Difference scheme on a uniform mesh of N+1 grid points.
%--------------------------------------------------------------------------
% Functions P,Q,& R are function handles to an arbitrary functions.
%--------------------------------------------------------------------------

x = linspace(a,b,N+1); % N sub-intervals, N+1 points
y = zeros(N+1,1); % pre-allocate solution (this is more efficient for MATLAB)
y(1) = alpha; % 1st Boundary Condition
y(N+1) = beta; % 2nd Boundary Condition
h = (b-a)/N; % size of uniform mesh

% For the internal points, it is possible to write the finite difference
% equations as Ly=b where L is a tridiagonal matrix
% if there are N+1 points in our mesh, then there are N-1 internal grid 
% for which we need to solve for.
L = diag(-(2-h*h*Q(x(2:N)))); % Populates the main diagonal
L = L + diag(1+0.5*h*P(x(2:N-1)),1); % Populates the right diagonal
L = L + diag(1-0.5*h*P(x(3:N)),-1); % Populates the left diagonal

b = h*h*R(x(2:N));
b(1) = b(1) - (1-.5*h*P(x(2)))*alpha;
b(end) = b(end) - (1+.5*h*P(x(N)))*beta;
b = b';

y_h = L\b;
y(2:N) = y_h;