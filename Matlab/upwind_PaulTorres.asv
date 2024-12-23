function [y,x] = upwind_PaulTorres(P,Q,R,a,b,alpha,beta,N)
%--------------------------------------------------------------------------
% This fuction returns the solution of the equation y''+P(x)y'+Q(x)y = R(x)
% with boundary values y(a) = alpha and y(b) = beta on the interval [a,b]
% using the Upwind Difference Scheme on a  uniform mesh of N grid points.
%--------------------------------------------------------------------------
% Functions P,Q,& R are function handles to an arbitrary functions.
%--------------------------------------------------------------------------

% Upwind difference scheme:
% y[i+1] - (2 - hP(xi)-h^2 Q(xi))y[i] + (1-hP(xi))y[i-1]=h^2 R(xi)
% Centered difference scheme:
% (1+h*P(xi)/2)y[i+1] - (2-h^2 Q(xi))y[i] + (1-h*P(xi)/2)y[i-1]=h^2 R(xi)

x = linspace(a,b,N+1); % N sub-intervals, N+1 points
y = zeros(N+1,1); % pre-allocate solution (this is more efficient for MATLAB)
y(1) = alpha; % initial condition
y(N+1) = beta; % 2nd Boundary Condition
h = (b-a)/N; % size of uniform mesh

main_diag = -(2 - h * P(x(2:N-1)) - h^2 * Q(x(2:N-1)));
super_diag = ones(size(x(2:N-2)));
sub_diag = 1 - h * P(x(3:N-1));

L = diag(main_diag) + diag(super_diag,1) + diag(sub_diag,-1);

b = h^2 * R(x(2:N));
b(1) = b(1) - (1 - h * P(x(2))) * alpha;
b(N-1) = b(N-1) - beta;

y_h = L\b'; %Solve for all values of y between a & b

y(2:N) = y_h;