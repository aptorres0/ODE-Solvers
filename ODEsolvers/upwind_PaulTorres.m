function [y,x] = upwind_PaulTorres(P,Q,R,a,b,alpha,beta,N)
%--------------------------------------------------------------------------
%This fuction returns the solution of the equation y''+P(x)y'+Q(x)y = R(x)
%with boundary values y(a) = alpha and y(b) = beta on the interval [a,b]
%using the Upwind Difference Scheme on a  uniform mesh of N grid points.
%--------------------------------------------------------------------------
%Functions P,Q,& R are function handles to an arbitrary functions.
%--------------------------------------------------------------------------

x = linspace(a,b,N+1); % N sub-intervals, N+1 points
y = zeros(N+1,1); % pre-allocate solution (this is more efficient for MATLAB)
y(1) = alpha; % initial condition
y(N+1) = beta; % 2nd Boundary Condition
h = (b-a)/N; % size of uniform mesh

L = zeros(N-1,N-1); % pre-allocate linear transformation matrix
L(1,1) = -(2-h*P(x(2))-h*h*Q(x(2)));
L(1,2) = 1;
L(N-1,N-2) = 1-h*P(x(N));
L(N-1,N-1) = -(2-h*P(x(N))-h*h*Q(x(N)));

b = zeros(N-1,1); % Setting up right hand side of linear equation
b(1,1) = h*h*R(x(2))-(1-h*P(x(2)))*alpha;
b(N-1,1) = h*h*R(x(N))-beta;

for i = 2:N-2
    L(i,i) = -(2-h*P(x(i+1))-h*h*Q(x(i+1)));
    L(i,i+1) = 1;
    L(i,i-1) = 1-h*P(x(i+1));
    
    b(i,1) = h*h*R(x(i+1));
    
end

y_unknown = L\b; %Solve for all values of y between a & b

for i = 2:N
    y(i,1) = y_unknown(i-1,1);
end



