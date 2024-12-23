function y = Eulers_Method(f,a,b,alpha,N)
%--------------------------------------------------------------------------
% This fuction returns the solution of the equation dY/dx=f(x,Y) with
% initial value Y(a)=alpha using Euler's Method.
%--------------------------------------------------------------------------
% Function f is a function handle to an arbitrary function f(x,y).
%--------------------------------------------------------------------------

x = linspace(a,b,N+1); % N sub-intervals, N+1 points
y = zeros(N+1,1); % pre-allocate solution (this is more efficient for MATLAB)
y(1) = alpha; % initial condition
h = (b-a)/N; % size of uniform mesh

% compute Euler's method
for i=1:N
    y(i+1) = y(i) + h*f(x(i),y(i)); % obtain numerical value for each term
end
