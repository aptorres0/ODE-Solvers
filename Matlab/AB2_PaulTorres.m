function [y,x] = AB2_PaulTorres(f,a,b,alpha,y1,N)

x = linspace(a,b,N+1); % N sub-intervals, N+1 points
y = zeros(N+1,1); % pre-allocate solution (this is more efficient for MATLAB)
y(1) = alpha; % initial condition
y(2) = y1;
h = (b-a)/N;

for i = 3:N+1
    y(i) = y(i-1) + h*(1.5*f(x(i-1),y(i-1))-.5*f(x(i-2),y(i-2)));
end