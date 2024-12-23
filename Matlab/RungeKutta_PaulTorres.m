function[y] = SecondTaylor_Torres(f,a,b,alpha,N)
%--------------------------------------------------------------------------
%This fuction returns the solution of the equation dY/dx=f(x,Y) with
%initial value Y(a)=alpha using the 4th Order RK Method.
%--------------------------------------------------------------------------
%Function f is a function handle to an arbitrary function f(x,y).
%--------------------------------------------------------------------------

x = linspace(a,b,N+1); % N sub-intervals, N+1 points
y = zeros(N+1,1); % pre-allocate solution (this is more efficient for MATLAB)
y(1) = alpha; % initial condition
h = (b-a)/N; % size of uniform mesh

% compute the second order taylor series method
for i=1:N
    k1 = f(x(i),y(i));
    k2 = f(x(i) +.5*h, y(i) +.5*k1*h);
    k3 = f(x(i) +.5*h, y(i) + .5*k2*h);
    k4 = f(x(i) + h, y(i) +k3*h);
    y(i+1) = y(i) + (1/6)*h*(k1 + 2*k2 + 2*k3 +k4 ); % obtain numerical value for each term
end
