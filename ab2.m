% Two-step Adams-Bashforth method
% currently used on the ODE: y' = 1 + sin(t)-y, y(2)=1
clear
clc
f = @(t, y) sin(t) - y;

h = 0.01;
N = 2000;
yvals = zeros(1,N);
tvals = zeros(1,N);
yvals(1) = 1; 
yvals(2) = 1 + h*f(0.5*h, 1-0.5*h);
tvals(1) = 0;
tvals(2) = h;
for n = 2:N
    yvals(n+1) = yvals(n) + h*(1.5*f(tvals(n), yvals(n)) - 0.5*f(tvals(n-1), yvals(n-1)));
    tvals(n+1) = tvals(n) + h;
end

plot(tvals,yvals)
