% Linear shooting method for solving BVPs
% of the form y'' = p(x)y' + q(x)y + r(x), y(a) = alpha, y(b) = beta
% currently with the system:
% p(x) = 2, q(x) = -1, r(x) = xe^x - x
% a = 0, b = 2, y(a) = 0, y(b) = -4.

p = @(x) 2;
q = @(x) -1;
r = @(x) x*exp(x)-x;

f = @(x,y1,y2) [y2, p(x)*y1 + q(x)*y2 + r(x)];
y = @(x) 1/6*x^3*exp(x) - 5/3*x*exp(x) + 2*exp(x) - x - 2;

a = 0; %boundary values
alpha = 0;
b = 2;
beta = -4;

h = 0.2;
N = (b-a)/h;

xvals = linspace(a,b,N+1);
w1vals = zeros(1,N+1);
w2vals = zeros(1,N+1);
v1vals = zeros(1,N+1);
v2vals = zeros(1,N+1);
yvals = zeros(1,N+1);
evals = zeros(1,N+1);

w1vals(1) = alpha; %initial conditions for w1,w2 and v1,v2
w2vals(1) = 0;
v1vals(1) = alpha;
v2vals(1) = 1;

for i = 1:N
    w1vals(i+1) = w1vals(i) + h*w2vals(i);
    w2vals(i+1) = w2vals(i) + h*(p(xvals(i))*w2vals(i) + q(xvals(i))*w1vals(i) + r(xvals(i)));
    v1vals(i+1) = v1vals(i) + h*v2vals(i);
    v2vals(i+1) = v2vals(i) + h*(p(xvals(i))*v2vals(i) + q(xvals(i))*v1vals(i) + r(xvals(i)));
end
lambda = (beta - v1vals(N+1))/(w1vals(N+1) - v1vals(N+1));
for i = 1:N+1
    yvals(i) = lambda*w1vals(i) + (1-lambda)*v1vals(i);
    evals(i) = abs(yvals(i) - y(xvals(i)));
end

figure(1)
plot(xvals,yvals)
figure(2)
plot(xvals, evals)