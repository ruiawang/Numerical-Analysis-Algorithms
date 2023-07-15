% 4th order Runge-Kutta method on a system of ODEs,
% currently used on the ODE problem
% y'' = -y, y(0) = 1, y'(0) = 0

% create mesh of equispaced points in time
N = 400; 
tvals = linspace(0, 4*pi, N); 
h = tvals(2)-tvals(1) ;  % step size

f = @(x,w) [w, -x];
y = @(t) cos(t);

xvals = zeros(N,1); % create blank arrays
wvals = zeros(N,1); % to be edited below 
evals = zeros(N,1);
xvals(1) = 1.0; % initial value for x
wvals(1) = 0.0;  % initial value for w

for i = 1:N-1
    % implement RK-4 method here
    xi = xvals(i);
    wi = wvals(i);
    k1 = f(xi, wi);
    k2 = f(xi + 0.5*h*k1(1), wi + 0.5*h*k1(2));
    k3 = f(xi + 0.5*h*k2(1), wi + 0.5*h*k2(2));
    k4 = f(xi + h*k3(1), wi + h*k3(2));
    xvals(i+1) = xi + h/6*(k1(1)+2*k2(1)+2*k3(1)+k4(1));
    wvals(i+1) = wi + h/6*(k1(2)+2*k2(2)+2*k3(2)+k4(2));
    evals(i+1) = abs(xvals(i+1) - y(tvals(i+1)));
end

% plot position x versus time
figure(1)
plot(tvals, xvals, 'r')

% make phase plot--plot velocity w versus x
figure(2)
plot(xvals, wvals, 'b')

% make error graph
figure(3)
plot(tvals, evals, 'g')