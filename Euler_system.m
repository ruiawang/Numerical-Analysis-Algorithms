% Euler's method on a system of ODEs,
% currently used on the ODE problem
% y'' = -sin(y), y(0) = -1, y'(0) = 0

% create mesh of equispaced points in time
N = 1000; 
tvals = linspace(0., 4*pi, N); 
h = tvals(2)-tvals(1) ;  % step size
xvals = zeros(N,1); % create blank arrays
wvals = zeros(N,1); % to be edited below 
xvals(1) = -1.0; % initial value for x
wvals(1) = 0.0;  % initial value for w
for i = 1:N-1
    % implement EM method here
    xi = xvals(i);
    wi = wvals(i);
    xvals(i+1) = xi + h*wi;
    wvals(i+1) = wi - h*sin(xi);
end
% plot position x versus time
figure(1)
plot(tvals, xvals, 'r')
% make phase plot--plot velocity w versus x
figure(2)
plot(xvals, wvals, 'b')