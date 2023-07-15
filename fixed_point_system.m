% Fixed point iteration for a system of functions
% currently used for the system of equations
% G(x,y) = [sqrt(1-y^2);
%           sqrt((9-5x^2)/21)]
% with initial point [0.5; 0.3]

clear
x0 = [0.5;0.3];
G = @(x,y) [sqrt(1-y^2); sqrt((9-5*x^2)/21)];
x = G(x0(1),x0(2));
tol = 1e-7;
N = 1000;

k = 1;

while k < N
    x = G(x0(1),x0(2));
    err = max([abs(x(1)-x0(1)), abs(x(2)-x0(2))]);
    if err < tol 
        break;
    end
    x0 = x;
    k = k + 1;
end
k