% Newton's method for a system of equations
% currently used for the system:
% F(x,y) = [x^2 + y - 11;
%           x + y^2 - 7]
% Jacobian is given below too
% given with initial guess [-0.164, 1]

clear
clc
F = @(x,y) [x^2 + y - 11; x + y^2 - 7]; 
J = @(x,y) [2*x, 1; 1, 2*y];

x0 = [-0.164;1];

N = 3000;
tol = 1e-7;

k = 1;
while k < N
    x = x0 - J(x0(1),x0(2))\F(x0(1),x0(2))
    err = max([abs(x(1)-x0(1)), abs(x(2)-x0(2))])
    if err <= tol
        break;
    end
    x0 = x;
    k = k + 1
end