% Broyden's Method for finding solutions to nonlinear systems of equations
% currently used for the system:
% F(x,y) = [x^2 + y - 11;
%           x + y^2 - 7]
% Jacobian is given below too
% given with initial guess [-0.164, 1]

clear
clc

F = @(x,y) [x^2 + y - 11; x + y^2 - 7];
J = @(x,y) [2*x, 1; 1, 2*y];

x0 = [-0.164;1]; % initial guess

N = 3000;
tol = 1e-7;

A0 = inv(J(x0(1),x0(2)));

k = 1;
while k < N
    x = x0 - A0*F(x0(1),x0(2))
    err =  max([abs(x(1)-x0(1)), abs(x(2)-x0(2))])
    if err <= tol
        break;
    end
    s = x - x0;
    y = F(x(1),x(2)) - F(x0(1),x0(2));
    A = A0 + (s - A0*y)*(transpose(s))*A0/(transpose(s)*A0*y);
    
    x0 = x;
    A0 = A;
    k = k + 1
end