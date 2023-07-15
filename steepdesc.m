% Implementation of steepest gradient descent
% for solving systems of nonlinear equations
% currently used for the system:
% F(x,y) = [x^2 + y - 11;
%           x + y^2 - 7]
% gradient for the function is given below,
% given with initial guess [-0.164, 1]

clear
clc

g = @(x,y) (x^2 + y - 11)^2 + (x + y^2 - 7)^2;
grad = @(x,y) [2*(x^2 + y - 11)*2*x + 2*(x + y^2 - 7);
               2*(x^2 + y - 11) + 2*(x + y^2 - 7)*2*y];
x0 = [-0.164;1];
N = 3000;
tol = 1e-7;

k = 1;
while k < N
    d0 = -1*grad(x0(1),x0(2));
    
    a = 1;
    s = 0.95;
    t = 0.45;
    while g(x0(1) +a*d0(1), x0(2) +a*d0(2)) > g(x0(1),x0(2)) - a*t*(d0(1)^2+d0(2)^2)
        a = s*a;
    end
    a0 = a;
    x = x0 + a0*d0;
    err = max([abs(x(1)-x0(1)), abs(x(2)-x0(2))]);
    gval = g(x(1),x(2));
    disp(strcat('k:',32,num2str(k),44,32,'x:',32,'[',num2str(x(1)),44, num2str(x(2)),'],',32,'err:',32,num2str(err),44,32,'gval:',32,num2str(gval)));
    if err <= tol
        break;
    end
    x0 = x;
    k = k + 1;
end
