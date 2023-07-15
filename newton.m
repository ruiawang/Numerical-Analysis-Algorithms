% Newton-Raphson method for finding roots
% currently set up with f(x) = e^x - 1 - x - 0.5x^2,
% tolerance of 10^-10, initial guess p0 = 1,
% and maximum 100 iterations

clear
clc

format long
%%% Edit %%%
f  = @(x) exp(x) - 1 - x - (x^2)/2;
F = @(x) exp(x) - 1 - x; 

tol = 1e-10;
p0 = 1;
N0 = 100;
%%%%%%%%%%%%


i = 1;
while i <= N0
    p = p0 - f(p0)/F(p0);
    if abs(p - p0) < tol
        break
    end
    p0 = p;
    i = i + 1;
end
fprintf(['Approximate root is ' num2str(p) ' after ' num2str(i) ' iterations']);