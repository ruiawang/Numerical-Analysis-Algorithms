% Modified Newton-Raphson method to avoid multiple roots
% This is similar to Halley's method, but instead of
% using f(x)/sqrt(abs(f'(x)) like in Halley's method,
% instead we use f(x)/f'(x).

clear
clc

format long

%%% Edit %%%
f  = @(x) exp(x) - 1 - x - (x^2)/2;
f1 = @(x) exp(x) - 1 - x; 
f2 = @(x) exp(x) - 1;   
tol = 1e-10;     
p0 = 1;          
N0 = 100;          
%%%%%%%%%%%%


n = 1;
while n <= N0
    p = p0 - f(p0)*f1(p0)/((f1(p0))^2-f(p0)*f2(p0));
    if abs(p-p0) < tol
        break
    end
    p0 = p;
    n = n+1;
end
fprintf(['Approximate root is ' num2str(p) ' after ' num2str(n) ' iterations']);