% Secant method of finding roots
% currently set up for f(x) = -x^3-cos(x)
% with tolerance 10^-20, and a maximum of 100 iterations
% and initial two guesses p0 = -1, p1 = 0

clear
clc

format long

%%% Edit %%%
f = @(x) -(x)^3 - cos(x);

tol = 1e-20;
p0 = -1;
p1 = 0;
N0 = 100;
%%%%%%%%%%%%

for n = 1:N0
    p = p1 - f(p1)*(p1-p0)/(f(p1)-f(p0));
    if abs(p-p1) <= tol
        break
    end
    p0 = p1;
    p1 = p;
end
fprintf(['Approximate root is ' num2str(p) ' after ' num2str(n) ' iterations']);