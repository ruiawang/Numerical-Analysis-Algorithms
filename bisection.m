% Bisection method for finding roots
% currently used with f(x) = x + cos(x) in the interval [-5,5]
% with tolerance 10^-10
clear
clc

format long
%%% Edit %%%
f = @(x) x + cos(x);
tol = 1e-10;
a = -5;
b = 5;
%%%%%%%%%%%%

n = 1;
while (b-a)/2 > tol
    m = a+(b-a)/2;
    if f(a)*f(m) < 0
        b = m;
    else
        a = m;
    end
    n = n+1;
end

fprintf(['approximation is ' num2str(m) ' after ' num2str(n) ' iterations']);