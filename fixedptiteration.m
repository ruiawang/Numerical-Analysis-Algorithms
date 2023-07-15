% Fixed point iteration method for finding roots
% currently used with the functions:
% g(x) = pi + 0.5*sin(0.5x), 
% g1(x) = (-2x^2 + x + 3)^0.25,
% g2(x) = ((-x^4 + x + 3)/2)^0.5,
% g3(x) = ((x+3)/(x^2+2))^0.5,
% g4(x) = ((3x^4+2x^2+3)/(4x^3+4x-1)),
% with initial guess p0 = pi, max of 100 iterations and tolerance of 0.01
clear
clc
format long

%%% EDIT %%% 
g1 = @(x) (3 + x - 2 * x^2)^(1/4);
g2 = @(x) ((x + 3 - x^4) / 2)^(1/2);
g3 = @(x) ((x + 3) / (x^2 + 2))^(1/2);
g4 = @(x) ((3 * x^4 + 2 * x^2 + 3) / (4 * x^3 + 4 * x - 1));


g = @(x) pi + (1/2) * sin(x/2);
TOL = 1e-2;
p0 = pi;
N0 = 100;
%%%%%%%%%%%%

i = 1;
while i <= N0
    p = g(p0);
    if abs(p-p0) < TOL
        break;
    end
    i = i + 1;
    p0 = p;
end


if i <= N0         % successful
    fprintf('\nFixed-Point Iteration approximated the the fixed-point p = %.9f after %d iterations.\n\n',p,i);
else                        % not successful 
    fprintf('\nFixed-Point Iteration did not converge within the tolerance in %d iterations.\n\n',N0)
end