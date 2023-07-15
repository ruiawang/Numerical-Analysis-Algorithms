% Aitken's Delta-Squared Process for Series Acceleration
% currently used on the sequence:
% {p_n} given by p_n = (1/3*e^p_{n-1})^0.5, p_1 = 0.75

clear
clc

format long

%%% Edit %%%
p(1) = 0.75;
g = @(x) (exp(x)/3)^(1/2);
for n = 1:7
p(n+1) = g(p(n));
end
%%%%%%
for n = 1:length(p)-2
q(n) = (p(n)*p(n+2)-p(n+1)^2)/(p(n)-2*p(n+1)+p(n+2))
end