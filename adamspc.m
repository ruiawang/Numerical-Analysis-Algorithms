% Adams 2nd order Predictor-Corrector Method
% using the explicit 2-step Adams-Bashforth method as the predictor
% and the implicit 2-step Adams-Moulton method as the corrector,
% currently used on the ODE y' = 1 + (t-y)^2, y(2) = 1.

clear
clc

f = @(t,y) 1 + (t-y)^2;

h1 = 0.5;
N1 = 8;

t1 = zeros(1, N1);
y_ab_1 = zeros(1,N1);
y_trap_1 = zeros(1,N1);

y_ab_1(1) = 1;
y_trap_1(1) = 1;
t1(1) = 2;

y_ab_1(2) = 1 + h1*f(2+0.5*h1, 1 + h1);
t1(2) = 2 + h1;

for n = 2:N1
    y_ab_1(n+1) = y_ab_1(n) + h1*(1.5*f(t1(n), y_ab_1(n)) - 0.5*f(t1(n-1), y_ab_1(n-1)));
    t1(n+1) = t1(n) + h1;
end

for n = 1:N1
    y_trap_1(n+1) = y_trap_1(n) + 0.5*h1*(f(t1(n), y_trap_1(n)) + f(t1(n+1), y_ab_1(n+1)));
end

h2 = 0.05;
N2 = 80;

t2 = zeros(1, N2);
y_ab_2 = zeros(1,N2);
y_trap_2 = zeros(1,N2);

y_ab_2(1) = 1;
y_trap_2(1) = 1;
t2(1) = 2;

y_ab_2(2) = 1 + h2*f(2+0.5*h2, 1 + h2);
t2(2) = 2 + h2;

for n = 2:N2
    y_ab_2(n+1) = y_ab_2(n) + h2*(1.5*f(t2(n), y_ab_2(n)) - 0.5*f(t2(n-1), y_ab_2(n-1)));
    t2(n+1) = t2(n) + h2;
end

for n = 1:N2
    y_trap_2(n+1) = y_trap_2(n) + 0.5*h2*(f(t2(n), y_trap_2(n)) + f(t2(n+1), y_ab_2(n+1)));
end

h3 = 0.005;
N3 = 800;

t3 = zeros(1, N3);
y_ab_3 = zeros(1,N3);
y_trap_3 = zeros(1,N3);

y_ab_3(1) = 1;
y_trap_3(1) = 1;
t3(1) = 2;

y_ab_3(2) = 1 + h3*f(2+0.5*h3, 1 + h3);
t3(2) = 2 + h3;

for n = 2:N3
    y_ab_3(n+1) = y_ab_3(n) + h3*(1.5*f(t3(n), y_ab_3(n)) - 0.5*f(t3(n-1), y_ab_3(n-1)));
    t3(n+1) = t3(n) + h3;
end
for n = 1:N3
    y_trap_3(n+1) = y_trap_3(n) + 0.5*h3*(f(t3(n), y_trap_3(n)) + f(t3(n+1), y_ab_3(n+1)));
end

g = @(t) t + 1/(1-t);
gvals = zeros(1, N3);
for n = 1:N3+1
    gvals(n) = g(t3(n));
end

figure(1)
plot(t1, y_trap_1,'r', t2, y_trap_2, 'm', t3, y_trap_3, 'k', t3, gvals, 'c--')