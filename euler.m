% Euler's Method for solving ODEs
% currently used for the ODE problem
% y' = 1 + (t-y)^2, y(2) = 1
% multiple implementations for different step sizes
% for h = 0.5, 0.05, 0.005

clear
clc
h = 0.5;
N = 8;
y1(1) = 1;
t1(1) = 2;
for n = 1:N
    y1(n+1) = y1(n) + h*(1+(t1(n) - y1(n))^2);
    t1(n+1) = t1(n) + h;
end

h = 0.05;
N = 80;
y2(1) = 1;
t2(1) = 2;
for n = 1:N
    y2(n+1) = y2(n) + h*(1+(t2(n) - y2(n))^2);
    t2(n+1) = t2(n) + h;
end

h = 0.005;
N = 800;
y3(1) = 1;
t3(1) = 2;
for n = 1:N
    y3(n+1) = y3(n) + h*(1+(t3(n) - y3(n))^2);
    t3(n+1) = t3(n) + h;
end
plot(t1,y1,'r',t2,y2,'m',t3,y3,'k')