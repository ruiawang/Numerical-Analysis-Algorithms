% Finite Difference method for solving BVPs
% of the form y'' = p(x)y' + q(x)y + r(x), y in [a,b], y(a) = alpha, y(b) = beta
% currently with the system:
% p(x) = 2, q(x) = -1, r(x) = xe^x - x
% a = 0, b = 2, y(a) = 0, y(b) = -4.

clear
clc

p = @(x) 2;
q = @(x) -1;
r = @(x) x*exp(x) - x;
y = @(x) 1/6*x^3*exp(x) - 5/3*x*exp(x) + 2*exp(x) - x - 2;

a = 0;
b = 2;
alpha = 0;
beta = -4;

h1 = 0.1;
N1 = (b-a)/h1;
A1 = zeros(N1-1, N1-1);
x1 = linspace(a,b,N1+1);
p1 = zeros(1,N1+1);
q1 = zeros(1,N1+1);
r1 = zeros(1,N1+1);
rhs1 = zeros(N1-1,1);
for i = 1:N1+1
    p1(i) = p(x1(i));
    q1(i) = q(x1(i));
    r1(i) = r(x1(i));
end
for i = 1:N1-1
    rhs1(i) = h1^2*r1(i);
    if i == 1
        rhs1(i) = h1^2*r1(i)-(1+0.5*h1*p1(i))*alpha;
    end
    if i == N1-1
        rhs1(i) = h1^2*r1(i)-(1-0.5*h1*p1(i))*beta;
    end
    A1(i,i) = -(2+h1^2*q1(i));
    if i >= 2
        A1(i,i-1) = (1+0.5*h1*p1(i));
    end
    if i <= N1-2
        A1(i,i+1) = (1-0.5*h1*p1(i));
    end
end
w1 = A1\rhs1;
w1vals = [alpha;w1;beta];
e1 = zeros(1,N1+1);
for i = 1:N1+1
    e1(i) = abs(w1vals(i) - y(x1(i)));
end

h2 = 0.05;
N2 = (b-a)/h2;
A2 = zeros(N2-1, N2-1);
x2 = linspace(a,b,N2+1);
p2 = zeros(1,N2+1);
q2 = zeros(1,N2+1);
r2 = zeros(1,N2+1);
rhs2 = zeros(N2-1,1);
for i = 1:N2+1
    p2(i) = p(x2(i));
    q2(i) = q(x2(i));
    r2(i) = r(x2(i));
end
for i = 1:N2-1
    rhs2(i) = h2^2*r2(i);
    if i == 1
        rhs2(i) = h2^2*r2(i)-(1+0.5*h2*p2(i))*alpha;
    end
    if i == N2-1
        rhs2(i) = h2^2*r2(i)-(1-0.5*h2*p2(i))*beta;
    end
    A2(i,i) = -(2+h2^2*q2(i));
    if i >= 2
        A2(i,i-1) = (1+0.5*h2*p2(i));
    end
    if i <= N2-2
        A2(i,i+1) = (1-0.5*h2*p2(i));
    end
end
w2 = A2\rhs2;
w2vals = [alpha;w2;beta];
e2 = zeros(1,N2+1);
for i = 1:N2+1
    e2(i) = abs(w2vals(i) - y(x2(i)));
end

h3 = 0.025;
N3 = (b-a)/h3;
A3 = zeros(N3-1, N3-1);
x3 = linspace(a,b,N3+1);
p3 = zeros(1,N3+1);
q3 = zeros(1,N3+1);
r3 = zeros(1,N3+1);
rhs3 = zeros(N3-1,1);
for i = 1:N3+1
    p3(i) = p(x3(i));
    q3(i) = q(x3(i));
    r3(i) = r(x3(i));
end
for i = 1:N3-1
    rhs3(i) = h3^2*r3(i);
    if i == 1
        rhs3(i) = h3^2*r3(i)-(1+0.5*h3*p3(i))*alpha;
    end
    if i == N3-1
        rhs3(i) = h3^2*r3(i)-(1-0.5*h3*p3(i))*beta;
    end
    A3(i,i) = -(2+h3^2*q3(i));
    if i >= 2
        A3(i,i-1) = (1+0.5*h3*p3(i));
    end
    if i <= N3-2
        A3(i,i+1) = (1-0.5*h3*p3(i));
    end
end
w3 = A3\rhs3;
w3vals = [alpha;w3;beta];
e3 = zeros(1,N3+1);
for i = 1:N3+1
    e3(i) = abs(w3vals(i) - y(x3(i)));
end

h4 = 0.0125;
N4 = (b-a)/h4;
A4 = zeros(N4-1, N4-1);
x4 = linspace(a,b,N4+1);
p4 = zeros(1,N4+1);
q4 = zeros(1,N4+1);
r4 = zeros(1,N4+1);
rhs4 = zeros(N4-1,1);
for i = 1:N4+1
    p4(i) = p(x4(i));
    q4(i) = q(x4(i));
    r4(i) = r(x4(i));
end
for i = 1:N4-1
    rhs4(i) = h4^2*r4(i);
    if i == 1
        rhs4(i) = h4^2*r4(i)-(1+0.5*h4*p4(i))*alpha;
    end
    if i == N4-1
        rhs4(i) = h4^2*r4(i)-(1-0.5*h4*p4(i))*beta;
    end
    A4(i,i) = -(2+h4^2*q4(i));
    if i >= 2
        A4(i,i-1) = (1+0.5*h4*p4(i));
    end
    if i <= N4-2
        A4(i,i+1) = (1-0.5*h4*p4(i));
    end
end
w4 = A4\rhs4;
w4vals = [alpha;w4;beta];
e4 = zeros(1,N4+1);
for i = 1:N4+1
    e4(i) = abs(w4vals(i) - y(x4(i)));
end

figure(1)
plot(x1,w1vals,'r',x2,w2vals,'g',x3,w3vals,'m',x4,w4vals,'k')

figure(2)
plot(x1,e1,'r',x2,e2,'g',x3,e3,'m',x4,e4,'k')

l1 = h1*sum(e1)
l2 = h2*sum(e2)
l3 = h3*sum(e3)
l4 = h4*sum(e4)

j1 = (h1^0.5)*sum(e1.^2)^0.5
j2 = (h2^0.5)*sum(e2.^2)^0.5
j3 = (h3^0.5)*sum(e3.^2)^0.5
j4 = (h4^0.5)*sum(e4.^2)^0.5

o2 = log(l1/l2)/log(2)
o3 = log(l2/l3)/log(2)
o4 = log(l3/l4)/log(2)

s2 = log(j1/j2)/log(2)
s3 = log(j2/j3)/log(2)
s4 = log(j3/j4)/log(2)


