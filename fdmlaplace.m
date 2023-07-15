% Finite difference for solving PDEs that satisfy Laplace's Equation
% currently used for the BVP u_xx + u_yy = 0,
% with (x,y) in [0,1] x [0,1], with boundaries
% u(0,y) = sin(2*pi*y), u(1,y) = 0, u(x,0) = 0, u(x,1) = 0

h = 0.05; % step size info
N = 1/h;

% boundary conditions
bl = @(y) sin(2*pi*y); %bl = u(0,y)
br = @(y) 0;           %br = u(1,y)
bb = @(x) 0;           %bb = u(x,0)
bt = @(x) 0;           %bt = u(x,1)

% initialize
xh = 0:h:1;
yh = xh;
M = N-1;
w = zeros(M);
v = zeros(M*M,1);
b = zeros(M*M,1);

D = zeros(M);
I = eye(M);

% building D
for i=1:M
    if i == 1
        D(i,i) = -4;
        D(i,i+1) = 1;
    elseif i == M
        D(i,i-1) = 1;
        D(i,i) = -4;
    else
        D(i,i-1) = 1;
        D(i,i) = -4;
        D(i,i+1) = 1;
    end
end

% building A
for i = 1:M
    if i == 1
        A((i-1)*M + 1:(i-1)*M + M,(i-1)*M + 1:(i-1)*M + M) = D;
        A((i-1)*M + 1:(i-1)*M + M,i*M + 1:i*M + M) = I;
    elseif i == M
        A((i-1)*M + 1:(i-1)*M + M,(i-1)*M + 1:(i-1)*M + M) = D;
        A((i-1)*M + 1:(i-1)*M + M,(i-2)*M + 1:(i-2)*M + M) = I;
    else
        A((i-1)*M + 1:(i-1)*M + M,(i-1)*M + 1:(i-1)*M + M) = D;
        A((i-1)*M + 1:(i-1)*M + M,i*M + 1:i*M + M)=I;
        A((i-1)*M + 1:(i-1)*M + M,(i-2)*M + 1:(i-2)*M + M) = I;
    end
end

% building b
for i = 1:M
    if i == 1
        b((i-1)*M + 1,1) = -bl(yh(2)) - bb(xh(2));
        b((i-1)*M + 2:(i-1)*M + M - 1, 1) = -bb(xh(3:M));
        b((i-1)*M + M,1) = -br(yh(2)) - bb(xh(M+1));
    elseif i == M
        b((i-1)*M + 1,1) = -bl(yh(M+1)) - bt(xh(2));
        b((i-1)*M + 2:(i-1)*M + M - 1,1) = -bt(xh(3:M));
        b((i-1)*M + M,1) = -br(yh(M+1)) - bt(xh(M+1));
    else
        b((i-1)*M + 1,1) = -bl(yh(i+1));
        b((i-1)*M + M,1) = -br(yh(i+1));
    end
end

% solving Av = b
v = A\b;

% turning v vector into w matrix. w(i,j) is approx for u(x_i, y_j)
for j = 1:M
    for i = 1:M
        w(i,j) = v((j-1)*M + i);
    end
end

% plot
surf(xh(2:N),yh(2:N)', w');