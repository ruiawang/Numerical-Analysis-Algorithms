%% Spline Interpolation
% Cubic splines with natural boundary
% Currently used for the table of values below
clear; clc

%%%%%%% Edit %%%%%%%
x = [0.1,0.2,0.3,0.4];  %List nodes

%%% Comment out %%%%
y = [-0.62049958,-0.28398668,0.00660095,0.24842440];  %Node outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(x) ~= length(y)
    disp('Error: x and y have different lengths.')
    return
end

plot(x,y,'o');
hold on

n = length(x);
A = zeros(4*(n-1),4*(n-1));
b = zeros(4*(n-1),1);

for j=1:n-1
    b(4*j-2) = y(j);
    b(4*j-1) = y(j+1);
end

for j=1:n-1
   %Left endpoint of spline
   A(4*j-2,4*j-3) = x(j)^3;
   A(4*j-2,4*j-2) = x(j)^2;
   A(4*j-2,4*j-1) = x(j);
   A(4*j-2,4*j)   = 1;
   
   %Right endpoint of spline
   A(4*j-1,4*j-3) = x(j+1)^3;
   A(4*j-1,4*j-2) = x(j+1)^2;
   A(4*j-1,4*j-1) = x(j+1);
   A(4*j-1,4*j)   = 1;
end

for j=1:n-2
   %Derivative match at internal nodes
   A(4*j,4*j-3) = 3*x(j+1)^2;
   A(4*j,4*j-2) = 2*x(j+1);
   A(4*j,4*j-1) = 1;
   A(4*j,4*j+1) = -3*x(j+1)^2;
   A(4*j,4*j+2) = -2*x(j+1);
   A(4*j,4*j+3) = -1;
   
  %2nd Derivative match at internal nodes
   A(4*j+1,4*j-3) = 6*x(j+1);
   A(4*j+1,4*j-2) = 2;
   A(4*j+1,4*j+1) = -6*x(j+1);
   A(4*j+1,4*j+2) = -2;
end

%Natural Boundary conditions
A(1,1) = 6*x(1);
A(1,2) = 2;
A(4*(n-1),4*(n-1)-3) = 6*x(n);
A(4*(n-1),4*(n-1)-2) = 2;

a = A\b;

%Pretty Plots
N = 1000;
X = linspace(x(1),x(end),N);


j = 1;
for i = 1:N
    if X(i) > x(j+1)
        j = j+1;
    end  
    Y(i) = a(4*j-3)*X(i)^3+a(4*j-2)*X(i)^2+a(4*j-1)*X(i)+a(4*j);
end
plot(X,Y)