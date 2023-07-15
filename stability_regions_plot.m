% Analysis of stability regions for RK-2, RK-3, and RK-4 methods

% specify range in x and y for plot below, and number of points
xrange = 3; 
yrange = 3; 
N = 301; 

% construct mesh
xv = linspace(-1.*xrange, xrange, N);
yv = linspace(-1.*yrange, yrange, N);
[x,y] = meshgrid(xv, yv);            % x and y each are NxN arrays to be
                                     % used in the call to 'contour' below

% calculate z 
z = x + 1i*y;         % notice the imaginary number i is written 1i

% second order RK     % define Q here for RK2
Q_rk2 = 1+z+z.^2;

% third order RK      % define Q here for RK3
Q_rk3 = 1+z+0.5*z.^2+1/24*z.^3;

% fourth order RK     % define Q here for RK4
Q_rk4 = 1+z+0.5*z.^2+1/6*z.^3+1/24*z.^4;

% compute complex modulus of Q
Q_rk2_mag = abs(Q_rk2);
Q_rk3_mag = abs(Q_rk3);
Q_rk4_mag = abs(Q_rk4); 

% plot contours--these plot the contour lines |Q(z)| = 1; since the region
%                of absolute stability is defined to be when |Q(z)| < 1, 
%                this means everything ~inside~ the curves is the stable
%                region, and everything outside is unstable. 
contour(x,y,Q_rk2_mag, [1 1], 'g-')
hold on; 
contour(x,y,Q_rk3_mag, [1 1], 'r-')
hold on; 
contour(x,y,Q_rk4_mag, [1 1], 'b-')
axis([-1.*xrange, xrange, -1.*yrange, yrange]);
axis('square') 
xlabel('Real \lambda h '); 
ylabel('Im \lambda h '); 
grid on; 
legend('RK2', 'RK3',  'RK4')