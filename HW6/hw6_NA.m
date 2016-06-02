% Luan Cong Doan - Numerical Analysis HW6
%% Problem 1d
close all; clear all; clc;
x = -3:0.01:0;
He = 0.5*x.^2 + x +1;
RK = (1/24)*x.^4 + (1/6)*x.^3 + (1/2)*x.^2 + x +1;
figure; plot(x,He,'r'); hold on;
plot(x,RK,'b--'); grid on;
xlabel('x'); ylabel('\lambda.h');
legend('Heun method', 'Runge-Kutta method');

%%
RKx = @(x) (1/24)*x.^3 + (1/6)*x.^2 + (1/2)*x + 1;
x0 = fzero(RKx,2.75);

%%
x = -1.5:0.001:1.5; y = -1.5:0.001:1.5;
fxy = x.^2 + y.^2;
plot(fxy);

%% Problem 1 - set of h\lambda
close all; clear all; clc;
 npts = 200;                        % the larger the number, the nicer the plot
 z_real = linspace(-3,1,npts);      % real parts of h*lambda
 z_imag = linspace(-3,3,npts);      % imaginary parts of h*lambda

 [Zr,Zi] = meshgrid(z_real,z_imag); % matrix of real, imag parts
 hLambda = Zr+1i*Zi;                % matrix of complex h*Lambda values

 figure(4), clf;
 contour(z_real,z_imag,abs(1+hLambda + (hLambda.^2)/2 )<1,[1 1]);
 axis equal; grid on; xlabel('real'); ylabel('imag');
 title('Set of h\lambda for Heun method');
 %axis([min(z_real) max(z_real) min(z_imag) max(z_imag)]);
 axis([-3 1 -2 2]);
 print('hw6_d1','-dpng');
 
 figure(6), clf;
 contour(z_real,z_imag,abs(1+hLambda + (hLambda.^2)/2 + (hLambda.^3)/6 + (hLambda.^4)/24 )<1,[1 1]);
 axis equal; grid on; xlabel('real'); ylabel('imag');
 title('Set of h\lambda for Rungae-Kutta method');
 axis([min(z_real) max(z_real) min(z_imag) max(z_imag)]);
 print('hw6_d2','-dpng');
 
 %% %% Problem 2 - set of h\lambda
 close all; clear all; clc;
 npts = 200;                        % the larger the number, the nicer the plot
 z_real = linspace(-0.1,0.1,npts);      % real parts of h*lambda01
 z_imag = linspace(-2,2,npts);      % imaginary parts of h*lambda

 [Zr,Zi] = meshgrid(z_real,z_imag); % matrix of real, imag parts
 hLambda = Zr+1i*Zi;                % matrix of complex h*Lambda values

 figure; clf; 
 contour(z_real,z_imag,abs(hLambda - sqrt(hLambda.^2 +1) )<1,[1 1],'k');
 grid on; xlabel('real'); ylabel('imag');
 title('Set of h\lambda for \gamma_1');
 print('hw6_2e1','-dpng');
 figure; clf; 
 contour(z_real,z_imag,abs(hLambda + sqrt(hLambda.^2 +1) )<1,[1 1],'k');
 grid on; xlabel('real'); ylabel('imag');
 title('Set of h\lambda for \gamma_2');
 print('hw6_2e2','-dpng');
 axis([min(z_real) max(z_real) min(z_imag) max(z_imag)])

 
%% Problem 3a
close all; clear all;
[t,x] = ode45(@(t,x) [0 1;0 10]*x, [0 1], [0;0.00045395]);
% c1 = 1/(exp(10)-1);
% ue = c1*exp(10*t) - c1;
figure; plot(t,x(:,1),'r'); grid on;
% hold on; plot(t,ue,'b--');
% legend('approx','true solution');
xlabel('x'); ylabel('u(x)');
title('Approximation for u(x) utilizing ode45, \epsilon = 10');
print('hw6_3a','-dpng');


%% Problem 3b
close all; clear all;
[t,x] = ode45(@(t,x) [0 1;0 50]*x, [0 1], [0;0.50395*1.0e-19]);
% c1 = 1/(exp(50)-1);
% ue = c1*exp(50*t) - c1;
figure; plot(t,x(:,1),'r'); grid on;
% hold on; plot(t,ue,'b--');
% legend('approx','true solution');
xlabel('x'); ylabel('u(x)');
title('Approximation for u(x) utilizing ode45, \epsilon = 50');
print('hw6_3b','-dpng');

%% Problem 4
close all; clear all; clc;
 npts = 200;                        % the larger the number, the nicer the plot
 z_real = linspace(-3,1,npts);      % real parts of h*lambda
 z_imag = linspace(-3,3,npts);      % imaginary parts of h*lambda

 [Zr,Zi] = meshgrid(z_real,z_imag); % matrix of real, imag parts
 hLambda = Zr+1i*Zi;                % matrix of complex h*Lambda values

% method one:  |1+h*lambda|<1 corresponds to yellow region
 figure(1), clf
 pcolor(z_real,z_imag,double(abs(1+hLambda)<1))
 axis equal
 axis([min(z_real) max(z_real) min(z_imag) max(z_imag)])
 colorbar

% method two:  draw |1+h*lambda|=1 using contour plotting
 figure(2), clf
 contour(z_real,z_imag,abs(1+hLambda )<1,[1 1])
 axis equal
 axis([min(z_real) max(z_real) min(z_imag) max(z_imag)])

 % method 3
 figure(3), clf
 pcolor(z_real,z_imag,double(abs(1+hLambda + (hLambda.^2)/2)<1))
 axis equal
 axis([min(z_real) max(z_real) min(z_imag) max(z_imag)])
 colorbar
 figure(4), clf
 contour(z_real,z_imag,abs(1+hLambda + (hLambda.^2)/2 )<1,[1 1])
 axis equal
 axis([min(z_real) max(z_real) min(z_imag) max(z_imag)])
 
 % method 3
 figure(5), clf
 pcolor(z_real,z_imag,double(abs(1+hLambda + (hLambda.^2)/2 + (hLambda.^3)/6 + (hLambda.^4)/24)<1))
 axis equal
 axis([min(z_real) max(z_real) min(z_imag) max(z_imag)])
 colorbar
 figure(6), clf
 contour(z_real,z_imag,abs(1+hLambda + (hLambda.^2)/2 + (hLambda.^2)/2 + (hLambda.^3)/6 + (hLambda.^4)/24 )<1,[1 1])
 axis equal
 axis([min(z_real) max(z_real) min(z_imag) max(z_imag)])
 
 
  
%%
close all; clear all; clc;
n = 50; Dx = 1/n;
x = 0:Dx:1; x(1) = [];

A = -1*eye(n) + [zeros(n-1,1),eye(n-1); zeros(1,n)];
A(n,1) = 1; A = A/Dx;

% For Dt is twice time than the maximum value of Dx
u1(:,1) = sin(2*pi*x);
Dt1 = 2*Dx; t1 = 0:Dt1:2; t1(1) = []; 
k1 = 2/Dt1;
[T1,X1] = meshgrid(t1,x);
for i =1:k1
    u1(:,i+1) = u1(:,i) + Dt1*A*u1(:,i);
end
u1(:,end) = [];
figure; plot(t1,u1(:,end)); grid on;
xlabel('n'); ylabel('u_k(2,x)');
title('u_k(2,x) for \Delta t = 2\Delta x');
print('hw6_4d1','-dpng');
figure; surf(X1,T1,u1); shading interp;
title('u_k(2,x) for \Delta t = 2\Delta x');
print('hw6_4d2','-dpng');

% For Dt is equal to the maximum value of Dx
u2(:,1) = sin(2*pi*x);
Dt2 = Dx; t2 = 0:Dt2:2; t2(1) = []; k2 = 2/Dt2;
[T2,X2] = meshgrid(t2,x);
for j =1:k2
    u2(:,j+1) = u2(:,j) + Dt2*A*u2(:,j);
end
u2(:,k2+1) = [];
figure; plot(t2,u2(end,:)); grid on;
xlabel('n'); ylabel('u_k(2,x)');
title('u_k(2,x) for \Delta t = \Delta x');
print('hw6_4d3','-dpng');
figure; surf(X2,T2,u2); shading interp;
title('u_k(2,x) for \Delta t = \Delta x');
print('hw6_4d4','-dpng');



 