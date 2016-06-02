% Luan Cong Doan - Numerical Analysis Homework 5

%% Problem 1 - check
n = 0.01:0.01:1;
for i = 1:length(n)
E(i) = -(0.99^3/(12*10))*(3*sin(n(i))/(4*n(i).^(5/2)) - cos(n(i))/n(i).^(3/2) - sin(n(i))/sqrt(n(i)));
E2(i) = 3*sin(n(i))/(4*n(i).^(5/2)) - cos(n(i))/n(i).^(3/2) - sin(n(i))/sqrt(n(i));
E3(i) = - cos(n(i))/n(i).^(3/2);
end
figure; plot(n,E); grid on;
figure; plot(n,E2); grid on;
figure; plot(n,E3); grid on;

%% Problem 2a
f5 = @(x) x^(4)*exp(-x);
f10 = @(x) x^(9)*exp(-x);
f05 = @(x) x^(-0.5)*exp(-x);
a = 0;
b = inf;
Gamma5 = gamma(5);
Gama5_s = simpson(f5,a,b,10);


%% Problem 2a
G5 = gamma(5);
G10 = gamma(10);
Gd5 = gamma(0.5);

%% Problem 2c
N = 1:1:25;
figure;
for i = 1:length(N)
    [x,~,~] = gauss_lag(i);
    plot(x, i*ones(i+1,1), '.','markersize',24); hold on; grid on;
end
xlabel('x'); ylabel('N'); title('Nodes of (n+1) points Gauss-Laguerre');
print('hw5_2c','-dpng');

%% Problem 2d


%% Problem 3d
x = 0:0.01:2;
px = (x.*x.*x + 12.*x)/(3.*x.*x +4);
figure; plot(x,px); grid on;

%% Problem 4b
close all; clear all; clc;
%x = 1:0.05:1.5;
x = -2.5:0.05:2;
%syms x;
f = sin(x).*exp(x);
fx = sin(1.25).*exp(1.25);
f1x = cos(x).*exp(x) + sin(x).*exp(x);
f1 = cos(1.25).*exp(1.25) + sin(1.25).*exp(1.25);
f2x = 2*cos(x).*exp(x);
f2 = 2*cos(1.25).*exp(1.25);
N2 = f1*(x-1.25) + fx;
N3 = f2*x.^2 -2*(f2*1.25-f1)*x + 2*(fx-f1*1.25);


figure; plot(x,f,'k','linewidth',1); grid on; hold on;
plot(x,N2,'r--'); plot(x,N3,'b.');
xlabel('x'); ylabel('f(x) = sin(x).e^x');
title('Newton method approximation of root at x_0 = 1.25');
legend('fx','tangent line','parabola');
%print('hw5_4b','-dpng');

%% 


%% Problem 4f
close all; clear all; clc;

fx = @(x) sin(x)*exp(x);
fprime = @(x) cos(x)*exp(x) + sin(x)*exp(x);
f2prime = @(x) 2*cos(x)*exp(x);

xstar1 = newton(fx,fprime,1.25);
xstar2 = newton2(fx,fprime,f2prime,1.25);

%%
x = -1:0.01:3;
fx = sin(x).*exp(x);
fprime = cos(x).*exp(x) + sin(x).*exp(x);
f2prime = 2.*cos(x).*exp(x);
figure; plot(x,fx); hold on; plot(x,fprime);
plot(x,f2prime); grid on;
legend('f','fprime','f2prime');

%% Problem 4h
f4h = @(x) x^2*exp(x);
f4h_prime = @(x) 2*x*exp(x) + x^2*exp(x);
f4h_2prime = @(x) 2*exp(x) + 4*x*exp(x) + x^2*exp(x);
x4_star1 = newton(f4h,f4h_prime,1);
x4_star2 = halley(f4h,f4h_prime,f4h_2prime,1);

%% Problem 5c
% c = poly(1:3);
c = poly(1:24);
n = length(c);
syms x;
px =0;
for i = 1:n
    px = px + c(i)*x.^(n-i);
end
C = zeros(n-1);
for i = 1:n-1
    for j=1:n-1
        if j-i == 1
            C(i,j) = 1;
        end
    end
end
for k = 1:n-1
    C(n-1,k) = -c(n-k+1)/c(1);
end
x = eig(C);

%px = @(x) c(1)*x^(n-1) + c(2)*x^(n-2) + c(3)*x^(n-3) + c(4)*x^(n-4);


%% Problem 5e
c = [1 -2 -1 2 -1/4];
n = length(c);
syms x;
px =0;
for i = 1:n
    px = px + c(i)*x.^(n-i);
end
C = zeros(n-1);
for i = 1:n-1
    for j=1:n-1
        if j-i == 1
            C(i,j) = 1;
        end
    end
end
for k = 1:n-1
    C(n-1,k) = -c(n-k+1)/c(1);
end
x = eig(C);

%% Problem 5e
a = [-1/8 1/2 0 -1/2 1/8];
n = length(a);
G = zeros(n-1);
G(1,2) = 1;
for i = 2:n-2
    for j = 1:n-1
        if (i-j ==1) 
            G(i,j) = 1/2;
        elseif (j-i==1)
            G(i,j) = 1/2;
        end
    end
end
for k = 1:n-1
    G(n-1,k) = -a(k)/(2*a(n));
end
x = eig(G);


        