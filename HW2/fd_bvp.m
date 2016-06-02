% solution of -u''(x) = g(x) with u(0)=u(1)=0 for g(x) = sin(pi*x).
%
% The exact solution is u(x) = (1/pi^2) sin(pi*x).
close all;
clear all;
clc;

% computed constants from 1b and 1c
A = -1/12; B = 4/3; C = -5/2;
D = 11/12; E = -5/3; F = 1/2; G = 1/3; H = -1/12;

 g = @(x) sin(pi*x);
 true_u  = @(x) (1/pi^2)*sin(pi*x);
 xx = linspace(0,1,500);

 n = 6;
 h = 1/n;
 x = [0:n]'*h;
 uexact = true_u(x);
 f = -g(x(2:n));

 % Quadratic interpolation
 A2 = (-2*eye(n-1)+diag(ones(n-2,1),1)+diag(ones(n-2,1),-1))/h^2;
 u2 = A2\f;
 u2 = [0;u2;0];            % add in Dirichlet values, u(0)=u(1)=0
 
 % Quartic interpolation
 A41 = [-E, -F,-G,-H, zeros(1,n-5)]/h^2;
 A42n2 = (-A*diag(ones(n-2,1),-1) - B*eye(n-1) -C*diag(ones(n-2,1),1) - B*diag(ones(n-3,1),2) -A*diag(ones(n-4,1),3))/h^2;
 A42n2(n-1,:) = [];
 A42n2(n-2,:) = [];
 A4n = [zeros(1,n-5), -H, -G,-F,-E]/h^2;
 A4 = [A41;A42n2;A4n];
 
 u4 = -A4\f;
 u4 = [0;u4;0];
 
 figure(1), clf
 plot(xx,true_u(xx), 'k-','linewidth',2)
 hold on
 plot(x,u2, 'r.-','linewidth',2,'markersize',24)
 plot(x,u4, 'b.--','linewidth',2,'markersize',24)
 xlabel('$x$','fontsize',18,'interpreter','latex')
 ylabel('u(x)','fontsize',18,'interpreter','latex')
 leg = legend('true solution', 'quadratic approx','quartic approx','location','northoutside','orientation','horizontal');
 set(leg,'interpreter','latex','fontsize',10)
 print('quartic_approx','-dpng');

 figure(2), clf
 plot(x,abs(true_u(x)-u2), 'r.--','linewidth',2,'markersize',24)
 hold on;
 plot(x,abs(true_u(x)-u4), 'b.--','linewidth',2,'markersize',24)
 xlabel('$x$','fontsize',18,'interpreter','latex')
 ylabel('$|{\rm error}|$','fontsize',18,'interpreter','latex')
 leg = legend('quadratic error','quartic error','location','northoutside','orientation','horizontal');
 set(leg,'interpreter','latex','fontsize',10)
 print('quartic_error','-dpng');
 
 