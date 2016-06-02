% LUAN CONG DOAN - Numerical Analysis - Homework 4
close all; clear all; clc;
%% 4.1a Clenshaw-Curtis quadrature
int1 = sqrt(pi)*erf(1);
int2 = 2*atan(5)/5;
int3 = 1;
N = 1:1:50;
int_approx1 = zeros(1,length(N));
int_approx2 = zeros(1,length(N));
int_approx3 = zeros(1,length(N));
error1 = zeros(1,length(N));
error2 = zeros(1,length(N));
error3 = zeros(1,length(N));

for N=1:length(N)
    [x,w] = clencurt(N);
    for j = 1:length(x)
        int_approx1(N) = int_approx1(N) + exp(-x(j).^2)*w(j);
        int_approx2(N) = int_approx2(N) + 1/(1+25*x(j).^2)*w(j);
        int_approx3(N) = int_approx3(N) + abs(x(j))*w(j);
    end
    error1(N) = abs(int_approx1(N) - int1);
    error2(N) = abs(int_approx2(N) - int2);
    error3(N) = abs(int_approx3(N) - int3);
end
figure; semilogy(error1); grid on;
xlabel('N'); ylabel('error'); 
title('Clenshaw-Curtis approximation error of exp(-x^2)');
print('hw4_1a1','-dpng');
figure; semilogy(error2); grid on;
xlabel('N'); ylabel('error'); 
title('Clenshaw-Curtis approximation error of (1+25*x^2)^(-1)');
print('hw4_1a2','-dpng');
figure; semilogy(error3); grid on;
xlabel('N'); ylabel('error'); 
title('Clenshaw-Curtis approximation error of |x|');
print('hw4_1a3','-dpng');

figure; semilogy(error1,'r--.'); hold on;
semilogy(error2,'b--'); semilogy(error3,'k'); grid on;
xlabel('N'); ylabel('error'); 
title('Clenshaw-Curtis approximation');
legend('f(x) = exp(-x^2)','(1+25x^2)^{-1}','|x|');
print('hw4_1a','-dpng');

%% tic and toc
tic
int_approx = 0;
[x,w] = clencurt(20);
    for j = 1:length(x)
        int_approx = int_approx + exp(-x(j).^2)*w(j);
    end
toc

%% tic toc quad function
tic
    f = @(x) exp(-x.^2);
    q = quad(f,-1,1);
toc
    
    
%% Problem 4.2f1
x = 0:0.01:2*pi;
Ifx1 = zeros(1,length(x));
for k =1:length(x)
    Ifx1(k) = Ifx1(k) + ((2*pi)/N)* exp(sin(2*pi*k/length(x)));
end
errorF1 = Ifx1 - 7.95492652101284527451322;
figure; loglog(errorF1); grid on;
xlabel('x'); ylabel('error');
title('loglog plot of error in trapezoid approximation of exp(sin(x))');
print('hw4_2f1','-dpng');

% error = p - 7.95492652101284527451322;
%% Problem 4.2f1 - way2
N = 1:1:50;
Ifx1 = zeros(1,length(N));
for i =1:length(N)
    for k=1:i
    Ifx1(i) = Ifx1(i) + ((2*pi)/i)* exp(sin((2*pi*k)/i));
    end
end
errorF1 = abs(Ifx1 - 7.9549265210128452745132) ;
figure; loglog(N,errorF1); grid on;
xlabel('N'); ylabel('error');
title('loglog plot of error in trapezoid approximation of exp(sin(x))');
print('hw4_2f1','-dpng');

%%
%N = 10^5;
	fx1 = @(x) exp(sin(x));
	fx2 = @(x) exp(sin(x/pi));
%	Ifx1 = zeros(1,length(N));
%	Ifx2 = zeros(1,length(N));
	for i =1:5
        N = 10^(i)
		Ifx1(i) = trapezoid(fx1,0,2*pi,N);
		Ifx2(i) = trapezoid(fx2,0,2*pi,N);
	end
	errorF1 = abs(Ifx1 - 7.9549265210128452745132) ;
	figure; loglog(errorF1); grid on;
	xlabel('N'); ylabel('error');
	title('loglog plot of error in trapezoid approximation of exp(sin(x))');
	print('hw4_2f1','-dpng');
		
	errorF2 = abs(Ifx2 - 13.3094551602297896414536) ;
	figure; loglog(errorF2); grid on;
	xlabel('N'); ylabel('error');
	title('loglog plot of error in trapezoid approximation of exp(sin(x/pi))');
	print('hw4_2f2','-dpng');
    
%% Problem 4.2f2 way 2
N = 1:1:50;
Ifx2 = zeros(1,length(N));
for i =1:length(N)
    h = 2*pi/i;
    Ifx2(i) = (h/2)*(exp(sin(0)) + exp(sin(2*pi)));
    for k = 1:length(N)-1
    Ifx2(i) = Ifx2(i) + h* exp(sin((2*pi*k)/(pi*i)));
    end
end
errorF2 = Ifx2 - 13.3094551602297896414536 ;
figure; loglog(errorF2); grid on;
xlabel('N'); ylabel('error');
title('loglog plot of error in trapezoid approximation of exp(sin(x/pi))');
print('hw4_2f2','-dpng');

%% Problem 4.2f
fx1 = @(x) exp(sin(x));
fx2 = @(x) exp(sin(x/pi));
it1 = 7.9549265210128452745132;
it2 = 13.3094551602297896414536;
for i =1:5 %length(N)
    N = 10^(i);
    Ifx1(i) = trapezoid(fx1,0,2*pi,N);
    
end
errorF1 = abs(Ifx1 - it1);

for i =1:5 %length(N)
    N = 10^(i);
    Ifx2(i) = trapezoid(fx2,0,2*pi,N);
end
errorF2 = abs(Ifx2 - it2);
figure; loglog(errorF1,'r'); hold on;
loglog(errorF2,'b--'); grid on;
xlabel('N'); ylabel('error');
title('loglog plot of error in trapezoid approximation of exp(sin(x))');
legend('Error 1','Error 2');
print('hw4_2f','-dpng');


%% Problem 4.4a
N = 10e4;
fx = @(x) sin(x);
int = integral(fx,0,2);

% Monte Carlo points
IfM = zeros(1,N);
for i = 1:N
    xi = rand(1,i)*2;
    for j = 1:i
    IfM(i) = IfM(i) + (2/i)*sin(xi(j));
    end
end
errorM = abs(IfM - int);
figure; loglog(errorM); grid on;
xlabel('N'); ylabel('error');
title('loglog plot of error by using Monte Carlo method for f(x) = sin(x)');
print('hw4_4a1','-dpng');
%figure; plot(N,ones(size(N))*int); hold on;
%plot(N,IfM); grid on;

%% trapezoid 
N = 10^5;
fx = @(x) sin(x);
int = integral(fx,0,2);
IfT = zeros(1,N);
for i = 1:N
IfT(i) = trapezoid(fx,0,2,i);
    %for j = 1:length(x)
    %    IfT(i) = IfT(i) + sin(x(j))*w(j);
    %end
end
errorT = abs(IfT - int);
figure; loglog(errorT); grid on;
xlabel('N'); ylabel('error');
title('loglog plot of error by using Trapezoid method for f(x) = sin(x)');
print('hw4_4a2','-dpng');

figure; loglog(errorM); hold on; 
loglog(errorT); grid on;
xlabel('N'); ylabel('error');
title('loglog plot of error from Monte Carlo and Trapezoid for f(x) = sin(x)');
legend('Monte Carlo','Trapezoid');
print('hw4_4a','-dpng');

%% test function
s = -1:0.01:1;
expx = exp(-x.^2);
fracx = 1/(1+25*x.^2);
absx = abs(x);
figure; plot(x,expx,'r'); hold on;
plot(x,fracx,'b--'); plot(x,absx,'g-o'); grid on;

%% Problem 4.4b
%symsf = @(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10); %sin(x1*x2*x3*x4*x5*x6*x7*x8*x9*x10);

%%
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 a b n;
h = (b-a)/n;
fsin = sin(x1*x2*x3*x4*x5*x6*x7*x8*x9*x10);
for i = 1:10
    x(i) = 
%fsin(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10);

%fsin = sin(x1*x2*x3);

%intfsin = sin(x1*x2*a) + sin(x1*x2*b) + sin(x1*x2*(a+b)/n);
end


%%
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 a b n;
x = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10];
f = f(x);
h = (b-a)/n;
for i = 1:10
    for j = 1:n
    x(i) = h/2*(f(a,x2,x3,x4,x5,x6,x7,x8,x9,x10) + f(b,x2,x3,x4,x5,x6,x7,x8,x9,x10) + ...
        + f(a+j*h,x2,x3,x4,x5,x6,x7,x8,x9,x10));
    end
end

    
%%
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10;
x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10];
xn = zeros(1,10);
fx = sin(prod(x(:)));
n = [3,4,5];
h(i) = 2*pi/n(i);

for i = 1:3
    for j = 1:10
        for k = 1:n(i)-1
            xn(j) = sin(prod(x(:))*0/x(j)) + sin((prod(x(:))*2*pi)/x(j)) + ...
                    + sin((prod(x(:))*k*h(i))/x(j));
            x(x==x(j)) = xn(j);
            fx = h*fx;
        end
    end
end

%%
tic
n = [3,4,5];
sum = [0,0,0];
e = [0,0,0];
for i = 1:length(n)
    N = n(i)^10;
    sum(i) = 0;
        for j = 1:N
            xi = rand(1,10)*2;
            fx = sin(prod(xi(:)));
            sum(i) = sum(i) + fx;
        end
        int(i) = ((2^10)*sum(i))/N;
end
toc


