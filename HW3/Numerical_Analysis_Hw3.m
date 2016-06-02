% Luan Cong Doan - Numerical Analysis Class - Homework 3
close all; clear all; clc;
%% Hw 3.1.a
%syms x;
x = 0:0.01:1;
%syms x;
fx = sqrt(x);
p1x = (sqrt(6)-sqrt(3))*x + (2-sqrt(2))/3;
error = fx - p1x;
figure; plot(x,error); grid on;
xlabel('x'); ylabel('error'); title('Error of sqrt(x) - p_1(x)');
print('hw3_11a','-dpng');

p2x = 1/8 + x;
error2 = fx - p2x;
figure; plot(x,error2); grid on;
xlabel('x'); ylabel('error from p*'); title('Error of sqrt(x) - p*(x)');
%print('hw3_11c','-dpng');

p3x = 4/15 + (4/5)*x;
error3 = fx - p3x;
figure; plot(x,error3); grid on;
xlabel('x'); ylabel('error from P*'); title('Error of sqrt(x) - P*(x)');
%print('hw3_21a','-dpng');

%xe = 1/x;
f2x = log(1) - log(x);
P0x = 1;
P1x = -3*x + 5/2;
P2x = 5*x.^2 - 8*x + 10/3;
P3x = -(35/3)*x.^3 + (45/2)*x.^2 - 15*x + 47/12;
figure; plot(x,f2x,'r','linewidth',2); hold on; grid on;
plot(x,ones(size(x))*1,'k','linewidth',1.5) %plot(x,P0x,'b'); 
plot(x,P1x,'g--'); plot(x,P2x,'b-'); plot(x,P3x,'m.');
legend('f(x)','P_0(x)','P_1(x)','P_2(x)','P_3(x)');
xlabel('x');ylabel('Least square approximation');
title('Four fisrt least square approximation');
print('hw3_34','-dpng');



    
