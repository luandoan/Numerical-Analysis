% hw3 Numerical analysis
close all; clear all; clc; 

x = 0:0.01:1;
%x = 1/(36-24*sqrt(2));      % test minimum error value
fx = sqrt(x);
p1x = (sqrt(6)-sqrt(3))*x + (2-sqrt(2))/sqrt(3);
error = fx - p1x;
figure; plot(x,fx,'r'); hold on;
plot(x,p1x,'b'); grid on; legend('f(x)','p_1(x)');

figure; plot(x,error); grid on;
title('Error of \sqrt(x) - p_1(x) over [0,1]');
xlabel('x'); ylabel('Error');print('hw3_11a','-dpng');

p2x = 1/8 + x;
error2 = fx - p2x;
figure; plot(x,error2); grid on;
figure; plot(x,fx,'r'); hold on;
plot(x,p2x,'b'); grid on; legend('f(x)','p_1(x)');


%%