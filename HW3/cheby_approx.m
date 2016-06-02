function pn = cheby_approx(k1,k2,lambda1,lambda2)
lambda = lambda1:0.01:lambda2;
x1 = (lambda1+lambda2)/(lambda2 - lambda1);
x2 = (lambda1+lambda2 - 2*lambda)/(lambda2 - lambda1);
T11 = cheby(x1,k1);
T21 = cheby(x2,k1);
pn1 = T21/T11;
T12 = cheby(x1,k2);
T22 = cheby(x2,k2);
pn2 = T22/T12;
pn = [pn1;pn2];
figure; plot(x2,pn1,'r'); hold on; grid on;
plot(x2,pn2,'b--'); xlabel('x'); ylabel('approximation');
title(['Minimizing polynomial for \lambda_1 = ',num2str(lambda1),', \lambda_N = ',num2str(lambda2)]);
legend('for k = 3',' for k = 5');
print(['hw3_4',num2str(lambda1),num2str(lambda2)],'-dpng');
end
