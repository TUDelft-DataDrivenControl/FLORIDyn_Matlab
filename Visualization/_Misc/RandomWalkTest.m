mu = [20 50];
sigma = [10 0; 0 10];
R = chol(sigma);
% z = repmat(mu,10,1) + randn(10,2)*R;

simL = 1000;
x1 = zeros(simL,2);
x2 = x1;
x3 = x1;
for i = 2:simL
    rn = randn(1,2);
    rStep1 = rn + (mu-x1(i-1,:))/R;
    x1(i,:) = x1(i-1,:) + rStep1*R;
    
    rStep2 = rn + 0.1*(mu-x2(i-1,:))*inv(R);
    x2(i,:) = x2(i-1,:) + 1*rStep2*R;
    
    x3(i,:) = x2(i-1,:) + rn*R + 1*(mu-x2(i-1,:));
end


plot(x1(:,1))
hold on
plot(x2(:,1))
plot(x3(:,1)+1)
hold off
legend('x1','x2','x3')