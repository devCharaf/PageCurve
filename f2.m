%% looping number
N = 1000;

%% Choose parameters
T = 2;

L1 = 100;
L2 = 1;
L_ = [L1 L2];

l1 = 1;
l2 = 10;
l_ = [l1 l2];

mu_lim = -10; %xlimit of \mu

%% Computed parameters (through formulas)
M = (2*pi*T)^2;

lambda_max = 1/l1 + 1/l2;
lambda_min = 1/l1 - 1/l2;
lambda_0 = sqrt(lambda_max * lambda_min);


lambda = lambda_min + 1/(100*l2); %lambda

A = (lambda_max^2-lambda^2)*(lambda^2-lambda_min^2);
B = 2*lambda^2*M;

rh_1 = sqrt(M*l1^2);
rh_2 = sqrt(M*l2^2);
rh_ = [rh_1 rh_2];

%%

mu = linspace(mu_lim,1,N);
s1 = (-lambda^2*(1 + mu) + lambda_0^2*(1 - mu) + ...
   2*lambda*sqrt((1 - mu)/l2^2 + (mu.^2 - mu)/l1^2 + mu*lambda^2))/A;
s2 = (-lambda^2*(1 + mu) + lambda_0^2*(1 - mu) - ...
   2*lambda*sqrt((1 - mu)/l2^2 + (mu.^2 - mu)/l1^2 + mu*lambda^2))/A;
ff2 = zeros(1,N);


for i=1:N
    
    int_f2 = @(x) l2/sqrt(A) * ((lambda^2 - lambda_0^2)*x - 1 + mu(i))./((x + mu(i)*l2^2).*sqrt(x.*(x - s1(i)).*(x - s2(i))));
    ff2(i) = integral(int_f2,s1(i),Inf);
    
end

figure(1)
plot(mu,ff2, 'r','LineWidth', 3)
xlim([mu_lim 1])
xlabel('$\mu$', 'Interpreter','latex')
ylabel('$f_2(\mu)$', 'Interpreter','latex');
