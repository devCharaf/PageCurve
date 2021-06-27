%% This MATLAB script plots the wall as a function of \sigma

%% Loop number
N = 100000;

%% Choose parameters
L1 = 100;
L2 = L1/100;
L_ = [L1 L2];

l1 = 1;
l2 = l1*10;
l_ = [l1 l2];

T = 2/L2;

%% Computed parameters
M = (2*pi*T)^2;
lambda_max = 1/l1 + 1/l2;
lambda_min = 1/l1 - 1/l2;
lambda_0 = sqrt(lambda_max * lambda_min);
lambda = lambda_min + 1/(100*l2);

A = (lambda_max^2-lambda^2)*(lambda^2-lambda_min^2);
B = 2*lambda^2*M;
sigma_1 = 0;
sigma_2 = -2*B/A;

rh_1 = sqrt(M*l1^2);
rh_2 = sqrt(M*l2^2);
rh_ = [rh_1 rh_2];

%% Constants used in \phi_sigma
alpha_1 = l1*(lambda^2+lambda_0^2)*rh_1/(2*lambda);
alpha_2 = l2*(lambda^2-lambda_0^2)*rh_2/(2*lambda);
alpha_ = [alpha_1 alpha_2];

a1 = A*l1^2/(4*lambda^2);
a2 = A*l2^2/(4*lambda^2);
a_ = [a1 a2];

%% plot
x = zeros(2,N);
sigma = linspace(0,10000000,N);
co = 'c';
alp = '';

for i=1:2
    
    if i== 2
        co = 'r';
    end
    l = l_(i);
    L = L_(i);
    alph = alpha_(i);
    a = a_(i);
    rh = rh_(i);
    
    x(i,:) = atanh(alph./sqrt(a*sigma+rh^2))/(2*pi*T)+L/2;
    
    
    figure(2)
    subplot(1,2,i)
    plot(sigma,x(i,:),co,'LineWidth', 3)
    hold on
    line([0,10000000],[L/2,L/2], 'LineStyle','--', 'color', [0.7,0.7,0.7], 'LineWidth',2)
    if i ==1
        ylim([L/2-x(1,1)/(100*L1) x(1,1)])
        alp = '$\alpha_1$';
    else
        ylim([x(2,1) L/2+x(2,1)/(100*L2)])
        alp = '$\alpha_2$';
    end
    xlabel('$\sigma$', 'Interpreter','latex')
    ylabel(alp, 'Interpreter','latex')
end
