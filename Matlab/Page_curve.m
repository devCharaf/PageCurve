%% Page curve

%% looping number
N = 100000;

%% Choose parameters
T = 2;

L1 = 100;
L2 = 1;
L_ = [L1 L2];

l1 = 1;
l2 = 10;
l_ = [l1 l2];

epsilon = 0.000000001;

%% Computed parameters (through formulas)
M = (2*pi*T)^2;

lambda_max = 1/l1 + 1/l2;
lambda_min = 1/l1 - 1/l2;
lambda_0 = sqrt(lambda_max * lambda_min);


lambda = lambda_min + 1/(100*l2); %lambda

A = (lambda_max^2-lambda^2)*(lambda^2-lambda_min^2);
B = 2*lambda^2*M;
sigma_1 = 0;
sigma_2 = -2*B/A;

rh_1 = sqrt(M*l1^2);
rh_2 = sqrt(M*l2^2);
rh_ = [rh_1 rh_2];

% Constants used in phi_sigma
alpha_1 = l1*(lambda^2+lambda_0^2)/(2*lambda);
alpha_2 = l2*(lambda^2-lambda_0^2)/(2*lambda);
alpha_ = [alpha_1 alpha_2];

a1 = A*l1^2/(4*lambda^2);
a2 = A*l2^2/(4*lambda^2);
a_ = [a1 a2];

%% q to t
t = linspace(0,N,N);
cosq = sqrt( (1-rh_1^2*epsilon^2*tanh(rh_1*t/l1).^2) ./ (1-rh_1^4*epsilon^4*tanh(rh_1*t/l1).^2));
invcosq = 1+epsilon^2*rh_2^2*tanh(rh_2*t/l2)/2;

%% curve plot (Do not change this part)

sigma = linspace(0,N,10);
q = linspace(0,pi/2,N);

bb = (L1-L2/10)/2; % 

x = 1/sqrt(M) * atanh(alpha_1) + L1/2;
b = tanh(sqrt(M)*(x-bb));
J = b*l1*rh_1;
S = 4*l1*(log(2/epsilon) - log(sqrt(rh_1^2-J^2/l1^2)))*q./q; %time dependant entropy

Sq = 2*l1*log(2./(epsilon*cos(q))); % constant entropy


figure(1)
h1 = plot(q,Sq, 'r', 'LineWidth', 3);
hold on
h2 = plot(q,S,'b', 'LineWidth', 3);
xlabel('$q$', 'Interpreter','latex')
ylabel('$S$', 'Interpreter','latex');
legend([h1,h2], {'Constant entropy $S_1$', 'Time dependant entropy $S_2$'}, 'Interpreter','latex','Location','northwest')

