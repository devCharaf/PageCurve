%% This MATLAB script plots 3D Penrose diagrams

N = 100;

%% Choose parameters
L1 = 100;
L2 = L1/100;
L_ = [L1 L2];

l1 = 1;
l2 = 10*l1;
l_ = [l1 l2];

epsilon = 0.000001;
tau2 = 2;
%% Computed parameters
T = tau2/L2;
M = (2*pi*T)^2;
lambda_max = 1/l1 + 1/l2;
lambda_min = 1/l1 - 1/l2;
lambda_0 = sqrt(lambda_max * lambda_min);
lambda = lambda_min + 1/(100*l2);

A = (lambda_max^2-lambda^2)*(lambda^2-lambda_min^2);
B = 2*lambda^2*M;
sigma_1 = 0;
sigma_2 = -2*B/A;
r_h1 = sqrt(l1^2*M);
r_h2 = sqrt(l2^2*M);
rh_ = [r_h1 r_h2];
%% Constants used in real(phi_sigma)
alpha_1 = l1*(lambda^2+lambda_0^2)*r_h1/(2*lambda);
alpha_2 = l2*(lambda^2-lambda_0^2)*r_h2/(2*lambda);
alpha_ = [alpha_1 alpha_2];

a1 = A*l1^2/(4*lambda^2);
a2 = A*l2^2/(4*lambda^2);
a_ = [a1 a2];

%%
qq = linspace(-pi/2,pi/2,N);
pp = qq;
[PP,QQ] = meshgrid(pp,qq);

co =  [0.8500 0.3250 0.0980];

for i=1:2

    l = l_(i);
    L = L_(i);
    alpha = alpha_(i);
    a = a_(i);
    rh = rh_(i);
    
    X = real(atanh(alpha./sqrt(a*rh^2*((cos(QQ)./cos(PP)).^2-1)+rh^2))/(2*pi*T))+L/2;
    
    x = linspace(-L,L,N);
    p = linspace(-pi/2,pi/2,N);

    q = p;

    [P,XX] = meshgrid(p,x);
    Q = P;
    
    figure(i)
    if i == 2
        co = [0.3010 0.7450 0.9330];
    end
    h1 = surf(X,PP,QQ,'EdgeColor', 'none');
    set(h1,'FaceColor',co)
    hold on
    h2 = surf(-X,PP,QQ,'EdgeColor', 'none');
    set(h2,'FaceColor',co)
    xlabel('$x$', 'Interpreter','latex')
    ylabel('$p$', 'Interpreter','latex')
    zlabel('$q$', 'Interpreter','latex')
    %set(gca,'XColor', 'none','YColor','none','ZColor','none')
end

