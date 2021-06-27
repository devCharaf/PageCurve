%% This MATLAB script plots the wall in the compactified coordinates (rho,x).
% We choose the origin of x to be in the middle of the boundary.

N = 1000;
%% Choose parameters
T = 2;
lambda = 901/1000;

L1 = 4;
L2 = 1;
L_ = [L1 L2];

l1 = 0.1;
l2 = 0.15;
l_ = [l1 l2];

%% Computed parameters
M = (2*pi*T)^2;
lambda_max = 1/l1 + 1/l2;
lambda_min = 1/l2 - 1/l1;
lambda_0 = sqrt(lambda_max * lambda_min);
A = (lambda_max^2-lambda^2)*(lambda^2-lambda_min^2);
B = 2*lambda^2*M;
sigma_1 = 0;
sigma_2 = -2*B/A;

rh_1 = M*l1^2;
rh_2 = M*l2^2;
rh_ = [rh_1 rh_2];

%% Constants used in real(phi_sigma)
alpha_1 = l1*(lambda^2+lambda_0^2)*rh_1/(2*lambda);
alpha_2 = l2*(lambda^2-lambda_0^2)*rh_2/(2*lambda);
alpha_ = [alpha_1 alpha_2];

a1 = A*l1^2/(4*lambda^2);
a2 = A*l2^2/(4*lambda^2);
a_ = [a1 a2];

%% Wall parametrization
rho_l = linspace(pi/2,pi/2,N);
phi_l = linspace(0,2*pi,N); 

xl = rho_l .* cos(real(phi_l));
yl = rho_l .* sin(real(phi_l));
eps = -1;

for i=1:2
    
    if i==2
        eps=1;
    end
    l = l_(i);
    L = L_(i);
    alpha = alpha_(i);
    a = a_(i);
    rh = rh_(i);
    
    R = L/2;
    E = rh/tanh(rh*R/2);
    %phi1 = linspace(-R/2,R/2,N);
    %r1 = rh * real(sqrt(( 1- (tanh(rh*phi1)).^2 ) ./ (1 - rh^2/E^2 * (tanh(rh*phi1)).^2)));
    rho1 = linspace(atan(E),pi/2,N);
    phi1 = atan(rh/E*sqrt((tan(rho1).^2-E^2)./(tan(rho1).^2-rh^2)))/rh;
    
    x1 = rho1 .* cos(real(phi1));
    y1 = rho1 .* sin(real(phi1));

    %x1 = x1(400:N);
    %y1 = y1(400:N);
    
    sigma = 500;
    phis = atanh(alpha./sqrt(a*sigma+rh^2))/(2*pi*T) + L/2;
    % phis = real(alpha * atan(beta*sqrt(sigma-sigma_2)) + L/2 - alpha * pi/2);
    
    phi_bh = linspace(0,2*pi,N); 

    xb = atan(rh) .* cos(real(phi_bh));
    yb = atan(rh) .* sin(real(phi_bh));
    
    rho_wall = linspace(atan(rh),pi/2,N);
    phi_wall = atanh(alpha./sqrt(a*(tan(rho_wall).^2-rh^2)+rh^2))/(2*pi*T) + L/2;
    %phi_wall = real(alpha * atan(beta*sqrt(tan(rho_wall).^2-rh^2-sigma_2)) + L/2 - alpha * pi/2);
    
    x_wall = rho_wall .* cos(phi_wall);
    y_wall = rho_wall .* sin(phi_wall);
    
    n = 800;
    E2 = rh*tan(rho_wall(n))/sqrt(tanh(rh*phis)^2*(tan(rho_wall(n))^2-rh^2)+rh^2);
    phi2 = linspace(0,phi_wall(n),N);
    r2 = rh * real(sqrt(( 1- (tanh(rh*phi2)).^2 ) ./ (rh^2/E2^2 - (tanh(rh*phi2)).^2))); 
    
    x2 = atan(r2) .* cos(real(phi2));
    y2 = atan(r2) .* sin(real(phi2));
    
    rho_l2 = linspace(pi/2,pi/2,N);
    phi_l2 = linspace(-phi_wall(N),phi_wall(N),N); 

    xl2 = rho_l2 .* cos(real(phi_l2));
    yl2 = rho_l2 .* sin(real(phi_l2));
    
    
    figure(i)
    fill(xb,yb,'k');
    hold on
    h3 = plot(xl,yl, '--', 'color', [0.7,0.7,0.7], 'LineWidth',2);
    hold on
    if i ==1
        h4 = plot(eps*x_wall,y_wall, 'color',[0.8500 0.3250 0.0980], 'LineWidth',3);
        hold on
        plot(eps*x_wall,-y_wall, 'color',[0.8500 0.3250 0.0980], 'LineWidth',3)
        hold on
        plot(eps*xl2,yl2, 'r-', 'LineWidth',3)
    else
        h4 = plot(eps*x_wall,y_wall, 'c-', 'LineWidth',3);
        hold on
        plot(eps*x_wall,-y_wall, 'c-', 'LineWidth',3)
        hold on
        plot(eps*xl2,yl2, 'b-', 'LineWidth',3)
    end
    axis equal
    legend([h3,h4], {'boundary $\rho=\frac{\pi}{2}$','wall'}, 'Interpreter','latex','Location','northwest')
    set(gca,'XColor', 'none','YColor','none')
   
end
