close all;
clear all;
delete all;
clc;
mu0=4*pi*1E-7;
Br=1; 
V=3.9270e-07;
m=Br*V/mu0;
x=0.17;
y=0;
z=0;
theta=pi/4;
phi=pi/2;

alpha = pi/20;
beta = pi/25;
gamma = pi/50;

% syms x y z theta phi 
% B1x=0.268639712609443;
% B1y=-1.61430149711401;
% B1z=-0.149575291787213;
% B2x=0.411181903491468;
% B2y=-1.72545709227836;
% B2z=-0.135157336695423;

% eq1 = mu0/(4*pi)*(3*x*((m*sin(theta)*cos(phi))*x+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt(x^2 + y^2+ z^2))^5-(m*sin(theta)*cos(phi))/(sqrt(x^2 + y^2+ z^2))^3) == B1x;
% eq2 = mu0/(4*pi)*(3*y*((m*sin(theta)*cos(phi))*x+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt(x^2 + y^2+ z^2))^5-(m*sin(theta)*sin(phi))/(sqrt(x^2 + y^2+ z^2))^3) == B1y;
% eq3 = mu0/(4*pi)*(3*z*((m*sin(theta)*cos(phi))*x+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt(x^2 + y^2+ z^2))^5-(m*cos(theta))/(sqrt(x^2 + y^2+ z^2))^3) == B1z;
% eq4 = mu0/(4*pi)*(3*(x-2)*((m*sin(theta)*cos(phi))*(x-2)+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt((x-2)^2 + y^2+ z^2))^5-(m*sin(theta)*cos(phi))/(sqrt((x-2)^2 + y^2+ z^2))^3) == B2x;
% eq5 = mu0/(4*pi)*(3*y*((m*sin(theta)*cos(phi))*(x-2)+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt((x-2)^2 + y^2+ z^2))^5-(m*sin(theta)*sin(phi))/(sqrt((x-2)^2 + y^2+ z^2))^3) == B2y;
% eq6 = mu0/(4*pi)*(3*z*((m*sin(theta)*cos(phi))*(x-2)+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt((x-2)^2 + y^2+ z^2))^5-(m*cos(theta))/(sqrt((x-2)^2 + y^2+ z^2))^3) == B2z;
% eqns=[eq1,eq2,eq3,eq4,eq5,eq6];
% 
% 
% S = solve(eqns,[x,y,z,theta,phi]);

[B1x,B1y,B1z,Bsum,ang1,ang2,Besum]=magnetic_field(m,x,y,z,theta,phi,alpha,beta,gamma);
[B2x,B2y,B2z,Bsum,ang1,ang2,Besum]=magnetic_field(m,(x-0.02),y,z,theta,phi,alpha,beta,gamma);

a0=[0.1 0.1 0.1 1 1 1 1 1];
opts = optimoptions(@fsolve,'Algorithm', 'trust-region','Display','none','PlotFcn',@optimplotfirstorderopt);%optimset('MaxIter',200000,'MaxFunEvals',2000000,'TolFun',1e-10,'TolX',1e-10);%optimoptions(@fsolve,'MaxFunctionEvaluations',500*numel(a0));
Inet = fsolve(@(a)solve_B(a,B1x,B1y,B1z,B2x,B2y,B2z,Besum,ang1,ang2),a0,opts);
% plot(V,Inet)
% Inet Function
function F = solve_B(a,B1x,B1y,B1z,B2x,B2y,B2z,Besum,ang1,ang2)
    mu0=4*pi*1E-7;
    Br=1; 
    V=3.9270e-07;
    m=Br*V/mu0;
    Bd=48.2862; %uT @ Gothernburg

    x = a(1);
    y = a(2);
    z = a(3);
    theta = a(4);
    phi = a(5);
    Bex = a(6);
    Bey = a(7);
    Bez = a(8);
    
    eq1 = mu0/(4*pi)*(3*x*((m*sin(theta)*cos(phi))*x+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt(x^2 + y^2+ z^2))^5-(m*sin(theta)*cos(phi))/(sqrt(x^2 + y^2+ z^2))^3) * 1e6 + Bex;
    eq2 = mu0/(4*pi)*(3*y*((m*sin(theta)*cos(phi))*x+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt(x^2 + y^2+ z^2))^5-(m*sin(theta)*sin(phi))/(sqrt(x^2 + y^2+ z^2))^3) * 1e6 + Bey;
    eq3 = mu0/(4*pi)*(3*z*((m*sin(theta)*cos(phi))*x+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt(x^2 + y^2+ z^2))^5-(m*cos(theta))/(sqrt(x^2 + y^2+ z^2))^3) * 1e6 + Bez;
    eq4 = mu0/(4*pi)*(3*(x-0.02)*((m*sin(theta)*cos(phi))*(x-0.02)+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt((x-0.02)^2 + y^2+ z^2))^5-(m*sin(theta)*cos(phi))/(sqrt((x-0.02)^2 + y^2+ z^2))^3) * 1e6 + Bex;
    eq5 = mu0/(4*pi)*(3*y*((m*sin(theta)*cos(phi))*(x-0.02)+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt((x-0.02)^2 + y^2+ z^2))^5-(m*sin(theta)*sin(phi))/(sqrt((x-0.02)^2 + y^2+ z^2))^3) * 1e6 + Bey;
    eq6 = mu0/(4*pi)*(3*z*((m*sin(theta)*cos(phi))*(x-0.02)+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt((x-0.02)^2 + y^2+ z^2))^5-(m*cos(theta))/(sqrt((x-0.02)^2 + y^2+ z^2))^3) * 1e6 + Bez;
    eq7 = sqrt(Bex^2+Bey^2+Bez^2);
    eq8 = Bex*sin(ang1)*sin(-1*ang2)+Bey*cos(ang2)*sin(-1*ang1)+Bez*cos(ang1);

    F = [eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]-[B1x;B1y;B1z;B2x;B2y;B2z;Besum;Bd];
end

% r = sqrt(x^2 + y^2+ z^2);
% Mx=m*sin(theta)*cos(phi);
% My=m*sin(theta)*sin(phi);
% Mz=m*cos(theta);
% MdotR=(m*sin(theta)*cos(phi))*x+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z;
% Bx=mu0/(4*pi)*(3*x*((m*sin(theta)*cos(phi))*x+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt(x^2 + y^2+ z^2))^5-(m*sin(theta)*cos(phi))/(sqrt(x^2 + y^2+ z^2))^3);
% By=mu0/(4*pi)*(3*y*((m*sin(theta)*cos(phi))*x+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt(x^2 + y^2+ z^2))^5-(m*sin(theta)*sin(phi))/(sqrt(x^2 + y^2+ z^2))^3);
% Bz=mu0/(4*pi)*(3*z*((m*sin(theta)*cos(phi))*x+(m*sin(theta)*sin(phi))*y+(m*cos(theta))*z)/(sqrt(x^2 + y^2+ z^2))^5-(m*cos(theta))/(sqrt(x^2 + y^2+ z^2))^3);
% Bsum = sqrt(Bx^2 + By^2+ Bz^2);










