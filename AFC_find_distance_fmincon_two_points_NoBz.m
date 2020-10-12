close all;
clear all;
delete all;
clc;
mu0=4*pi*1E-7;
Br=1; 
V=3.9270e-07;
m=Br*V/mu0;
xarray = [0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];

for i = 1 : length (xarray)

    x=xarray(i);
    y1=0.05;
    y2=-0.05;
    z=0;
    theta=pi/2;
    phi1=pi/7;
    phi2=phi1+pi/6;

    alpha = pi/20;
    beta = pi/20;
    gamma = 0;

    [B1x,B1y,B1z,Bsum,ang1,ang2,Besum]=magnetic_field(m,x,y1,z,theta,phi1,alpha,beta,gamma);
    [B2x,B2y,B2z,Bsum,ang1,ang2,Besum]=magnetic_field(m,(x-0.02),y1,z,theta,phi1,alpha,beta,gamma);

    [A1x,A1y,A1z,Bsum,ang3,ang4,Besum]=magnetic_field(m,x,y2,z,theta,phi2,alpha,beta,gamma);
    [A2x,A2y,A2z,Bsum,ang3,ang4,Besum]=magnetic_field(m,(x-0.02),y2,z,theta,phi2,alpha,beta,gamma);

    %% Find the distance

    %    1     2     3     4     5        6   7   8    9    10    11      12  13  14  15    
    a0=[0.1   0.1   0.1   0.1    0.1      1   1   1   0.1  0.1    pi/20   1   1   1   1 ];
    lb=[0.05 -0.15 -0.1   0     -pi/2   -60 -60 -60 -0.15 -0.1   -pi/4   -60 -60 -60 -60];
    ub=[0.6   0.15  0.1   pi     pi/2    60  60  60  0.15  0.1    pi/4    60  60  60  60];

    opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt','Display',...
        'none','PlotFcn',@optimplotfirstorderopt,'MaxIterations',200000,...
        'MaxFunctionEvaluations',2000000,'FunctionTolerance',1e-10);%optimset('MaxIter',200000,'MaxFunEvals',2000000,'TolFun',1e-10,'TolX',1e-10);%optimoptions(@fsolve,'MaxFunctionEvaluations',500*numel(a0));
    Inet = fmincon(@(a)solve_B(a,B1x,B1y,B1z,B2x,B2y,B2z,Besum,ang1,ang2,A1x,A1y,A1z,A2x,A2y,A2z,ang3,ang4),a0,[],[],[],[],lb,ub);
    dist(i) = Inet(1);
    
end

figure;plot(dist,'o-');hold on
plot(xarray,'*--')
legend('Estimated distance','Distance')
xlabel('#No')
ylabel('Distance (m)')

%% Functions
function F = solve_B(a,B1x,B1y,B1z,B2x,B2y,B2z,Besum,ang1,ang2,A1x,A1y,A1z,A2x,A2y,A2z,ang3,ang4)
    mu0=4*pi*1E-7;
    Br=1; 
    V=3.9270e-07;
    m=Br*V/mu0;
    
    x = a(1);
    y1 = a(2);
    z1 = a(3);
    theta = a(4);
    phi1 = a(5);
    Bex1 = a(6);
    Bey1 = a(7);
    Bez1 = a(8);
    
    y2 = a(9);
    z2 = a(10);
    DeltaPhi = a(11);
    Bex2 = a(12);
    Bey2 = a(13);
    Bez2 = a(14);
    Bdd = a(15);
    
    eq1 = mu0/(4*pi)*(3*x*((m*sin(theta)*cos(phi1))*x+(m*sin(theta)*sin(phi1))*y1+(m*cos(theta))*z1)/(sqrt(x^2 + y1^2+ z1^2))^5-(m*sin(theta)*cos(phi1))/(sqrt(x^2 + y1^2+ z1^2))^3) * 1e6 + Bex1;
    eq2 = mu0/(4*pi)*(3*y1*((m*sin(theta)*cos(phi1))*x+(m*sin(theta)*sin(phi1))*y1+(m*cos(theta))*z1)/(sqrt(x^2 + y1^2+ z1^2))^5-(m*sin(theta)*sin(phi1))/(sqrt(x^2 + y1^2+ z1^2))^3) * 1e6 + Bey1;
    eq3 = mu0/(4*pi)*(3*z1*((m*sin(theta)*cos(phi1))*x+(m*sin(theta)*sin(phi1))*y1+(m*cos(theta))*z1)/(sqrt(x^2 + y1^2+ z1^2))^5-(m*cos(theta))/(sqrt(x^2 + y1^2+ z1^2))^3) * 1e6 + Bez1;
    eq4 = mu0/(4*pi)*(3*(x-0.02)*((m*sin(theta)*cos(phi1))*(x-0.02)+(m*sin(theta)*sin(phi1))*y1+(m*cos(theta))*z1)/(sqrt((x-0.02)^2 + y1^2+ z1^2))^5-(m*sin(theta)*cos(phi1))/(sqrt((x-0.02)^2 + y1^2+ z1^2))^3) * 1e6 + Bex1;
    eq5 = mu0/(4*pi)*(3*y1*((m*sin(theta)*cos(phi1))*(x-0.02)+(m*sin(theta)*sin(phi1))*y1+(m*cos(theta))*z1)/(sqrt((x-0.02)^2 + y1^2+ z1^2))^5-(m*sin(theta)*sin(phi1))/(sqrt((x-0.02)^2 + y1^2+ z1^2))^3) * 1e6 + Bey1;
    eq6 = mu0/(4*pi)*(3*z1*((m*sin(theta)*cos(phi1))*(x-0.02)+(m*sin(theta)*sin(phi1))*y1+(m*cos(theta))*z1)/(sqrt((x-0.02)^2 + y1^2+ z1^2))^5-(m*cos(theta))/(sqrt((x-0.02)^2 + y1^2+ z1^2))^3) * 1e6 + Bez1;
    eq7 = sqrt(Bex1^2+Bey1^2+Bez1^2);
    eq8 = Bex1*sin(ang1)*sin(-1*ang2)+Bey1*cos(ang2)*sin(-1*ang1)+Bez1*cos(ang1) - a(15);
    
    eq9 = mu0/(4*pi)*(3*x*((m*sin(theta)*cos((phi1+DeltaPhi)))*x+(m*sin(theta)*sin((phi1+DeltaPhi)))*y2+(m*cos(theta))*z2)/(sqrt(x^2 + y2^2+ z2^2))^5-(m*sin(theta)*cos((phi1+DeltaPhi)))/(sqrt(x^2 + y2^2+ z2^2))^3) * 1e6 + Bex2;
    eq10 = mu0/(4*pi)*(3*y2*((m*sin(theta)*cos((phi1+DeltaPhi)))*x+(m*sin(theta)*sin((phi1+DeltaPhi)))*y2+(m*cos(theta))*z2)/(sqrt(x^2 + y2^2+ z2^2))^5-(m*sin(theta)*sin((phi1+DeltaPhi)))/(sqrt(x^2 + y2^2+ z2^2))^3) * 1e6 + Bey2;
    eq11 = mu0/(4*pi)*(3*z2*((m*sin(theta)*cos((phi1+DeltaPhi)))*x+(m*sin(theta)*sin((phi1+DeltaPhi)))*y2+(m*cos(theta))*z2)/(sqrt(x^2 + y2^2+ z2^2))^5-(m*cos(theta))/(sqrt(x^2 + y2^2+ z2^2))^3) * 1e6 + Bez2;
    eq12 = mu0/(4*pi)*(3*(x-0.02)*((m*sin(theta)*cos((phi1+DeltaPhi)))*(x-0.02)+(m*sin(theta)*sin((phi1+DeltaPhi)))*y2+(m*cos(theta))*z2)/(sqrt((x-0.02)^2 + y2^2+ z2^2))^5-(m*sin(theta)*cos((phi1+DeltaPhi)))/(sqrt((x-0.02)^2 + y2^2+ z2^2))^3) * 1e6 + Bex2;
    eq13 = mu0/(4*pi)*(3*y2*((m*sin(theta)*cos((phi1+DeltaPhi)))*(x-0.02)+(m*sin(theta)*sin((phi1+DeltaPhi)))*y2+(m*cos(theta))*z2)/(sqrt((x-0.02)^2 + y2^2+ z2^2))^5-(m*sin(theta)*sin((phi1+DeltaPhi)))/(sqrt((x-0.02)^2 + y2^2+ z2^2))^3) * 1e6 + Bey2;
    eq14 = mu0/(4*pi)*(3*z2*((m*sin(theta)*cos((phi1+DeltaPhi)))*(x-0.02)+(m*sin(theta)*sin((phi1+DeltaPhi)))*y2+(m*cos(theta))*z2)/(sqrt((x-0.02)^2 + y2^2+ z2^2))^5-(m*cos(theta))/(sqrt((x-0.02)^2 + y2^2+ z2^2))^3) * 1e6 + Bez2;
    eq15 = sqrt(Bex2^2+Bey2^2+Bez2^2);
    eq16 = Bex2*sin(ang3)*sin(-1*ang4)+Bey2*cos(ang4)*sin(-1*ang3)+Bez2*cos(ang3) - a(15);

    f = [eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8;eq9;eq10;eq11;eq12;eq13;eq14;eq15;eq16]-[B1x;B1y;B1z;B2x;B2y;B2z;Besum;0;A1x;A1y;A1z;A2x;A2y;A2z;Besum;0];
    F = sum(abs(f));
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










