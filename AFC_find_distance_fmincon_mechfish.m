close all;
clear all;
delete all;
clc;
mu0=4*pi*1E-7;
Br=1; 
V=3.9270e-07;
m=Br*V/mu0;
load('testset.mat');

points = [2369 3903 5448 6923 8381 9624 11120];
xarray = [0.15 0.2  0.25 0.3  0.35 0.4  0.45 ];

for i = 1 : length(points)
    
    params = testset(points(i),:);

    gx = params(1)/10000;
    gy = params(2)/10000;
    gz = params(3)/10000;
    B1x = params(4);
    B1y = params(5);
    B1z = params(6);
    B2x = params(7);
    B2y = params(8);
    B2z = params(9);
    % Besum = 50.9162266438706; %uT @ Gothenburg
    Besum = 57.26; %uT @ Lab

    %% Find angles mapping Earth magnetic field in the direction of gravity (Bz)
    opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt','Display','none','PlotFcn',@optimplotfirstorderopt);%optimset('MaxIter',200000,'MaxFunEvals',2000000,'TolFun',1e-10,'TolX',1e-10);%optimoptions(@fsolve,'MaxFunctionEvaluations',500*numel(a0));
    b0=[0.1 0.1];
    angles = fsolve(@(b)solve_g(b, gx, gy, gz),b0,opts);
    ang1=angles(1);
    ang2=angles(2);

    %%

    a0=[0.1   0.1  0.1  0.1    0.1     1   1   1 ];
    lb=[0.05 -0.1 -0.1  0      0      -60 -60 -60];
    ub=[0.6   0.1  0.1  pi     pi      60  60  60];
    
    if i > 1
        a0(4)=Inet0(4);
        lb(4)=Inet0(4)-pi/10;
        ub(4)=Inet0(4)+pi/10;
%         lb(5)=Inet0(5)-pi/5;
%         ub(5)=Inet0(5)+pi/5;
    end
    
    opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt','Display',...
        'none','PlotFcn',@optimplotfirstorderopt,'MaxIterations',200000,...
        'MaxFunctionEvaluations',2000000,'FunctionTolerance',1e-10);%optimset('MaxIter',200000,'MaxFunEvals',2000000,'TolFun',1e-10,'TolX',1e-10);%optimoptions(@fsolve,'MaxFunctionEvaluations',500*numel(a0));
    Inet = fmincon(@(a)solve_B(a,B1x,B1y,B1z,B2x,B2y,B2z,Besum,ang1,ang2),a0,[],[],[],[],lb,ub);
    
    if i==1
        Inet0 = Inet;
    end
    
    dist(i) = Inet(1);
    thetas(i) = Inet(4);
    phis(i) = Inet(5);
end

figure;plot (dist,'o--'); hold on
plot (thetas,'+--');
plot(phis,'*--')
figure;plot(dist,'o-');hold on
plot(xarray,'*--')
legend('Estimated distance','Distance')
xlabel('#No')
ylabel('Distance (m)')

%%  Functions
function F = solve_B(a,B1x,B1y,B1z,B2x,B2y,B2z,Besum,ang1,ang2)
    mu0=4*pi*1E-7;
    Br=1; 
    V=3.9270e-07;
    m=Br*V/mu0;
    %Bd=48.2862; %uT @ Gothernburg
    Bd=56; %uT @ Lab

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

    f = [eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]-[B1x;B1y;B1z;B2x;B2y;B2z;Besum;Bd];
    F = sum(abs(f));
end

function ge = solve_g(b, gx, gy, gz)
    ge = [cos(b(1)) sin(b(1))*cos(b(2)) sin(b(1))*sin(b(2))] - [gz gy gx];
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










