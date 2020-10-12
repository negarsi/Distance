close all;
clear all;
delete all;
clc;
mu0=4*pi*1E-7;
Br=1; 
V=3.9270e-07;
m=Br*V/mu0;

load('testset.mat');

% points1 = [2369 3903 5448 6923 8381 9624 11120];
% points2 = [2616 4203 5555 7045 8471 9890 11320];
points1 =  [1647 3186 4749 6150 7602 9078 10580 12000];
points2 =  [1616 3235 4689 6132 7762 9102 10560 11970];
realdist = [0.15 0.2  0.25 0.3  0.35 0.4  0.45  0.5  ];

for i = 1 : length (points1)

    params1 = testset(points1(i),:);
    params2 = testset(points2(i),:);

    g1x = params1(1)/10000;
    g1y = params1(2)/10000;
    g1z = params1(3)/10000;
    B1x = params1(4);
    B1y = params1(5);
    B1z = params1(6);
    B2x = params1(7);
    B2y = params1(8);
    B2z = params1(9);

    g2x = params2(1)/10000;
    g2y = params2(2)/10000;
    g2z = params2(3)/10000;
    A1x = params2(4);
    A1y = params2(5);
    A1z = params2(6);
    A2x = params2(7);
    A2y = params2(8);
    A2z = params2(9);

    % Besum = 50.9162266438706; %uT @ Gothenburg
    Besum = 57.26; %uT @ Lab
    
    %% Find angles mapping Earth magnetic field in the direction of gravity (Bz)
    opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt','Display','none','PlotFcn',@optimplotfirstorderopt);%optimset('MaxIter',200000,'MaxFunEvals',2000000,'TolFun',1e-10,'TolX',1e-10);%optimoptions(@fsolve,'MaxFunctionEvaluations',500*numel(a0));
    b0=[0.1 0.1];
    angles = fsolve(@(b)solve_g(b, g1x, g1y, g1z),b0,opts);
    ang1=angles(1);
    ang2=angles(2);
    angles = fsolve(@(b)solve_g(b, g2x, g2y, g2z),b0,opts);
    ang3=angles(1);
    ang4=angles(2);


    %% Distance estimation
    %    1     2     3     4     5        6   7    8    9     10    11    12  13  14     
    a0=[0.1    0.1   0.1   pi/2   2*pi      1   1    1   0.1   0.1   pi/20  1   1   1 ];
    lb=[0.05  -0.15 -0.1   0     1.8*pi    -60 -60  -60 -0.15 -0.1  -pi/2  -60 -60 -60];
    ub=[0.6    0.15  0.1   pi    3*pi     60  60   60  0.15  0.1   pi/2   60  60  60];

    opts = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','Display',...
        'none','PlotFcn',@optimplotfirstorderopt,'MaxIterations',200000,...
        'MaxFunctionEvaluations',2000000,'FunctionTolerance',1e-10);%optimset('MaxIter',200000,'MaxFunEvals',2000000,'TolFun',1e-10,'TolX',1e-10);%optimoptions(@fsolve,'MaxFunctionEvaluations',500*numel(a0));
    for j = 1 : 5
        Inet = fmincon(@(a)solve_B(a,B1x,B1y,B1z,B2x,B2y,B2z,Besum,ang1,ang2,A1x,A1y,A1z,A2x,A2y,A2z,ang3,ang4),a0,[],[],[],[],lb,ub);
        a0 = Inet
    end
    dist(i) = Inet(1);
    thetas(i) = Inet(4);
    phis(i) = Inet(5);
    
end

figure;plot (dist,'o--'); hold on
plot (thetas,'+--');
plot(phis,'*--')
figure;plot(dist,'o-');hold on
plot(realdist,'*--')
legend('Estimated distance','Distance')
xlabel('#No')
ylabel('Distance (m)')

%% Functions
function F = solve_B(a,B1x,B1y,B1z,B2x,B2y,B2z,Besum,ang1,ang2,A1x,A1y,A1z,A2x,A2y,A2z,ang3,ang4)
    mu0=4*pi*1E-7;
    Br=1; 
    V=3.9270e-07;
    m=Br*V/mu0;
    %Bd=48.2862; %uT @ Gothernburg
    Bd=56; %uT @ Lab

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
    
    eq1 = mu0/(4*pi)*(3*x*((m*sin(theta)*cos(phi1))*x+(m*sin(theta)*sin(phi1))*y1+(m*cos(theta))*z1)/(sqrt(x^2 + y1^2+ z1^2))^5-(m*sin(theta)*cos(phi1))/(sqrt(x^2 + y1^2+ z1^2))^3) * 1e6 + Bex1;
    eq2 = mu0/(4*pi)*(3*y1*((m*sin(theta)*cos(phi1))*x+(m*sin(theta)*sin(phi1))*y1+(m*cos(theta))*z1)/(sqrt(x^2 + y1^2+ z1^2))^5-(m*sin(theta)*sin(phi1))/(sqrt(x^2 + y1^2+ z1^2))^3) * 1e6 + Bey1;
    eq3 = mu0/(4*pi)*(3*z1*((m*sin(theta)*cos(phi1))*x+(m*sin(theta)*sin(phi1))*y1+(m*cos(theta))*z1)/(sqrt(x^2 + y1^2+ z1^2))^5-(m*cos(theta))/(sqrt(x^2 + y1^2+ z1^2))^3) * 1e6 + Bez1;
    eq4 = mu0/(4*pi)*(3*(x-0.02)*((m*sin(theta)*cos(phi1))*(x-0.02)+(m*sin(theta)*sin(phi1))*y1+(m*cos(theta))*z1)/(sqrt((x-0.02)^2 + y1^2+ z1^2))^5-(m*sin(theta)*cos(phi1))/(sqrt((x-0.02)^2 + y1^2+ z1^2))^3) * 1e6 + Bex1;
    eq5 = mu0/(4*pi)*(3*y1*((m*sin(theta)*cos(phi1))*(x-0.02)+(m*sin(theta)*sin(phi1))*y1+(m*cos(theta))*z1)/(sqrt((x-0.02)^2 + y1^2+ z1^2))^5-(m*sin(theta)*sin(phi1))/(sqrt((x-0.02)^2 + y1^2+ z1^2))^3) * 1e6 + Bey1;
    eq6 = mu0/(4*pi)*(3*z1*((m*sin(theta)*cos(phi1))*(x-0.02)+(m*sin(theta)*sin(phi1))*y1+(m*cos(theta))*z1)/(sqrt((x-0.02)^2 + y1^2+ z1^2))^5-(m*cos(theta))/(sqrt((x-0.02)^2 + y1^2+ z1^2))^3) * 1e6 + Bez1;
    eq7 = sqrt(Bex1^2+Bey1^2+Bez1^2);
    eq8 = Bex1*sin(-1*ang1)*sin(-1*ang2)+Bey1*cos(ang2)*sin(-1*ang1)+Bez1*cos(ang1);
    
    eq9 = mu0/(4*pi)*(3*x*((m*sin(theta)*cos((phi1+DeltaPhi)))*x+(m*sin(theta)*sin((phi1+DeltaPhi)))*y2+(m*cos(theta))*z2)/(sqrt(x^2 + y2^2+ z2^2))^5-(m*sin(theta)*cos((phi1+DeltaPhi)))/(sqrt(x^2 + y2^2+ z2^2))^3) * 1e6 + Bex2;
    eq10 = mu0/(4*pi)*(3*y2*((m*sin(theta)*cos((phi1+DeltaPhi)))*x+(m*sin(theta)*sin((phi1+DeltaPhi)))*y2+(m*cos(theta))*z2)/(sqrt(x^2 + y2^2+ z2^2))^5-(m*sin(theta)*sin((phi1+DeltaPhi)))/(sqrt(x^2 + y2^2+ z2^2))^3) * 1e6 + Bey2;
    eq11 = mu0/(4*pi)*(3*z2*((m*sin(theta)*cos((phi1+DeltaPhi)))*x+(m*sin(theta)*sin((phi1+DeltaPhi)))*y2+(m*cos(theta))*z2)/(sqrt(x^2 + y2^2+ z2^2))^5-(m*cos(theta))/(sqrt(x^2 + y2^2+ z2^2))^3) * 1e6 + Bez2;
    eq12 = mu0/(4*pi)*(3*(x-0.02)*((m*sin(theta)*cos((phi1+DeltaPhi)))*(x-0.02)+(m*sin(theta)*sin((phi1+DeltaPhi)))*y2+(m*cos(theta))*z2)/(sqrt((x-0.02)^2 + y2^2+ z2^2))^5-(m*sin(theta)*cos((phi1+DeltaPhi)))/(sqrt((x-0.02)^2 + y2^2+ z2^2))^3) * 1e6 + Bex2;
    eq13 = mu0/(4*pi)*(3*y2*((m*sin(theta)*cos((phi1+DeltaPhi)))*(x-0.02)+(m*sin(theta)*sin((phi1+DeltaPhi)))*y2+(m*cos(theta))*z2)/(sqrt((x-0.02)^2 + y2^2+ z2^2))^5-(m*sin(theta)*sin((phi1+DeltaPhi)))/(sqrt((x-0.02)^2 + y2^2+ z2^2))^3) * 1e6 + Bey2;
    eq14 = mu0/(4*pi)*(3*z2*((m*sin(theta)*cos((phi1+DeltaPhi)))*(x-0.02)+(m*sin(theta)*sin((phi1+DeltaPhi)))*y2+(m*cos(theta))*z2)/(sqrt((x-0.02)^2 + y2^2+ z2^2))^5-(m*cos(theta))/(sqrt((x-0.02)^2 + y2^2+ z2^2))^3) * 1e6 + Bez2;
    eq15 = sqrt(Bex2^2+Bey2^2+Bez2^2);
    eq16 = Bex2*sin(-1*ang3)*sin(-1*ang4)+Bey2*cos(ang4)*sin(-1*ang3)+Bez2*cos(ang3);

    f = [eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8;eq9;eq10;eq11;eq12;eq13;eq14;eq15;eq16]-[B1x;B1y;B1z;B2x;B2y;B2z;Besum;Bd;A1x;A1y;A1z;A2x;A2y;A2z;Besum;Bd];
    F = sum(abs(f));
end

function ge = solve_g(b, gx, gy, gz)
    ge = [cos(b(1)) sin(b(1))*cos(b(2)) sin(b(1))*sin(b(2))] - [gz gy gx];
end










