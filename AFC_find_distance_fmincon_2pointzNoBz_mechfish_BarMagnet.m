close all;
clear all;
delete all;
clc;
mu0=4*pi*1E-7;
Br=1; 
V=3.9270e-07;
m=Br*V/mu0;
load('testsetBarMagnet_MovingCap.mat');

% points1 = [2369 3903 5448 6923 8381 9624 11120];
% points2 = [2616 4203 5555 7045 8471 9890 11320];
% points1 =  [1647 3186 4749 6150 7602 9078 10580 12000];
% points2 =  [1616 3235 4689 6132 7762 9102 10560 11970];
% points1  = [2688 4749 6554 8220 10020 11580 13550 15740];% file testsetBarMagnet.txt
points1  = [1860 3463 9000 12830];
realdist = [0.15 0.2  0.25 0.3  0.35  0.4   0.45  0.5  ];

for k = 1 : length(points1)
    for i = 1 : 15 
        params1 = testset(points1(k) + i*50,:);
        params2 = testset(points1(k) + i*50 + 5,:);

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
        Besum = 53.2894;%57.26; %uT @ Lab


        %% Find angles mapping Earth magnetic field in the direction of gravity (Bz)
        
        %optimset('MaxIter',200000,'MaxFunEvals',2000000,'TolFun',1e-10,'TolX',1e-10);
        %optimoptions(@fsolve,'MaxFunctionEvaluations',500*numel(a0));
        opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt','Display','none');%,'PlotFcn',@optimplotfirstorderopt
        b0=[0.1 0.1];
        angles = fsolve(@(b)solve_g(b, g1x, g1y, g1z),b0,opts);
        ang1=angles(1);
        ang2=angles(2);
        angles = fsolve(@(b)solve_g(b, g2x, g2y, g2z),b0,opts);
        ang3=angles(1);
        ang4=angles(2);


        %% Find the distance
        %    1     2     3     4      5       6   7   8    9     10    11      12  13  14  15    
        a0=[0.3   0.1   0.1   2.5*pi  2.5*pi  1   1   1   0.1   0.1    pi/20   1   1   1   1 ];
        lb=[0.1 -0.15 -0.1   1*pi    1*pi   -60 -60 -60 -0.15 -0.1   -pi/4   -60 -60 -60 -60];
        ub=[0.6   0.15  0.1   3.2*pi  3.2*pi  60  60  60  0.15  0.1    pi/4    60  60  60  60];

%         if i > 1
%             a0 = Inet;
%             lb=[Inet(1)-0.1 Inet(2)-0.1  Inet(3)-0.1  Inet(4)-pi/6  Inet(5)-pi/4   Inet(6:8)-10  Inet(9:10)-0.1 -pi/6  Inet(12:end)-10];
%             ub=[Inet(1)+0.1 Inet(2)+0.1  Inet(3)+0.1  Inet(4)+pi/6  Inet(5)+pi/4   Inet(6:8)+10  Inet(9:10)+0.1 +pi/6  Inet(12:end)+10];
%             lb(5)= Inet(5)-pi/6;lb(11)= Inet(11)-pi/8;
%             ub(5)= Inet(5)+pi/6;ub(11)= Inet(11)+pi/8;
%         end

        opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt','Display',...
            'none','PlotFcn',@optimplotfirstorderopt,'MaxIterations',200000,...
            'MaxFunctionEvaluations',2000000,'FunctionTolerance',1e-10);%optimset('MaxIter',200000,'MaxFunEvals',2000000,'TolFun',1e-10,'TolX',1e-10);%optimoptions(@fsolve,'MaxFunctionEvaluations',500*numel(a0));
        for j = 1 : 10
            a0(1)=0.6*rand+0.1;
            Inet = fmincon(@(a)solve_B(a,B1x,B1y,B1z,B2x,B2y,B2z,Besum,ang1,ang2,A1x,A1y,A1z,A2x,A2y,A2z,ang3,ang4),a0,[],[],[],[],lb,ub);
            res(:,j) = Inet;
%             lb=[Inet(1)-0.1 Inet(2:3)-0.1  Inet(4)-pi/4  Inet(5)-pi/4   Inet(6:8)-10  Inet(9:10)-0.1 -pi/6  Inet(12:end)-10];
%             ub=[Inet(1)+0.3 Inet(2:3)+0.1  Inet(4)+pi/4  Inet(5)+pi/4   Inet(6:8)+10  Inet(9:10)+0.1 +pi/6  Inet(12:end)+10];
        end

        results(i,:) = median(res')';
        try 
            bbb = results(i,:) - results(i-1,:);
        end
        dist(i) = Inet(1);
        thetas(i) = Inet(4);
        phis(i) = Inet(5);

    end
    figure;plot(results(:,1:5),'*-');hold on
    plot(results(:,11),'*-')
    legend;
    meandist(k)=median(dist);
    meanthetas(k) = (mean(thetas) - 2*pi)/pi;
    meanphis(k) = (mean(phis) - pi)/pi;
end

figure;plot (meandist,'o--'); hold on
plot (meanthetas,'+--');
plot(meanphis,'*--')
figure;plot(meandist,'o-');hold on
plot(realdist,'*--')
legend('Estimated distance','Distance')
xlabel('#No')
ylabel('Distance (m)')

%% Functions
function F = solve_B(a,B1x,B1y,B1z,B2x,B2y,B2z,Besum,ang1,ang2,A1x,A1y,A1z,A2x,A2y,A2z,ang3,ang4)
    mu0=4*pi*1E-7;
    Br=1.2; 
    V=pi*(1.5e-3)^2*19e-3;
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

function ge = solve_g(b, gx, gy, gz)
    ge = [cos(b(1)) sin(b(1))*cos(b(2)) sin(b(1))*sin(b(2))] - [gz gy gx];
end










