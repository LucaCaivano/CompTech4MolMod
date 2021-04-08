clear
clc
close all
%% Parameters
m1 = 1;
m2 = 0.2;
k1 = 100;
k2 = 10;
l1 = 1;
l2 = 2;
g = [0  -9.81];

T0 = 0;
Tf = 4;
dt = 1e-3;

theta = -pi/6;
phi = -pi/12;

scheme = "PV";

%% Initial conditions
K = @(p1, p2) 0.5*(norm(p1)^2 / m1 + norm(p2)^2 / m2);
T = @(q1, q2) -m1*g*q1' - m2*g*q2' + 0.5*k1*(norm(q1) - l1)^2 + 0.5*k2*(norm(q1-q2) - l2)^2;

dK1 = @(p1, p2) p1 ./ m1;
dK2 = @(p1, p2) p2 ./ m2;

dT1 = @(q1, q2) - m1*g + k1*(norm(q1) - l1).*q1 / norm(q1) + k2*(norm(q1-q2) - l2)*(q1-q2) / norm(q1-q2);
dT2 = @(q1, q2) - m2*g - k2*(norm(q1-q2) - l2) *(q1-q2) / norm(q1-q2);

N = (Tf - T0) / dt;

q1(1,:) = [l1*cos(theta)  l1*sin(theta)];
q2(1,:) = [q1(1) + l2*cos(phi)  q1(2) + l2*sin(phi)];

p1(1,:) = [0  0];
p2(1,:) = [0  0];

% PV
if (scheme == "PV")
    qq1 = q1(1,:) + 0.5 * dt * p1(1,:) / m1;
    qq2 = q2(1,:) + 0.5 * dt * p2(1,:) / m2;
    
    pp1 = p1(1,:) + 0.5 * dt * dT1(q1(1,:), q2(1,:));
    pp2 = p2(1,:) + 0.5 * dt * dT2(q1(1,:), q2(1,:));
    
    q1(2,:) = q1(1,:) + dt * pp1 / m1;
    q2(2,:) = q2(1,:) + dt * pp2 / m2;
end

KK(1) = K(p1(1,:), p2(1,:));
TT(1) = T(q1(1,:), q2(1,:));
HH(1) = KK(1) + TT(1);

%% Solution
% PV
if (scheme == "PV")
    disp ("Solution with PV");
    for ii = 2:N    
        q1(ii + 1,:) = 2*q1(ii,:) - q1(ii - 1,:) - dt^2 * dT1(q1(ii,:), q2(ii,:)) / m1;
        q2(ii + 1,:) = 2*q2(ii,:) - q2(ii - 1,:) - dt^2 * dT2(q1(ii,:), q2(ii,:)) / m2;

        p1(ii,:) = 0.5 * m1 * (q1(ii + 1,:) - q1(ii - 1, :)) / dt;
        p2(ii,:) = 0.5 * m2 * (q2(ii + 1,:) - q2(ii - 1, :)) / dt;

        KK(ii) = K(p1(ii,:), p2(ii,:));
        TT(ii) = T(q1(ii,:), q2(ii,:));
        HH(ii) = KK(ii) + TT(ii);
    end
end

% VV
if (scheme == "VV")
    disp ("Solution with VV");
    for ii = 1:N
    pp1 = p1(ii,:) - dt * dT1(q1(ii,:), q2(ii,:)) / 2;
    pp2 = p2(ii,:) - dt * dT2(q1(ii,:), q2(ii,:)) / 2;
    
    q1(ii + 1,:) = q1(ii,:) + dt * pp1 / m1;
    q2(ii + 1,:) = q2(ii,:) + dt * pp2 / m2;
    
    p1(ii + 1,:) = pp1 - dt * dT1(q1(ii + 1,:), q2(ii + 1,:)) / 2;
    p2(ii + 1,:) = pp2 - dt * dT2(q1(ii + 1,:), q2(ii + 1,:)) / 2;
    
    KK(ii + 1) = K(p1(ii + 1,:), p2(ii + 1,:));
    TT(ii + 1) = T(q1(ii + 1,:), q2(ii + 1,:));
    HH(ii + 1) = KK(ii + 1) + TT(ii + 1);
    end
end
    
%% Plotting
figure('units', 'pixel', 'outerposition', [0, 0, 600, 600], 'Name', 'Double pendulum - ' + scheme)
for ii = 1:N
    if mod(ii,100)==0
        x1 = q1(ii,1);
        y1 = q1(ii,2);
        x2 = q2(ii,1);
        y2 = q2(ii,2);
        
        % Pendulum
        subplot(2,3,1);
        plot(x1,y1,'ob');
        hold on
        plot(x2,y2,'or');
        title(['t=',num2str(ii*dt)])
        xlim([min(min(q1(:,1),q2(:,1)))-1 max(max(q1(:,1),q2(:,1)))+1]);
        ylim([min(min(q1(:,2),q2(:,2)))-1 max(max(q1(:,2),q2(:,2)))+1]);
        quiver(0,0,x1,y1,1,'b');
        quiver(x1,y1,x2-x1,y2-y1,1,'r');
        hold off
        
        % Energy
        subplot(2,3,4)
        start = 1 + (scheme == "PV");
        plot(KK(start:ii));
        hold on
        plot(TT(start:ii));
        plot(HH(start:ii));
        title(['Energy = ',num2str(HH(1))])
        xlim([0 length(HH)]);
        ylim([min(min(KK),min(TT))-1 max(max(KK),max(TT))+1]);
        legend(["K","T","H"]);
        hold off
        
        % Trajectories 1
        subplot(2,3,2)
        plot(q1(ii,1),q1(ii,2),'.k');
        hold on
        title(['Trajectories 1'])
        xlim([min(q1(:,1))-1 max(q1(:,1))+1]);
        ylim([min(q1(:,2))-1 max(q1(:,2))+1]);
        
        % Trajectories 2
        subplot(2,3,5)
        plot(q2(ii,1),q2(ii,2),'.k');
        hold on
        title(['Trajectories 2'])
        xlim([min(q2(:,1))-1 max(q2(:,1))+1]);
        ylim([min(q2(:,2))-1 max(q2(:,2))+1]);
        
%         % Trajectories 1
%         subplot(2,3,3)
%         plot(norm(q1(ii,:)),norm(p1(ii,:)),'.b');
%         hold on
% 
%         % Trajectories 1
%         subplot(2,3,6)
%         plot(norm(q2(ii,:)),norm(p2(ii,:)),'.b');
%         hold on
        
        pause(2e-3);
    end
end