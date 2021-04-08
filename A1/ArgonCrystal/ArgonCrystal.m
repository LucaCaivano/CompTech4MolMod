clear
clc
close all
%% Parameters
N = 7;
m = 66.34e-27; % [kg]
M = m * ones(1, N);
kB = 1.380658e-23; % [J/K]
epsilon = 119.8 * kB; % [J]
Epsilon = epsilon * ones(N, N);
sigma = 0.341; % [nm]
Sigma = sigma * ones(N, N);
r = 0.4;
theta = pi/3;

T0 = 0; % [ns]
Tf = 0.2; % [ns]
dt = 1e-6; % [ns]

%% Initial conditions
N_step = (Tf - T0) / dt;

% q(1,:,1) = [0.0 0.0];
% q(2,:,1) = [.02 .39];
% q(3,:,1) = [.34 .17];
% q(4,:,1) = [.36 -.21];
% q(5,:,1) = [-.02 -.40];
% q(6,:,1) = [-.35 -.16];
% q(7,:,1) = [-.31 .21];
% 
% p(1,:,1) = M(1)*[-30 -20];
% p(2,:,1) = M(2)*[50 -90];
% p(3,:,1) = M(3)*[-70 -60];
% p(4,:,1) = M(4)*[90 40];
% p(5,:,1) = M(5)*[80 90];
% p(6,:,1) = M(6)*[-40 100];
% p(7,:,1) = M(7)*[-80 -60];

q(1,:,1) = [0.0 0.0];
q(2,:,1) = [r*cos(theta) r*sin(theta)];
q(3,:,1) = [r*cos(2*theta) r*sin(2*theta)];
q(4,:,1) = [r*cos(3*theta) r*sin(3*theta)];
q(5,:,1) = [r*cos(4*theta) r*sin(4*theta)];
q(6,:,1) = [r*cos(5*theta) r*sin(5*theta)];
q(7,:,1) = [r*cos(6*theta) r*sin(6*theta)];

p(1,:,1) = M(1)*[0 0];
p(2,:,1) = M(2)*[0 0];
p(3,:,1) = M(3)*[0 0];
p(4,:,1) = M(4)*[0 0];
p(5,:,1) = M(5)*[0 0];
p(6,:,1) = M(6)*[0 0];
p(7,:,1) = M(7)*[0 0];

qq = zeros(N,2);
for jj = 1:N
    qq(jj,:) = q(jj,:,1) + 0.5 * dt * p(jj,:,1) / M(jj);
    pp = p(jj,:,1) + 0.5 * dt * dT_pot_dqk(q(:,:,1),Sigma,Epsilon,jj);
    q(jj,:,2) = q(jj,:,1) + dt * pp / M(jj); 
end

% for jj = 1:N 
%     p(jj,:,2) = p(jj,:,1) + dt * dT_pot_dqk(qq,Sigma,Epsilon,jj);
% end

KK(1) = K_cin(p(:,:,1),M);
TT(1) = T_pot(q(:,:,1),Sigma,Epsilon);
HH(1) = KK(1) + TT(1);
% 
% KK(2) = K_cin(p(:,:,2),M);
% TT(2) = T_pot(q(:,:,2),Sigma,Epsilon);
% HH(2) = KK(2) + TT(2);

%% Solution
% POSITION VERLET (multistep)
for ii = 2:N_step
    for jj = 1:N
        q(jj,:,ii+1) = 2*q(jj,:,ii) - q(jj,:,ii-1) - dt^2 * dT_pot_dqk(q(:,:,ii),Sigma,Epsilon,jj) / M(jj);
        p(jj,:,ii) = 0.5 * M(jj) * (q(jj,:,ii+1) - q(jj,:,ii-1)) / dt;
    end
    KK(ii) = K_cin(p(:,:,ii),M);
    TT(ii) = T_pot(q(:,:,ii),Sigma,Epsilon);
    HH(ii) = KK(ii) + TT(ii);
end

%% Plotting
figure('units', 'pixel', 'outerposition', [0, 0, 1000, 500], 'Name', 'Argon crystal - PV')
for ii = 1:N_step
    if mod(ii,1000)==0    
        % Crystal
        subplot(1,2,1);
            for jj = 1:N
                plot(q(jj,1,ii), q(jj,2,ii),'-o');
                hold on
            end
        title(['t=',num2str(ii*dt)])
        xlim([min(q,[],"all")-1 max(q,[],"all")+1]);
        ylim([min(q,[],"all")-1 max(q,[],"all")+1]);
        hold off
        
        % Energy
        subplot(1,2,2)
        plot(KK(1:ii));
        hold on
        plot(TT(1:ii));
        plot(HH(1:ii));
        title(['Energy = ',num2str(HH(1))])
        legend(["K","T","H"]);
        hold off        
        
        pause(2e-4);
    end
end