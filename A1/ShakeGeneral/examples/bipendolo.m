%
%
%	(-1,0)		(1,0)
%
%    5	*		*  6
%	|		|
%	|		|
%	|		|
%	4		5
%	|		|
%	|		|
%    3	O------3--------O  4
%	|		|
%	|		|
%	1		2
%	|		|
%	|		|
%	|		|
%	O 1		O 2
%
%
%
%

addpath("../code");

NumeroParticelle = 4;


n = 3;
m = 2;

L = [sqrt(2),sqrt(2),2,1,1];
M = [1,1,1,1]';
g = 10;

V = [L(1),0,1,1,0;	 ... 	%1
     0,L(2),1,0,1;	 ...	%2
     1,1,L(3),1,1;	 ...	%3
     1,0,1,L(4),0;	 ...	%4
     0,1,1,0,L(5)];		%5

VP = [3,1; ...		%1
      4,2; ...		%2
      4,3; ...		%3
      5,3; ...		%4
      6,4];		%5
			

th0 = 0;
psi0 = 0;
phi0 = 5*pi/12;

q0 = [L(4)*sin(th0)-1+L(1)*sin(psi0);-L(4)*cos(th0)-L(1)*cos(psi0); ... 	%1
      L(5)*sin(th0)+1+L(2)*sin(phi0);-L(5)*cos(th0)-L(2)*cos(phi0); ... 	%2
      L(4)*sin(th0)-1;-L(4)*cos(th0); ... 					%3
      L(5)*sin(th0)+1;-L(5)*cos(th0); ... 					%4
      -1;0; ... 								%5
      1;0];	   								%6


%q0 = [-(L(3)/2)+L(1)*sin(psi0);-sqrt(1-(1-(L(3)/2))^2)-L(1)*cos(psi0); ...      %1
%      (L(3)/2)+L(2)*sin(phi0);-sqrt(1-(1-(L(3)/2))^2)-L(2)*cos(phi0); ...       %2
%      -(L(3)/2);-sqrt(L(4)^2-(1-(L(3)/2))^2); ...                               %3
%      (L(3)/2);-sqrt(L(5)^2-(1-(L(3)/2))^2); ...                                %4
%      -1;0; ...                                                                 %5
%      1;0];                                                                     %6


DimensioneQ = numel(q0);

v0 = zeros(DimensioneQ,1);

%no repulsione coulombiana
%U = @(qt) M(1)*g*qt(2) + M(2)*g*qt(4) + M(3)*g*qt(6) + M(4)*g*qt(8);
%Force = @(qt) [0;-g*M(1);0;-g*M(2);0;-g*M(3);0;-g*M(4)];

%repulsione coulombiana tra particelle 1 e 2
kQQ = 1;
U = @(qt) M(1)*g*qt(2) + M(2)*g*qt(4) + M(3)*g*qt(6) + M(4)*g*qt(8) + kQQ/norm(qt([1,2])-qt([3,4]));;
Force = @(qt) [0;-g*M(1);0;-g*M(2);0;-g*M(3);0;-g*M(4)];
Force = @(qt) Force(qt) + kQQ*[(qt([1,2])-qt([3,4]))/(norm(qt([1,2])-qt([3,4]))^3);(qt([3,4])-qt([1,2]))/(norm(qt([1,2])-qt([3,4]))^3);zeros(4,1)];

T = 50;
dt = 0.01;
err = 1e-10;


timer = tic();

[q,v] = shake(NumeroParticelle,n,m,M,V,VP,q0,v0,Force,T,dt,err);

toc(timer)



N = T/dt + 1;

TT = zeros(1,N);
UU = zeros(1,N);
MM = zeros(1,N);

for i = 1:1:N
  for k = 1:1:NumeroParticelle  
    TT(i) = TT(i) + 0.5*M(k)*norm(v([2*k-1,2*k],i))^2;    
  endfor

  UU(i) = UU(i) + U(q(:,i));
  MM(i) = TT(i) + UU(i);
endfor

figure
plot(0:dt:T, TT)
grid on
hold on
plot(0:dt:T, UU)
plot(0:dt:T, MM)
legend("Energia cinetica","Energia potenziale","Energia totale")

figure, plot(0:dt:T, MM)

figure, mostraStruttura(VP,q,NumeroParticelle,[-3,3,-3,1],dt,2)
