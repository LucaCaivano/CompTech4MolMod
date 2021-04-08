%
%
%	(-L5/2,0)	(L5/2,0)
%
%    8	*		*  7
%	|		|
%	|		|
%	|		|
%	6		7
%	|		|
%	|		|
%    5	O------5--------O  6
%	|		|
%	|		|
%	3		4
%	|		|
%	|		|
%	|		|
%	O 4		O 3
%	|		|
%	|		|
%	1		2
%	|		|
%	|		|
%	|		|
%	O 1		O 2
%
%

addpath("../code");

NumeroParticelle = 6;


n = 5;
m = 2;

L = [1,1,1,1,2,1,1];
M = [1,1,1,1,1,1]';
g = 10;

V = [L(1),0,1,0,0,0,0; ... 	%1
     0,L(2),0,1,0,0,0; ...	%2
     1,0,L(3),0,1,1,0; ...	%3
     0,1,0,L(4),1,0,1; ...	%4
     0,0,1,1,L(5),1,1; ...	%5
     0,0,1,0,1,L(6),0; ...	%6
     0,0,0,1,1,0,L(7)];	   	%7

VP = [4,1; ...		%1
      3,2; ...		%2
      5,4; ...		%3
      6,3; ...		%4
      6,5; ...		%5
      8,5; ...		%6
      7,6];		%7
			

th0 = 0;
phi0 = pi/3;
psi0 = -pi/3;

q0 = [L(6)*sin(th0)-L(5)/2+L(3)*sin(phi0);-L(6)*cos(th0)-L(3)*cos(phi0)-L(1); ... 	%1
      L(7)*sin(th0)+L(5)/2+L(4)*sin(psi0);-L(7)*cos(th0)-L(4)*cos(psi0)-L(2); ... 	%2
      L(7)*sin(th0)+L(5)/2+L(4)*sin(psi0);-L(7)*cos(th0)-L(4)*cos(psi0); ... 	%3
      L(6)*sin(th0)-L(5)/2+L(3)*sin(phi0);-L(6)*cos(th0)-L(3)*cos(phi0); ... 	%4
      L(6)*sin(th0)-L(5)/2;-L(6)*cos(th0); ... 					%5
      L(7)*sin(th0)+L(5)/2;-L(7)*cos(th0); ... 					%6
      L(5)/2;0; ... 									%7
      -L(5)/2;0];	   								%8

DimensioneQ = numel(q0);

v0 = zeros(DimensioneQ,1);

U = @(qt) M(1)*g*qt(2) + M(2)*g*qt(4) + M(3)*g*qt(6) + M(4)*g*qt(8) + M(5)*g*qt(10) + M(6)*g*qt(12);
Force = @(qt) [0;-g*M(1);0;-g*M(2);0;-g*M(3);0;-g*M(4);0;-g*M(5);0;-g*M(6)];

T = 10;
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

figure, mostraStruttura(VP,q,NumeroParticelle,[-3,3,-3,1],dt,1)


