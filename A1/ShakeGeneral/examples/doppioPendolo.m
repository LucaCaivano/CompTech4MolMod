%
%
%	(0,0)
%
%    3	*
%	|
%	|
%	|
%	2
%	|
%	|
%    2	O
%	|
%	|
%	1
%	|
%	|
%	|
%	O 1
%
%
%
%

addpath("../code");

NumeroParticelle = 2;


n = 1;
m = 1;

L = [1,1];
M = [1,1]';
g = 10;

V = [L(1),1; ...	%1
     1,L(2)];		%2

VP = [2,1;		%1
      3,2];		%2

th0 = 0;
phi0 = 5*pi/12;

q0 = [L(2)*sin(th0) + L(1)*sin(phi0); -L(2)*cos(th0) - L(1)*cos(phi0); ...	%1
      L(2)*sin(th0); -L(2)*cos(th0); ...					%2
      0; 0];									%3

DimensioneQ = numel(q0);

v0 = zeros(DimensioneQ,1);

U = @(qt) M(1)*g*qt(2) + M(2)*g*qt(4);
Force = @(qt) [0;-g*M(1);0;-g*M(2)];

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

figure, mostraStruttura(VP,q,NumeroParticelle,[-1,1,-2,0],dt,1)
