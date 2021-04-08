%
%
%	(-l,0)		(l,0)
%
%    5	*		*  6
%	|		|
%	|		|
%	|		|
%	4		5
%	|		|
%	|		|
%    3	O		O  4
%	|		|
%	|		|
%	2		3
%	|		|
%	|		|
%	|		|
%     1 O--------1------O 2
%
%
%



addpath("../code");

NumeroParticelle = 4;


n = 3;
m = 2;

l = 1;

L = [2*l,l,l,l,l];
M = [1,1,1,1]';
g = 10;

V = [L(1),1,1,0,0; ...	%1
     1,L(2),0,1,0; ...	%2
     1,0,L(3),0,1; ...	%3
     0,1,0,L(4),0; ...	%4
     0,0,1,0,L(5)];	%5

VP = [2,1; ...	%1
      3,1; ...	%2
      4,2; ...	%3
      5,3; ...	%4
      6,4];	%5		

th0 = pi/6;

q0 = [-l;-2*l*cos(th0); ...	
      l;-2*l*cos(th0); ...	
      -l-l*sin(th0);-l*cos(th0); ...	
      l+l*sin(th0);-l*cos(th0); ...	
      -l;0; ...	
      l;0];	

DimensioneQ = numel(q0);

v0 = zeros(DimensioneQ,1);

U = @(qt) M(1)*g*qt(2) + M(2)*g*qt(4) + M(3)*g*qt(6) + M(4)*g*qt(8);
Force = @(qt) [0;-g*M(1);0;-g*M(2);0;-g*M(3);0;-g*M(4)];


T = 20;
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

figure('position',[1800,0,1000,1000])
pause(1)
%figure
mostraStruttura(VP,q,NumeroParticelle,[-2,2,-3,1],dt,1)


