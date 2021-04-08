
timer = tic();

L1 = 1;
L2 = 1;
m1 = 1;
m2 = 1;
g = 10;

th0 = 0;
phi0 = 5*pi/12;
q0 = [L1*sin(th0),-L1*cos(th0), L1*sin(th0)+L2*sin(phi0), -L1*cos(th0)-L2*cos(phi0)]';
v0 = [0,0,0,0]';

F1 = m1*[0,-g]';
F2 = m2*[0,-g]';

T = 10;
dt = 0.01;
err = 1e-10;

N = T/dt + 1;
q = zeros(4,N);
v = zeros(4,N);
q(:,1) = q0;
v(:,1) = v0;
TT = zeros(1,N);
UU = zeros(1,N);
MM = zeros(1,N);


q(:,2) = q0 + v0*dt;


for i = 3:1:N

  q1t = q([1,2],i-1);
  q2t = q([3,4],i-1);
  q1tmh = q([1,2],i-2);
  q2tmh = q([3,4],i-2);

  a1 = F1/m1;
  a2 = F2/m2;
  
  q1tilda = 2*q1t - q1tmh + a1*(dt^2);
  q2tilda = 2*q2t - q2tmh + a2*(dt^2);
  
  A = zeros(2,2);

  A(1,1) = (q1tilda'*q1t)/m1;
  A(1,2) = (q1tilda'*(q1t-q2t))/m1;
  A(2,1) = ((q1tilda-q2tilda)'*q1t)/m1;
  A(2,2) = ((q1tilda-q2tilda)'*(q1t-q2t))*(1/m1 + 1/m2);
  A = 4*(dt^2)*A;

  b = [L1^2 - norm(q1tilda)^2, L2^2 - norm(q1tilda- q2tilda)^2]';

  sol = A\b;
  
  
  F = @(gamma) [4*(dt^4)*((norm(q1t)^2)*(gamma(1)^2) + 2*(norm(q1t)^2 - q1t'*q2t)*(gamma(1)*gamma(2)) + (norm(q1t-q2t)^2)*(gamma(2)^2))/(m1^2) + ...
		4*(dt^2)*(( (q1tilda'*q1t)/m1 )*gamma(1) + ( (q1tilda'*(q1t - q2t))/m1 )*gamma(2)) + ...
		norm(q1tilda)^2 - L1^2; ...
		4*(dt^4)*((1/(m1^2))*(norm(q1t)^2)*(gamma(1)^2) + 2*((1/m1 + 1/m2)/m1)*(norm(q1t)^2 - q1t'*q2t)*(gamma(1)*gamma(2)) + ((1/m1 + 1/m2)^2)*(norm(q1t - q2t)^2)*(gamma(2)^2)) + ...
		4*(dt^2)*((1/m1)*( (q1tilda - q2tilda)'*q1t )*gamma(1) + (1/m1 + 1/m2)*( (q1tilda - q2tilda)'*(q1t - q2t) )*gamma(2)) + ...
		norm(q2tilda - q1tilda)^2 - L2^2];

  [sol, fval, info] = fsolve(F, sol, optimset("TolFun", err));
  
  gamma1 = sol(1);
  gamma2 = sol(2);
  q([1,2],i) = q1tilda + 2*(gamma1*q1t + gamma2*(q1t-q2t))*(dt^2)/m1;
  q([3,4],i) = q2tilda + 2*(gamma2*(q2t-q1t))*(dt^2)/m2;

  v(:,i-1) = (q(:,i) - q(:,i-2))/(2*dt);  
  
endfor

v(:,N) = v(:,N-1); %brutto, da aggiustare

toc(timer)

for i = 1:1:N


  TT(i) = 0.5 * m1 * norm(v([1,2],i))^2 + 0.5 * m2 * norm(v([3,4],i))^2;
  UU(i) = m1*g*(q(2,i)) + m2*g*(q(4,i));
  MM(i) = TT(i) + UU(i);
  
endfor

figure
subplot(2,1,2);
plot(0:dt:T, TT)
grid on
hold on
plot(0:dt:T, UU)
plot(0:dt:T, MM)
legend("Energia cinetica","Energia potenziale","Energia totale")

subplot(2,2,1)
plot(0:dt:T, sqrt((q(1,:)').^2 + (q(2,:)').^2) )
hold on
grid on
plot(0:dt:T, sqrt(((q(3,:) - q(1,:))').^2 +  ((q(4,:) - q(2,:))').^2))

max(abs(sqrt((q(1,:)').^2 + (q(2,:)').^2) - L1))
max(abs(sqrt(((q(3,:) - q(1,:))').^2 +  ((q(4,:) - q(2,:))').^2) - L2))

subplot(2,2,2)
mostraDoppioPendolo(q(1,:),q(2,:),q(3,:),q(4,:),[-2,2,-2,2],dt,1)



figure, plot(0:dt:T, MM)
