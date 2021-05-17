function [x,fk,gk,alphak]= backtrack(fun,xk,fk0,gk0,dk,varargin)
%BACKTRACK Metodo backtracking per line search.
%  [X,ALPHAK] = BACKTRACK(FUN,XK,FK,GK,DK) calcola
%  x_{k+1}=x_k+alpha_k d_k del metodo di discesa,
%  in cui alpha_k e' costruito con la tecnica di
%  backtracking, con sigma=1.e-4 e rho=1/4.
%  [X,ALPHAK] = BACKTRACK(FUN,XK,FK,GK,DK,SIGMA,RHO)
%  permette di precisare i valori dei parametri
%  sigma e rho. Tipicamente 1.e-4<sigma<0.1 e
%  1/10< rho <1/2. FUN e' un function handle
%  associato alla funzione obiettivo.
%  XK contiene l'elemento x_k della successione,
%  GK il gradiente di FUN in XK e DK la direzione d_k.
if nargin==5
  sigma=1.e-4; rho=1/4;
else
  sigma=varargin{1}; rho=varargin{2};
end
alphamin = 1.e-3; % valore minimo per il passo alpha
alphak = 1; 
k=0; x=xk+alphak*dk;
[fk, gk] = fun(x);
while fk > fk0+sigma*alphak*gk0'*dk && alphak>alphamin
    alphak = alphak*rho;
    x = xk+alphak*dk; k = k+1;
    [fk, gk] = fun(x);
end
