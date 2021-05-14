function [x,alphak]= backtrack(fun,xk,gk,dk,varargin)
%BACKTRACK Metodo backtracking per line search.
%  [X,ALPHAK] = BACKTRACK(FUN,XK,GK,DK) calcola
%  x_{k+1}=x_k+alpha_k d_k del metodo di discesa,
%  in cui alpha_k e' costruito con la tecnica di
%  backtracking, con sigma=1.e-4 e rho=1/4.
%  [X,ALPHAK] = BACKTRACK(FUN,XK,GK,DK,SIGMA,RHO)
%  permette di precisare i valori dei parametri
%  sigma e rho. Tipicamente 1.e-4<sigma<0.1 e
%  1/10< rho <1/2. FUN e' un function handle
%  associato alla funzione obiettivo.
%  XK contiene l'elemento x_k della successione,
%  GK il gradiente di FUN in XK e DK la direzione d_k.
if nargin==4
    sigma=1.e-4; rho=1/4;
else
    sigma=varargin{1}; rho=varargin{2};
end
alphamin=1.e-1; % valore minimo per il passo alpha
alphak = 1; fk = fun(xk);
k=0; x=xk+alphak*dk;
while fun(x)>fk+sigma*alphak*gk'*dk && alphak>alphamin
    alphak = alphak*rho;
    x = xk+alphak*dk; k = k+1;
end
