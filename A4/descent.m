function [x,err,iter,trace]= descent (fun, x0, tol, kmax, ...
                                meth, varargin)
%DESCENT Metodo di discesa per il calcolo di minimi
%  [X,ERR,ITER]=DESCENT(FUN,GRAD,X0,TOL,KMAX,METH,HESS)
%  approssima un punto di minimo della funzione FUN
%  mediante il metodo di discesa con direzioni di
%  Newton (METH=1), BFGS (METH=2), del gradiente
%  (METH=3) o del gradiente coniugato con
%  beta_k di Fletcher and Reeves (METH=41),
%  beta_k di Polak and Ribiere (METH=42),
%  beta_k di Hestenes and Stiefel (METH=43).
%  Il passo e' costruito con la tecnica di back-
%  tracking. FUN, GRAD ed HESS (quest'ultima usata
%  solo se METH=1) sono function handle
%  associati alla funzione obiettivo, al suo gradiente
%  ed alla matrice Hessiana. Se METH=2, HESS e' una
%  matrice approssimante l'Hessiana nel punto iniziale
%  X0 della successione. TOL e' la tolleranza per il
%  test d'arresto e KMAX e' il numero massimo di
%  iterazioni. Si richiama le function backtrack.m
if nargin>=6
if meth==1, hess=varargin{1};
elseif meth==2, H=varargin{1}; end
end
err=tol+1; k=0; xk=x0(:); [~,gk]=fun(xk); dk=-gk;
eps2=sqrt(eps);
trace.param (:, 1) = xk;
trace.err (:, 1) = err;
while err>tol && k< kmax
  if meth==1;        H=hess(xk); dk=-H\gk; % Newton
  elseif meth==2     dk=-H\gk;             % BFGS
  elseif meth==3     dk=-gk;               % gradient
  end
  [xk1,alphak]= backtrack(fun,xk,gk,dk);
  [res,gk1]=fun(xk1);
  if meth==2 % BFGS update
    yk=gk1-gk; sk=xk1-xk; yks=yk'*sk;
    if yks> eps2*norm(sk)*norm(yk)
      Hs=H*sk;
      H=H+(yk*yk')/yks-(Hs*Hs')/(sk'*Hs);
    end
  elseif meth>=40 % CG upgrade
    if meth == 41
      betak=(gk1'*gk1)/(gk'*gk); % FR
    elseif meth == 42
      betak=(gk1'*(gk1-gk))/(gk'*gk); % PR
    elseif meth == 43
      betak=(gk1'*(gk1-gk))/(dk'*(gk1-gk)); % HS
    end
    dk=-gk1+betak*dk;
  end
  xk=xk1; gk=gk1; k=k+1; xkt=xk1;

  for i=1:length(xk1); xkt(i)=max([abs(xk1(i)),1]); end
  err=norm((gk1.*xkt)/max([abs(res),1]),inf);

  trace.param (:, k+1) = xk;
  trace.err (:, k+1) = err;
  printf("Errore = %g\n", err);
  printf("xk = %g\n", xk);

end

x=xk; iter=k;
if (k==kmax && err > tol)
  fprintf(['descent si e'' arrestato senza aver ',...
   'soddisfatto l''accuratezza richiesta, avendo\n',...
   'raggiunto il massimo numero di iterazioni\n']);
end
