function s = Urrextrapolation(X)
  
%  s = Urrextrapolation(X)
%  RRE vector extrapolation see 
%  Smith, Ford & Sidi SIREV 29 II 06/1987
  
  if (Ucolumns(X)>Urows(X))
    X=X';
  end
  
  % compute first and second variations
  U = X(:,2:end) - X(:,1:end-1);
  V = U(:,2:end) - U(:,1:end-1);
  
  % eliminate unused u_k column
  U(:,end) = [];
  
  s = X(:,1) - U * pinv(V) * U(:,1);
  
