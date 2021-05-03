function [phi, L, U, P, Q, R] = poisson (msh, epsilon, rho, L, U, P, Q, R)

  tnodes = boundary_nodes (msh.nx, msh.ny, 1);
  bnodes = boundary_nodes (msh.nx, msh.ny, 2);
  lnodes = boundary_nodes (msh.nx, msh.ny, 3);
  rnodes = boundary_nodes (msh.nx, msh.ny, 4);

  inodes = setdiff (1:msh.nx*msh.ny, tnodes);
  inodes = setdiff (inodes, rnodes);

  if (nargin == 3)
    A = structured_laplacian (msh, epsilon);
    %% periodic BC
    A(:, bnodes) += +A(:, tnodes);
    A(bnodes, :) += +A(tnodes, :);
    A(:, lnodes) += +A(:, rnodes);
    A(lnodes, :) += +A(rnodes, :);

    M =  A(inodes(2:end),inodes(2:end));
  endif

  f = structured_rhs (msh, rho(:));
  f(bnodes)    += f(tnodes);
  f(lnodes)    += f(rnodes);

  %% solve setting phi=0 at one point
  phi = zeros (msh.nx*msh.ny, 1);
  if (nargin > 3)
    phi(inodes(2:end)) = Q * (U \ (L \ (P * (R \ f(inodes(2:end))))));
  else
    if (nargout > 1)
      [L, U, P, Q, R] = lu (M);
      phi(inodes(2:end)) = Q * (U \ (L \ (P * (R \ f(inodes(2:end))))));
    else
      phi(inodes(2:end)) = M \ f(inodes(2:end));
    endif
  endif

  phi(tnodes) = phi(bnodes);
  phi(rnodes) = phi(lnodes);

  %% normalize so that integral mean of solution is 0
  phi -= sum (structured_rhs (msh, phi(:))) / ((max(msh.x)-min(msh.x))*(max(msh.y)-min(msh.y)));
endfunction
