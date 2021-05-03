function A = structured_laplacian (msh, D)

  N  = msh.nx*msh.ny;

  [ix, iy, ii, hastop, hasbot, haslft, hasrgt] = ...
  indices (msh.nx, msh.ny);

  fltop = flbot = flrgt = fllft = zeros (size (ii));

  hx_over_hy = (msh.hx(:) ./ msh.hy(:).').'(:);
  fltop(hastop & hasrgt) += D .* hx_over_hy / 2;
  fltop(hastop & haslft) += D .* hx_over_hy / 2;

  flbot(hasbot & hasrgt) += D .* hx_over_hy / 2;
  flbot(hasbot & haslft) += D .* hx_over_hy / 2;

  flrgt(hasrgt & hastop) += D ./ hx_over_hy / 2;
  flrgt(hasrgt & hasbot) += D ./ hx_over_hy / 2;

  fllft(haslft & hastop) += D ./ hx_over_hy / 2;
  fllft(haslft & hasbot) += D ./ hx_over_hy / 2;

  A = flux_assembly (N, ii,
                     msh.nx, msh.ny,
                     hastop, hasbot,
                     haslft, hasrgt,
                     fltop, flbot,
                     fllft, flrgt);

endfunction
