function rhs = structured_rhs (msh, ndrhs)

  [ix, iy, ii, hastop, hasbot, haslft, hasrgt] = ...
  indices (msh.nx, msh.ny);

  rhs = zeros (size (ii));
  hxhy = (msh.hx(:) .* msh.hy(:).').'(:);
  rhs(hastop & hasrgt) += 1*(hxhy / 4);
  rhs(hastop & haslft) += 1*(hxhy / 4);
  rhs(hasbot & hasrgt) += 1*(hxhy / 4);
  rhs(hasbot & haslft) += 1*(hxhy / 4);

  rhs .*= ndrhs;

endfunction
