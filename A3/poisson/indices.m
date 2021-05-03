function [ix, iy, ii, hastop, hasbot, haslft, hasrgt] = ...
         indices (nx, ny)

  [ix, iy] = meshgrid (1:nx, 1:ny);
  ii = sub2ind (size(ix), iy, ix);
  ix = ix(:); iy = iy(:); ii = ii(:);

  if (nargout > 3)
    hastop = iy > 1;
    hasbot = iy < ny;

    hasrgt = ix < nx;
    haslft = ix > 1;
  endif

endfunction
