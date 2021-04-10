function grd_to_ptcl = init_ptcl_mesh (grd, ptcls)
  ridx = lookup (grd.y, ptcls.x(2, :));
  cidx = lookup (grd.x, ptcls.x(1, :));

  ptcl_to_grd = sub2ind ([grd.ncy, grd.ncx], ridx, cidx);
  grd_to_ptcl = cell (grd.ncy,grd.ncx);
  for ii = 1 : numel (ptcl_to_grd)
    grd_to_ptcl{ptcl_to_grd(ii)}(end+1) = ii;
  end
end
