/*

Copyright (C) 2021 Carlo de Falco

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Author: Carlo de Falco <carlo@guglielmo.local>
Created: 2021-04-19

*/

#include <octave/oct.h>
#include <iostream>

DEFUN_DLD(phigradient, args, nargout,
          "-*- texinfo -*-\n\
@deftypefn {} {[@var{dphidpx}, @var{dphidpy}] =} \
phigradient (@var{ptcls}, @var{grd_to_ptcl}, @var{msh}, @var{phi})\n\
Compute the gradient of the field @var{phi} at the locations of the particles in @var{ptcls}.\n\
@var{grd_to_ptcl} is a cell array that defines which particles belong to a given grid cell.\n\
The struct @var{msh} describes the cartesian grid.\n\
@seealso{}\n\
@end deftypefn")
{
  octave_value_list retval;
  int nargin = args.length ();

  octave_scalar_map pctls = args(0).scalar_map_value ();
  Matrix px = pctls.contents ("x").matrix_value ();
  octave_idx_type nptcls = px.cols ();

  Cell grd_to_ptcl = args(1).cell_value ();
  octave_scalar_map msh = args(2).scalar_map_value ();
  octave_idx_type ncx = msh.contents ("ncx").idx_type_value ();
  octave_idx_type ncy = msh.contents ("ncy").idx_type_value ();
  octave_idx_type nx = msh.contents ("nx").idx_type_value ();
  octave_idx_type ny = msh.contents ("ny").idx_type_value ();

  ColumnVector phi = args(3).column_vector_value ();

  RowVector dphidpx(nptcls, 0.0);
  RowVector dphidpy(nptcls, 0.0);

  for (octave_idx_type icell = 0; icell < ncy; ++icell)
    for (octave_idx_type jcell = 0; jcell < ncx; ++jcell)
    {
      Array<octave_idx_type> plist = grd_to_ptcl(icell, jcell).octave_idx_type_vector_value ();

      if (! plist.isempty ())
        {
          const double x0 = msh.contents ("x").column_vector_value () (jcell);
          const double x1 = msh.contents ("x").column_vector_value () (jcell+1);
          const double hx = x1-x0;
          const double y0 = msh.contents ("y").column_vector_value () (icell);
          const double y1 = msh.contents ("y").column_vector_value () (icell+1);
          const double hy = y1-y0;

          const double phi00 = phi(jcell*ny + icell);
          const double phi10 = phi((jcell+1)*ny + icell);
          const double phi01 = phi(jcell*ny + icell+1);
          const double phi11 = phi((jcell+1)*ny + icell+1);

          for (octave_idx_type iptc = 0; iptc < plist.numel (); ++iptc)
            {
              const octave_idx_type nptc = plist(iptc)-1;
              const double pclx = px(0,nptc);
              const double pcly = px(1,nptc);


              dphidpx(nptc) = (phi00 * (-1) * (y1 - pcly) -
                               phi01 * (-1) * (y0 - pcly) -
                               phi10 * (-1) * (y1 - pcly) +
                               phi11 * (-1) * (y0 - pcly))/(hy*hx);

              dphidpy(nptc) = (phi00 * (x1 - pclx) * (-1.) -
                               phi01 * (x1 - pclx) * (-1.) -
                               phi10 * (x0 - pclx) * (-1.) +
                               phi11 * (x0 - pclx) * (-1.))/(hy*hx);

            }
        }
    }

  retval(1) = dphidpy;
  retval(0) = dphidpx;

  return retval;
}
