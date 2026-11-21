#!/usr/bin/env python3
"""Applies to a copy of the nec2c 1.3.1 sources the fixes that nec2-rs takes
from the original NEC-2 Fortran (see docs/nec2c-deviations.md).

Usage: patch-nec2c.py <directory with the nec2c sources>

The result is the "fixed" reference that the deviating test decks are
compared with. Whitespace in the patterns is matched loosely because the
nec2c sources mix tabs and spaces.
"""
import pathlib
import re
import sys

root = pathlib.Path(sys.argv[1])


def patch(name, pairs):
    path = root / name
    # newline="" keeps the line endings: somnec.c has CRLF.
    with open(path, newline="") as source:
        text = source.read()
    for old, new in pairs:
        pattern = r"\s+".join(re.escape(part) for part in old.split())
        found = re.findall(pattern, text)
        if len(found) != 1:
            sys.exit(f"{name}: pattern found {len(found)} times: {old[:70]!r}")
        text = re.sub(pattern, lambda _: new.strip("\n"), text)
    with open(path, "w", newline="") as target:
        target.write(text)


# somnec.c: the contour selection of evlua and the convergence test of gshank
# follow the Fortran SOMNEC; the interpolation is told about a new grid.
patch("somnec.c", [
    ("""
  jump = FALSE;
  if( (rmis >= 2.*ck2) && (rho >= 1.e-10) )
  {
	if(zph >= 1.e-10)
	{
	  bk=cmplx(-zph,rho)*(ck1-cp3);
	  rmis=-creal(bk)/fabs(cimag(bk));
	  if(rmis > 4.*rho/zph)
		jump = TRUE;
	}

	if( ! jump )
	{
""", """
  /* jump: integrate up between the branch cuts (Fortran labels 6-7),
   * otherwise below the branch points (Fortran label 8) */
  jump = FALSE;
  if( (rmis >= 2.*ck2) && (rho >= 1.e-10) )
  {
	jump = TRUE;
	if(zph >= 1.e-10)
	{
	  bk=cmplx(-zph,rho)*(ck1-cp3);
	  rmis=-creal(bk)/fabs(cimag(bk));
	  if(rmis > 4.*rho/zph)
		jump = FALSE;
	}

	if( jump )
	{
"""),
    ("""
	jump = TRUE;

  } /* if( (rmis >= 2.*ck2) || (rho >= 1.e-10) ) */
  else
	jump = FALSE;

  if( ! jump )
""", """
  } /* if( (rmis >= 2.*ck2) || (rho >= 1.e-10) ) */

  if( ! jump )
"""),
    ("amg=fabs(creal(a1)+fabs(cimag(a1)));", "amg=fabs(creal(a1))+fabs(cimag(a1));"),
    ("""
  /* No convergence */
  abort_on_error(-6);
""", """
  /* No convergence: return the last estimate, as the Fortran SOMNEC does */
  for( i = 0; i < nans; i++ )
	sum[i]=.5*(q1[i][MAXH-1]+q2[i][MAXH-1]);
"""),
    ("""
  secnds(&tim);
  tim -= tst;
""", """
  secnds(&tim);
  tim -= tst;

  /* the interpolation must not reuse coefficients of the previous grid */
  intrp_new_grid = 1;
"""),
    ("""
/* common /cntour/ */
static complex double a, b;
""", """
/* common /cntour/ */
static complex double a, b;

extern int intrp_new_grid;
"""),
])

# calculations.c: zint selects one of its three formulas (the returns were
# lost and cc5 had the wrong exponent); intrp forgets its previous region
# when the grid has been recomputed.
patch("calculations.c", [
    ("#define cc5		( 0.         - I*9.765e4)", "#define cc5		( 0.         - I*9.765e-4)"),
    ("""
	  br1= br1/ br2;
	  *zint= CPLX_01* sqrt( cmotp/sigl )* br1/ rolam;
""", """
	  br1= br1/ br2;
	  *zint= CPLX_01* sqrt( cmotp/sigl )* br1/ rolam;
	  return;
"""),
    ("""
	br1= br1/ br2;
	*zint= CPLX_01* sqrt( cmotp/ sigl)* br1/ rolam;
""", """
	br1= br1/ br2;
	*zint= CPLX_01* sqrt( cmotp/ sigl)* br1/ rolam;
	return;
"""),
    ("""
  jump = FALSE;
  if( (x < xs) || (y < ys) )
	jump = TRUE;
""", """
  if( intrp_new_grid )
  {
	intrp_new_grid = 0;
	ixs=-10; iys=-10; igrs=-10; ixeg=0; iyeg=0;
	dx = 1.; dy = 1.; xs = 0.; ys = 0.;
  }

  jump = FALSE;
  if( (x < xs) || (y < ys) )
	jump = TRUE;
"""),
    ("/* intrp uses bivariate cubic interpolation to obtain */", """/* set by somnec() when the interpolation grid has been recomputed */
int intrp_new_grid = 0;

/* intrp uses bivariate cubic interpolation to obtain */"""),
])

# input.c: qdsrc sets ind2 (not ind1) for end 2 and indexes the connected
# segment consistently, as cmww does.
patch("input.c", [
    ("""
		ipr--;
		if( -data.icon1[ipr-1] != jp1 )
""", """
		ipr--;
		if( -data.icon1[ipr] != jp1 )
"""),
    ("""
		if( -data.icon2[ipr] != jp1 )
		  dataj.ind1=2;
		else
		{
		  xi= fabs( dataj.cabj* data.cab[ipr]+ dataj.sabj*
			  data.sab[ipr]+ dataj.salpj* data.salp[ipr]);
		  if( (xi < 0.999999) || (fabs(data.bi[ipr]/dataj.b-1.) > 1.0e-6) )
			dataj.ind1=2;
		  else
			dataj.ind1=0;
		}
""", """
		if( -data.icon2[ipr] != jp1 )
		  dataj.ind2=2;
		else
		{
		  xi= fabs( dataj.cabj* data.cab[ipr]+ dataj.sabj*
			  data.sab[ipr]+ dataj.salpj* data.salp[ipr]);
		  if( (xi < 0.999999) || (fabs(data.bi[ipr]/dataj.b-1.) > 1.0e-6) )
			dataj.ind2=2;
		  else
			dataj.ind2=0;
		}
"""),
    ("""
			  if( dataj.cabj* dataj.cabj+ dataj.sabj* dataj.sabj > 1.0e-8)
				dataj.ind1=2;
			  else
				dataj.ind1=0;
			}
		  } /* else */

	  } /* if( dataj.iexk != 0) */
""", """
			  if( dataj.cabj* dataj.cabj+ dataj.sabj* dataj.sabj > 1.0e-8)
				dataj.ind2=2;
			  else
				dataj.ind2=0;
			}
		  } /* else */

	  } /* if( dataj.iexk != 0) */
"""),
])

# ground.c: rom2 limits its halving by the constant 65536 of the Fortran,
# not by the number of segments.
patch("ground.c", [
    ("ep= s/(1.e4* data.npm);", "ep= s/(1.e4* 65536);"),
    ("if( ns < data.npm )", "if( ns < 65536 )"),
])

# matrix.c: a current source on the centre of a segment or of a patch gives
# zero entries of the right-hand side instead of leaving the old values, and
# the count of the patches goes on. In cmsw the component of a segment's own
# basis function is taken after trio has been called for that segment, as in
# the Fortran (AX(JSNO)).
patch("matrix.c", [
    ("""
	rs= pxl* pxl+ pyl* pyl+ pzl* pzl;
	if( rs < 1.0e-30)
	  continue;
""", """
	rs= pxl* pxl+ pyl* pyl+ pzl* pzl;
	if( rs < 1.0e-30)
	{
	  if( i < data.n ) e[i]=CPLX_00;
	  else { e[i1]=CPLX_00; e[i2]=CPLX_00; is++; }
	  continue;
	}
"""),
    ("""
  jsnox = segj.jsno-1;

  if( itrp >= 0)
""", """
  if( itrp >= 0)
"""),
    ("""
			  trio(i+1);

			  il= i-ncw;
""", """
			  trio(i+1);
			  jsnox = segj.jsno-1;

			  il= i-ncw;
"""),
])

# geometry.c: a strip of patches goes on with SC cards. The test for the SC
# card used the index of another card (the Fortran tests for 'SC').
patch("geometry.c", [
    ("""
	if( gm_num != 10 )
	  isct=0;
""", """
	if( gm_num != 9 )
	  isct=0;
"""),
])

# fields.c: nefld sets ind2 (not ind1) for end 2, as cmww does.
patch("fields.c", [
    ("""
		  if( -data.icon2[iprx] != ix )
			dataj.ind1=2;
		  else
		  {
			xi= fabs( dataj.cabj* data.cab[iprx]+ dataj.sabj*
				data.sab[iprx]+ dataj.salpj* data.salp[iprx]);
			if( (xi < 0.999999) || (fabs(data.bi[iprx]/dataj.b-1.) > 1.0e-6) )
			  dataj.ind1=2;
			else
			  dataj.ind1=0;
		  }
""", """
		  if( -data.icon2[iprx] != ix )
			dataj.ind2=2;
		  else
		  {
			xi= fabs( dataj.cabj* data.cab[iprx]+ dataj.sabj*
				data.sab[iprx]+ dataj.salpj* data.salp[iprx]);
			if( (xi < 0.999999) || (fabs(data.bi[iprx]/dataj.b-1.) > 1.0e-6) )
			  dataj.ind2=2;
			else
			  dataj.ind2=0;
		  }
"""),
    ("""
			  if( dataj.cabj* dataj.cabj+ dataj.sabj* dataj.sabj > 1.0e-8)
				dataj.ind1=2;
			  else
				dataj.ind1=0;
			}

		  } /* else */

	  } /* if( dataj.iexk != 0) */

	  efld( xob, yob, zob, ax,1);
""", """
			  if( dataj.cabj* dataj.cabj+ dataj.sabj* dataj.sabj > 1.0e-8)
				dataj.ind2=2;
			  else
				dataj.ind2=0;
			}

		  } /* else */

	  } /* if( dataj.iexk != 0) */

	  efld( xob, yob, zob, ax,1);
"""),
])

# network.c: a source that is not on a network port is not counted as
# network loss; the asymmetry check gets its work arrays without networks.
patch("network.c", [
    ("""
		irow1= ndimn- (j+1);
		cux= rhnx[irow1];
		for( j = 0; j < nteq; j++ )
		  cux -= cmn[j+irow1*ndimn]*rhnt[j];
		cux=(einc[isc1]+ cux)* data.wlam;
		irow1++;
""", """
		if( j == ntsc )
		{
		  cux= einc[isc1]* data.wlam;
		  irow1=0;
		}
		else
		{
		  irow1= ndimn- (j+1);
		  cux= rhnx[irow1];
		  for( j = 0; j < nteq; j++ )
			cux -= cmn[j+irow1*ndimn]*rhnt[j];
		  cux=(einc[isc1]+ cux)* data.wlam;
		  irow1++;
		}
"""),
    ("""
	  mreq = (size_t)j;
	  mreq *= sizeof(int);
	  mem_alloc( (void *)&ipnt, mreq );
	}
""", """
	  mreq = (size_t)j;
	  mreq *= sizeof(int);
	  mem_alloc( (void *)&ipnt, mreq );
	  mreq = (size_t)data.np3m;
	  mreq *= sizeof(complex double);
	  mem_alloc( (void *)&rhs, mreq );
	  mreq = (size_t)(j*j);
	  mreq *= sizeof(complex double);
	  mem_alloc( (void *)&cmn, mreq );
	}
"""),
])
print("patched", root)
