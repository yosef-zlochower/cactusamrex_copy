#include "derivs.hxx"
#include "physics.hxx"
#include "z4c_vars.hxx"

#include <loop_device.hxx>
#include <mat.hxx>
#include <mempool.hxx>
#include <simd.hxx>
#include <vec.hxx>

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace Z4c {
using namespace Arith;
using namespace Loop;
using namespace std;

extern "C" void Z4c_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_RHS;
  DECLARE_CCTK_PARAMETERS;
  int i, j, k;
  CCTK_LOOP3_ALL(Z4c, cctkGH, i,j,k) {
    const int ind = CCTK_GFINDEX3D(cctkGH, i,j,k);
    chi_rhs[ind] = 0.0;
    gammatxx_rhs[ind] = 0.0;
    gammatxy_rhs[ind] = 0.0;
    gammatxz_rhs[ind] = 0.0;
    gammatyy_rhs[ind] = 0.0;
    gammatyz_rhs[ind] = 0.0;
    gammatzz_rhs[ind] = 0.0;
    Kh_rhs[ind] = 0.0;
    Atxx_rhs[ind] = 0.0;
    Atxy_rhs[ind] = 0.0;
    Atxz_rhs[ind] = 0.0;
    Atyy_rhs[ind] = 0.0;
    Atyz_rhs[ind] = 0.0;
    Atzz_rhs[ind] = 0.0;
    Gamtx_rhs[ind] = 0.0;
    Gamty_rhs[ind] = 0.0;
    Gamtz_rhs[ind] = 0.0;
    Theta_rhs[ind] = 0.0;
    alphaG_rhs[ind] = 0.0;
    betaGx_rhs[ind] = 0.0;
    betaGy_rhs[ind] = 0.0;
    betaGz_rhs[ind] = 0.0;
  } CCTK_ENDLOOP3_ALL(Z4c);
}

} // namespace Z4c
