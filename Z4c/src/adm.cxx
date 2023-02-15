#include "z4c_vars.hxx"

#include <loop_device.hxx>
#include <mat.hxx>
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

extern "C" void Z4c_ADM(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_ADM;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const GF3D2<const CCTK_REAL> gf_chi1(layout1, chi);

  const smat<GF3D2<const CCTK_REAL>, 3, DN, DN> gf_gammat1{
      GF3D2<const CCTK_REAL>(layout1, gammatxx),
      GF3D2<const CCTK_REAL>(layout1, gammatxy),
      GF3D2<const CCTK_REAL>(layout1, gammatxz),
      GF3D2<const CCTK_REAL>(layout1, gammatyy),
      GF3D2<const CCTK_REAL>(layout1, gammatyz),
      GF3D2<const CCTK_REAL>(layout1, gammatzz)};

  const GF3D2<const CCTK_REAL> gf_Kh1(layout1, Kh);

  const smat<GF3D2<const CCTK_REAL>, 3, DN, DN> gf_At1{
      GF3D2<const CCTK_REAL>(layout1, Atxx),
      GF3D2<const CCTK_REAL>(layout1, Atxy),
      GF3D2<const CCTK_REAL>(layout1, Atxz),
      GF3D2<const CCTK_REAL>(layout1, Atyy),
      GF3D2<const CCTK_REAL>(layout1, Atyz),
      GF3D2<const CCTK_REAL>(layout1, Atzz)};

  const vec<GF3D2<const CCTK_REAL>, 3, UP> gf_Gamt1{
      GF3D2<const CCTK_REAL>(layout1, Gamtx),
      GF3D2<const CCTK_REAL>(layout1, Gamty),
      GF3D2<const CCTK_REAL>(layout1, Gamtz)};

  const GF3D2<const CCTK_REAL> gf_Theta1(layout1, Theta);

  const GF3D2<const CCTK_REAL> gf_alphaG1(layout1, alphaG);

  const vec<GF3D2<const CCTK_REAL>, 3, UP> gf_betaG1{
      GF3D2<const CCTK_REAL>(layout1, betaGx),
      GF3D2<const CCTK_REAL>(layout1, betaGy),
      GF3D2<const CCTK_REAL>(layout1, betaGz)};

  const GF3D2<const CCTK_REAL> gf_eTtt1(layout1, eTtt);

  const vec<GF3D2<const CCTK_REAL>, 3, DN> gf_eTti1{
      GF3D2<const CCTK_REAL>(layout1, eTtx),
      GF3D2<const CCTK_REAL>(layout1, eTty),
      GF3D2<const CCTK_REAL>(layout1, eTtz)};

  const smat<GF3D2<const CCTK_REAL>, 3, DN, DN> gf_eTij1{
      GF3D2<const CCTK_REAL>(layout1, eTxx),
      GF3D2<const CCTK_REAL>(layout1, eTxy),
      GF3D2<const CCTK_REAL>(layout1, eTxz),
      GF3D2<const CCTK_REAL>(layout1, eTyy),
      GF3D2<const CCTK_REAL>(layout1, eTyz),
      GF3D2<const CCTK_REAL>(layout1, eTzz)};

  const smat<GF3D2<CCTK_REAL>, 3, DN, DN> gf_g1{
      GF3D2<CCTK_REAL>(layout1, gxx), GF3D2<CCTK_REAL>(layout1, gxy),
      GF3D2<CCTK_REAL>(layout1, gxz), GF3D2<CCTK_REAL>(layout1, gyy),
      GF3D2<CCTK_REAL>(layout1, gyz), GF3D2<CCTK_REAL>(layout1, gzz)};

  const smat<GF3D2<CCTK_REAL>, 3, DN, DN> gf_K1{
      GF3D2<CCTK_REAL>(layout1, kxx), GF3D2<CCTK_REAL>(layout1, kxy),
      GF3D2<CCTK_REAL>(layout1, kxz), GF3D2<CCTK_REAL>(layout1, kyy),
      GF3D2<CCTK_REAL>(layout1, kyz), GF3D2<CCTK_REAL>(layout1, kzz)};

  const GF3D2<CCTK_REAL> gf_alp1(layout1, alp);

  const GF3D2<CCTK_REAL> gf_dtalp1(layout1, dtalp);

  const vec<GF3D2<CCTK_REAL>, 3, UP> gf_beta1{GF3D2<CCTK_REAL>(layout1, betax),
                                              GF3D2<CCTK_REAL>(layout1, betay),
                                              GF3D2<CCTK_REAL>(layout1, betaz)};

  const vec<GF3D2<CCTK_REAL>, 3, UP> gf_dtbeta1{
      GF3D2<CCTK_REAL>(layout1, dtbetax), GF3D2<CCTK_REAL>(layout1, dtbetay),
      GF3D2<CCTK_REAL>(layout1, dtbetaz)};

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_all_device<0, 0, 0, vsize>(
      grid.nghostzones,
      [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D2index index1(layout1, p.I);

        // Load and calculate
        const z4c_vars_noderivs<vreal> vars(
            kappa1, kappa2, f_mu_L, f_mu_S, eta, //
            gf_chi1(mask, index1, 1),
            gf_gammat1(mask, index1, one<smat<int, 3, DN, DN> >()()),
            gf_Kh1(mask, index1), gf_At1(mask, index1), gf_Gamt1(mask, index1),
            gf_Theta1(mask, index1), gf_alphaG1(mask, index1, 1),
            gf_betaG1(mask, index1), //
            gf_eTtt1(mask, index1), gf_eTti1(mask, index1),
            gf_eTij1(mask, index1));

        // Store
        gf_g1.store(mask, index1, vars.g);
        gf_K1.store(mask, index1, vars.K);
        gf_alp1.store(mask, index1, vars.alpha);
        gf_dtalp1.store(mask, index1, vars.dtalpha);
        gf_beta1.store(mask, index1, vars.beta);
        gf_dtbeta1.store(mask, index1, vars.dtbeta);
      });
}

} // namespace Z4c
