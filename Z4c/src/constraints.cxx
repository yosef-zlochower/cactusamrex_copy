#include "derivs.hxx"
#include "physics.hxx"
#include "z4c_vars.hxx"

#include <loop_device.hxx>
#include <mempool.hxx>
#include <simd.hxx>

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace Z4c {
using namespace Arith;
using namespace Loop;
using namespace std;

extern "C" void Z4c_Constraints(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Constraints;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < deriv_order / 2 + 1)
      CCTK_VERROR("Need at least %d ghost zones", deriv_order / 2 + 1);

  //

  const array<int, dim> indextype = {0, 0, 0};
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};
  vect<int, dim> imin, imax;
  GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  // Suffix 1: with ghost zones, suffix 0: without ghost zones
  const GF3D2layout layout1(cctkGH, indextype);
  const GF3D5layout layout0(imin, imax);

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

  //

  const size_t mempool_id = GetCallFunctionCount();
  mempool_t &restrict mempool = mempools.get_mempool(mempool_id);

  const auto make_gf = [&]() { return GF3D5<CCTK_REAL>(layout0, mempool); };
  const auto make_vec = [&](const auto &f) {
    return vec<result_of_t<decltype(f)()>, 3, DN>([&](int) { return f(); });
  };
  const auto make_mat = [&](const auto &f) {
    return smat<result_of_t<decltype(f)()>, 3, DN, DN>(
        [&](int, int) { return f(); });
  };
  const auto make_vec_gf = [&]() { return make_vec(make_gf); };
  const auto make_mat_gf = [&]() { return make_mat(make_gf); };
  const auto make_vec_vec_gf = [&]() { return make_vec(make_vec_gf); };
  // const auto make_vec_mat_gf = [&]() { return make_vec(make_mat_gf); };
  const auto make_mat_vec_gf = [&]() { return make_mat(make_vec_gf); };
  const auto make_mat_mat_gf = [&]() { return make_mat(make_mat_gf); };

  const GF3D5<CCTK_REAL> gf_chi0(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3, DN> gf_dchi0(make_vec_gf());
  const smat<GF3D5<CCTK_REAL>, 3, DN, DN> gf_ddchi0(make_mat_gf());
  calc_derivs2(cctkGH, gf_chi1, gf_chi0, gf_dchi0, gf_ddchi0, layout0);

  const smat<GF3D5<CCTK_REAL>, 3, DN, DN> gf_gammat0(make_mat_gf());
  const smat<vec<GF3D5<CCTK_REAL>, 3, DN>, 3, DN, DN> gf_dgammat0(
      make_mat_vec_gf());
  const smat<smat<GF3D5<CCTK_REAL>, 3, DN, DN>, 3, DN, DN> gf_ddgammat0(
      make_mat_mat_gf());
  calc_derivs2(cctkGH, gf_gammat1, gf_gammat0, gf_dgammat0, gf_ddgammat0,
               layout0);

  const GF3D5<CCTK_REAL> gf_Kh0(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3, DN> gf_dKh0(make_vec_gf());
  calc_derivs(cctkGH, gf_Kh1, gf_Kh0, gf_dKh0, layout0);

  const smat<GF3D5<CCTK_REAL>, 3, DN, DN> gf_At0(make_mat_gf());
  const smat<vec<GF3D5<CCTK_REAL>, 3, DN>, 3, DN, DN> gf_dAt0(
      make_mat_vec_gf());
  calc_derivs(cctkGH, gf_At1, gf_At0, gf_dAt0, layout0);

  const vec<GF3D5<CCTK_REAL>, 3, UP> gf_Gamt0(make_vec_gf());
  const vec<vec<GF3D5<CCTK_REAL>, 3, DN>, 3, UP> gf_dGamt0(make_vec_vec_gf());
  calc_derivs(cctkGH, gf_Gamt1, gf_Gamt0, gf_dGamt0, layout0);

  const GF3D5<CCTK_REAL> gf_Theta0(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3, DN> gf_dTheta0(make_vec_gf());
  calc_derivs(cctkGH, gf_Theta1, gf_Theta0, gf_dTheta0, layout0);

  //

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

  //

  const vec<GF3D2<CCTK_REAL>, 3, UP> gf_ZtC1{GF3D2<CCTK_REAL>(layout1, ZtCx),
                                             GF3D2<CCTK_REAL>(layout1, ZtCy),
                                             GF3D2<CCTK_REAL>(layout1, ZtCz)};

  const GF3D2<CCTK_REAL> gf_HC1(layout1, HC);

  const vec<GF3D2<CCTK_REAL>, 3, UP> gf_MtC1{GF3D2<CCTK_REAL>(layout1, MtCx),
                                             GF3D2<CCTK_REAL>(layout1, MtCy),
                                             GF3D2<CCTK_REAL>(layout1, MtCz)};

  const GF3D2<CCTK_REAL> gf_allC1(layout1, allC);

  //

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones,
      [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D2index index1(layout1, p.I);
        const GF3D5index index0(layout0, p.I);

        // Load and calculate

        const auto alphaG = one<CCTK_REAL>()();
        const auto dalphaG = zero<vec<CCTK_REAL, 3, DN> >()();
        const auto ddalphaG = zero<smat<CCTK_REAL, 3, DN, DN> >()();

        const auto betaG = zero<vec<CCTK_REAL, 3, UP> >()();
        const auto dbetaG = zero<vec<vec<CCTK_REAL, 3, DN>, 3, UP> >()();
        const auto ddbetaG = zero<vec<smat<CCTK_REAL, 3, DN, DN>, 3, UP> >()();

        const z4c_vars<vreal> vars(
            kappa1, kappa2, f_mu_L, f_mu_S, eta, //
            gf_chi0(mask, index0, 1), gf_dchi0(mask, index0),
            gf_ddchi0(mask, index0), //
            gf_gammat0(mask, index0, 1), gf_dgammat0(mask, index0),
            gf_ddgammat0(mask, index0),                           //
            gf_Kh0(mask, index0), gf_dKh0(mask, index0),          //
            gf_At0(mask, index0), gf_dAt0(mask, index0),          //
            gf_Gamt0(mask, index0), gf_dGamt0(mask, index0),      //
            gf_Theta0(mask, index0, 1), gf_dTheta0(mask, index0), //
            alphaG, dalphaG, ddalphaG,                            //
            betaG, dbetaG, ddbetaG,                               //
            gf_eTtt1(mask, index1), gf_eTti1(mask, index1),
            gf_eTij1(mask, index1));

        // Store
        gf_ZtC1.store(mask, index1, vars.ZtC);
        gf_HC1.store(mask, index1, vars.HC);
        gf_MtC1.store(mask, index1, vars.MtC);
        gf_allC1.store(mask, index1, vars.allC);
      });
}

} // namespace Z4c
