#include "derivs.hxx"
#include "physics.hxx"
#include "tensor.hxx"
#include "z4c_vars.hxx"

#include <loop_device.hxx>
#include <mempool.hxx>
#include <simd.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
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
  const GF3D2layout layout0(imin, imax);

  const GF3D2<const CCTK_REAL> gf_chi1(layout1, chi);

  const mat3<GF3D2<const CCTK_REAL>, DN, DN> gf_gammat1(
      GF3D2<const CCTK_REAL>(layout1, gammatxx),
      GF3D2<const CCTK_REAL>(layout1, gammatxy),
      GF3D2<const CCTK_REAL>(layout1, gammatxz),
      GF3D2<const CCTK_REAL>(layout1, gammatyy),
      GF3D2<const CCTK_REAL>(layout1, gammatyz),
      GF3D2<const CCTK_REAL>(layout1, gammatzz));

  const GF3D2<const CCTK_REAL> gf_Kh1(layout1, Kh);

  const mat3<GF3D2<const CCTK_REAL>, DN, DN> gf_At1(
      GF3D2<const CCTK_REAL>(layout1, Atxx),
      GF3D2<const CCTK_REAL>(layout1, Atxy),
      GF3D2<const CCTK_REAL>(layout1, Atxz),
      GF3D2<const CCTK_REAL>(layout1, Atyy),
      GF3D2<const CCTK_REAL>(layout1, Atyz),
      GF3D2<const CCTK_REAL>(layout1, Atzz));

  const vec3<GF3D2<const CCTK_REAL>, UP> gf_Gamt1(
      GF3D2<const CCTK_REAL>(layout1, Gamtx),
      GF3D2<const CCTK_REAL>(layout1, Gamty),
      GF3D2<const CCTK_REAL>(layout1, Gamtz));

  const GF3D2<const CCTK_REAL> gf_Theta1(layout1, Theta);

  //

  const size_t mempool_id = GetCallFunctionCount();
  mempool_t &restrict mempool = mempools.get_mempool(mempool_id);

  const auto make_gf = [&]() { return GF3D5<CCTK_REAL>(layout0, mempool); };
  const auto make_vec_gf = [&](int) { return make_gf(); };
  const auto make_mat_gf = [&](int, int) { return make_gf(); };
  const auto make_vec_vec_gf = [&](int) { return make_vec_gf; };
  // const auto make_vec_mat_gf = [&](int) { return make_mat_gf; };
  const auto make_mat_vec_gf = [&](int, int) { return make_vec_gf; };
  const auto make_mat_mat_gf = [&](int, int) { return make_mat_gf; };

  const GF3D5<CCTK_REAL> gf_chi0(make_gf());
  const vec3<GF3D5<CCTK_REAL>, DN> gf_dchi0(make_vec_gf);
  const mat3<GF3D5<CCTK_REAL>, DN, DN> gf_ddchi0(make_mat_gf);
  calc_derivs2(cctkGH, gf_chi1, gf_chi0, gf_dchi0, gf_ddchi0, layout0);

  const mat3<GF3D5<CCTK_REAL>, DN, DN> gf_gammat0(make_mat_gf);
  const mat3<vec3<GF3D5<CCTK_REAL>, DN>, DN, DN> gf_dgammat0(make_mat_vec_gf);
  const mat3<mat3<GF3D5<CCTK_REAL>, DN, DN>, DN, DN> gf_ddgammat0(
      make_mat_mat_gf);
  calc_derivs2(cctkGH, gf_gammat1, gf_gammat0, gf_dgammat0, gf_ddgammat0,
               layout0);

  const GF3D5<CCTK_REAL> gf_Kh0(make_gf());
  const vec3<GF3D5<CCTK_REAL>, DN> gf_dKh0(make_vec_gf);
  calc_derivs(cctkGH, gf_Kh1, gf_Kh0, gf_dKh0, layout0);

  const mat3<GF3D5<CCTK_REAL>, DN, DN> gf_At0(make_mat_gf);
  const mat3<vec3<GF3D5<CCTK_REAL>, DN>, DN, DN> gf_dAt0(make_mat_vec_gf);
  calc_derivs(cctkGH, gf_At1, gf_At0, gf_dAt0, layout0);

  const vec3<GF3D5<CCTK_REAL>, UP> gf_Gamt0(make_vec_gf);
  const vec3<vec3<GF3D5<CCTK_REAL>, DN>, UP> gf_dGamt0(make_vec_vec_gf);
  calc_derivs(cctkGH, gf_Gamt1, gf_Gamt0, gf_dGamt0, layout0);

  const GF3D5<CCTK_REAL> gf_Theta0(make_gf());
  const vec3<GF3D5<CCTK_REAL>, DN> gf_dTheta0(make_vec_gf);
  calc_derivs(cctkGH, gf_Theta1, gf_Theta0, gf_dTheta0, layout0);

  //

  const GF3D2<const CCTK_REAL> gf_eTtt1(layout1, eTtt);

  const vec3<GF3D2<const CCTK_REAL>, DN> gf_eTti1(
      GF3D2<const CCTK_REAL>(layout1, eTtx),
      GF3D2<const CCTK_REAL>(layout1, eTty),
      GF3D2<const CCTK_REAL>(layout1, eTtz));

  const mat3<GF3D2<const CCTK_REAL>, DN, DN> gf_eTij1(
      GF3D2<const CCTK_REAL>(layout1, eTxx),
      GF3D2<const CCTK_REAL>(layout1, eTxy),
      GF3D2<const CCTK_REAL>(layout1, eTxz),
      GF3D2<const CCTK_REAL>(layout1, eTyy),
      GF3D2<const CCTK_REAL>(layout1, eTyz),
      GF3D2<const CCTK_REAL>(layout1, eTzz));

  //

  const vec3<GF3D2<CCTK_REAL>, UP> gf_ZtC1(GF3D2<CCTK_REAL>(layout1, ZtCx),
                                           GF3D2<CCTK_REAL>(layout1, ZtCy),
                                           GF3D2<CCTK_REAL>(layout1, ZtCz));

  const GF3D2<CCTK_REAL> gf_HC1(layout1, HC);

  const vec3<GF3D2<CCTK_REAL>, UP> gf_MtC1(GF3D2<CCTK_REAL>(layout1, MtCx),
                                           GF3D2<CCTK_REAL>(layout1, MtCy),
                                           GF3D2<CCTK_REAL>(layout1, MtCz));

  const GF3D2<CCTK_REAL> gf_allC1(layout1, allC);

  //

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones, [=] Z4C_INLINE Z4C_GPU(const PointDesc &p) {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);

        // Load and calculate

        const vreal alphaG{1};
        const vec3<vreal, DN> dalphaG{0, 0, 0};
        const mat3<vreal, DN, DN> ddalphaG{0, 0, 0, 0, 0, 0};

        const vec3<vreal, UP> betaG{0, 0, 0};
        const vec3<vec3<vreal, DN>, UP> dbetaG{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        const vec3<mat3<vreal, DN, DN>, UP> ddbetaG{
            {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};

        const z4c_vars<vreal> vars(
            kappa1, kappa2, f_mu_L, f_mu_S, eta, //
            gf_chi0(mask, layout0, p.I, 1), gf_dchi0(mask, layout0, p.I),
            gf_ddchi0(mask, layout0, p.I), //
            gf_gammat0(mask, layout0, p.I, 1), gf_dgammat0(mask, layout0, p.I),
            gf_ddgammat0(mask, layout0, p.I),                                 //
            gf_Kh0(mask, layout0, p.I), gf_dKh0(mask, layout0, p.I),          //
            gf_At0(mask, layout0, p.I), gf_dAt0(mask, layout0, p.I),          //
            gf_Gamt0(mask, layout0, p.I), gf_dGamt0(mask, layout0, p.I),      //
            gf_Theta0(mask, layout0, p.I, 1), gf_dTheta0(mask, layout0, p.I), //
            alphaG, dalphaG, ddalphaG,                                        //
            betaG, dbetaG, ddbetaG,                                           //
            gf_eTtt1(mask, p.I), gf_eTti1(mask, p.I), gf_eTij1(mask, p.I));

        // Store
        gf_ZtC1.store(mask, p.I, vars.ZtC);
        gf_HC1.store(mask, p.I, vars.HC);
        gf_MtC1.store(mask, p.I, vars.MtC);
        gf_allC1.store(mask, p.I, vars.allC);
      });
}

} // namespace Z4c
