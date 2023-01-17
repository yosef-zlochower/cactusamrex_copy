#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include "utils.hxx"
#include "reconstruct.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;

enum class flux_t { LxF, HLLE };

// Calculate the fluxes in direction `dir`. This function is more
// complex because it has to handle any direction, but as reward,
// there is only one function, not three.
template <int dir> void CalcFlux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  static_assert(dir >= 0 && dir < 3, "");

  reconstruction_t reconstruction;
  if (CCTK_EQUALS(reconstruction_method, "Godunov"))
    reconstruction = reconstruction_t::Godunov;
  else if (CCTK_EQUALS(reconstruction_method, "minmod"))
    reconstruction = reconstruction_t::minmod;
  else if (CCTK_EQUALS(reconstruction_method, "monocentral"))
    reconstruction = reconstruction_t::monocentral;
  else if (CCTK_EQUALS(reconstruction_method, "ppm"))
    reconstruction = reconstruction_t::ppm;
  else
    CCTK_ERROR("Unknown value for parameter \"reconstruction_method\"");

  flux_t fluxtype;
  if (CCTK_EQUALS(flux_type, "LxF")) {
    fluxtype = flux_t::LxF;
  } else if (CCTK_EQUALS(flux_type, "HLLE")) {
    fluxtype = flux_t::HLLE;
  } else {
    CCTK_ERROR("Unknown value for parameter \"flux_type\"");
  }

  switch (reconstruction) {
  case reconstruction_t::Godunov:
    assert(cctk_nghostzones[dir] >= 1);
    break;
  case reconstruction_t::minmod:
    assert(cctk_nghostzones[dir] >= 2);
    break;
  case reconstruction_t::monocentral:
    assert(cctk_nghostzones[dir] >= 2);
    break;
  case reconstruction_t::ppm:
    assert(cctk_nghostzones[dir] >= 3);
    break;
  }

  const auto reconstruct_pt =
      [=] CCTK_DEVICE(const GF3D2<const CCTK_REAL> &var, const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return reconstruct(var, p, reconstruction, dir);
          };

  const auto eigenvalues =
      [=] CCTK_DEVICE(CCTK_REAL alp_avg, CCTK_REAL beta_avg, CCTK_REAL u_avg,
                      vec<CCTK_REAL, 2> vel, vec<CCTK_REAL, 2> rho,
                      vec<CCTK_REAL, 2> cs2, vec<CCTK_REAL, 2> w_lor,
                      vec<CCTK_REAL, 2> h,
                      vec<CCTK_REAL, 2> bsq) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // computing characteristics for the minus side
        // See Eq. (28) of Giacomazzo & Rezzolla (2007) with b^i=0
        vec<CCTK_REAL, 3> a_m{
            (bsq(0) + cs2(0) * h(0) * rho(0)) *
                    (pow2(beta_avg) - pow2(alp_avg) * u_avg) -
                (-1 + cs2(0)) * h(0) * rho(0) *
                    pow2(beta_avg - alp_avg * vel(0)) * pow2(w_lor(0)),

            2 * beta_avg * (bsq(0) + cs2(0) * h(0) * rho(0)) -
                2 * (-1 + cs2(0)) * h(0) * rho(0) *
                    (beta_avg - alp_avg * vel(0)) * pow2(w_lor(0)),

            bsq(0) + h(0) * rho(0) *
                         (cs2(0) + pow2(w_lor(0)) - cs2(0) * pow2(w_lor(0)))};

        CCTK_REAL det_m = pow2(a_m(1)) - 4.0 * a_m(2) * a_m(0);
        if (det_m < 0.0)
          det_m = 0.0;

        vec<CCTK_REAL, 4> lambda_m{
            ((-a_m(1) + sqrt(det_m)) / (2.0 * a_m(2))) / alp_avg,
            ((-a_m(1) + sqrt(det_m)) / (2.0 * a_m(2))) / alp_avg,
            ((-a_m(1) - sqrt(det_m)) / (2.0 * a_m(2))) / alp_avg,
            ((-a_m(1) - sqrt(det_m)) / (2.0 * a_m(2))) / alp_avg};

        // computing characteristics for the plus side

        vec<CCTK_REAL, 3> a_p{
            (bsq(1) + cs2(1) * h(1) * rho(1)) *
                    (pow2(beta_avg) - pow2(alp_avg) * u_avg) -
                (-1 + cs2(1)) * h(1) * rho(1) *
                    pow2(beta_avg - alp_avg * vel(1)) * pow2(w_lor(1)),

            2 * beta_avg * (bsq(1) + cs2(1) * h(1) * rho(1)) -
                2 * (-1 + cs2(1)) * h(1) * rho(1) *
                    (beta_avg - alp_avg * vel(1)) * pow2(w_lor(1)),

            bsq(1) + h(1) * rho(1) *
                         (cs2(1) + pow2(w_lor(1)) - cs2(1) * pow2(w_lor(1)))};

        CCTK_REAL det_p = pow2(a_p(1)) - 4.0 * a_p(2) * a_p(0);
        if (det_p < 0.0)
          det_p = 0.0;

        vec<CCTK_REAL, 4> lambda_p{
            ((-a_p(1) + sqrt(det_p)) / (2.0 * a_p(2))) / alp_avg,
            ((-a_p(1) + sqrt(det_p)) / (2.0 * a_p(2))) / alp_avg,
            ((-a_p(1) - sqrt(det_p)) / (2.0 * a_p(2))) / alp_avg,
            ((-a_p(1) - sqrt(det_p)) / (2.0 * a_p(2))) / alp_avg};

        // 2D array containing characteristics for left (minus) and right (plus)
        // sides
        vec<vec<CCTK_REAL, 4>, 2> lambda{lambda_m, lambda_p};
        return lambda;
      };

  const auto calcflux =
      [=] CCTK_DEVICE(vec<vec<CCTK_REAL, 4>, 2> lam, vec<CCTK_REAL, 2> var,
                      vec<CCTK_REAL, 2> flux) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL flx;
        switch (fluxtype) {
        case flux_t::LxF: {
          const CCTK_REAL charmax =
              max({0.0, fabs(lam(0)(0)), fabs(lam(0)(1)), fabs(lam(0)(2)),
                   fabs(lam(0)(3)), fabs(lam(1)(0)), fabs(lam(1)(1)),
                   fabs(lam(1)(2)), fabs(lam(1)(3))});

          flx = 0.5 * ((flux(0) + flux(1)) - charmax * (var(1) - var(0)));
          break;
        }

        case flux_t::HLLE: {
          const CCTK_REAL charmax =
              max({0.0, lam(0)(0), lam(0)(1), lam(0)(2), lam(0)(3), lam(1)(0),
                   lam(1)(1), lam(1)(2), lam(1)(3)});

          const CCTK_REAL charmin =
              min({0.0, lam(0)(0), lam(0)(1), lam(0)(2), lam(0)(3), lam(1)(0),
                   lam(1)(1), lam(1)(2), lam(1)(3)});

          const CCTK_REAL charpm = charmax - charmin;

          flx = (charmax * flux(1) - charmin * flux(0) +
                 charmax * charmin * (var(1) - var(0))) /
                charpm;
          break;
        }

        default:
          assert(0);
        }

        return flx;
      };

  /* grid functions for fluxes */
  const vec<GF3D2<CCTK_REAL>, dim> fluxdenss{fxdens, fydens, fzdens};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomxs{fxmomx, fymomx, fzmomx};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomys{fxmomy, fymomy, fzmomy};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomzs{fxmomz, fymomz, fzmomz};
  const vec<GF3D2<CCTK_REAL>, dim> fluxtaus{fxtau, fytau, fztau};
  const vec<GF3D2<CCTK_REAL>, dim> fluxBxs{fxBx, fyBx, fzBx};
  const vec<GF3D2<CCTK_REAL>, dim> fluxBys{fxBy, fyBy, fzBy};
  const vec<GF3D2<CCTK_REAL>, dim> fluxBzs{fxBz, fyBz, fzBz};
  /* grid functions */
  const vec<GF3D2<const CCTK_REAL>, dim> gf_vels{velx, vely, velz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_Bvecs{Bvecx, Bvecy, Bvecz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_beta{betax, betay, betaz};
  const smat<GF3D2<const CCTK_REAL>, dim> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};


  // Cell-centered layout with one ghost cell in each direction
  const array<int, dim> ngh_m1 = {
    cctk_nghostzones[0] - 1,
    cctk_nghostzones[1] - 1,
    cctk_nghostzones[2] - 1
  };

  vect<int, dim> imin, imax;
  GridDescBase(cctkGH).box_intp1<1, 1, 1>(ngh_m1, imin, imax);
  const GF3D5layout layout(imin, imax);


  /* Two temporary tile-sized local arrays for each primitive to store the two
   * reconstructed states from the center to the faces of each cell             */
  // TODO: use the staggered magnetic field instead of reconstructing it
  const int ntmps = 16;
  GF3D5vector<CCTK_REAL> tmps(layout, ntmps);

  const vec<GF3D5<CCTK_REAL>, 2> tile_rho_rc_cell   {tmps(0),  tmps(1)};
  const vec<GF3D5<CCTK_REAL>, 2> tile_eps_rc_cell   {tmps(2),  tmps(3)};
  const vec<GF3D5<CCTK_REAL>, 2> tile_velx_rc_cell  {tmps(4),  tmps(5)};
  const vec<GF3D5<CCTK_REAL>, 2> tile_vely_rc_cell  {tmps(6),  tmps(7)};
  const vec<GF3D5<CCTK_REAL>, 2> tile_velz_rc_cell  {tmps(8),  tmps(9)};
  const vec<GF3D5<CCTK_REAL>, 2> tile_Bvecx_rc_cell {tmps(10), tmps(11)};
  const vec<GF3D5<CCTK_REAL>, 2> tile_Bvecy_rc_cell {tmps(12), tmps(13)};
  const vec<GF3D5<CCTK_REAL>, 2> tile_Bvecz_rc_cell {tmps(14), tmps(15)};


  // Loop over all internal cells plus one ghost cell on each side
  grid.loop_intp1_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE { 
        // Reconstruct the primitives on both faces of each cell 
        const vec<CCTK_REAL, 2> rho_rc_cell   {reconstruct_pt(rho,   p)};
        const vec<CCTK_REAL, 2> eps_rc_cell   {reconstruct_pt(eps,   p)};
        const vec<CCTK_REAL, 2> velx_rc_cell  {reconstruct_pt(velx,  p)};
        const vec<CCTK_REAL, 2> vely_rc_cell  {reconstruct_pt(vely,  p)};
        const vec<CCTK_REAL, 2> velz_rc_cell  {reconstruct_pt(velz,  p)};
        const vec<CCTK_REAL, 2> Bvecx_rc_cell {reconstruct_pt(Bvecx, p)};
        const vec<CCTK_REAL, 2> Bvecy_rc_cell {reconstruct_pt(Bvecy, p)};
        const vec<CCTK_REAL, 2> Bvecz_rc_cell {reconstruct_pt(Bvecz, p)};

        // Get the GF3D5 index corresponding to the current point
        const GF3D5index index(layout, p.I);

        /* Save the reconstructed values in the tile-sized local arrays built
         * above                                                                */
        tile_rho_rc_cell(0).store(index,   rho_rc_cell(0));
        tile_eps_rc_cell(0).store(index,   eps_rc_cell(0));
        tile_velx_rc_cell(0).store(index,  velx_rc_cell(0));
        tile_vely_rc_cell(0).store(index,  vely_rc_cell(0));
        tile_velz_rc_cell(0).store(index,  velz_rc_cell(0));
        tile_Bvecx_rc_cell(0).store(index, Bvecx_rc_cell(0));
        tile_Bvecy_rc_cell(0).store(index, Bvecy_rc_cell(0));
        tile_Bvecz_rc_cell(0).store(index, Bvecz_rc_cell(0));

        tile_rho_rc_cell(1).store(index,   rho_rc_cell(1));
        tile_eps_rc_cell(1).store(index,   eps_rc_cell(1));
        tile_velx_rc_cell(1).store(index,  velx_rc_cell(1));
        tile_vely_rc_cell(1).store(index,  vely_rc_cell(1));
        tile_velz_rc_cell(1).store(index,  velz_rc_cell(1));
        tile_Bvecx_rc_cell(1).store(index, Bvecx_rc_cell(1));
        tile_Bvecy_rc_cell(1).store(index, Bvecy_rc_cell(1));
        tile_Bvecz_rc_cell(1).store(index, Bvecz_rc_cell(1));
    });



  /* Loop over all interior faces (including the ones delimiting a ghost cell on
   * each side)                                                                 */
  constexpr auto DI = PointDesc::DI;

  constexpr array<int, dim> face_centred = {
    not (dir == 0), not (dir == 1), not (dir == 2)
  };

  grid.loop_int_device<face_centred[0], face_centred[1], face_centred[2]>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* Get the GF3D5 indices corresponding to the current and previous
         * points                                                               */
        const GF3D5index index   (layout, p.I);
        const GF3D5index index_m1(layout, p.I - DI[dir]);

        /* Retrieve the values of the primitives reconstructed on each side of
         * the current face from the tile-sized local arrays                    */
        const vec<CCTK_REAL, 2> rho_rc_face {
          tile_rho_rc_cell(1)(index_m1),
          tile_rho_rc_cell(0)(index)
        };

        const vec<CCTK_REAL, 2> eps_rc_face {
          tile_eps_rc_cell(1)(index_m1),
          tile_eps_rc_cell(0)(index)
        };

        const vec<CCTK_REAL, 2> velx_rc_face {
          tile_velx_rc_cell(1)(index_m1),
          tile_velx_rc_cell(0)(index)
        };

        const vec<CCTK_REAL, 2> vely_rc_face {
          tile_vely_rc_cell(1)(index_m1),
          tile_vely_rc_cell(0)(index)
        };

        const vec<CCTK_REAL, 2> velz_rc_face {
          tile_velz_rc_cell(1)(index_m1),
          tile_velz_rc_cell(0)(index)
        };

        const vec<CCTK_REAL, 2> Bvecx_rc_face {
          tile_Bvecx_rc_cell(1)(index_m1),
          tile_Bvecx_rc_cell(0)(index)
        };

        const vec<CCTK_REAL, 2> Bvecy_rc_face {
          tile_Bvecy_rc_cell(1)(index_m1),
          tile_Bvecy_rc_cell(0)(index)
        };

        const vec<CCTK_REAL, 2> Bvecz_rc_face {
          tile_Bvecz_rc_cell(1)(index_m1),
          tile_Bvecz_rc_cell(0)(index)
       };


        // Auxiliary vectors
        const vec<vec<CCTK_REAL, 2>, 3> vels_rc_face {
          velx_rc_face, vely_rc_face, velz_rc_face
        };

        const vec<vec<CCTK_REAL, 2>, 3> Bs_rc_face {
          Bvecx_rc_face, Bvecy_rc_face, Bvecz_rc_face
        };


        // Interpolate spacetime variables from the vertices to the current face
        const CCTK_REAL alp_avg = calc_avg_v2f(alp, p, dir);
        const vec<CCTK_REAL, 3> betas_avg([&](int i) ARITH_INLINE {
          return calc_avg_v2f(gf_beta(i), p, dir);
        });
        const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2f(gf_g(i, j), p, dir);
        });

        // Determinant of the spatial metric
        const CCTK_REAL detg_avg = calc_det(g_avg);
        const CCTK_REAL sqrtg    = sqrt(detg_avg);

        // Co-velocity measured by the Eulerian observer: v_j
        const vec<vec<CCTK_REAL, 2>, 3> vlows_rc_face =
            calc_contraction(g_avg, vels_rc_face);

        // vtilde^i = alpha * v^i - beta^i
        const vec<vec<CCTK_REAL, 2>, 3> vtildes_rc_face([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return alp_avg * vels_rc_face(i)(f) - betas_avg(i);
          });
        });

        // Lorentz factor
        const vec<CCTK_REAL, 2> w_lorentz_rc_face([&](int f) ARITH_INLINE {
          return 1.0 / sqrt(1.0 - calc_contraction(vlows_rc_face, vels_rc_face)(f));
        });

        // alpha * b0 = W * B^i * v_i
        const vec<CCTK_REAL, 2> alp_b0_rc_face([&](int f) ARITH_INLINE {
          return w_lorentz_rc_face(f) * calc_contraction(Bs_rc_face, vlows_rc_face)(f);
        });

        // Covariant magnetic field measured by the Eulerian observer
        const vec<vec<CCTK_REAL, 2>, 3> Blows_rc_face =
            calc_contraction(g_avg, Bs_rc_face);

        // B^2
        const vec<CCTK_REAL, 2> B2_rc_face = calc_contraction(Bs_rc_face, Blows_rc_face);

        /* Covariant magnetic field measured by the comoving observer:
         *  b_i = B_i / W + alpha * b^0 * v_i                                   */
        const vec<vec<CCTK_REAL, 2>, 3> blows_rc_face([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return Blows_rc_face(i)(f) / w_lorentz_rc_face(f) +
                   alp_b0_rc_face(f) * vlows_rc_face(i)(f);
          });
        });

        // b^2 = b^\mu * b_\mu
        const vec<CCTK_REAL, 2> b2_rc_face([&](int f) ARITH_INLINE {
          return (B2_rc_face(f) + pow2(alp_b0_rc_face(f))) / pow2(w_lorentz_rc_face(f));
        });

        // Components of some vectors along the current direction
        const CCTK_REAL beta_avg = betas_avg(dir);
        const vec<CCTK_REAL, 2> vel_rc_face{vels_rc_face(dir)};
        const vec<CCTK_REAL, 2> B_rc_face{Bs_rc_face(dir)};
        const vec<CCTK_REAL, 2> vtilde_rc_face{vtildes_rc_face(dir)};

        // Pressure (ideal gas EOS)
        // TODO: compute this from a user-specified EOS
        const vec<CCTK_REAL, 2> press_rc_face([&](int f) ARITH_INLINE {
          return eps_rc_face(f) * rho_rc_face(f) * (gamma - 1.0);
        });

        // Sound speed squared (ideal gas EOS )
        // TODO: compute this from a user-specified EOS
        const vec<CCTK_REAL, 2> cs2_rc_face([&](int f) ARITH_INLINE {
          return (gamma - 1.0) * eps_rc_face(f) / (eps_rc_face(f) + 1.0 / gamma);
        });

        // Enthalpy for ideal gas EOS
        // TODO: compute this from a user-specified EOS
        const vec<CCTK_REAL, 2> h_rc_face([&](int f) ARITH_INLINE {
          return 1.0 + eps_rc_face(f) + press_rc_face(f) / rho_rc_face(f);
        });


        /* Compute the conservatives from the primitives *
         * --------------------------------------------- */
        // dens = sqrt(g) * D = sqrt(g) * rho * W
        const vec<CCTK_REAL, 2> dens_rc_face([&](int f) ARITH_INLINE {
          return sqrtg * rho_rc_face(f) * w_lorentz_rc_face(f);
        });

        // Auxiliary: dens * h * W = sqrt(g) * rho * h * W^2
        const vec<CCTK_REAL, 2> dens_h_W_rc_face([&](int f) ARITH_INLINE {
          return dens_rc_face(f) * h_rc_face(f) * w_lorentz_rc_face(f);
        });

        // Auxiliary: sqrt(g) * (rho*h + b^2) * W^2
        const vec<CCTK_REAL, 2> dens_h_W_plus_sqrtg_W2b2_rc_face =
          dens_h_W_rc_face + sqrtg * (pow2(alp_b0_rc_face) + B2_rc_face);

        // Auxiliary: ptot = pgas + pmag (pmag = b^2/2)
        const vec<CCTK_REAL, 2> ptot_rc_face = press_rc_face + 0.5 * b2_rc_face;

        // mom_i = sqrt(g) * S_i = sqrt(g)*((rho*h + b^2)*W^2*v_i - alpha*b^0*b_i)
        const vec<vec<CCTK_REAL, 2>, 3> moms_rc_face([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return dens_h_W_plus_sqrtg_W2b2_rc_face(f) * vlows_rc_face(i)(f) -
                   sqrtg * alp_b0_rc_face(f) * blows_rc_face(i)(f);
          });
        });

        /* tau = sqrt(g)*t =
         *   sqrt(g)((rho*h + b^2)*W^2 - (pgas+pmag) - (alpha*b^0)^2 - D)       */
        const vec<CCTK_REAL, 2> tau_rc_face =
          dens_h_W_rc_face - dens_rc_face + sqrtg * (B2_rc_face - ptot_rc_face);

        // Btildes^i = sqrt(g) * B^i
        const vec<vec<CCTK_REAL, 2>, 3> Btildes_rc_face(
            [&](int i) ARITH_INLINE { return sqrtg * Bs_rc_face(i); });


        /* Compute the fluxes of the conservatives *
         * --------------------------------------- */
        // Auxiliary: unit vector along the current direction
        const vec<CCTK_REAL, 3> unit_dir{vec<int, 3>::unit(dir)};

        // Auxiliary: alpha * sqrt(g)
        const CCTK_REAL alp_sqrtg = alp_avg * sqrtg;

        // Auxiliary: B^i / W
        const vec<CCTK_REAL, 2> B_over_w_lorentz_rc_face(
            [&](int f) ARITH_INLINE { return B_rc_face(f) / w_lorentz_rc_face(f); });

        // flux(dens) = sqrt(g) * D * vtilde^i = sqrt(g) * rho * W * vtilde^i */
        const vec<CCTK_REAL, 2> flux_dens(
            [&](int f) ARITH_INLINE { return dens_rc_face(f) * vtilde_rc_face(f); });

        /* flux(mom_j)^i = sqrt(g)*(
         *   S_j * vtilde^i + alpha*((pgas+pmag)*delta^i_j - b_j*B^i/W) )       */
        const vec<vec<CCTK_REAL, 2>, 3> flux_moms([&](int j) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return moms_rc_face(j)(f) * vtilde_rc_face(f) +
                   alp_sqrtg * (ptot_rc_face(f) * unit_dir(j) -
                                blows_rc_face(j)(f) * B_over_w_lorentz_rc_face(f));
          });
        });

        /* flux(tau) = sqrt(g)*(
         *   t*vtilde^i + alpha*((pgas+pmag)*v^i-alpha*b0*B^i/W) )              */
        const vec<CCTK_REAL, 2> flux_tau([&](int f) ARITH_INLINE {
          return tau_rc_face(f) * vtilde_rc_face(f) +
                 alp_sqrtg * (ptot_rc_face(f) * vel_rc_face(f) -
                              alp_b0_rc_face(f) * B_over_w_lorentz_rc_face(f));
        });

        // Electric field E_i = \epsilon_{ijk} Btilde_j * vtilde_k */
        const vec<vec<CCTK_REAL, 2>, 3> Es_rc_face =
            calc_cross_product(Btildes_rc_face, vtildes_rc_face);

        // flux(Btildes) = {{0, Ez, -Ey}, {-Ez, 0, Ex}, {Ey, -Ex, 0}}
        const vec<vec<CCTK_REAL, 2>, 3> flux_Btildes =
            calc_cross_product(unit_dir, Es_rc_face);


        /* Calculate eigenvalues *
         * --------------------- */
        // Alias for g^xx, g^yy or g^zz depending on the direction
        const CCTK_REAL u_avg = calc_inv(g_avg, detg_avg)(dir, dir);

        // Eigenvalues
        vec<vec<CCTK_REAL, 4>, 2> lambda =
            eigenvalues(alp_avg, beta_avg, u_avg,
                        vel_rc_face, rho_rc_face, cs2_rc_face,
                        w_lorentz_rc_face, h_rc_face, b2_rc_face);

        // Numerical fluxes
        fluxdenss(dir)(p.I) = calcflux(lambda, dens_rc_face,    flux_dens);
        fluxmomxs(dir)(p.I) = calcflux(lambda, moms_rc_face(0), flux_moms(0));
        fluxmomys(dir)(p.I) = calcflux(lambda, moms_rc_face(1), flux_moms(1));
        fluxmomzs(dir)(p.I) = calcflux(lambda, moms_rc_face(2), flux_moms(2));
        fluxtaus(dir)(p.I)  = calcflux(lambda, tau_rc_face,     flux_tau);
        fluxBxs(dir)(p.I)   = (dir != 0) ? calcflux(lambda, Btildes_rc_face(0), flux_Btildes(0)) : 0.0;
        fluxBys(dir)(p.I)   = (dir != 1) ? calcflux(lambda, Btildes_rc_face(1), flux_Btildes(1)) : 0.0;
        fluxBzs(dir)(p.I)   = (dir != 2) ? calcflux(lambda, Btildes_rc_face(2), flux_Btildes(2)) : 0.0;
      });
}





void CalcAuxForAvecPsi(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  const vec<GF3D2<const CCTK_REAL>, dim> gf_Avecs{Avec_x, Avec_y, Avec_z};
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* interpolate A to vertices */
        const vec<CCTK_REAL, 3> A_vert([&](int i) ARITH_INLINE {
          return calc_avg_e2v(gf_Avecs(i), p, i);
        });
        const smat<CCTK_REAL, 3> g{gxx(p.I), gxy(p.I), gxz(p.I),
                                   gyy(p.I), gyz(p.I), gzz(p.I)};
        const vec<CCTK_REAL, 3> betas{betax(p.I), betay(p.I), betaz(p.I)};
        const CCTK_REAL detg = calc_det(g);
        const CCTK_REAL sqrtg = sqrt(detg);
        const smat<CCTK_REAL, 3> ug = calc_inv(g, detg);
        const vec<CCTK_REAL, 3> Aup = calc_contraction(ug, A_vert);

        Fx(p.I) = alp(p.I) * sqrtg * Aup(0);
        Fy(p.I) = alp(p.I) * sqrtg * Aup(1);
        Fz(p.I) = alp(p.I) * sqrtg * Aup(2);
        Fbetax(p.I) = betas(0) * Psi(p.I);
        Fbetay(p.I) = betas(1) * Psi(p.I);
        Fbetaz(p.I) = betas(2) * Psi(p.I);
        G(p.I) = alp(p.I) * Psi(p.I) / sqrtg - calc_contraction(betas, A_vert);
      });
}





extern "C" void AsterX_Fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  CalcFlux<0>(cctkGH);
  CalcFlux<1>(cctkGH);
  CalcFlux<2>(cctkGH);

  /* Set auxiliary variables for the rhs of A and Psi  */
  CalcAuxForAvecPsi(cctkGH);
}

} // namespace AsterX
