#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include "utils.hxx"
#include "fluxes.hxx"
#include "reconstruct.hxx"
#include <eos.hxx>
#include <eos_idealgas.hxx>

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;
using namespace EOSX;

enum class flux_t { LxF, HLLE };
enum class eos_t { IdealGas, Hybrid, Tabulated };

// Calculate the fluxes in direction `dir`. This function is more
// complex because it has to handle any direction, but as reward,
// there is only one function, not three.
template <int dir, typename EOSType>
void CalcFlux(CCTK_ARGUMENTS, EOSType &eos_th) {
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
  const auto calcflux =
      [=] CCTK_DEVICE(vec<vec<CCTK_REAL, 4>, 2> lam, vec<CCTK_REAL, 2> var,
                      vec<CCTK_REAL, 2> flux) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL flx;
        switch (fluxtype) {

        case flux_t::LxF: {
          flx = laxf(lam, var, flux);
          break;
        }

        case flux_t::HLLE: {
          flx = hlle(lam, var, flux);
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

  // Face-centred grid functions (in direction `dir`)
  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};

  grid.loop_int_device<face_centred[0], face_centred[1], face_centred[2]>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* Reconstruct primitives from the cells on left (indice 0) and right
         * (indice 1) side of this face rc = reconstructed variables or
         * computed from reconstructed variables */
        const vec<CCTK_REAL, 2> rho_rc{reconstruct_pt(rho, p)};
        const vec<vec<CCTK_REAL, 2>, 3> vels_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>{reconstruct_pt(gf_vels(i), p)};
        });
        const vec<CCTK_REAL, 2> eps_rc{reconstruct_pt(eps, p)};
        const vec<vec<CCTK_REAL, 2>, 3> Bs_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>{reconstruct_pt(gf_Bvecs(i), p)};
        });

        /* Interpolate metric components from vertices to faces */
        const CCTK_REAL alp_avg = calc_avg_v2f(alp, p, dir);
        const vec<CCTK_REAL, 3> betas_avg([&](int i) ARITH_INLINE {
          return calc_avg_v2f(gf_beta(i), p, dir);
        });
        const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2f(gf_g(i, j), p, dir);
        });

        /* determinant of spatial metric */
        const CCTK_REAL detg_avg = calc_det(g_avg);
        const CCTK_REAL sqrtg = sqrt(detg_avg);
        /* co-velocity measured by Euleian observer: v_j */
        const vec<vec<CCTK_REAL, 2>, 3> vlows_rc =
            calc_contraction(g_avg, vels_rc);
        /* vtilde^i = alpha * v^i - beta^i */
        const vec<vec<CCTK_REAL, 2>, 3> vtildes_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return alp_avg * vels_rc(i)(f) - betas_avg(i);
          });
        });
        /* Lorentz factor: W = 1 / sqrt(1 - v^2) */
        const vec<CCTK_REAL, 2> w_lorentz_rc([&](int f) ARITH_INLINE {
          return 1 / sqrt(1 - calc_contraction(vlows_rc, vels_rc)(f));
        });

        /* alpha * b0 = W * B^i * v_i */
        const vec<CCTK_REAL, 2> alp_b0_rc([&](int f) ARITH_INLINE {
          return w_lorentz_rc(f) * calc_contraction(Bs_rc, vlows_rc)(f);
        });
        /* covariant magnetic field measured by the Eulerian observer */
        const vec<vec<CCTK_REAL, 2>, 3> Blows_rc =
            calc_contraction(g_avg, Bs_rc);
        /* B^2 = B^i * B_i */
        const vec<CCTK_REAL, 2> B2_rc = calc_contraction(Bs_rc, Blows_rc);
        /* covariant magnetic field measured by the comoving observer:
         *  b_i = B_i/W + alpha*b^0*v_i */
        const vec<vec<CCTK_REAL, 2>, 3> blows_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return Blows_rc(i)(f) / w_lorentz_rc(f) +
                   alp_b0_rc(f) * vlows_rc(i)(f);
          });
        });
        /* b^2 = b^{\mu} * b_{\mu} */
        const vec<CCTK_REAL, 2> bsq_rc([&](int f) ARITH_INLINE {
          return (B2_rc(f) + pow2(alp_b0_rc(f))) / pow2(w_lorentz_rc(f));
        });

        /* componets correspond to the dir we are considering */
        const CCTK_REAL beta_avg = betas_avg(dir);
        const vec<CCTK_REAL, 2> vel_rc{vels_rc(dir)};
        const vec<CCTK_REAL, 2> B_rc{Bs_rc(dir)};
        const vec<CCTK_REAL, 2> vtilde_rc{vtildes_rc(dir)};

        // TODO: Compute pressure based on user-specified EOS.
        // Currently, computing press for classical ideal gas from reconstructed
        // vars

        // TODO: Correctly reconstruct Ye
        const vec<CCTK_REAL, 2> ye_rc{ye_min, ye_max};

        // Ideal gas case {
        /* pressure for ideal gas EOS */
        const vec<CCTK_REAL, 2> press_rc([&](int f) ARITH_INLINE {
          return eos_th.press_from_valid_rho_eps_ye(rho_rc(f), eps_rc(f),
                                                    ye_rc(f));
        });
        /* cs2 for ideal gas EOS */
        const vec<CCTK_REAL, 2> cs2_rc([&](int f) ARITH_INLINE {
          return eos_th.csnd_from_valid_rho_eps_ye(rho_rc(f), eps_rc(f),
                                                   ye_rc(f)) *
                 eos_th.csnd_from_valid_rho_eps_ye(rho_rc(f), eps_rc(f),
                                                   ye_rc(f));
        });
        /* enthalpy h for ideal gas EOS */
        const vec<CCTK_REAL, 2> h_rc([&](int f) ARITH_INLINE {
          return 1 + eps_rc(f) + press_rc(f) / rho_rc(f);
        });
        // } Ideal gas case

        /* Computing conservatives from primitives: */

        /* dens = sqrt(g) * D = sqrt(g) * (rho * W) */
        const vec<CCTK_REAL, 2> dens_rc([&](int f) ARITH_INLINE {
          return sqrtg * rho_rc(f) * w_lorentz_rc(f);
        });

        /* auxiliary: dens * h * W = sqrt(g) * rho * h * W^2 */
        const vec<CCTK_REAL, 2> dens_h_W_rc([&](int f) ARITH_INLINE {
          return dens_rc(f) * h_rc(f) * w_lorentz_rc(f);
        });
        /* auxiliary: sqrt(g) * (rho*h + b^2)*W^2 */
        const vec<CCTK_REAL, 2> dens_h_W_plus_sqrtg_W2b2_rc =
            dens_h_W_rc + sqrtg * (pow2(alp_b0_rc) + B2_rc);
        /* auxiliary: (pgas + pmag) */
        const vec<CCTK_REAL, 2> press_plus_pmag_rc = press_rc + 0.5 * bsq_rc;

        /* mom_i = sqrt(g)*S_i = sqrt(g)((rho*h+b^2)*W^2*v_i - alpha*b^0*b_i) */
        const vec<vec<CCTK_REAL, 2>, 3> moms_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return dens_h_W_plus_sqrtg_W2b2_rc(f) * vlows_rc(i)(f) -
                   sqrtg * alp_b0_rc(f) * blows_rc(i)(f);
          });
        });

        /* tau = sqrt(g)*t =
         *  sqrt(g)((rho*h + b^2)*W^2 - (pgas+pmag) - (alpha*b^0)^2 - D) */
        const vec<CCTK_REAL, 2> tau_rc =
            dens_h_W_rc - dens_rc + sqrtg * (B2_rc - press_plus_pmag_rc);

        /* Btildes^i = sqrt(g) * B^i */
        const vec<vec<CCTK_REAL, 2>, 3> Btildes_rc(
            [&](int i) ARITH_INLINE { return sqrtg * Bs_rc(i); });

        /* Computing fluxes of conserved variables: */

        /* auxiliary: unit in 'dir' */
        const vec<CCTK_REAL, 3> unit_dir{vec<int, 3>::unit(dir)};
        /* auxiliary: alpha * sqrt(g) */
        const CCTK_REAL alp_sqrtg = alp_avg * sqrtg;
        /* auxiliary: B^i / W */
        const vec<CCTK_REAL, 2> B_over_w_lorentz_rc(
            [&](int f) ARITH_INLINE { return B_rc(f) / w_lorentz_rc(f); });

        /* flux(dens) = sqrt(g) * D * vtilde^i = sqrt(g) * rho * W * vtilde^i */
        const vec<CCTK_REAL, 2> flux_dens(
            [&](int f) ARITH_INLINE { return dens_rc(f) * vtilde_rc(f); });

        /* flux(mom_j)^i = sqrt(g)*(
         *  S_j*vtilde^i + alpha*((pgas+pmag)*delta^i_j - b_jB^i/W) ) */
        const vec<vec<CCTK_REAL, 2>, 3> flux_moms([&](int j) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return moms_rc(j)(f) * vtilde_rc(f) +
                   alp_sqrtg * (press_plus_pmag_rc(f) * unit_dir(j) -
                                blows_rc(j)(f) * B_over_w_lorentz_rc(f));
          });
        });

        /* flux(tau) = sqrt(g)*(
         *  t*vtilde^i + alpha*((pgas+pmag)*v^i-alpha*b0*B^i/W) ) */
        const vec<CCTK_REAL, 2> flux_tau([&](int f) ARITH_INLINE {
          return tau_rc(f) * vtilde_rc(f) +
                 alp_sqrtg * (press_plus_pmag_rc(f) * vel_rc(f) -
                              alp_b0_rc(f) * B_over_w_lorentz_rc(f));
        });

        /* electric field E_i = \tilde\epsilon_{ijk} Btilde_j * vtilde_k */
        const vec<vec<CCTK_REAL, 2>, 3> Es_rc =
            calc_cross_product(Btildes_rc, vtildes_rc);
        /* flux(Btildes) = {{0, Ez, -Ey}, {-Ez, 0, Ex}, {Ey, -Ex, 0}} */
        const vec<vec<CCTK_REAL, 2>, 3> flux_Btildes =
            calc_cross_product(unit_dir, Es_rc);

        /* Calculate eigenvalues: */

        /* variable for either g^xx, g^yy or g^zz depending on the direction */
        const CCTK_REAL u_avg = calc_inv(g_avg, detg_avg)(dir, dir);
        /* eigenvalues */
        vec<vec<CCTK_REAL, 4>, 2> lambda =
            eigenvalues(alp_avg, beta_avg, u_avg, vel_rc, rho_rc, cs2_rc,
                        w_lorentz_rc, h_rc, bsq_rc);

        /* Calculate numerical fluxes */
        fluxdenss(dir)(p.I) = calcflux(lambda, dens_rc, flux_dens);
        fluxmomxs(dir)(p.I) = calcflux(lambda, moms_rc(0), flux_moms(0));
        fluxmomys(dir)(p.I) = calcflux(lambda, moms_rc(1), flux_moms(1));
        fluxmomzs(dir)(p.I) = calcflux(lambda, moms_rc(2), flux_moms(2));
        fluxtaus(dir)(p.I) = calcflux(lambda, tau_rc, flux_tau);
        fluxBxs(dir)(p.I) =
            (dir != 0) * calcflux(lambda, Btildes_rc(0), flux_Btildes(0));
        fluxBys(dir)(p.I) =
            (dir != 1) * calcflux(lambda, Btildes_rc(1), flux_Btildes(1));
        fluxBzs(dir)(p.I) =
            (dir != 2) * calcflux(lambda, Btildes_rc(2), flux_Btildes(2));
        if (isnan(flux_dens(0)) || isnan(flux_dens(1)) || isnan(dens_rc(0)) ||
            isnan(dens_rc(1)) || isnan(fluxdenss(dir)(p.I)) ||
            isnan(fluxdenss(dir)(p.I))) {
          printf("dens_rc = %f, %f, vtilde_rc = %f, %f, flux_den = %f, %f",
                 dens_rc(0), dens_rc(1), vtilde_rc(0), vtilde_rc(1),
                 flux_dens(0), flux_dens(1));
        }
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

  eos_t eostype;
  eos::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
      rgye(ye_min, ye_max);

  if (CCTK_EQUALS(evolution_eos, "IdealGas")) {
    eostype = eos_t::IdealGas;
  } else if (CCTK_EQUALS(evolution_eos, "Hybrid")) {
    eostype = eos_t::Hybrid;
  } else if (CCTK_EQUALS(evolution_eos, "Tabulated")) {
    eostype = eos_t::Tabulated;
  } else {
    CCTK_ERROR("Unknown value for parameter \"evolution_eos\"");
  }

  switch (eostype) {
  case eos_t::IdealGas: {
    eos_idealgas eos_th(gl_gamma, particle_mass, rgeps, rgrho, rgye);
    CalcFlux<0>(cctkGH, eos_th);
    CalcFlux<1>(cctkGH, eos_th);
    CalcFlux<2>(cctkGH, eos_th);
    break;
  }
  case eos_t::Hybrid: {
    CCTK_ERROR("Hybrid EOS is not yet supported");
    break;
  }
  case eos_t::Tabulated: {
    CCTK_ERROR("Tabulated EOS is not yet supported");
    break;
  }
  default:
    assert(0);
  }

  /* Set auxiliary variables for the rhs of A and Psi  */
  CalcAuxForAvecPsi(cctkGH);
}

} // namespace AsterX
