//#include <loop.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

namespace HydroToyCarpetX {
using namespace std;

constexpr int dim = 3;

extern "C" void HydroToyCarpetX_Evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_Evolve;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  // Transport
  // dt rho + d_i (rho vel^i) = 0
  // dt mom_j + d_i (mom_j vel^i) = 0
  // dt etot + d_i (etot vel^i) = 0

  const CCTK_REAL dt_dx = dt / dx;
  const CCTK_REAL dt_dy = dt / dy;
  const CCTK_REAL dt_dz = dt / dz;

#if 0
    // CPU

    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      auto calcupdate{[&](auto &fx, auto &fy, auto &fz) {
        return dt_dx * (fx(p.I + p.DI(0)) - fx(p.I)) +
               dt_dy * (fy(p.I + p.DI(1)) - fy(p.I)) +
               dt_dz * (fz(p.I + p.DI(2)) - fz(p.I));
      }};

      rho_(p.I) = rho_p_(p.I) - calcupdate(fxrho_, fyrho_, fzrho_);

      // momx_(p.I) = momx_p_(p.I) - calcupdate(fxmomx_, fymomx_, fzmomx_);
      // momy_(p.I) = momy_p_(p.I) - calcupdate(fxmomy_, fymomy_, fzmomy_);
      // momz_(p.I) = momz_p_(p.I) - calcupdate(fxmomz_, fymomz_, fzmomz_);

      // etot_(p.I) = etot_p_(p.I) - calcupdate(fxetot_, fyetot_, fzetot_);
    });
#else
   // CPU or GPU
   // TODO: &p? 
  const auto calcupdate =
      [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE(
          CCTK_REAL fx_p, CCTK_REAL fy_p, CCTK_REAL fz_p, CCTK_REAL fx_m, CCTK_REAL fy_m, CCTK_REAL fz_m) {  
	  return dt_dx * (fx_p - fx_m) +
             dt_dy * (fy_p - fy_m) +
             dt_dz * (fz_p - fz_m);
      };

  printf("Before loop\n");

  // Determine loop extent
  const array<int, dim> nghostzones{cctkGH->cctk_nghostzones[0],
                                    cctkGH->cctk_nghostzones[1],
                                    cctkGH->cctk_nghostzones[2]};
  const Loop::GridDescBaseDevice griddesc(cctkGH);
  griddesc.loop_int_device<1, 1, 1>(
      nghostzones, [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE(
                       const Loop::PointDesc &p) {

        const auto px = griddesc.point_desc<0, 1, 1>(p);
        const auto py = griddesc.point_desc<1, 0, 1>(p);
        const auto pz = griddesc.point_desc<1, 1, 0>(p);

        rho[p.idx] = rho_p[p.idx];
        momx[p.idx] = momx_p[p.idx];
        momy[p.idx] = momy_p[p.idx];
        momz[p.idx] = momz_p[p.idx];
        etot[p.idx] = etot_p[p.idx];
      
//        printf("p.idx:= %d\n", p.idx);
//        printf("px.idx:= %d\n", px.idx);

//        if (p.i == 3 && p.j == 4 && p.k == 5) {
       
        rho[p.idx] = rho_p[p.idx] - calcupdate(fxrho[px.idx+px.di], fyrho[py.idx+py.dj], fzrho[pz.idx+pz.dk], fxrho[px.idx], fyrho[py.idx], fzrho[pz.idx]);

        printf("p.idx = %d, px.idx: = %d, p.di=%d, px.di = %d\n", p.idx, px.idx, p.di, px.di);
        printf("fxrho_p = %g, fxrho_m = %g\n", fxrho[px.idx+px.di], fxrho[px.idx]);
//        printf("rho = %g\n", rho[p.idx]);
//        if (p.idx == 972 || p.idx==971) {
//          printf("p: x = %g, y = %g, z = %g\n", p.x, p.y, p.z);
//          printf("px: x = %g, y = %g, z = %g\n", px.x, px.y, px.z);
//          printf("py: x = %g, y = %g, z = %g\n", py.x, py.y, py.z);
//          printf("pz: x = %g, y = %g, z = %g\n", pz.x, pz.y, pz.z);
//          printf("p.idx = %g, px.idx: = %d, px.di=%d, px.dx = %d\n", p.idx, px.idx, px.dx,p.di);
//          printf("fxrho_p = %g, fxrho_m = %g\n", fxrho[px.idx+px.di], fxrho[px.idx]);
//          printf("rho = %g\n", rho[p.idx]);
//        }
 
//        rho[p.idx] = 1.0;
//        momx[p.idx] = momx_p[p.idx] - calcupdate(fxmomx[px.idx+px.di], fymomx[py.idx+py.dj], fzmomx[pz.idx+pz.dk], fxmomx[px.idx], fymomx[py.idx], fzmomx[pz.idx]);
//        momy[p.idx] = momy_p[p.idx] - calcupdate(fxmomy[px.idx+px.di], fymomy[py.idx+py.dj], fzmomy[pz.idx+pz.dk], fxmomy[px.idx], fymomy[py.idx], fzmomy[pz.idx]);
//        momz[p.idx] = momz_p[p.idx] - calcupdate(fxmomz[px.idx+px.di], fymomz[py.idx+py.dj], fzmomz[pz.idx+pz.dk], fxmomz[px.idx], fymomz[py.idx], fzmomz[pz.idx]);
//        etot[p.idx] = etot_p[p.idx] - calcupdate(fxetot[px.idx+px.di], fyetot[py.idx+py.dj], fzetot[pz.idx+pz.dk], fxetot[px.idx], fyetot[py.idx], fzetot[pz.idx]);

//        psi[p.idx] =
//            psi_p[p.idx] +
//            dt * (ddx_phi + ddy_phi + ddz_phi - pow(mass, 2) * phi_p[p.idx] +
//                  4 * M_PI * central_potential(t, p.x, p.y, p.z));
//        phi[p.idx] = phi_p[p.idx] + dt * psi[p.idx];
      });
#endif

      printf("After loop\n");
   
}

} // namespace HydroToyCarpetX
