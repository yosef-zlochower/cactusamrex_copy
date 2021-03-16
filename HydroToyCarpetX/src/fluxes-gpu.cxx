//#include <loop.hxx>

//#include <vectors.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>

namespace HydroToyCarpetX {
using namespace std;

constexpr int dim = 3;

extern "C" void HydroToyCarpetX_Fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dAx = dt * dy * dz;
  const CCTK_REAL dAy = dt * dx * dz;
  const CCTK_REAL dAz = dt * dx * dy;

  // frho^i = rho vel^i
  // fmom^i_j = mom_j vel^i + delta^i_j press
  // fetot^i = (etot + press) vel^i

  const auto calcflux =
          [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE(
              CCTK_REAL var_p, CCTK_REAL var_m, CCTK_REAL flux_p, CCTK_REAL flux_m) {

            CCTK_REAL lambda_m = 1.0;
            CCTK_REAL lambda_p = -1.0;
//            CCTK_REAL var_m = u[p.idx-p.di];
//            CCTK_REAL var_p = u[p.idx];
//            CCTK_REAL flux_m = f[p1.idx-p1.di];
//            CCTK_REAL flux_p = f[p1.idx];
            CCTK_REAL llf = 0.5 *((flux_m + flux_p) - fmax(fabs(lambda_m), fabs(lambda_p)) * (var_p - var_m));
//            printf("llf = %g, dAx = %g\n",llf,dAx);
            return dAx * llf;
          };

//  printf("Before loop\n");

  // Determine loop extent
  const array<int, dim> nghostzones{cctkGH->cctk_nghostzones[0],
                                    cctkGH->cctk_nghostzones[1],
                                    cctkGH->cctk_nghostzones[2]};

  const Loop::GridDescBaseDevice griddesc(cctkGH);
  printf("nghost = %g %g %g\n", nghostzones[0],nghostzones[1],nghostzones[2]);
  griddesc.loop_all_device<0, 1, 1>(
      nghostzones, [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE(
                       const Loop::PointDesc &p) {
      const auto p1 = griddesc.point_desc<1, 1, 1>(p);

      if (p.i == 3 && p.j == 3 && p.k == 3) {
        printf("p: x = %g, y = %g, z = %g\n", p.x, p.y, p.z);
        printf("p1: x = %g, y = %g, z = %g\n", p1.x, p1.y, p1.z);
        printf("p: i = %d, j = %d, k = %d\n", p.i, p.j, p.k);
        printf("p1: i = %d, j = %d, k = %d\n", p1.i, p1.j, p1.k);
        printf("rho_p = %g, rho_m = %g\n", rho[p1.idx], rho[p1.idx-p1.di]);
        printf("momx_p = %g, momx_m = %g\n", momx[p1.idx], momx[p1.idx-p1.di]);
        printf("momy_p = %g, momy_m = %g\n", momy[p1.idx], momy[p1.idx-p1.di]);
        printf("momz_p = %g, momz_m = %g\n", momz[p1.idx], momz[p1.idx-p1.di]);
        printf("etot_p = %g, etot_m = %g\n", etot[p1.idx], etot[p1.idx-p1.di]);
        printf("fxrho_p = %f, fxrho_m =%f\n",(rho[p1.idx] * velx[p1.idx]),(rho[p1.idx-p1.di] * velx[p1.idx-p1.di]));
        printf("fxmomx_p = %f, fxmomx_m =%f\n",(momx[p1.idx] * velx[p1.idx] + press[p1.idx],momx[p1.idx-p1.di] * velx[p1.idx-p1.di] + press[p1.idx-p1.di]));
        printf("fxmomy_p = %f, fxmomy_m =%f\n",(momy[p1.idx] * velx[p1.idx] + press[p1.idx],momx[p1.idx-p1.di] * velx[p1.idx-p1.di] + press[p1.idx-p1.di]));
        printf("fxmomz_p = %f, fxmomz_m =%f\n",(momz[p1.idx] * velx[p1.idx] + press[p1.idx],momx[p1.idx-p1.di] * velx[p1.idx-p1.di] + press[p1.idx-p1.di]));
        printf("fxetot_p = %f, fxetot_m =%f\n",(etot[p1.idx] + press[p1.idx]) * velx[p1.idx], (etot[p1.idx-p1.di] + press[p1.idx-p1.di]) * velx[p1.idx-p1.di]);
      }

      fxrho[p.idx] = calcflux(rho[p1.idx], rho[p1.idx-p1.di], rho[p1.idx] * velx[p1.idx], rho[p1.idx-p1.di] * velx[p1.idx-p1.di]);
      fxmomx[p.idx] = calcflux(momx[p1.idx], momx[p1.idx-p1.di],momx[p1.idx] * velx[p1.idx] + press[p1.idx],momx[p1.idx-p1.di] * velx[p1.idx-p1.di] + press[p1.idx-p1.di]);
      fxmomy[p.idx] = calcflux(momy[p1.idx], momy[p1.idx-p1.di],momy[p1.idx] * velx[p1.idx] + press[p1.idx],momy[p1.idx-p1.di] * velx[p1.idx-p1.di] + press[p1.idx-p1.di]);
      fxmomz[p.idx] = calcflux(momz[p1.idx], momz[p1.idx-p1.di],momz[p1.idx] * velx[p1.idx] + press[p1.idx],momz[p1.idx-p1.di] * velx[p1.idx-p1.di] + press[p1.idx-p1.di]);
      fxetot[p.idx] = calcflux(etot[p1.idx], etot[p1.idx-p1.di], (etot[p1.idx] + press[p1.idx]) * velx[p1.idx], (etot[p1.idx-p1.di] + press[p1.idx-p1.di]) * velx[p1.idx-p1.di]);
//      fxrho[p.idx] = 0.0;

      if (p.i == 3 && p.j == 3 && p.k == 3) {
        printf("fxrho = %f\n",fxrho[p.idx]);
        printf("fxmomx = %f\n",fxmomx[p.idx]);
        printf("fxmomy = %f\n",fxmomy[p.idx]);
        printf("fxmomz = %f\n",fxmomz[p.idx]);
        printf("fxetot = %f\n",fxetot[p.idx]);
      }

  });

  griddesc.loop_all_device<1, 0, 1>(
      nghostzones, [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE(
                       const Loop::PointDesc &p) {
      const auto p1 = griddesc.point_desc<1, 1, 1>(p);
      if (p1.i == 3 && p1.j == 4 && p1.k == 5) {
        printf("p.idx = %d, p1.idx: = %d, p.di=%d, p1.di = %d\n", p.idx, p1.idx, p.di, p1.di);
        printf("p: x = %g, y = %g, z = %g\n", p.x, p.y, p.z);
        printf("p1: x = %g, y = %g, z = %g\n", p1.x, p1.y, p1.z);
        printf("p: i = %d, j = %d, k = %d\n", p.i, p.j, p.k);
        printf("p1: i = %d, j = %d, k = %d\n", p1.i, p1.j, p1.k);
        printf("rho_p = %g, rho_m = %g\n", rho[p1.idx], rho[p1.idx-p1.dj]);
        printf("momx_p = %g, momx_m = %g\n", momx[p1.idx], momx[p1.idx-p1.di]);
        printf("momy_p = %g, momy_m = %g\n", momy[p1.idx], momy[p1.idx-p1.di]);
        printf("momz_p = %g, momz_m = %g\n", momz[p1.idx], momz[p1.idx-p1.di]);
        printf("etot_p = %g, etot_m = %g\n", etot[p1.idx], etot[p1.idx-p1.di]);
        printf("fyrho_p = %f, fyrho_m =%f\n",(rho[p1.idx] * vely[p1.idx]),(rho[p1.idx-p1.dj] * vely[p1.idx-p1.dj]));
        printf("fymomx_p = %f, fymomx_m =%f\n",(momx[p1.idx] * vely[p1.idx] + press[p1.idx],momx[p1.idx-p1.dj] * vely[p1.idx-p1.dj] + press[p1.idx-p1.dj]));
        printf("fymomy_p = %f, fymomy_m =%f\n",(momy[p1.idx] * vely[p1.idx] + press[p1.idx],momx[p1.idx-p1.dj] * vely[p1.idx-p1.dj] + press[p1.idx-p1.dj]));
        printf("fymomz_p = %f, fymomz_m =%f\n",(momz[p1.idx] * vely[p1.idx] + press[p1.idx],momx[p1.idx-p1.dj] * vely[p1.idx-p1.dj] + press[p1.idx-p1.dj]));
        printf("fyetot_p = %f, fyetot_m =%f\n",(etot[p1.idx] + press[p1.idx]) * vely[p1.idx], (etot[p1.idx-p1.dj] + press[p1.idx-p1.dj]) * vely[p1.idx-p1.dj]);
      }

      fyrho[p.idx] = calcflux(rho[p1.idx], rho[p1.idx-p1.dj], rho[p1.idx] * vely[p1.idx], rho[p1.idx-p1.dj] * vely[p1.idx-p1.dj]);
      fymomx[p.idx] = calcflux(momx[p1.idx], momx[p1.idx-p1.dj],momx[p1.idx] * vely[p1.idx] + press[p1.idx],momx[p1.idx-p1.dj] * vely[p1.idx-p1.dj] + press[p1.idx-p1.dj]);
      fymomy[p.idx] = calcflux(momy[p1.idx], momy[p1.idx-p1.dj],momy[p1.idx] * vely[p1.idx] + press[p1.idx],momy[p1.idx-p1.dj] * vely[p1.idx-p1.dj] + press[p1.idx-p1.dj]);
      fymomz[p.idx] = calcflux(momz[p1.idx], momz[p1.idx-p1.dj],momz[p1.idx] * vely[p1.idx] + press[p1.idx],momz[p1.idx-p1.dj] * vely[p1.idx-p1.dj] + press[p1.idx-p1.dj]);
      fyetot[p.idx] = calcflux(etot[p1.idx], etot[p1.idx-p1.dj], (etot[p1.idx] + press[p1.idx]) * vely[p1.idx], (etot[p1.idx-p1.dj] + press[p1.idx-p1.dj]) * vely[p1.idx-p1.dj]);

      if (p1.i == 3 && p1.j == 4 && p1.k == 5) {
        printf("fyrho = %f\n",fyrho[p.idx]);
        printf("fymomx = %f\n",fymomx[p.idx]);
        printf("fymomy = %f\n",fymomy[p.idx]);
        printf("fymomz = %f\n",fymomz[p.idx]);
        printf("fyetot = %f\n",fyetot[p.idx]);
      }

  });

  griddesc.loop_all_device<1, 1, 0>(
      nghostzones, [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE(
                       const Loop::PointDesc &p) {
      const auto p1 = griddesc.point_desc<1, 1, 1>(p);
      if (p1.i == 3 && p1.j == 4 && p1.k == 5) {
        printf("p: x = %g, y = %g, z = %g\n", p.x, p.y, p.z);
        printf("p1: x = %g, y = %g, z = %g\n", p1.x, p1.y, p1.z);
        printf("p: i = %d, j = %d, k = %d\n", p.i, p.j, p.k);
        printf("p1: i = %d, j = %d, k = %d\n", p1.i, p1.j, p1.k);
        printf("rho_p = %g, rho_m = %g\n", rho[p1.idx], rho[p1.idx-p1.dk]);
        printf("momx_p = %g, momx_m = %g\n", momx[p1.idx], momx[p1.idx-p1.dk]);
        printf("momy_p = %g, momy_m = %g\n", momy[p1.idx], momy[p1.idx-p1.dk]);
        printf("momz_p = %g, momz_m = %g\n", momz[p1.idx], momz[p1.idx-p1.dk]);
        printf("etot_p = %g, etot_m = %g\n", etot[p1.idx], etot[p1.idx-p1.dk]);
        printf("fzrho_p = %f, fzrho_m =%f\n",(rho[p1.idx] * velz[p1.idx]),(rho[p1.idx-p1.dk] * velz[p1.idx-p1.dk]));
        printf("fzmomx_p = %f, fzmomx_m =%f\n",(momx[p1.idx] * velz[p1.idx] + press[p1.idx],momx[p1.idx-p1.dk] * velz[p1.idx-p1.dk] + press[p1.idx-p1.dk]));
        printf("fzmomy_p = %f, fzmomy_m =%f\n",(momy[p1.idx] * velz[p1.idx] + press[p1.idx],momx[p1.idx-p1.dk] * velz[p1.idx-p1.dk] + press[p1.idx-p1.dk]));
        printf("fzmomz_p = %f, fzmomz_m =%f\n",(momz[p1.idx] * velz[p1.idx] + press[p1.idx],momx[p1.idx-p1.dk] * velz[p1.idx-p1.dk] + press[p1.idx-p1.dk]));
        printf("fzetot_p = %f, fzetot_m =%f\n",(etot[p1.idx] + press[p1.idx]) * velz[p1.idx], (etot[p1.idx-p1.dk] + press[p1.idx-p1.dk]) * velz[p1.idx-p1.dk]);
      }

      fzrho[p.idx] = calcflux(rho[p1.idx], rho[p1.idx-p1.dk], rho[p1.idx] * velz[p1.idx], rho[p1.idx-p1.dk] * velz[p1.idx-p1.dk]);
      fzmomx[p.idx] = calcflux(momx[p1.idx], momx[p1.idx-p1.dk],momx[p1.idx] * velz[p1.idx] + press[p1.idx],momx[p1.idx-p1.dk] * velz[p1.idx-p1.dk] + press[p1.idx-p1.dk]);
      fzmomy[p.idx] = calcflux(momy[p1.idx], momy[p1.idx-p1.dk],momy[p1.idx] * velz[p1.idx] + press[p1.idx],momy[p1.idx-p1.dk] * velz[p1.idx-p1.dk] + press[p1.idx-p1.dk]);
      fzmomz[p.idx] = calcflux(momz[p1.idx], momz[p1.idx-p1.dk],momz[p1.idx] * velz[p1.idx] + press[p1.idx],momz[p1.idx-p1.dk] * velz[p1.idx-p1.dk] + press[p1.idx-p1.dk]);
      fzetot[p.idx] = calcflux(etot[p1.idx], etot[p1.idx-p1.dk], (etot[p1.idx] + press[p1.idx]) * velz[p1.idx], (etot[p1.idx-p1.dk] + press[p1.idx-p1.dk]) * velz[p1.idx-p1.dk]);

      if (p1.i == 3 && p1.j == 4 && p1.k == 5) {
        printf("fzrho = %f\n",fzrho[p.idx]);
        printf("fzmomx = %f\n",fzmomx[p.idx]);
        printf("fzmomy = %f\n",fzmomy[p.idx]);
        printf("fzmomz = %f\n",fzmomz[p.idx]);
        printf("fzetot = %f\n",fzetot[p.idx]);
      }

  });
 }
} // namespace HydroToyCarpetX
