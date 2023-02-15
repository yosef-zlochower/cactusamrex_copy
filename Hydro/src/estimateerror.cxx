#include "defs.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>

namespace Hydro {
using namespace std;

////////////////////////////////////////////////////////////////////////////////

extern "C" void Hydro_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Hydro_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  const Loop::vect<int, dim> lsh{{cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]}};
  const Loop::vect<int, dim> nghostzones{
      {cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]}};

  array<CCTK_INT, dim> tmin_, tmax_;
  GetTileExtent(cctkGH, tmin_.data(), tmax_.data());
  const Loop::vect<int, dim> tmin(tmin_);
  const Loop::vect<int, dim> tmax(tmax_);

  const Loop::vect<int, dim> imin = max(tmin, nghostzones);
  const Loop::vect<int, dim> imax = min(tmax, lsh - nghostzones);

  const Loop::vect<int, dim> ash{{cctk_ash[0], cctk_ash[1], cctk_ash[2]}};
  constexpr ptrdiff_t di = 1;
  const ptrdiff_t dj = di * cctk_ash[0];
  const ptrdiff_t dk = dj * cctk_ash[1];

  for (int k = imin[2]; k < imax[2]; ++k) {
    for (int j = imin[1]; j < imax[1]; ++j) {
      for (int i = imin[0]; i < imax[0]; i += vsize) {
        CCTK_BOOLVEC mask = mask_for_loop_tail<CCTK_BOOLVEC>(i, imax[0]);
        ptrdiff_t ind = i + dj * j + dk * k;

        auto calcerr = [&](const CCTK_REAL *restrict var) {
          CCTK_REALVEC varxx = maskz_loadu(mask, &var[ind - 1]) -
                               CCTK_REAL(2.0) * maskz_loadu(mask, &var[ind]) +
                               maskz_loadu(mask, &var[ind + 1]);
          CCTK_REALVEC varyy = maskz_loadu(mask, &var[ind - dj]) -
                               CCTK_REAL(2.0) * maskz_loadu(mask, &var[ind]) +
                               maskz_loadu(mask, &var[ind + dj]);
          CCTK_REALVEC varzz = maskz_loadu(mask, &var[ind - dk]) -
                               CCTK_REAL(2.0) * maskz_loadu(mask, &var[ind]) +
                               maskz_loadu(mask, &var[ind + dk]);
          return fmax3(varxx, varyy, varzz);
        };

        CCTK_REALVEC regrid_error1 =
            fmax5(calcerr(dens), calcerr(momx), calcerr(momy), calcerr(momz),
                  calcerr(etot));
        mask_storeu(mask, &regrid_error[ind], regrid_error1);
      }
    }
  }
}

} // namespace Hydro
