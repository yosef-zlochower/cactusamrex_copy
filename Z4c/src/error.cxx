#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace Z4c {
using namespace Loop;
using namespace std;

extern "C" void Z4c_Estimate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Estimate;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_chi_(cctkGH, chi);
  const GF3D<CCTK_REAL, 1, 1, 1> gf_regrid_error_(cctkGH, regrid_error);

  loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    // Load
    const CCTK_REAL chi = gf_chi_(p.I);

    // assume chi = psi^4 = (1+M/r)^4
    const CCTK_REAL err = dx*fabs(pow(chi, 0.25)-1.);

    // Store
    gf_regrid_error_(p.I) = err;
  });
}

} // namespace Z4c
