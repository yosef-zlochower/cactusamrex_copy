#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <iostream>
#include <vector>

#include "loop.hxx"
#include "loopcontrol.h"

namespace TestLoops {
using namespace std;

CCTK_REAL fun1d(const CCTK_REAL x, const CCTK_REAL dx, const bool avgx,
                const int order) {
  if (!avgx)
    return pow(x, order);
  else
    return (pow(x + dx / 2, order + 1) - pow(x - dx / 2, order + 1)) /
           ((order + 1) * dx);
}

// The grid stores the average values of the underlying function (x*y*z)**n
// which amounts to storing differences of the anti-derivative 1/(n+1)**3 *
// (x*y*z)**(n+1)
CCTK_REAL fun(const Loop::PointDesc &p, const bool avgx, const bool avgy,
              const bool avgz, const int order) {
  return fun1d(p.x, p.dx, avgx, order) * fun1d(p.y, p.dy, avgy, order) *
         fun1d(p.z, p.dz, avgz, order);
}

extern "C" void TestLoops_Set(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  std::vector<ptrdiff_t> inds(cctk_ash[0] * cctk_ash[1] * cctk_ash[2]);

  const int order = 3; // an order that is higher than linear so that the value
                       // at the cell center is not the same as the average
                       // over the cell

  // check that template and macro loop in the same way
  int ind_template = 0;
  Loop::loop_all<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    inds[ind_template++] = p.idx;
    gf000[p.idx] = fun(p, 0, 0, 0, order);
  });

  *template_and_macro_diff = 0;
  int ind_macro = 0;
  CCTK_LOOP3_ALL(TestLoops000, cctkGH, i, j, k) {
    int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
    if (ind_macro >= ind_template) {
      *template_and_macro_diff = 1;
    } else if (inds[ind_macro] != idx) {
      *template_and_macro_diff -= 1;
    }
    ind_macro += 1;
  }
  CCTK_ENDLOOP3_ALL(TestLoops000);
  if (ind_macro != ind_template) {
    *template_and_macro_diff = 2;
  }

  // check
}

} // namespace TestLoops
