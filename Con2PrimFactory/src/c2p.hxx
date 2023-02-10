/*! \file c2p.hxx
\brief Defines a c2p
\author Jay Kalinani

c2p is effectively an interface to be used by different c2p implementations.

*/

#ifndef C2P_HXX
#define C2P_HXX
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <boost/math/tools/roots.hpp>
#include "c2p_utils.hxx"
#include "eos.hxx"
#include "eos_idealgas.hxx"

constexpr CCTK_INT X = 0;
constexpr CCTK_INT Y = 1;
constexpr CCTK_INT Z = 2;

#include <math.h>
#include "prims.hxx"
#include "cons.hxx"
#include "atmo.hxx"

namespace Con2PrimFactory {

/* Abstract class c2p */
class c2p {
public:
  /* The constructor must initialize the following variables */

  CCTK_INT maxIterations;
  CCTK_REAL tolerance;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Ssq_Exact(vec<CCTK_REAL, 3> &mom, const smat<CCTK_REAL, 3> &gup) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Bsq_Exact(vec<CCTK_REAL, 3> &B_up, const smat<CCTK_REAL, 3> &glo) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_BiSi_Exact(vec<CCTK_REAL, 3> &Bvec, vec<CCTK_REAL, 3> &mom) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 2>
  get_WLorentz_bsq_Seeds(vec<CCTK_REAL, 3> &B_up, vec<CCTK_REAL, 3> &v_up,
                         const smat<CCTK_REAL, 3> &glo) const;
};

} // namespace Con2PrimFactory

#endif
