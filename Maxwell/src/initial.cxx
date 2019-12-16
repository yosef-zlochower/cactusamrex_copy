#include "dual.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>

namespace Maxwell {
using namespace std;

namespace {
inline int bitsign(bool s) { return s ? -1 : +1; }

template <typename T> T pow2(T x) { return x * x; }
template <typename T> T sinc(T x) { return x == T(0) ? T(1) : sin(x) / x; }
} // namespace

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct potential {
  // Electric scalar potential
  T phi;
  // Electric vector potential (to ensure div E = 0)
  T cyz, czx, cxy;
  // Magnetic vector potential
  T ax, ay, az;
};

// Continuous derivative
template <typename F, typename T>
potential<T> calc_dt(const F &f, T t, T x, T y, T z) {
  auto fd = f(dual<T>(t, 1), dual<T>(x), dual<T>(y), dual<T>(z));
  return {
      fd.phi.eps, fd.cyz.eps, fd.czx.eps, fd.cxy.eps,
      fd.ax.eps,  fd.ay.eps,  fd.az.eps,
  };
}

// Discrete derivatives (centred)
template <typename F, typename T>
potential<T> calc_dxc(const F &f, T t, T x, T y, T z, T dx) {
  auto fm = f(t, x - dx / 2, y, z);
  auto fp = f(t, x + dx / 2, y, z);
  return {
      .phi = (fp.phi - fm.phi) / dx,
      .cyz = (fp.cyz - fm.cyz) / dx,
      .czx = (fp.czx - fm.czx) / dx,
      .cxy = (fp.cxy - fm.cxy) / dx,
      .ax = (fp.ax - fm.ax) / dx,
      .ay = (fp.ay - fm.ay) / dx,
      .az = (fp.az - fm.az) / dx,
  };
}

template <typename F, typename T>
potential<T> calc_dyc(const F &f, T t, T x, T y, T z, T dy) {
  auto fm = f(t, x, y - dy / 2, z);
  auto fp = f(t, x, y + dy / 2, z);
  return {
      .phi = (fp.phi - fm.phi) / dy,
      .cyz = (fp.cyz - fm.cyz) / dy,
      .czx = (fp.czx - fm.czx) / dy,
      .cxy = (fp.cxy - fm.cxy) / dy,
      .ax = (fp.ax - fm.ax) / dy,
      .ay = (fp.ay - fm.ay) / dy,
      .az = (fp.az - fm.az) / dy,
  };
}

template <typename F, typename T>
potential<T> calc_dzc(const F &f, T t, T x, T y, T z, T dz) {
  auto fm = f(t, x, y, z - dz / 2);
  auto fp = f(t, x, y, z + dz / 2);
  return {
      .phi = (fp.phi - fm.phi) / dz,
      .cyz = (fp.cyz - fm.cyz) / dz,
      .czx = (fp.czx - fm.czx) / dz,
      .cxy = (fp.cxy - fm.cxy) / dz,
      .ax = (fp.ax - fm.ax) / dz,
      .ay = (fp.ay - fm.ay) / dz,
      .az = (fp.az - fm.az) / dz,
  };
}

////////////////////////////////////////////////////////////////////////////////

// Plane wave
template <typename T>
potential<T> plane_wave(const T t, const T x, const T y, const T z) {
  DECLARE_CCTK_PARAMETERS;
  // wave number
  T kx = M_PI * spatial_frequency_x;
  T ky = M_PI * spatial_frequency_y;
  T kz = M_PI * spatial_frequency_z;
  assert(kx == 0);
  assert(ky == 0);
  T omega = sqrt(pow2(kx) + pow2(ky) + pow2(kz));
  // amplitude
  T hx = amplitude_x;
  T hy = amplitude_y;
  T hz = amplitude_z;
  assert(hy == 0);
  assert(hz == 0);
  // solution
  assert(t == 0);
  T u = sin(omega * t - kz * z);
  return {
      .phi = 0,
      .cyz = 0,
      .czx = hx / kz * u,
      .cxy = 0,
      .ax = -hx / kz * u,
      .ay = 0,
      .az = 0,
  };
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void Maxwell_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_Initial;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> ax_(cctkGH, ax);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> ay_(cctkGH, ay);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> az_(cctkGH, az);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> ex_(cctkGH, ex);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> ey_(cctkGH, ey);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> ez_(cctkGH, ez);

  const Loop::GF3D<CCTK_REAL, 0, 1, 1> byz_(cctkGH, byz);
  const Loop::GF3D<CCTK_REAL, 1, 0, 1> bzx_(cctkGH, bzx);
  const Loop::GF3D<CCTK_REAL, 1, 1, 0> bxy_(cctkGH, bxy);

  const auto loop_setup{[&](const auto &f4) {
    const auto f{[&](const auto &p) { return f4(t, p.x, p.y, p.z); }};
    const auto dtf{[&](const auto &p) {
      return calc_dt(
          [&](auto t, auto x, auto y, auto z) { return f4(t, x, y, z); }, t,
          p.x, p.y, p.z);
    }};
    const auto dxf{[&](const auto &p) {
      return calc_dxc(
          [&](auto t, auto x, auto y, auto z) { return f4(t, x, y, z); }, t,
          p.x, p.y, p.z, dx);
    }};
    const auto dyf{[&](const auto &p) {
      return calc_dyc(
          [&](auto t, auto x, auto y, auto z) { return f4(t, x, y, z); }, t,
          p.x, p.y, p.z, dy);
    }};
    const auto dzf{[&](const auto &p) {
      return calc_dzc(
          [&](auto t, auto x, auto y, auto z) { return f4(t, x, y, z); }, t,
          p.x, p.y, p.z, dz);
    }};

    Loop::loop_int<0, 0, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { phi_(p.I) = f(p).phi; });

    Loop::loop_int<1, 0, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { ax_(p.I) = f(p).ax; });
    Loop::loop_int<0, 1, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { ay_(p.I) = f(p).ay; });
    Loop::loop_int<0, 0, 1>(
        cctkGH, [&](const Loop::PointDesc &p) { az_(p.I) = f(p).az; });

    Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ex_(p.I) = -dxf(p).phi + dyf(p).cxy - dzf(p).czx;
    });
    Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ey_(p.I) = -dyf(p).phi + dzf(p).cyz - dxf(p).cxy;
    });
    Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      ez_(p.I) = -dzf(p).phi + dxf(p).czx - dyf(p).cyz;
    });
    // Loop::loop_int<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    //   dive_(p.I) = dxm(ex_, p) /*+ dym(ey_, p) + dzm(ez_, p)*/;
    // });

    Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      byz_(p.I) = dyf(p).az - dzf(p).ay;
    });
    Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      bzx_(p.I) = dzf(p).ax - dxf(p).az;
    });
    Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      bxy_(p.I) = dxf(p).ay - dyf(p).ax;
    });
  }};

  if (CCTK_EQUALS(setup, "plane wave")) {
    loop_setup(
        [&](auto t, auto x, auto y, auto z) { return plane_wave(t, x, y, z); });
  } else {
    assert(0);
  }
}

} // namespace Maxwell
