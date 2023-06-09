#include "derivs.hxx"
#include "physics.hxx"
#include "tensor.hxx"

#include <cctk.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <random>

namespace Weyl {
using namespace std;

// TODO: Use GoogleTest instead of assert

extern "C" void Weyl_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  // Test tensors
  {
    mt19937 engine(42);
    uniform_int_distribution<int> dist(-10, 10);
    const auto rand10{[&]() { return double(dist(engine)); }};
    const auto randmat10{[&]() {
      array<array<double, 3>, 3> arr;
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          arr[a][b] = rand10();
      return mat3<double, DN, DN>(
          [&](int a, int b) { return arr[min(a, b)][max(a, b)]; });
    }};

    const mat3<double, DN, DN> Z([&](int a, int b) { return double(0); });
    const mat3<double, DN, DN> I([&](int a, int b) { return double(a == b); });
    assert(I != Z);
    const mat3<double, UP, UP> Zup([&](int a, int b) { return double(0); });
    const mat3<double, UP, UP> Iup(
        [&](int a, int b) { return double(a == b); });
    assert(Iup != Zup);

    for (int n = 0; n < 100; ++n) {
      const mat3<double, DN, DN> A = randmat10();
      const mat3<double, DN, DN> B = randmat10();
      const mat3<double, DN, DN> C = randmat10();
      const double a = rand10();
      const double b = rand10();

      assert((A + B) + C == A + (B + C));
      assert(Z + A == A);
      assert(A + Z == A);
      assert(A + (-A) == Z);
      assert((-A) + A == Z);
      assert(A - B == A + (-B));
      assert(A + B == B + A);

      assert(1 * A == A);
      assert(0 * A == Z);
      assert(-1 * A == -A);
      // assert(mul(a * A, B) == a * mul(A, B));
      assert((a * b) * A == a * (b * A));
      assert(a * (A + B) == a * A + a * B);
      assert((a + b) * A == a * A + b * A);

      // assert(mul(mul(A, B), C) == mul(A, mul(B, C)));
      // DNUP  assert(mul(I, A) == A);
      // DNUP  assert(mul(A, I) == A);
      // DNUP  assert(mul(Z, A) == Z);
      // DNUP  assert(mul(A, Z) == Z);

      assert(Z.det() == 0);
      assert(I.det() == 1);
      assert((a * A).det() == pow3(a) * A.det());

      assert(Z.inv(1) == Zup);
      assert(I.inv(1) == Iup);

      // DNUP assert(mul(A.inv(1), A) == A.det() * Iup);
      // DNUP assert(mul(A, A.inv(1)) == A.det() * Iup);

      assert((a * A).inv(1) == pow2(a) * A.inv(1));
    }
  }

  // Test 4-tensors
  {
    mt19937 engine(42);
    uniform_int_distribution<int> dist(-10, 10);
    const auto rand10{[&]() { return double(dist(engine)); }};
    const auto randmat10{[&]() {
      array<array<double, 4>, 4> arr;
      for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
          arr[a][b] = rand10();
      return mat4<double, DN, DN>(
          [&](int a, int b) { return arr[min(a, b)][max(a, b)]; });
    }};

    const mat4<double, DN, DN> Z([&](int a, int b) { return double(0); });
    const mat4<double, DN, DN> I([&](int a, int b) { return double(a == b); });
    assert(I != Z);
    const mat4<double, UP, UP> Zup([&](int a, int b) { return double(0); });
    const mat4<double, UP, UP> Iup(
        [&](int a, int b) { return double(a == b); });
    assert(Iup != Zup);

    for (int n = 0; n < 100; ++n) {
      const mat4<double, DN, DN> A = randmat10();
      const mat4<double, DN, DN> B = randmat10();
      const mat4<double, DN, DN> C = randmat10();
      const double a = rand10();
      const double b = rand10();

      assert((A + B) + C == A + (B + C));
      assert(Z + A == A);
      assert(A + Z == A);
      assert(A + (-A) == Z);
      assert((-A) + A == Z);
      assert(A - B == A + (-B));
      assert(A + B == B + A);

      assert(1 * A == A);
      assert(0 * A == Z);
      assert(-1 * A == -A);
      // assert(mul(a * A, B) == a * mul(A, B));
      assert((a * b) * A == a * (b * A));
      assert(a * (A + B) == a * A + a * B);
      assert((a + b) * A == a * A + b * A);

      // assert(mul(mul(A, B), C) == mul(A, mul(B, C)));
      // DNUP  assert(mul(I, A) == A);
      // DNUP  assert(mul(A, I) == A);
      // DNUP  assert(mul(Z, A) == Z);
      // DNUP  assert(mul(A, Z) == Z);

      assert(Z.det() == 0);
      assert(I.det() == 1);
      assert((a * A).det() == pow4(a) * A.det());

      assert(Z.inv(1) == Zup);
      assert(I.inv(1) == Iup);

      // DNUP assert(mul(A.inv(1), A) == A.det() * Iup);
      // DNUP assert(mul(A, A.inv(1)) == A.det() * Iup);

      assert((a * A).inv(1) == pow3(a) * A.inv(1));
    }
  }

  // Test antisymmetric 4-tensors
  {
    mt19937 engine(42);
    uniform_int_distribution<int> dist(-10, 10);
    const auto rand10{[&]() { return double(dist(engine)); }};
    const auto randamat10{[&]() {
      array<array<double, 4>, 4> arr;
      for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
          arr[a][b] = rand10();
      return amat4<double, DN, DN>(
          [&](int a, int b) { return arr[min(a, b)][max(a, b)]; });
    }};

    const amat4<double, DN, DN> Z([&](int a, int b) { return double(0); });

    for (int n = 0; n < 100; ++n) {
      const amat4<double, DN, DN> A = randamat10();
      const amat4<double, DN, DN> B = randamat10();
      const amat4<double, DN, DN> C = randamat10();
      const double a = rand10();
      const double b = rand10();

      assert((A + B) + C == A + (B + C));
      assert(Z + A == A);
      assert(A + Z == A);
      assert(A + (-A) == Z);
      assert((-A) + A == Z);
      assert(A - B == A + (-B));
      assert(A + B == B + A);

      assert(1 * A == A);
      assert(0 * A == Z);
      assert(-1 * A == -A);
      assert((a * b) * A == a * (b * A));
      assert(a * (A + B) == a * A + a * B);
      assert((a + b) * A == a * A + b * A);
    }
  }

  // Test Riemann tensors
  {
    mt19937 engine(42);
    uniform_int_distribution<int> dist(-10, 10);
    const auto rand10{[&]() { return double(dist(engine)); }};
    const auto randrten10{[&]() {
      array<array<array<array<double, 4>, 4>, 4>, 4> arr;
      for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
          for (int c = 0; c < 4; ++c)
            for (int d = 0; d < 4; ++d)
              arr[a][b][c][d] = rand10();
      return rten4<double, DN, DN, DN, DN>(
          [&](int a, int b, int c, int d) { return arr[a][b][c][d]; });
    }};

    const rten4<double, DN, DN, DN, DN> Z(
        [&](int a, int b, int c, int d) { return double(0); });

    for (int n = 0; n < 100; ++n) {
      const rten4<double, DN, DN, DN, DN> A = randrten10();
      const rten4<double, DN, DN, DN, DN> B = randrten10();
      const rten4<double, DN, DN, DN, DN> C = randrten10();
      const double a = rand10();
      const double b = rand10();

      assert((A + B) + C == A + (B + C));
      assert(Z + A == A);
      assert(A + Z == A);
      assert(A + (-A) == Z);
      assert((-A) + A == Z);
      assert(A - B == A + (-B));
      assert(A + B == B + A);

      assert(1 * A == A);
      assert(0 * A == Z);
      assert(-1 * A == -A);
      assert((a * b) * A == a * (b * A));
      assert(a * (A + B) == a * A + a * B);
      assert((a + b) * A == a * A + b * A);
    }
  }

  // Test derivatives

  static_assert(deriv_order % 2 == 0, "");
  constexpr int required_ghosts = deriv_order / 2 + 1;
  constexpr int fences = 3;
  constexpr int vsize = tuple_size_v<simd<double> >;
  const double eps = 1.0e-12;

  // deriv
  for (int npoints = 1; npoints <= vsize; ++npoints) {
    for (int order = 0; order <= deriv_order; ++order) {
      // CCTK_VINFO("Testing deriv order=%d", order);
      array<double, 2 * (fences + required_ghosts) + vsize> arr;
      for (size_t i = 0; i < arr.size(); ++i)
        arr[i] = NAN;
      double *const var = &arr[fences + required_ghosts];
      for (int i = -deriv_order / 2; i < vsize + deriv_order / 2; ++i)
        var[i] = pown(i, order);
      const simd<double> expected =
          order == 0 ? 0 : order * pown(iota<simd<double> >(), order - 1);
      const simdl<double> mask = mask_for_loop_tail<simdl<double> >(0, npoints);
      const simd<double> found = deriv1d(mask, var, 1, 1.0);
      if (!(all(fabs(found - expected) <= eps || !mask)))
        cout << "deriv:\n"
             << "  npoints: " << npoints << "\n"
             << "  order: " << order << "\n"
             << "  expected: " << expected << "\n"
             << "  found: " << found << "\n";
      assert(all(fabs(found - expected) <= eps || !mask));
    }
  }

  // deriv (upwind)
  for (int npoints = 1; npoints <= vsize; ++npoints) {
    for (int order = 0; order <= deriv_order; ++order) {
      for (int sign = 0; sign <= 1; ++sign) {
        // CCTK_VINFO("Testing deriv (upwind) order=%d sign=%d", order, sign);
        array<double, 2 * (fences + required_ghosts) + vsize> arr;
        for (size_t i = 0; i < arr.size(); ++i)
          arr[i] = NAN;
        double *const var = &arr[fences + required_ghosts];
        for (int i = -deriv_order / 2 - 1; i < vsize + deriv_order / 2 + 1; ++i)
          var[i] = pown(i, order);
        const simd<double> vel = sign ? -1 : +1;
        const simd<double> expected =
            vel *
            (order == 0 ? 0 : order * pown(iota<simd<double> >(), order - 1));
        const simdl<double> mask =
            mask_for_loop_tail<simdl<double> >(0, npoints);
        const simd<double> found = deriv1d_upwind(mask, var, 1, vel, 1.0);
        if (!(all(fabs(found - expected) <= eps || !mask)))
          cout << "deriv_upwind:\n"
               << "  npoints: " << npoints << "\n"
               << "  order: " << order << "\n"
               << "  sign: " << sign << "\n"
               << "  expected: " << expected << "\n"
               << "  found: " << found << "\n";
        assert(all(fabs(found - expected) <= eps || !mask));
      }
    }
  }

  // deriv2
  for (int npoints = 1; npoints <= vsize; ++npoints) {
    for (int order = 0; order <= deriv_order; ++order) {
      // CCTK_VINFO("Testing deriv2 order=%d", order);
      array<double, 2 * (fences + required_ghosts) + vsize> arr;
      for (size_t i = 0; i < arr.size(); ++i)
        arr[i] = NAN;
      double *const var = &arr[fences + required_ghosts];
      for (int i = -deriv_order / 2; i < vsize + deriv_order / 2; ++i)
        var[i] = pown(i, order);
      const simd<double> expected =
          order < 2
              ? 0
              : order * (order - 1) * pown(iota<simd<double> >(), order - 2);
      const simdl<double> mask = mask_for_loop_tail<simdl<double> >(0, npoints);
      const simd<double> found = deriv2_1d(mask, var, 1, 1.0);
      if (!(all(fabs(found - expected) <= eps || !mask)))
        cout << "deriv_upwind:\n"
             << "  npoints: " << npoints << "\n"
             << "  order: " << order << "\n"
             << "  expected: " << expected << "\n"
             << "  found: " << found << "\n";
      assert(all(fabs(found - expected) <= eps || !mask));
    }
  }

  // deriv2 (mixed)
  for (int npoints = 1; npoints <= vsize; ++npoints) {
    for (int orderj = 0; orderj <= deriv_order; ++orderj) {
      for (int orderi = 0; orderi <= deriv_order; ++orderi) {
        // CCTK_VINFO("Testing deriv2 (mixed) order=%d,%d", orderi, orderj);
        array<array<double, 2 * (fences + required_ghosts) + vsize>,
              2 * (fences + required_ghosts) + 1>
            arr;
        for (size_t j = 0; j < arr.size(); ++j)
          for (size_t i = 0; i < arr[0].size(); ++i)
            arr[j][i] = NAN;
        const int di = 1;
        const int dj = arr[0].size();
        double *const var =
            &arr[fences + required_ghosts][fences + required_ghosts];
        for (int j = -deriv_order / 2; j < 1 + deriv_order / 2; ++j)
          for (int i = -deriv_order / 2; i < vsize + deriv_order / 2; ++i)
            var[j * dj + i * di] = pown(i, orderi) * pown(j, orderj);
        const simd<double> expected =
            (orderi == 0 ? 0
                         : orderi * pown(iota<simd<double> >(), orderi - 1)) *
            (orderj == 0 ? 0 : orderj * pown(0.0, orderj - 1));
        const simdl<double> mask =
            mask_for_loop_tail<simdl<double> >(0, npoints);
        const simd<double> found =
            deriv2_2d(npoints, mask, var, di, dj, 1.0, 1.0);
        if (!(all(fabs(found - expected) <= eps || !mask)))
          cout << "deriv2_mixed:\n"
               << "  npoints: " << npoints << "\n"
               << "  orderi: " << orderi << "\n"
               << "  orderj: " << orderj << "\n"
               << "  expected: " << expected << "\n"
               << "  found: " << found << "\n";
        assert(all(fabs(found - expected) <= eps || !mask));
      }
    }
  }
}

} // namespace Weyl
