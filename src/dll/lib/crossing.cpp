// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "dll/sim/crossing.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <utility>

namespace dll {

namespace {

/// Squared distance below which two segments touch (lattice coordinates; a bond has length² 2).
constexpr double kContact = 1e-12;

using Poly = std::array<std::int64_t, 4>;  // c0 + c1 t + c2 t^2 + c3 t^3

/// A vector whose components are linear in time: v0 + t v1.
struct Linear {
  Vec3i v0, v1;
};
/// A vector whose components are quadratic in time.
using Quadratic = std::array<Poly, 3>;

Linear difference(Path from, Path to) { return {to.start - from.start, to.step - from.step}; }

std::int64_t at(const Vec3i& v, int k) { return k == 0 ? v.x : (k == 1 ? v.y : v.z); }

Quadratic cross(const Linear& a, const Linear& b) {
  Quadratic c{};
  for (int k = 0; k < 3; ++k) {
    const int i = (k + 1) % 3, j = (k + 2) % 3;
    c[k][0] = at(a.v0, i) * at(b.v0, j) - at(a.v0, j) * at(b.v0, i);
    c[k][1] = at(a.v0, i) * at(b.v1, j) - at(a.v0, j) * at(b.v1, i) + at(a.v1, i) * at(b.v0, j) - at(a.v1, j) * at(b.v0, i);
    c[k][2] = at(a.v1, i) * at(b.v1, j) - at(a.v1, j) * at(b.v1, i);
  }
  return c;
}

Poly dot(const Quadratic& q, const Linear& l) {
  Poly p{};
  for (int k = 0; k < 3; ++k) {
    for (int d = 0; d < 3; ++d) {
      p[d] += q[k][d] * at(l.v0, k);
      p[d + 1] += q[k][d] * at(l.v1, k);
    }
  }
  return p;
}

bool is_zero(const Poly& p) { return p[0] == 0 && p[1] == 0 && p[2] == 0 && p[3] == 0; }

double eval(const Poly& p, double t) {
  return ((static_cast<double>(p[3]) * t + static_cast<double>(p[2])) * t + static_cast<double>(p[1])) * t + static_cast<double>(p[0]);
}

/// Real roots in [0, 1] of a polynomial of degree at most 3 that is not identically zero.
/// Returns their number; double roots are found at the critical points.
int unit_roots(const Poly& p, std::array<double, 8>& roots) {
  int n = 0;
  const auto add = [&](double t) {
    for (int i = 0; i < n; ++i) {
      if (std::fabs(roots[i] - t) < 1e-12) return;
    }
    roots[n++] = t;
  };
  const double scale = static_cast<double>(std::abs(p[0]) + std::abs(p[1]) + std::abs(p[2]) + std::abs(p[3]));
  const double tiny = 1e-12 * scale;
  // Critical points: roots of p' = p1 + 2 p2 t + 3 p3 t^2 inside (0, 1); with t = 0 and t = 1
  // they split [0, 1] into intervals on which p is monotonic.
  std::array<double, 2> crit{};
  int nc = 0;
  const double a = 3.0 * static_cast<double>(p[3]), b = 2.0 * static_cast<double>(p[2]), c = static_cast<double>(p[1]);
  const auto critical = [&](double t) {
    if (t > 0.0 && t < 1.0) crit[nc++] = t;
  };
  if (a != 0.0) {
    const double disc = b * b - 4 * a * c;
    if (disc >= 0) {
      const double q = -0.5 * (b + std::copysign(std::sqrt(disc), b));
      if (q != 0.0) {  // q = 0 only if b = c = 0: the critical point is t = 0, already a break
        critical(q / a);
        critical(c / q);
      }
    }
  } else if (b != 0.0) {
    critical(-c / b);
  }
  if (nc == 2 && crit[0] > crit[1]) std::swap(crit[0], crit[1]);
  std::array<double, 4> breaks{0.0, 1.0, 1.0, 1.0};
  int nb = 1;
  for (int i = 0; i < nc; ++i) breaks[nb++] = crit[i];
  breaks[nb++] = 1.0;
  for (int i = 0; i < nb; ++i) {
    if (std::fabs(eval(p, breaks[i])) <= tiny) add(breaks[i]);
  }
  for (int i = 0; i + 1 < nb; ++i) {
    double lo = breaks[i], hi = breaks[i + 1];
    double flo = eval(p, lo), fhi = eval(p, hi);
    if (std::fabs(flo) <= tiny || std::fabs(fhi) <= tiny || (flo < 0) == (fhi < 0)) continue;
    for (int it = 0; it < 200 && hi - lo > 1e-16; ++it) {
      const double mid = 0.5 * (lo + hi);
      const double fm = eval(p, mid);
      if ((fm < 0) == (flo < 0)) {
        lo = mid;
        flo = fm;
      } else {
        hi = mid;
      }
    }
    add(0.5 * (lo + hi));
  }
  return n;
}

struct Point {
  double x, y, z;
};
Point position(Path p, double t) {
  return {p.start.x + t * p.step.x, p.start.y + t * p.step.y, p.start.z + t * p.step.z};
}
Point sub(Point a, Point b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
double dot(Point a, Point b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

/// Squared distance between the segments p1–q1 and p2–q2 (either may be a single point).
double segment_distance2(Point p1, Point q1, Point p2, Point q2) {
  // Closest points of two segments (C. Ericson, Real-Time Collision Detection, 5.1.9).
  const Point d1 = sub(q1, p1), d2 = sub(q2, p2), r = sub(p1, p2);
  const double a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
  double s = 0.0, t = 0.0;
  if (a <= 0.0 && e <= 0.0) return dot(r, r);
  if (a <= 0.0) {
    t = std::clamp(f / e, 0.0, 1.0);
  } else {
    const double c = dot(d1, r);
    if (e <= 0.0) {
      s = std::clamp(-c / a, 0.0, 1.0);
    } else {
      const double b = dot(d1, d2);
      const double denom = a * e - b * b;
      s = denom > 0.0 ? std::clamp((b * f - c * e) / denom, 0.0, 1.0) : 0.0;
      t = (b * s + f) / e;
      if (t < 0.0) {
        t = 0.0;
        s = std::clamp(-c / a, 0.0, 1.0);
      } else if (t > 1.0) {
        t = 1.0;
        s = std::clamp((b - c) / a, 0.0, 1.0);
      }
    }
  }
  const Point c1{p1.x + d1.x * s, p1.y + d1.y * s, p1.z + d1.z * s};
  const Point c2{p2.x + d2.x * t, p2.y + d2.y * t, p2.z + d2.z * t};
  const Point d = sub(c1, c2);
  return dot(d, d);
}

bool touch_at(Path a, Path b, Path c, Path d, double t) {
  return segment_distance2(position(a, t), position(b, t), position(c, t), position(d, t)) <= kContact;
}

/// Whether the bead x meets the segment y–z: x must then be collinear with y and z.
bool point_meets_segment(Path x, Path y, Path z) {
  const Quadratic w = cross(difference(x, y), difference(x, z));
  std::array<double, 8> roots{};
  for (const Poly& component : w) {
    if (is_zero(component)) continue;
    const int n = unit_roots(component, roots);
    for (int i = 0; i < n; ++i) {
      if (touch_at(x, x, y, z, roots[i])) return true;
    }
    return false;  // a common root of all components is a root of this one
  }
  return false;  // always collinear: x could only reach the segment through one of its ends
}

}  // namespace

bool segments_touch(Path a, Path b, Path c, Path d) noexcept {
  // The segments can meet only when a, b, c and d are coplanar: det[b - a, d - c, c - a] = 0.
  const Poly coplanar = dot(cross(difference(a, b), difference(c, d)), difference(a, c));
  if (!is_zero(coplanar)) {
    std::array<double, 8> roots{};
    const int n = unit_roots(coplanar, roots);
    for (int i = 0; i < n; ++i) {
      if (touch_at(a, b, c, d, roots[i])) return true;
    }
    return false;
  }
  // Coplanar throughout: the first contact of two segments in a plane puts an end of one on
  // the other.
  if (point_meets_segment(c, a, b)) return true;
  if (d == c) return false;
  return point_meets_segment(d, a, b) || point_meets_segment(a, c, d) || point_meets_segment(b, c, d);
}

}  // namespace dll
