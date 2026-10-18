#include "euler.hpp"
#include <algorithm>
namespace euler {
// Conservative state ordering is (rho, rho*vx, rho*vy, total energy).
// This accuracy case fixes gamma=3.
double pressure(const State &u) {
  return 2 * (u[3] - .5 * (u[1] * u[1] + u[2] * u[2]) / u[0]);
}
bool positive(const State &u) {
  for (int c = 0; c < 4; ++c) { if (!std::isfinite(u[c])) { return false; } }
  return u[0] > 1e-13 && pressure(u) > 1e-13;
}
State flux(const State &u, int axis) {
  double vn = u[axis + 1] / u[0], p = pressure(u);
  State f = {u[axis + 1], u[1] * vn, u[2] * vn, (u[3] + p) * vn};
  f[axis + 1] += p;
  return f;
}
void jacobian(const State &u, int axis, Matrix4 &j) {
  double vx = u[1] / u[0], vy = u[2] / u[0], vn = u[axis + 1] / u[0];
  State dp = {vx * vx + vy * vy, -2 * vx, -2 * vy, 2}, dv = {-vn / u[0], 0, 0, 0};
  dv[axis + 1] = 1 / u[0];
  j.Set(0.);
  j(0, axis + 1) = 1;
  for (int d = 0; d < 2; d++) {
    for (int c = 0; c < 4; c++) {
      j(d + 1, c) = u[d + 1] * dv[c] + (axis == d ? dp[c] : 0);
    }
    j(d + 1, d + 1) += vn;
  }
  for (int c = 0; c < 4; c++) {
    j(3, c) = vn * dp[c] + (u[3] + pressure(u)) * dv[c];
  }
  j(3, 3) += vn;
}
std::array<double, 2> speeds(const State &l, const State &r, int axis, double safety) {
  require(positive(l) && positive(r), "Nonpositive old face state; refine mesh");
  double sl = std::sqrt(l[0]), sr = std::sqrt(r[0]), v[2];
  for (int d = 0; d < 2; d++) {
    v[d] = (sl * l[d + 1] / l[0] + sr * r[d + 1] / r[0]) / (sl + sr);
  }
  double h =
      (sl * (l[3] + pressure(l)) / l[0] + sr * (r[3] + pressure(r)) / r[0]) / (sl + sr);
  double c2 = 2 * (h - .5 * (v[0] * v[0] + v[1] * v[1]));
  require(c2 > 0, "Invalid Roe sound speed");
  double ct = std::sqrt(c2), cl = std::sqrt(3 * pressure(l) / l[0]),
         cr = std::sqrt(3 * pressure(r) / r[0]);
  return {safety * std::max({0., cl - l[axis + 1] / l[0], ct - v[axis]}),
          safety * std::max({0., cr + r[axis + 1] / r[0], ct + v[axis]})};
}
State exact(double x, double y, double time, double alpha) {
  require(time >= 0 && time * pi < 1, "Exact solution requires T<1/pi");
  // Solve xi + t*sin(pi*xi) = phase using bracketed Newton.
  double phase = std::cos(alpha) * x + std::sin(alpha) * y;
  phase -= 2 * std::floor((phase + 1) / 2);
  double lo = phase - time, hi = phase + time, z = phase;
  for (int i = 0; i < 80; i++) {
    double r = z + time * std::sin(pi * z) - phase;
    if (std::abs(r) < 2e-15) {
      break;
    }
    if (r > 0) {
      hi = z;
    } else {
      lo = z;
    }
    double q = z - r / (1 + time * pi * std::cos(pi * z));
    z = (q >= lo && q <= hi) ? q : .5 * (lo + hi);
  }
  double jm = std::sin(pi * z), v = 1 + .5 * jm, rho = (1 - .5 * jm) / std::sqrt(3.);
  return {rho, rho * v * std::cos(alpha), rho * v * std::sin(alpha),
          .5 * rho * rho * rho + .5 * rho * v * v};
}
} // namespace euler
