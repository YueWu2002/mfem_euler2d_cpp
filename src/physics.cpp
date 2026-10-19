#include "euler.hpp"
#include <algorithm>
namespace euler {
// Conservative state ordering is (rho, rho*vx, rho*vy, total energy).
double EulerGas::pressure(const State &u) const {
  return gamma_minus_one * (u[3] - .5 * (u[1] * u[1] + u[2] * u[2]) / u[0]);
}
bool EulerGas::positive(const State &u) const {
  for (int c = 0; c < 4; ++c) { if (!std::isfinite(u[c])) { return false; } }
  return u[0] > 1e-13 && pressure(u) > 1e-13;
}
State EulerGas::flux(const State &u, int axis) const {
  double vn = u[axis + 1] / u[0], p = pressure(u);
  State f = {u[axis + 1], u[1] * vn, u[2] * vn, (u[3] + p) * vn};
  f[axis + 1] += p;
  return f;
}
void EulerGas::jacobian(const State &u, int axis, Matrix4 &j) const {
  const double inverse_density = 1. / u[0];
  const double vx = u[1] * inverse_density, vy = u[2] * inverse_density;
  const double vn = u[axis + 1] * inverse_density;
  const double enthalpy_density = u[3] + pressure(u);
  State dp = {.5 * gamma_minus_one * (vx * vx + vy * vy),
              -gamma_minus_one * vx, -gamma_minus_one * vy, gamma_minus_one},
        dv = {-vn * inverse_density, 0, 0, 0};
  dv[axis + 1] = inverse_density;
  j.Set(0.);
  j(0, axis + 1) = 1;
  for (int d = 0; d < 2; d++) {
    for (int c = 0; c < 4; c++) {
      j(d + 1, c) = u[d + 1] * dv[c] + (axis == d ? dp[c] : 0);
    }
    j(d + 1, d + 1) += vn;
  }
  for (int c = 0; c < 4; c++) {
    j(3, c) = vn * dp[c] + enthalpy_density * dv[c];
  }
  j(3, 3) += vn;
}
std::array<double, 2> EulerGas::speeds(const State &l, const State &r, int axis,
                                      double safety) const {
  bool finite = true;
  for (int c = 0; c < 4; ++c) {
    finite = finite && std::isfinite(l[c]) && std::isfinite(r[c]);
  }
  require(finite && l[0] > 1e-13 && r[0] > 1e-13,
          "Nonpositive old face state; refine mesh");
  const double pl = pressure(l), pr = pressure(r);
  require(pl > 1e-13 && pr > 1e-13, "Nonpositive old face state; refine mesh");
  const double inverse_l = 1. / l[0], inverse_r = 1. / r[0];
  double sl = std::sqrt(l[0]), sr = std::sqrt(r[0]), v[2];
  for (int d = 0; d < 2; d++) {
    v[d] = (sl * l[d + 1] * inverse_l + sr * r[d + 1] * inverse_r) / (sl + sr);
  }
  double h =
      (sl * (l[3] + pl) * inverse_l + sr * (r[3] + pr) * inverse_r) / (sl + sr);
  double c2 = gamma_minus_one * (h - .5 * (v[0] * v[0] + v[1] * v[1]));
  require(c2 > 0, "Invalid Roe sound speed");
  double ct = std::sqrt(c2), cl = std::sqrt(gamma * pl * inverse_l),
         cr = std::sqrt(gamma * pr * inverse_r);
  return {safety * std::max({0., cl - l[axis + 1] * inverse_l, ct - v[axis]}),
          safety * std::max({0., cr + r[axis + 1] * inverse_r, ct + v[axis]})};
}
State EulerGas::exact(double x, double y, double time, double alpha) const {
  const double drift = (3. - gamma) / 2., amplitude = (gamma + 1.) / 4.;
  require(time >= 0 && time * pi * amplitude < 1,
          "Exact solution requires T < 4/((gamma+1)*pi)");
  // J_+=2, J_-=sin(pi*xi): solve the left-going characteristic map.
  double phase = std::cos(alpha) * x + std::sin(alpha) * y;
  phase -= 2 * std::floor((phase + 1) / 2);
  phase -= time * drift;
  phase -= 2 * std::floor((phase + 1) / 2);
  double lo = phase - time * amplitude, hi = phase + time * amplitude, z = phase;
  for (int i = 0; i < 80; i++) {
    double r = z + time * amplitude * std::sin(pi * z) - phase;
    if (std::abs(r) < 2e-15) {
      break;
    }
    if (r > 0) {
      hi = z;
    } else {
      lo = z;
    }
    double q = z - r / (1 + time * amplitude * pi * std::cos(pi * z));
    z = (q >= lo && q <= hi) ? q : .5 * (lo + hi);
  }
  const double jm = std::sin(pi * z), v = 1 + .5 * jm;
  const double sound = gamma_minus_one * (2. - jm) / 4.;
  const double rho = std::pow(sound * sound / gamma, 1. / gamma_minus_one);
  const double pressure = std::pow(rho, gamma);
  return {rho, rho * v * std::cos(alpha), rho * v * std::sin(alpha),
          pressure / gamma_minus_one + .5 * rho * v * v};
}
} // namespace euler
