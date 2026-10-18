#include "euler.hpp"
#include <algorithm>
#include <iostream>
namespace euler {
void VerifyPhysics();
void VerifyWeighted(const Options &base);
void Verify(const Options &base) {
  VerifyPhysics();
  VerifyWeighted(base);
  std::ofstream report(base.output + "_verification.csv");
  report << "k,jacobian_fd_relative,dg_hdg_increment_relative,assembled_mf_"
            "relative,mass_off_diagonal,partial_dg_relative,partial_schur_"
            "relative,off_manifold_increment_relative\n"
         << std::setprecision(17);
  for (int k : {0, 1, 2, 3, 4, 5, 6, 7, 8}) {
    Options o = base;
    o.k = k;
    o.nx = 3; o.ny = 2;
    o.precond = "block-jacobi";
    o.matrix_free = o.solver == "gmres";
    o.linear_tol = 1e-12;
    Problem p(o);
    Vector u = p.Initial();
    auto background = exact(0, 0, 0, p.alpha);
    for (int i = 0; i < p.nu; i++) {
      u[i] = .05 * u[i] + .95 * background[i % 4];
    }
    p.Freeze(u, .001);
    auto dg = p.DG(u);
    auto a = dg->Assemble();
    Vector v(p.nu), z(p.nu), az(p.nu);
    for (int i = 0; i < p.nu; i++) {
      v[i] = std::sin(.37 * (i + 1));
    }
    dg->Mult(v, z);
    a->Mult(v, az);
    az -= z;
    double mfrel = az.Norml2() / z.Norml2();
    double eps = 1e-7;
    Vector up = u, um = u;
    up.Add(eps, v);
    um.Add(-eps, v);
    Vector fd = p.Residual(up);
    fd -= p.Residual(um);
    fd *= .5 / eps;
    fd -= z;
    double fdrel = fd.Norml2() / z.Norml2();
    Vector r = p.Residual(u), rhs = r;
    rhs *= -1;
    auto d = Solve(*dg, rhs, p.m, o, p.cell_metric);
    Vector lam = p.Trace(u);
    auto hr = p.Residual(u, &lam);
    auto h = p.HDG(u);
    auto s = p.Condense(h, hr);
    auto l = Solve(*s.op, s.rhs, p.t, o, p.trace_metric);
    Vector hd = p.Recover(s, l.x);
    hd -= d.x;
    double equiv = hd.Norml2() / d.x.Norml2();
    auto pa = p.PartialDG(u);
    Vector pz(p.nu);
    pa->Mult(v, pz);
    pz -= z;
    double parel = pz.Norml2() / z.Norml2();
    Options po = o;
    po.assembly = "partial";
    Problem pp(po);
    pp.Freeze(u, .001);
    auto ph = pp.HDG(u);
    auto ps = pp.Condense(ph, hr);
    Vector sz(ps.op->Height()), ez(s.op->Height());
    for (int i = 0; i < sz.Size(); i++) {
      sz[i] = std::sin(.27 * i);
    }
    ps.op->Mult(sz, ez);
    Vector ref(ez.Size());
    s.op->Mult(sz, ref);
    ez -= ref;
    double schurpa = ez.Norml2() / ref.Norml2();
    double off = 0;
    if (k <= 3) {
      Vector offlam = lam;
      for (int i = 0; i < p.nl; i++) {
        offlam[i] += 1e-3 * std::sin(.31 * i);
      }
      auto offr = p.Residual(u, &offlam), offg = p.Skeleton(u, offlam);
      auto os = p.Condense(h, offr, &offg);
      auto ol = Solve(*os.op, os.rhs, p.t, o, p.trace_metric);
      auto od = p.Recover(os, ol.x);
      od -= d.x;
      off = od.Norml2() / d.x.Norml2();
    }
    require(parel < 1e-12 && schurpa < 1e-12 && off < 1e-9,
            "Partial/off-manifold verification failed");
    double massoff = 0;
    for (int i = 0; i < p.p; i++) {
      for (int j = 0; j < p.p; j++) {
        if (i != j) {
          massoff = std::max(massoff, std::abs(p.mass(i, j)));
        }
      }
    }
    require(fdrel < 1e-7 && equiv < 1e-8 && mfrel < 1e-13 && massoff < 1e-12,
            "Algebra verification failed");
    report << k << ',' << fdrel << ',' << equiv << ',' << mfrel << ',' << massoff << ','
           << parel << ',' << schurpa << ',' << (k <= 3 ? off : -1.) << std::endl;
    std::cout << "verified k=" << k << " FD=" << fdrel << " equivalence=" << equiv
              << '\n';
  }
}

// A dense, nonorthogonal change within each block requires the FULL Gram
// matrix. Whiten it by Cholesky, then compare the physical correction and the
// actual restarted-GMRES iteration count with the original diagonal metric.
void VerifyWeighted(const Options &base) {
  constexpr int n = 6;
  using Mat = mfem::TMatrix<n,n>;
  Mat k{}, t{}, kt{}, ht{}, l{}, li{}, normalized{};
  mfem::TVector<n> mass{.2,.7,1.1,.3,.8,1.4};
  for (int i = 0; i < n; ++i) {
    k(i, i) = 4 + .3*i;
    t(i, i) = .7 + .2*i;
    if (i+1 < n) { k(i, i+1) = .3; k(i+1, i) = -.2; }
  }
  t(0, 1) = .25; t(1, 2) = -.15; t(3, 5) = .2;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int a = 0; a < n; ++a) {
        ht(i, j) += t(a, i)*mass[a]*t(a, j);
        for (int b = 0; b < n; ++b) { kt(i, j) += t(a, i)*k(a, b)*t(b, j); }
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double value = ht(i, j);
      for (int a = 0; a < j; ++a) { value -= l(i, a)*l(j, a); }
      if (i == j) { require(value > 0, "Test Gram must be SPD"); l(i, j) = std::sqrt(value); }
      else { l(i, j) = value/l(j, j); }
    }
  }
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      double value = i == j ? 1. : 0.;
      for (int a = 0; a < i; ++a) { value -= l(i, a)*li(a, j); }
      li(i, j) = value/l(i, i);
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int a = 0; a < n; ++a) {
        for (int b = 0; b < n; ++b) { normalized(i, j) += li(i, a)*kt(a, b)*li(j, b); }
      }
    }
  }
  Vector rhs(n), transformed_rhs(n), metric(n), ones(n);
  ones = 1.; transformed_rhs = 0.;
  std::vector<int> ids{0,1,2,3,4,5};
  DenseMatrix original(n), changed(n);
  for (int i = 0; i < n; ++i) {
    rhs[i] = std::sin(i+1.); metric[i] = 1/mass[i];
    for (int j = 0; j < n; ++j) { original(i,j) = k(i, j); changed(i,j) = normalized(i, j); }
  }
  for (int i = 0; i < n; ++i) {
    for (int a = 0; a < n; ++a) {
      for (int b = 0; b < n; ++b) { transformed_rhs[i] += li(i, a)*t(b, a)*rhs[b]; }
    }
  }
  ElementOperator a(n), b(n);
  a.Add(ids, ids, original); b.Add(ids, ids, changed);
  std::ofstream report(base.output + "_weighted_gmres.csv");
  report << "preconditioner,physical_relative_difference,original_iterations,transformed_iterations\n"
         << std::setprecision(17);
  for (const std::string pre : {"block-jacobi", "block-ilu0"}) {
    Options o = base;
    o.solver = "gmres"; o.precond = pre; o.ilu_ordering = "none";
    o.matrix_free = true; o.restart = 2; o.linear_tol = 1e-10;
    auto reference = Solve(a, rhs, 3, o, metric);
    auto transformed = Solve(b, transformed_rhs, 3, o, ones);
    Vector physical(n); physical = 0.;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        for (int h = 0; h < n; ++h) { physical[i] += t(i, j)*li(h, j)*transformed.x[h]; }
      }
    }
    physical -= reference.x;
    double error = physical.Norml2()/reference.x.Norml2();
    require(error < 1e-8 && reference.count == transformed.count,
            "Weighted GMRES dense basis invariance failed");
    report << pre << ',' << error << ',' << reference.count << ',' << transformed.count << '\n';
    Vector zero(n); zero = 0.;
    auto z = Solve(a, zero, 3, o, metric);
    require(z.count == 0 && z.x.Norml2() == 0, "Zero RHS solve failed");
  }
}

// One-sided Einfeldt/HLL branches and the exact homogeneous Euler PDE.
void VerifyPhysics() {
  State l{1, 5, -5, 25.5}, r{1.1, 5.5, -5.5, 28.};
  require(speeds(l, r, 0, 1.1)[0] == 0 && speeds(l, r, 1, 1.1)[1] == 0,
          "One-sided HLL branch failed");
  double eps = 1e-6, maxerr = 0;
  for (int i = 0; i < 20; i++) {
    double x = .07 * i, y = .13 * i, t = .07, a = pi / 6;
    auto ut1 = exact(x, y, t + eps, a), ut0 = exact(x, y, t - eps, a);
    auto fx1 = flux(exact(x + eps, y, t, a), 0), fx0 = flux(exact(x - eps, y, t, a), 0),
         fy1 = flux(exact(x, y + eps, t, a), 1), fy0 = flux(exact(x, y - eps, t, a), 1);
    for (int c = 0; c < 4; c++) {
      maxerr = std::max(
          maxerr,
          std::abs((ut1[c] - ut0[c] + fx1[c] - fx0[c] + fy1[c] - fy0[c]) / (2 * eps)));
    }
  }
  require(maxerr < 2e-8, "Exact Euler PDE verification failed");
  std::cout << "One-sided HLL and exact PDE verified; FD residual=" << maxerr
            << std::endl;
}
} // namespace euler
