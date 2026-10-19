#include "euler.hpp"
#include <algorithm>
#include <iostream>
namespace euler {
void VerifyPhysics(const EulerGas &gas);
void VerifyWeighted(const Options &base);
void Verify(const Options &base, const EulerGas &gas) {
  VerifyPhysics(gas);
  VerifyWeighted(base);
  std::ofstream storage(base.output + "_storage.csv");
  storage << "k,dg_bsr_ilu_relative,hdg_bsr_ilu_relative,streamed_schur_relative,"
             "schur_diagonal_relative\n" << std::setprecision(17);
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
    Problem p(o, gas);
    Vector u = p.Initial();
    auto background = gas.exact(0, 0, 0, p.alpha);
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
    auto streamed = pa->Assemble();
    streamed->Mult(v,az);
    az -= z;
    require(az.Norml2()/z.Norml2()<1e-12,
            "Streaming tensor DG CSR differs from reference");
    for (int e=0; e<p.n; ++e) {
      for (const auto *coupling : {&h.B[e], &h.C[e]}) {
        DenseMatrix expanded;
        coupling->Dense(expanded);
        Vector input(coupling->Width()), actual, expected(coupling->Height());
        for (int i=0; i<input.Size(); ++i) { input[i]=std::sin(.19*(i+1)); }
        coupling->Mult(input,actual);
        expanded.Mult(input,expected);
        actual-=expected;
        require(actual.Norml2()<1e-12*std::max(1.,expected.Norml2()),
                "Tensor face action differs from expanded matrix");
      }
    }
    Vector pz(p.nu);
    pa->Mult(v, pz);
    pz -= z;
    double parel = pz.Norml2() / z.Norml2();
    Options po = o;
    po.assembly = "partial";
    Problem pp(po, gas);
    pp.Freeze(u, .001);
    auto ph = pp.HDG(u);
    auto ps = pp.Condense(ph, hr);
    auto check_ilu = [&](const ElementOperator &op,int block) {
      FixedBlockILU test(op.BlockGraph(block),block,"none");
      auto csr=op.Assemble();
      Vector input(op.Height()), actual(op.Height()), expected(op.Height());
      for(int i=0;i<input.Size();++i) { input[i]=std::sin(.17*(i+1)); }
      double maximum=0.;
      for(int repeat=0;repeat<2;++repeat) {
        // Refill overwritten factors: verifies reuse without graph reconstruction.
        op.FillBlocks(test); test.Factor(); test.Mult(input,actual);
        mfem::BlockILU native(*csr,block,mfem::BlockILU::Reordering::NONE,0);
        native.Mult(input,expected); actual-=expected;
        maximum=std::max(maximum,actual.Norml2()/expected.Norml2());
        input*=.7;
      }
      require(maximum<1e-10,"Direct BSR ILU differs from native fixed-order ILU");
      return maximum;
    };
    const double dg_ilu=check_ilu(*pa,p.m), hdg_ilu=check_ilu(*ps.op,pp.t);
    Hybrid streamed_h; Condensed streamed_s;
    pp.PrepareHybrid(streamed_h,streamed_s);
    pp.StreamCondense(u,hr,streamed_h,streamed_s);
    require(streamed_h.A.empty() && streamed_s.cell_matrix.Height()==pp.m,
            "Production HDG retained all cell matrices");
    const auto sd=streamed_s.op->DiagonalBlocks(pp.t);
    double diag_error=0.,diag_scale=0.;
    for(int f=0;f<pp.nf;++f) {
      auto reference=ps.op->DiagonalBlock(f*pp.t,pp.t);
      for(int j=0;j<pp.t;++j) { for(int i=0;i<pp.t;++i) {
        diag_error=std::max(diag_error,std::abs(sd[f](i,j)-reference(i,j)));
        diag_scale=std::max(diag_scale,std::abs(reference(i,j)));
      } }
    }
    require(diag_error<1e-12*std::max(1.,diag_scale),"Panel-only Schur diagonal differs");
    Vector sz(ps.op->Height()), ez(s.op->Height());
    for (int i = 0; i < sz.Size(); i++) {
      sz[i] = std::sin(.27 * i);
    }
    ps.op->Mult(sz, ez);
    Vector streamed_value(ez.Size());
    streamed_s.op->Mult(sz,streamed_value); streamed_value-=ez;
    const double streamed_error=streamed_value.Norml2()/ez.Norml2();
    require(streamed_error<1e-12,"Streaming Schur differs from reference");
    storage << k << ',' << dg_ilu << ',' << hdg_ilu << ',' << streamed_error
            << ',' << diag_error/std::max(1.,diag_scale) << std::endl;
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
void VerifyPhysics(const EulerGas &gas) {
  const double configured_gamma = gas.Gamma();
  for (double value : {1.4, 1.67, configured_gamma}) {
    const EulerGas test_gas(value);
    State u{1.2, .3, -.2, 3.};
    const double expected_pressure = (value-1.)*(3.-.5*(.3*.3+.2*.2)/1.2);
    require(std::abs(test_gas.pressure(u)-expected_pressure) < 1e-14,
            "Parameterized pressure failed");
    for (int axis = 0; axis < 2; ++axis) {
      Matrix4 j;
      test_gas.jacobian(u, axis, j);
      for (int c = 0; c < 4; ++c) {
        State up = u, um = u;
        up[c] += 1e-6; um[c] -= 1e-6;
        auto fp = test_gas.flux(up, axis), fm = test_gas.flux(um, axis);
        for (int d = 0; d < 4; ++d) {
          require(std::abs((fp[d]-fm[d])/2e-6-j(d,c)) < 1e-8,
                  "Parameterized flux Jacobian failed");
        }
      }
      const double sound = std::sqrt(value*expected_pressure/u[0]);
      const double velocity = u[axis+1]/u[0];
      auto bounds = test_gas.speeds(u,u,axis,1.1);
      require(std::abs(bounds[0]-1.1*(sound-velocity)) < 1e-12 &&
                  std::abs(bounds[1]-1.1*(sound+velocity)) < 1e-12,
              "Parameterized HLL sound speed failed");
    }
  }
  State l{1, 5, -5, 25.5}, r{1.1, 5.5, -5.5, 28.};
  require(gas.speeds(l, r, 0, 1.1)[0] == 0 && gas.speeds(l, r, 1, 1.1)[1] == 0,
          "One-sided HLL branch failed");
  double eps = 1e-6, maxerr = 0;
  for (int i = 0; i < 20; i++) {
    double x = .07 * i, y = .13 * i, t = .07, a = pi / 6;
    auto ut1 = gas.exact(x, y, t + eps, a), ut0 = gas.exact(x, y, t - eps, a);
    auto fx1 = gas.flux(gas.exact(x + eps, y, t, a), 0), fx0 = gas.flux(gas.exact(x - eps, y, t, a), 0),
         fy1 = gas.flux(gas.exact(x, y + eps, t, a), 1), fy0 = gas.flux(gas.exact(x, y - eps, t, a), 1);
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
