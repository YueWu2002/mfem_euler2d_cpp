#include "euler.hpp"
#include <algorithm>
namespace euler {
namespace {
DenseMatrix zeros(int r, int c) {
  DenseMatrix a(r, c);
  a = 0.;
  return a;
}
// MFEM mutable operator() writes into a; cppcheck lacks its external definition.
// cppcheck-suppress constParameterReference
void add_face(DenseMatrix &a, const DenseMatrix &row, const DenseMatrix &col, int q,
              const Matrix4 &j, double scale, int ro = 0, int co = 0) {
  for (int i = 0; i < row.Width(); i++) {
    for (int h = 0; h < col.Width(); h++) {
      if (row(q, i) != 0 && col(q, h) != 0) {
        for (int c = 0; c < 4; c++) {
          for (int d = 0; d < 4; d++) {
            a(ro + 4 * i + c, co + 4 * h + d) +=
                scale * row(q, i) * col(q, h) * j(c, d);
          }
        }
      }
    }
  }
}
} // namespace
Problem::Problem(const Options &o)
    : opt(o), alpha(o.alpha * pi / 180),
      lx(2 / std::cos(alpha)), ly(2 / std::sin(alpha)),
      hx(lx / o.nx), hy(ly / o.ny), area(hx * hy), n(o.nx * o.ny),
      p((o.k + 1) * (o.k + 1)), m(4 * p), t(4 * (o.k + 1)), nf(2 * n), nu(n * m),
      nl(nf * t), mesh(mfem::Mesh::MakeCartesian2D(
                      o.nx, o.ny, mfem::Element::QUADRILATERAL, true, lx, ly, false)),
      fec(o.k, 2, mfem::BasisType::GaussLegendre),
      trace_fec(o.k, 1, mfem::BasisType::GaussLegendre),
      fes(&mesh, &fec), fe(fes.GetFE(0)) {
  old.SetSize(nu);
  require(fe->GetDof() == p, "Expected native Qk element");
  require(o.alpha > 0 && o.alpha < 90, "Angle must be strictly between 0 and 90 degrees");
  // Shift the actual MFEM mesh, including exported geometry, to the requested box.
  for (int v = 0; v < mesh.GetNV(); ++v) {
    mesh.GetVertex(v)[0] -= .5 * lx;
    mesh.GetVertex(v)[1] -= .5 * ly;
  }
  // Geometry is obtained from MFEM, not assumed to follow its element ordering.
  origx.resize(n);
  origy.resize(n);
  std::vector<int> grid(n);
  mfem::IntegrationPoint center;
  center.Set2(.5, .5);
  Vector pt(2);
  for (int e = 0; e < n; e++) {
    auto *tr = mesh.GetElementTransformation(e);
    tr->SetIntPoint(&center);
    tr->Transform(center, pt);
    const auto &j = tr->Jacobian();
    require(std::abs(j(0, 0) - hx) < 1e-12 * std::max(1., hx) &&
                std::abs(j(1, 1) - hy) < 1e-12 * std::max(1., hy) &&
                std::abs(j(0, 1)) + std::abs(j(1, 0)) < 1e-12 * std::max({1., hx, hy}),
            "Affine axis-aligned rectangles required");
    origx[e] = pt[0] - .5 * hx;
    origy[e] = pt[1] - .5 * hy;
    int ix = std::lround((pt[0] + .5 * lx) / hx - .5), iy = std::lround((pt[1] + .5 * ly) / hy - .5);
    grid[iy * o.nx + ix] = e;
  }
  cell_faces.resize(n);
  for (auto &a : cell_faces) {
    a.fill(-1);
  }
  for (int iy = 0; iy < o.ny; iy++) {
    for (int ix = 0; ix < o.nx; ix++) {
      int l = grid[iy * o.nx + ix];
      Face fx{l, grid[iy * o.nx + (ix + 1) % o.nx], 1, 0, 0, hy, {}, {}},
          fy{l, grid[((iy + 1) % o.ny) * o.nx + ix], 3, 2, 1, hx, {}, {}};
      for (Face f : {fx, fy}) {
        int id = faces.size();
        f.a.resize(o.k+1); f.b.resize(o.k+1);
        faces.push_back(f);
        cell_faces[f.l][f.ls] = id;
        cell_faces[f.r][f.rs] = id;
      }
    }
  }
  int order = o.volume_order < 0 ? 2 * o.k + 1 : o.volume_order;
  require(order >= 2 * o.k, "Volume quadrature exactness must be at least 2k");
  const auto &ir = mfem::IntRules.Get(mfem::Geometry::SQUARE, order);
  int nq = ir.GetNPoints();
  if (o.volume_order < 0) { require(nq == p, "Exactly (k+1)^2 GL volume points required"); }
  V.SetSize(nq, p);
  Gx.SetSize(nq, p);
  Gy.SetSize(nq, p);
  mass = zeros(p, p);
  Vector sh(p);
  DenseMatrix ds(p, 2);
  for (int q = 0; q < nq; q++) {
    auto &ip = ir.IntPoint(q);
    require(ip.weight > 0, "Positive volume weights required");
    w.push_back(ip.weight);
    fe->CalcShape(ip, sh);
    fe->CalcDShape(ip, ds);
    for (int i = 0; i < p; i++) {
      V(q, i) = (o.volume_order < 0 && std::abs(sh[i]) < 1e-14) ? 0 : sh[i];
      Gx(q, i) = (o.volume_order < 0 && std::abs(ds(i, 0)) < 1e-14) ? 0 : ds(i, 0) / hx;
      Gy(q, i) = (o.volume_order < 0 && std::abs(ds(i, 1)) < 1e-14) ? 0 : ds(i, 1) / hy;
      for (int j = 0; j < p; j++) {
        mass(i, j) += ip.weight * sh[i] * sh[j];
      }
    }
  }
  mfem::DenseMatrixInverse mi(mass);
  mi.GetInverseMatrix(inv_mass);
  if (o.volume_order < 0) {
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        if (i != j) {
          require(std::abs(mass(i, j)) < 1e-12, "Collocation mass not diagonal");
          mass(i, j) = 0;
          inv_mass(i, j) = 0;
        }
      }
    }
  }
  const auto &edge = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2 * o.k + 1);
  require(edge.GetNPoints() == o.k + 1, "Exactly k+1 GL edge points required");
  const auto *trace_fe = trace_fec.FiniteElementForGeometry(mfem::Geometry::SEGMENT);
  require(trace_fe->GetDof() == o.k+1, "Expected native Pk trace element");
  Vector trace_shape(o.k+1);
  psi.SetSize(o.k + 1, o.k + 1);
  for (auto &a : fv) {
    a.SetSize(o.k + 1, p);
  }
  for (int q = 0; q <= o.k; q++) {
    const auto &ip = edge.IntPoint(q);
    ew.push_back(ip.weight);
    ex.push_back(ip.x);
    trace_fe->CalcShape(ip, trace_shape);
    for (int i = 0; i <= o.k; ++i) {
      const double nodal_value = q == i ? 1. : 0.;
      require(std::abs(trace_shape[i]-nodal_value) < 1e-12,
              "Native Pk trace/GL quadrature ordering mismatch");
      psi(q,i) = nodal_value;
    }
    for (int s = 0; s < 4; s++) {
      mfem::IntegrationPoint fp;
      if (s < 2) {
        fp.Set2(s, ip.x);
      } else {
        fp.Set2(ip.x, s - 2);
      }
      fe->CalcShape(fp, sh);
      for (int i = 0; i < p; i++) {
        fv[s](q, i) = std::abs(sh[i]) < 1e-14 ? 0 : sh[i];
      }
    }
  }
  volume_tests.resize(nq);
  for (int q = 0; q < nq; ++q) {
    for (int i = 0; i < p; ++i) {
      if (Gx(q,i) != 0 || Gy(q,i) != 0) {
        volume_tests[q].push_back({i, Gx(q,i), Gy(q,i)});
      }
    }
  }
  collocated = o.volume_order < 0;
  volume_node.resize(nq, -1);
  if (collocated) {
    for (int q = 0; q < nq; ++q) {
      for (int i = 0; i < p; ++i) {
        if (std::abs(V(q, i) - 1.) < 1e-12) { volume_node[q] = i; }
      }
      require(volume_node[q] >= 0, "GL quadrature/nodal ordering mismatch");
    }
  }
  for (int side = 0; side < 4; ++side) {
    face_basis[side].resize(o.k + 1);
    for (int q = 0; q <= o.k; ++q) {
      for (int i = 0; i < p; ++i) {
        if (fv[side](q, i) != 0) { face_basis[side][q].push_back({i, fv[side](q, i)}); }
      }
    }
  }
  cell_metric.SetSize(nu);
  trace_metric.SetSize(nl);
  for (int e = 0; e < n; ++e) {
    for (int i = 0; i < p; ++i) {
      for (int c = 0; c < 4; ++c) { cell_metric[e*m+4*i+c] = area * mass(i,i); }
    }
  }
  // Stored cell equations are mass-inverted; stored skeleton equations have
  // physical length divided out. These metrics undo those row scalings.
  for (int f = 0; f < nf; ++f) {
    for (int q = 0; q <= o.k; ++q) {
      for (int c = 0; c < 4; ++c) { trace_metric[f*t+4*q+c] = faces[f].length / ew[q]; }
    }
  }
}
std::vector<int> Problem::Uids(int e) const {
  std::vector<int> r(m);
  for (int j = 0; j < m; j++) {
    r[j] = e * m + j;
  }
  return r;
}
std::vector<int> Problem::Lids(int f) const {
  std::vector<int> r(t);
  for (int j = 0; j < t; j++) {
    r[j] = f * t + j;
  }
  return r;
}
State Problem::Eval(const Vector &u, int e, const DenseMatrix &v, int q) const {
  State a{};
  if (collocated && &v == &V) {
    for (int c = 0; c < 4; ++c) { a[c] = u[e*m+4*volume_node[q]+c]; }
    return a;
  }
  for (int side = 0; side < 4; ++side) {
    if (&v == &fv[side]) {
      for (const auto &entry : face_basis[side][q]) {
        for (int c = 0; c < 4; ++c) { a[c] += entry.second * u[e*m+4*entry.first+c]; }
      }
      return a;
    }
  }
  for (int i = 0; i < p; i++) {
    for (int c = 0; c < 4; c++) {
      a[c] += v(q, i) * u[e * m + 4 * i + c];
    }
  }
  return a;
}
bool Problem::Admissible(const Vector &u) const {
  for (int e = 0; e < n; e++) {
    for (int q = 0; q < V.Height(); q++) {
      if (!positive(Eval(u, e, V, q))) {
        return false;
      }
    }
    for (int s = 0; s < 4; s++) {
      for (int q = 0; q <= opt.k; q++) {
        if (!positive(Eval(u, e, fv[s], q))) {
          return false;
        }
      }
    }
  }
  return true;
}
void Problem::MassSolve(DenseMatrix &a) const {
  if (collocated) {
    for (int i = 0; i < p; ++i) {
      for (int c = 0; c < 4; ++c) {
        for (int j = 0; j < a.Width(); ++j) { a(4*i+c,j) *= inv_mass(i,i); }
      }
    }
    return;
  }
  DenseMatrix z(a.Height(), a.Width());
  z = 0.;
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++) {
      if (inv_mass(i, j) != 0) {
        for (int c = 0; c < 4; c++) {
          for (int h = 0; h < a.Width(); h++) {
            z(4 * i + c, h) += inv_mass(i, j) * a(4 * j + c, h);
          }
        }
      }
    }
  }
  a = z;
}
void Problem::MassSolve(Vector &r) const {
  if (collocated) {
    for (int e = 0; e < n; ++e) {
      for (int i = 0; i < p; ++i) {
        for (int c = 0; c < 4; ++c) { r[e*m+4*i+c] *= inv_mass(i,i); }
      }
    }
    return;
  }
  Vector z(r.Size());
  z = 0.;
  for (int e = 0; e < n; e++) {
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        if (inv_mass(i, j) != 0) {
          for (int c = 0; c < 4; c++) {
            z[e * m + 4 * i + c] += inv_mass(i, j) * r[e * m + 4 * j + c];
          }
        }
      }
    }
  }
  r = z;
}
Vector Problem::Initial() const {
  Vector u(nu);
  u = 0.;
  const auto &ir =
      mfem::IntRules.Get(mfem::Geometry::SQUARE, std::max(24, 2 * opt.k + 8));
  Vector sh(p);
  for (int e = 0; e < n; e++) {
    for (int q = 0; q < ir.GetNPoints(); q++) {
      auto &ip = ir.IntPoint(q);
      fe->CalcShape(ip, sh);
      auto v = exact(origx[e] + hx * ip.x, origy[e] + hy * ip.y, 0, alpha);
      for (int i = 0; i < p; i++) {
        for (int c = 0; c < 4; c++) {
          u[e * m + 4 * i + c] += ip.weight * sh[i] * v[c];
        }
      }
    }
  }
  MassSolve(u);
  return u;
}
void Problem::Freeze(const Vector &u, double step) {
  ++freeze_generation;
  require(Admissible(u), "Old state is inadmissible; refine mesh");
  old = u;
  dt = step;
  for (auto &f : faces) {
    f.a.resize(opt.k + 1);
    f.b.resize(opt.k + 1);
    for (int q = 0; q <= opt.k; q++) {
      auto s = speeds(Eval(u, f.l, fv[f.ls], q), Eval(u, f.r, fv[f.rs], q), f.axis,
                      opt.safety);
      f.a[q] = s[0];
      f.b[q] = s[1];
      require(s[0] + s[1] > 0, "Degenerate HLL speeds");
    }
  }
}
double Problem::WaveRate() const {
  double speed = 0.;
  for (const auto &f : faces) {
    for (int q = 0; q <= opt.k; ++q) { speed = std::max({speed, f.a[q], f.b[q]}); }
  }
  // CFL = max_{faces,q}(|S_L|,|S_R|) dt / min(hx,hy), including safety.
  return speed / std::min(hx, hy);
}
double Problem::CellNorm(const Vector &r) const {
  require(r.Size() == nu, "Cell norm size mismatch");
  double norm = 0.;
  for (int i = 0; i < nu; ++i) { norm += cell_metric[i] * r[i] * r[i]; }
  return std::sqrt(norm);
}
Vector Problem::Trace(const Vector &u) const {
  Vector lam(nl);
  lam = 0.;
  for (int f = 0; f < nf; f++) {
    for (int q = 0; q <= opt.k; q++) {
      const auto &fc = faces[f];
      auto l = Eval(u, fc.l, fv[fc.ls], q), r = Eval(u, fc.r, fv[fc.rs], q),
           fl = flux(l, fc.axis), fr = flux(r, fc.axis);
      for (int i = 0; i <= opt.k; i++) {
        for (int c = 0; c < 4; c++) {
          lam[f * t + 4 * i + c] += psi(q, i) *
                                    (fc.a[q] * l[c] + fc.b[q] * r[c] + fl[c] - fr[c]) /
                                    (fc.a[q] + fc.b[q]);
        }
      }
    }
  }
  return lam;
}
Vector Problem::Skeleton(const Vector &u, const Vector &lam) const {
  Vector g(nl);
  g = 0.;
  for (int f = 0; f < nf; f++) {
    for (int q = 0; q <= opt.k; q++) {
      const auto &fc = faces[f];
      auto l = Eval(u, fc.l, fv[fc.ls], q), r = Eval(u, fc.r, fv[fc.rs], q),
           fl = flux(l, fc.axis), fr = flux(r, fc.axis);
      State lv{};
      for (int i = 0; i <= opt.k; i++) {
        for (int c = 0; c < 4; c++) {
          lv[c] += psi(q, i) * lam[f * t + 4 * i + c];
        }
      }
      for (int i = 0; i <= opt.k; i++) {
        for (int c = 0; c < 4; c++) {
          g[f * t + 4 * i + c] += ew[q] * psi(q, i) *
                                  ((fc.a[q] + fc.b[q]) * lv[c] - fc.a[q] * l[c] -
                                   fc.b[q] * r[c] - fl[c] + fr[c]);
        }
      }
    }
  }
  return g;
}
Vector Problem::Residual(const Vector &u, const Vector *lam) const {
  Vector r(nu);
  Residual(u, r, lam);
  return r;
}
void Problem::Residual(const Vector &u, Vector &r, const Vector *lam) const {
  r.SetSize(nu);
  r = 0.;
  for (int e = 0; e < n; e++) {
    for (int q = 0; q < V.Height(); q++) {
      auto v = Eval(u, e, V, q), fx = flux(v, 0), fy = flux(v, 1);
      for (const auto &test : volume_tests[q]) {
        for (int c = 0; c < 4; c++) {
          r[e*m+4*test.node+c] -= dt*w[q]*(test.dx*fx[c] + test.dy*fy[c]);
        }
      }
    }
  }
  for (int f = 0; f < nf; f++) {
    for (int q = 0; q <= opt.k; q++) {
      const auto &fc = faces[f];
      auto l = Eval(u, fc.l, fv[fc.ls], q), v = Eval(u, fc.r, fv[fc.rs], q),
           fl = flux(l, fc.axis), fr = flux(v, fc.axis);
      State hl{}, hr{}, lv{};
      if (lam) {
        for (int i = 0; i <= opt.k; i++) {
          for (int c = 0; c < 4; c++) {
            lv[c] += psi(q, i) * (*lam)[f * t + 4 * i + c];
          }
        }
      }
      for (int c = 0; c < 4; c++) {
        if (lam) {
          hl[c] = fl[c] + fc.a[q] * (l[c] - lv[c]);
          hr[c] = -fr[c] + fc.b[q] * (v[c] - lv[c]);
        } else {
          hl[c] =
              (fc.b[q] * fl[c] + fc.a[q] * fr[c] + fc.a[q] * fc.b[q] * (l[c] - v[c])) /
              (fc.a[q] + fc.b[q]);
          hr[c] = -hl[c];
        }
      }
      double scale = dt * fc.length / area * ew[q];
      for (const auto &entry : face_basis[fc.ls][q]) {
        for (int c = 0; c < 4; ++c) { r[fc.l*m+4*entry.first+c] += scale*entry.second*hl[c]; }
      }
      for (const auto &entry : face_basis[fc.rs][q]) {
        for (int c = 0; c < 4; ++c) { r[fc.r*m+4*entry.first+c] += scale*entry.second*hr[c]; }
      }
    }
  }
  MassSolve(r);
  r += u;
  r -= old;
}
std::vector<DenseMatrix> Problem::VolumeJacobian(const Vector &u) const {
  std::vector<DenseMatrix> a;
  VolumeJacobian(u, a);
  return a;
}
void Problem::VolumeJacobian(const Vector &u, std::vector<DenseMatrix> &a,
                             const std::vector<Matrix4> *cache_x,
                             const std::vector<Matrix4> *cache_y) const {
  a.resize(n);
  for (int e = 0; e < n; e++) {
    a[e].SetSize(m); a[e] = 0.;
    for (int q = 0; q < V.Height(); q++) {
      if (!cache_x) {
        auto v = Eval(u,e,V,q);
        jacobian(v,0,jac_x); jacobian(v,1,jac_y);
      }
      const auto &jx = cache_x ? (*cache_x)[e*V.Height()+q] : jac_x;
      const auto &jy = cache_y ? (*cache_y)[e*V.Height()+q] : jac_y;
      for (const auto &test : volume_tests[q]) {
        const int i = test.node;
        const int begin = collocated ? volume_node[q] : 0;
        const int end = collocated ? begin + 1 : p;
        for (int j = begin; j < end; ++j) {
          if (V(q,j) == 0) { continue; }
          for (int c = 0; c < 4; ++c) {
            for (int d = 0; d < 4; ++d) {
              a[e](4*i+c,4*j+d) -= dt*w[q]*V(q,j)*
                  (test.dx*jx(c, d) + test.dy*jy(c, d));
            }
          }
        }
      }
    }
    MassSolve(a[e]);
    for (int i = 0; i < m; i++) {
      a[e](i, i) += 1;
    }
  }
}
std::unique_ptr<ElementOperator> Problem::DG(const Vector &u) const {
  auto op = std::make_unique<ElementOperator>(nu);
  auto vol = VolumeJacobian(u);
  for (int e = 0; e < n; e++) {
    op->Add(Uids(e), Uids(e), vol[e]);
  }
  for (const auto &f : faces) {
    DenseMatrix ll = zeros(m, m), lr = zeros(m, m), rl = zeros(m, m), rr = zeros(m, m);
    for (int q = 0; q <= opt.k; q++) {
      auto &jl = face_left; auto &jr = face_right;
      jacobian(Eval(u, f.l, fv[f.ls], q), f.axis, jl);
      jacobian(Eval(u, f.r, fv[f.rs], q), f.axis, jr);
      double a = f.a[q], b = f.b[q];
      for (int c = 0; c < 4; c++) {
        for (int d = 0; d < 4; d++) {
          jl(c, d) = b / (a + b) * jl(c, d) + (c == d ? a * b / (a + b) : 0);
          jr(c, d) = a / (a + b) * jr(c, d) - (c == d ? a * b / (a + b) : 0);
        }
      }
      double s = dt * f.length / area * ew[q];
      add_face(ll, fv[f.ls], fv[f.ls], q, jl, s);
      add_face(lr, fv[f.ls], fv[f.rs], q, jr, s);
      add_face(rl, fv[f.rs], fv[f.ls], q, jl, -s);
      add_face(rr, fv[f.rs], fv[f.rs], q, jr, -s);
    }
    MassSolve(ll);
    MassSolve(lr);
    MassSolve(rl);
    MassSolve(rr);
    op->Add(Uids(f.l), Uids(f.l), ll);
    op->Add(Uids(f.l), Uids(f.r), lr);
    op->Add(Uids(f.r), Uids(f.l), rl);
    op->Add(Uids(f.r), Uids(f.r), rr);
  }
  return op;
}
// Build cell-local blocks at the current Newton state, with frozen face speeds.
// The four face restrictions in ids[e] follow the cell side order.
Hybrid Problem::HDG(const Vector &u) const {
  Hybrid h;
  HDG(u, h);
  return h;
}
void Problem::HDG(const Vector &u, Hybrid &h) const {
  require(collocated, "HDG requires the specified GL collocation");
  VolumeJacobian(u, h.A);
  h.B.resize(n);
  h.C.resize(n);
  h.D.resize(nf);
  h.ids.resize(n);

  // B and D are constant throughout this implicit time step. Only A and C
  // contain current-state physical flux derivatives and need every-Newton work.
  const bool refresh_static = h.freeze_generation != freeze_generation;
  if (refresh_static) {
    for (int f = 0; f < nf; ++f) {
      h.D[f].SetSize(t);
      h.D[f] = 0.;
      for (int q = 0; q <= opt.k; ++q) {
        const double diagonal = ew[q] * (faces[f].a[q] + faces[f].b[q]);
        for (int c = 0; c < 4; ++c) {
          h.D[f](4*q+c, 4*q+c) = diagonal;
        }
      }
    }
  }

  for (int e = 0; e < n; ++e) {
    h.B[e].SetSize(m, 4*t);
    h.C[e].SetSize(4*t, m);
    if (refresh_static) {
      h.B[e] = 0.;
    }
    h.C[e] = 0.;
    if (h.ids[e].empty()) {
      h.ids[e].reserve(4*t);
      for (int side = 0; side < 4; ++side) {
        const int f = cell_faces[e][side];
        for (int i = 0; i < t; ++i) {
          h.ids[e].push_back(f*t+i);
        }
      }
    }

    for (int side = 0; side < 4; ++side) {
      const auto &face = faces[cell_faces[e][side]];
      const bool left = face.l == e && face.ls == side;
      const double normal_sign = left ? 1. : -1.;
      for (int q = 0; q <= opt.k; ++q) {
        auto &derivative = face_left;
        jacobian(Eval(u, e, fv[side], q), face.axis, derivative);
        const double stabilization = left ? face.a[q] : face.b[q];
        for (int d = 0; d < 4; ++d) {
          for (int c = 0; c < 4; ++c) {
            derivative(c,d) = normal_sign*derivative(c,d) +
                              (c == d ? stabilization : 0.);
          }
        }

        // Trace basis is nodal: only column/node q contributes. Accumulate
        // inverse-mass-scaled A/B rows directly, without dense face scratch.
        const double scale = dt * face.length / area * ew[q];
        const int trace_node = side*t + 4*q;
        const auto &restriction = face_basis[side][q];
        for (const auto &test : restriction) {
          const int i = test.first;
          const double row_weight = scale*test.second*inv_mass(i,i);
          if (refresh_static) {
            for (int c = 0; c < 4; ++c) {
              h.B[e](4*i+c, trace_node+c) -= row_weight*stabilization;
            }
          }
          for (const auto &trial : restriction) {
            const int j = trial.first;
            const double value_weight = row_weight*trial.second;
            for (int d = 0; d < 4; ++d) {
              for (int c = 0; c < 4; ++c) {
                h.A[e](4*i+c, 4*j+d) += value_weight*derivative(c,d);
              }
            }
          }
        }
        // Skeleton rows have neither dt nor physical face length. The sign
        // matches g = g_R - g_L, so the trace penalty D is positive.
        for (const auto &trial : restriction) {
          const int j = trial.first;
          const double weight = -ew[q]*trial.second;
          for (int d = 0; d < 4; ++d) {
            for (int c = 0; c < 4; ++c) {
              h.C[e](trace_node+c, 4*j+d) += weight*derivative(c,d);
            }
          }
        }
      }
    }
  }
  h.freeze_generation = freeze_generation;
}
std::array<double, 4> Problem::Errors(const Vector &u, double time) const {
  std::array<double, 4> a{};
  const auto &ir =
      mfem::IntRules.Get(mfem::Geometry::SQUARE, std::max(24, 2 * opt.k + 8));
  Vector sh(p);
  for (int e = 0; e < n; e++) {
    for (int q = 0; q < ir.GetNPoints(); q++) {
      auto &ip = ir.IntPoint(q);
      fe->CalcShape(ip, sh);
      auto v = exact(origx[e] + hx * ip.x, origy[e] + hy * ip.y, time, alpha);
      for (int c = 0; c < 4; c++) {
        double z = 0;
        for (int i = 0; i < p; i++) {
          z += sh[i] * u[e * m + 4 * i + c];
        }
        a[c] += area * ip.weight * (z - v[c]) * (z - v[c]);
      }
    }
  }
  std::transform(a.begin(), a.end(), a.begin(),
                 [](double value) { return std::sqrt(value); });
  return a;
}
std::array<double, 4> Problem::Integral(const Vector &u) const {
  std::array<double, 4> a{};
  for (int e = 0; e < n; e++) {
    for (int q = 0; q < V.Height(); q++) {
      auto v = Eval(u, e, V, q);
      for (int c = 0; c < 4; c++) {
        a[c] += area * w[q] * v[c];
      }
    }
  }
  return a;
}
double Problem::Difference(const Vector &u, const Vector &v) const {
  double a = 0;
  for (int e = 0; e < n; e++) {
    for (int q = 0; q < V.Height(); q++) {
      auto x = Eval(u, e, V, q), y = Eval(v, e, V, q);
      a += area * w[q] * (x[0] - y[0]) * (x[0] - y[0]);
    }
  }
  return std::sqrt(a);
}
void Problem::Save(const Vector &u, const std::string &prefix) const {
  std::ofstream ms(prefix + ".mesh");
  mesh.Print(ms);
  const char *names[] = {"rho", "rhou", "rhov", "energy"};
  for (int c = 0; c < 4; c++) {
    mfem::GridFunction gf(const_cast<mfem::FiniteElementSpace *>(&fes));
    mfem::Array<int> ids;
    for (int e = 0; e < n; e++) {
      fes.GetElementDofs(e, ids);
      for (int i = 0; i < p; i++) {
        gf[ids[i]] = u[e * m + 4 * i + c];
      }
    }
    std::ofstream fs(prefix + "_" + names[c] + ".gf");
    fs << std::setprecision(17);
    gf.Save(fs);
  }
  std::ofstream cs(prefix + "_coefficients.csv");
  cs << "element,node,rho,rhou,rhov,energy\n" << std::setprecision(17);
  for (int e = 0; e < n; e++) {
    for (int i = 0; i < p; i++) {
      cs << e << ',' << i;
      for (int c = 0; c < 4; c++) {
        cs << ',' << u[e * m + 4 * i + c];
      }
      cs << '\n';
    }
  }
}
} // namespace euler
