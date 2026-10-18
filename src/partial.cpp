#include "euler.hpp"
namespace euler {
namespace {
// Quadrature-level partial assembly: cache flux Jacobians instead of the dense
// cell and face matrices. Shape/derivative tables are shared across rectangles.
class PartialDGOperator final : public ElementOperator {
  const Problem &p;
  Vector state;
  std::vector<Matrix4> jx, jy, jl, jr;
  mutable std::unique_ptr<ElementOperator> assembled_blocks;
  const ElementOperator &Blocks() const {
    if (!assembled_blocks) {
      assembled_blocks = p.DG(state);
    }
    return *assembled_blocks;
  }

public:
  PartialDGOperator(const Problem &problem, const Vector &u)
      : ElementOperator(problem.nu), p(problem), state(u) {
    const int nq = p.V.Height();
    jx.resize(p.n*nq); jy.resize(p.n*nq);
    jl.resize(p.nf*(p.opt.k+1)); jr.resize(p.nf*(p.opt.k+1));
  }
  void Refresh(const Vector &u) override {
    state = u;
    assembled_blocks.reset();
    const int nq = p.V.Height();
    for (int e = 0; e < p.n; ++e) {
      for (int q = 0; q < nq; ++q) {
        auto v = p.Eval(u, e, p.V, q);
        jacobian(v, 0, jx[e*nq+q]); jacobian(v, 1, jy[e*nq+q]);
      }
    }
    int id = 0;
    for (const auto &f : p.faces) {
      for (int q = 0; q <= p.opt.k; q++, ++id) {
        auto &l = jl[id]; auto &r = jr[id];
        jacobian(p.Eval(u, f.l, p.fv[f.ls], q), f.axis, l);
        jacobian(p.Eval(u, f.r, p.fv[f.rs], q), f.axis, r);
        double a = f.a[q], b = f.b[q];
        for (int c = 0; c < 4; c++) {
          for (int d = 0; d < 4; d++) {
            l(c, d) = b / (a + b) * l(c, d) + (c == d ? a * b / (a + b) : 0);
            r(c, d) = a / (a + b) * r(c, d) - (c == d ? a * b / (a + b) : 0);
          }
        }
      }
    }
  }
  // Analytic Jacobian action at fixed Newton state; no residual finite difference.
  // Apply the mass inverse after accumulating weak volume/face contributions.
  void Mult(const Vector &x, Vector &y) const override {
    y.SetSize(height);
    y = 0.;
    int nq = p.V.Height();
    for (int e = 0; e < p.n; e++) {
      for (int q = 0; q < nq; q++) {
        auto v = p.Eval(x, e, p.V, q);
        State fx{}, fy{};
        int id = e * nq + q;
        for (int c = 0; c < 4; c++) {
          for (int d = 0; d < 4; d++) {
            fx[c] += jx[id](c, d) * v[d];
            fy[c] += jy[id](c, d) * v[d];
          }
        }
        for (const auto &test : p.volume_tests[q]) {
          for (int c = 0; c < 4; ++c) {
            y[e*p.m+4*test.node+c] -= p.dt*p.w[q]*(test.dx*fx[c] + test.dy*fy[c]);
          }
        }
      }
    }
    int id = 0;
    for (const auto &f : p.faces) {
      for (int q = 0; q <= p.opt.k; q++, id++) {
        auto l = p.Eval(x, f.l, p.fv[f.ls], q), r = p.Eval(x, f.r, p.fv[f.rs], q);
        State z{};
        for (int c = 0; c < 4; c++) {
          for (int d = 0; d < 4; d++) {
            z[c] += jl[id](c, d) * l[d] + jr[id](c, d) * r[d];
          }
        }
        double scale = p.dt * f.length / p.area * p.ew[q];
        for (const auto &entry : p.face_basis[f.ls][q]) {
          for (int c = 0; c < 4; ++c) { y[f.l*p.m+4*entry.first+c] += scale*entry.second*z[c]; }
        }
        for (const auto &entry : p.face_basis[f.rs][q]) {
          for (int c = 0; c < 4; ++c) { y[f.r*p.m+4*entry.first+c] -= scale*entry.second*z[c]; }
        }
      }
    }
    p.MassSolve(y);
    y += x;
  }
  std::unique_ptr<mfem::SparseMatrix> Assemble() const override {
    return Blocks().Assemble();
  }
  std::vector<DenseMatrix> DiagonalBlocks(int size) const override {
    // Build only the DG diagonal blocks; never build off-diagonal cell blocks.
    std::vector<DenseMatrix> diag;
    DiagonalBlocks(size, diag);
    return diag;
  }
  void DiagonalBlocks(int size, std::vector<DenseMatrix> &diag) const override {
    require(size == p.m && p.collocated, "Partial DG expects collocated cell blocks");
    p.VolumeJacobian(state,diag,&jx,&jy);
    int id = 0;
    for (const auto &f : p.faces) {
      for (int q = 0; q <= p.opt.k; q++, id++) {
        double scale = p.dt * f.length / p.area * p.ew[q];
        for (const auto &i : p.face_basis[f.ls][q]) {
          for (const auto &j : p.face_basis[f.ls][q]) {
            for (int c = 0; c < 4; ++c) {
              for (int d = 0; d < 4; ++d) {
                diag[f.l](4*i.first+c,4*j.first+d) +=
                    scale*p.inv_mass(i.first,i.first)*i.second*j.second*jl[id](c,d);
              }
            }
          }
        }
        for (const auto &i : p.face_basis[f.rs][q]) {
          for (const auto &j : p.face_basis[f.rs][q]) {
            for (int c = 0; c < 4; ++c) {
              for (int d = 0; d < 4; ++d) {
                diag[f.r](4*i.first+c,4*j.first+d) -=
                    scale*p.inv_mass(i.first,i.first)*i.second*j.second*jr[id](c,d);
              }
            }
          }
        }
      }
    }
  }
  DenseMatrix DiagonalBlock(int start, int size) const override {
    return Blocks().DiagonalBlock(start, size);
  }
};
} // namespace
std::unique_ptr<ElementOperator> Problem::PartialDG(const Vector &u, bool refresh) const {
  auto op = std::make_unique<PartialDGOperator>(*this, u);
  if (refresh) { op->Refresh(u); }
  return op;
}
} // namespace euler
