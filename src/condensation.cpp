#include "euler.hpp"
#include <algorithm>

namespace euler {
namespace {
// Homogeneous trace-to-cell response: v -> -A_e^{-1} B_e R_e v.
// Borrow pointers to vector elements (stable across Condensed moves), not to
// the movable vector objects themselves. Hybrid and Condensed must outlive it.
class RecoveryOperator final : public mfem::Operator {
  struct CellResponse {
    const std::vector<int> *trace_ids;
    const DenseMatrix *response;
  };
  std::vector<CellResponse> cells;
  int cell_size, trace_size;
  mutable Vector local_trace, local_cell;

public:
  RecoveryOperator(const Problem &p, const Hybrid &h, const Condensed &s)
      : Operator(p.nu, p.nl), cell_size(p.m), trace_size(4 * p.t),
        local_trace(trace_size), local_cell(cell_size) {
    cells.reserve(p.n);
    for (int e = 0; e < p.n; ++e) {
      cells.push_back({h.ids.data() + e, s.X.data() + e});
    }
  }

  void Mult(const Vector &trace_increment, Vector &cell_increment) const override {
    cell_increment.SetSize(height);
    for (int e = 0; e < static_cast<int>(cells.size()); ++e) {
      const auto &cell = cells[e];
      for (int i = 0; i < trace_size; ++i) {
        local_trace[i] = trace_increment[(*cell.trace_ids)[i]];
      }
      cell.response->Mult(local_trace, local_cell);
#pragma omp simd
      for (int i = 0; i < cell_size; ++i) {
        cell_increment[e * cell_size + i] = -local_cell[i];
      }
    }
  }
};
} // namespace

// Allocate topology, cached responses and one cell factor workspace before
// the time loop. Streaming numerical assembly/condensation stays in the timer.
void Problem::PrepareHybrid(Hybrid &h, Condensed &s) const {
  h.B.resize(n); h.C.resize(n); h.D.SetSize(nl); h.ids.resize(n);
  s.X.resize(n); s.y.resize(n);
  s.rhs.SetSize(nl); s.cell_scratch.SetSize(m); s.trace_scratch.SetSize(4*t);
  s.local_scratch.SetSize(4*t); s.scaled_scratch.SetSize(nl);
  s.product_scratch.SetSize(4*t);
  s.cell_matrix.SetSize(m); s.cell_pivots.resize(m);
  for (int e = 0; e < n; ++e) {
    for (auto *coupling : {&h.B[e], &h.C[e]}) {
      coupling->basis = &face_basis;
      coupling->inverse_mass = &inv_mass;
      coupling->nodes = opt.k+1;
      coupling->cell_size = m;
      coupling->flux.resize(4*(opt.k+1));
    }
    h.B[e].lift = true; h.C[e].lift = false;
    s.X[e].SetSize(m,4*t); s.y[e].SetSize(m);
    h.ids[e].reserve(4*t);
    for (int side = 0; side < 4; ++side) {
      int f = cell_faces[e][side];
      for (int i = 0; i < t; ++i) { h.ids[e].push_back(f*t+i); }
    }
  }
  s.op = std::make_unique<ElementOperator>(nl);
  s.op->entries.reserve(opt.assembly == "partial" ? 0 : n);
  s.op->factors.reserve(opt.assembly == "partial" ? n : 0);
  h.D=0.; s.op->diagonal.SetSize(nl); s.op->diagonal=0.;
  for (int e = 0; e < n; ++e) {
    if (opt.assembly == "partial") {
      s.op->factors.push_back({h.ids[e], &h.C[e], &s.X[e], -1});
    } else {
      s.product_scratch = 0.;
      s.op->Add(h.ids[e], h.ids[e], s.product_scratch);
    }
  }
  s.op->PrepareScratch(4*t, 4*t, m);
  s.recovery = std::make_unique<RecoveryOperator>(*this, h, s);
}
Condensed Problem::Condense(const Hybrid &h, const Vector &r, const Vector *g) const {
  // Convenience path for isolated verification calls. Production uses the
  // preallocated overload below and never replaces its factor/trace storage.
  Condensed s;
  s.X.resize(n); s.y.resize(n); s.rhs.SetSize(nl);
  s.cell_scratch.SetSize(m); s.trace_scratch.SetSize(4*t);
  s.product_scratch.SetSize(4*t);
  s.op = std::make_unique<ElementOperator>(nl);
  s.op->diagonal=h.D;
  for (int e = 0; e < n; ++e) {
    s.X[e].SetSize(m,4*t); s.y[e].SetSize(m);
    s.pivots.emplace_back(m);
    if (opt.assembly == "partial") {
      s.op->factors.push_back({h.ids[e], &h.C[e], &s.X[e], -1});
    } else {
      s.product_scratch = 0.; s.op->Add(h.ids[e],h.ids[e],s.product_scratch);
    }
  }
  s.op->PrepareScratch(4*t,4*t,m);
  s.recovery = std::make_unique<RecoveryOperator>(*this,h,s);
  Hybrid work = h; // Reference path preserves the caller's unfactored matrices.
  Condense(work,r,s,g);
  return s;
}
void Problem::Condense(Hybrid &h, const Vector &r, Condensed &s,
                       const Vector *g) const {
  s.rhs = 0.;
  if (g) { s.rhs = *g; s.rhs *= -1.; }
  s.op->diagonal=h.D;
  for (int e = 0; e < n; ++e) {
    mfem::LUFactors lu(h.A[e].GetData(),s.pivots[e].data());
    require(lu.Factor(m), "Singular HDG cell factor");
    h.B[e].Dense(s.X[e]);
    lu.Solve(m,4*t,s.X[e].GetData());
    Vector re(const_cast<double *>(r.GetData()) + e*m,m);
    s.y[e] = re;
    lu.Solve(m,1,s.y[e].GetData());
    h.C[e].Mult(s.y[e],s.trace_scratch);
    for (int i = 0; i < 4*t; ++i) { s.rhs[h.ids[e][i]] += s.trace_scratch[i]; }
    if (opt.assembly != "partial") {
      auto &cx = s.op->entries[e].block;
      h.C[e].Mult(s.X[e],cx); cx *= -1.;
    }
  }
}
void Problem::StreamCondense(const Vector &u,const Vector &r,Hybrid &h,Condensed &s) const {
  HybridDiagonal(h);
  s.op->diagonal=h.D; s.rhs=0.;
  for(int e=0;e<n;++e) {
    HybridCell(u,e,h,s.cell_matrix);
    mfem::LUFactors lu(s.cell_matrix.GetData(),s.cell_pivots.data());
    require(lu.Factor(m),"Singular streaming HDG cell factor");
    h.B[e].Dense(s.X[e]); lu.Solve(m,4*t,s.X[e].GetData());
    std::copy_n(r.GetData()+e*m,m,s.y[e].GetData());
    lu.Solve(m,1,s.y[e].GetData());
    h.C[e].Mult(s.y[e],s.trace_scratch);
    for(int i=0;i<4*t;++i) { s.rhs[h.ids[e][i]]+=s.trace_scratch[i]; }
    if(opt.assembly!="partial") {
      auto &cx=s.op->entries[e].block;
      h.C[e].Mult(s.X[e],cx); cx*=-1.;
    }
  }
}
Vector Problem::Recover(const Condensed &s, const Vector &dl) const {
  Vector du(nu);
  Recover(s,dl,du);
  return du;
}
void Problem::Recover(const Condensed &s, const Vector &dl, Vector &du) const {
  s.recovery->Mult(dl, du);
  // Recovery is affine: add the particular response -A_e^{-1} r_e.
  for (int e = 0; e < n; ++e) {
#pragma omp simd
    for (int i = 0; i < m; ++i) {
      du[e * m + i] -= s.y[e][i];
    }
  }
}
} // namespace euler
