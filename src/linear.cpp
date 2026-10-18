#include "euler.hpp"
#include <algorithm>
#include <limits>

namespace euler {
// The stored cell equation is M^{-1}r, whereas stored skeleton rows have
// face length divided out. 'metric' accounts for these row scalings, so dot()
// always measures the PHYSICAL weak residual's inverse-mass (dual) norm.
// With right preconditioning this is equivalent to the paired Gram transform
// L^{-1} K L^{-T} in note.pdf, section 3.2.1.
struct LinearSolver::Impl {
  int n, block, restart;
  Vector metric, residual, action, preconditioned, combination, correction;
  DenseMatrix basis, h;
  Vector cosine, sine, rotated_rhs, coefficients;
  std::vector<DenseMatrix> diagonal;
  std::vector<std::unique_ptr<mfem::DenseMatrixInverse>> diagonal_lu;
  std::unique_ptr<mfem::BlockILU> ilu;
  std::unique_ptr<mfem::SparseMatrix> assembled;
  DenseMatrix dense_matrix;
  std::unique_ptr<mfem::DenseMatrixInverse> dense_lu;
#ifdef MFEM_USE_SUITESPARSE
  std::unique_ptr<mfem::UMFPackSolver> direct;
#endif

  Impl(int size, int b, const Options &o, const Vector &weights)
      : n(size), block(b), restart(std::min(size,o.restart)), metric(weights),
        residual(n), action(n), preconditioned(n), combination(n), correction(n),
        basis(n,restart+1), h(restart+1,restart), cosine(restart), sine(restart),
        rotated_rhs(restart+1), coefficients(restart) {
    require(block > 0 && n % block == 0 && metric.Size() == n,
            "Invalid linear workspace dimensions");
    for (int i = 0; i < n; ++i) {
      require(std::isfinite(metric[i]) && metric[i] > 0, "Invalid dual residual metric");
    }
    if (o.solver == "gmres" && o.precond == "block-jacobi") {
      diagonal.resize(n/block);
      diagonal_lu.reserve(diagonal.size());
      for (auto &a : diagonal) {
        a.SetSize(block);
        // Allocate factor buffers, without factoring uninitialized storage.
        diagonal_lu.push_back(std::make_unique<mfem::DenseMatrixInverse>(&a));
      }
    }
    if (o.solver == "gmres" && o.precond == "block-ilu0") {
      auto order = o.ilu_ordering == "mdf"
          ? mfem::BlockILU::Reordering::MINIMUM_DISCARDED_FILL
          : mfem::BlockILU::Reordering::NONE;
      ilu = std::make_unique<mfem::BlockILU>(block,order,0);
    }
    if (o.solver == "dense") {
      require(n <= 5000, "Dense direct solve limited to 5000 unknowns");
      dense_matrix.SetSize(n);
      dense_lu = std::make_unique<mfem::DenseMatrixInverse>(&dense_matrix);
    }
#ifdef MFEM_USE_SUITESPARSE
    if (o.solver == "direct") { direct = std::make_unique<mfem::UMFPackSolver>(); }
#else
    require(o.solver != "direct", "--solver direct requires MFEM with SuiteSparse");
#endif
  }
  double Dot(const Vector &a, const Vector &b) const {
    double sum = 0.;
    const double *ap = a.GetData(), *bp = b.GetData(), *mp = metric.GetData();
    for (int i = 0; i < n; ++i) { sum += mp[i]*ap[i]*bp[i]; }
    return sum;
  }
  double Norm(const Vector &a) const { return std::sqrt(Dot(a,a)); }
  void Precondition(const Vector &b, Vector &x) const {
    if (ilu) { ilu->Mult(b,x); return; }
    // Non-owning block views: no allocation, gather or scatter in block Jacobi.
    for (int e = 0; e < n/block; ++e) {
      diagonal_lu[e]->Mult(b.GetData()+e*block,x.GetData()+e*block);
    }
  }
  void GMRES(const mfem::Operator &a, const Vector &rhs, double tolerance,
             int limit, LinearResult &out, bool zero_guess) {
    if (zero_guess) {
      out.x = 0.;
      residual = rhs;
    } else {
      a.Mult(out.x, action);
      residual = rhs;
      residual -= action;
    }
    const double rhs_norm = Norm(rhs);
    if (rhs_norm == 0.) { return; }
    int iterations = 0;
    while (iterations < limit) {
      const double beta = Norm(residual);
      if (beta <= tolerance*rhs_norm) { break; }
      Vector first(basis.GetData(),n);
      first = residual; first *= 1/beta;
      h = 0.; rotated_rhs = 0.; rotated_rhs[0] = beta;
      int count = 0;
      for (int j = 0; j < restart && iterations < limit; ++j) {
        Vector v(basis.GetData()+j*n,n);
        Precondition(v,preconditioned);
        a.Mult(preconditioned,action);
        // Twice-modified Gram-Schmidt controls loss of orthogonality without
        // allocating any Arnoldi vectors inside the Newton or Krylov loop.
        for (int pass = 0; pass < 2; ++pass) {
          for (int i = 0; i <= j; ++i) {
            Vector vi(basis.GetData()+i*n,n);
            const double value = Dot(vi,action);
            h(i,j) += value;
            action.Add(-value,vi);
          }
        }
        const double next = Norm(action);
        h(j+1,j) = next;
        if (next > 0) {
          Vector vnext(basis.GetData()+(j+1)*n,n);
          vnext = action; vnext *= 1/next;
        }
        for (int i = 0; i < j; ++i) {
          const double u = h(i,j), vrot = h(i+1,j);
          h(i,j) = cosine[i]*u + sine[i]*vrot;
          h(i+1,j) = -sine[i]*u + cosine[i]*vrot;
        }
        const double d = std::hypot(h(j,j),h(j+1,j));
        require(std::isfinite(d) && d > 0, "GMRES Arnoldi breakdown");
        cosine[j] = h(j,j)/d; sine[j] = h(j+1,j)/d;
        h(j,j) = d; h(j+1,j) = 0.;
        rotated_rhs[j+1] = -sine[j]*rotated_rhs[j];
        rotated_rhs[j] *= cosine[j];
        ++iterations; count = j+1;
        if (std::abs(rotated_rhs[j+1]) <= tolerance*rhs_norm ||
            next <= std::numeric_limits<double>::min()) { break; }
      }
      for (int i = count-1; i >= 0; --i) {
        double value = rotated_rhs[i];
        for (int j = i+1; j < count; ++j) { value -= h(i,j)*coefficients[j]; }
        coefficients[i] = value/h(i,i);
      }
      // P is fixed and linear during one solve: apply P^{-1} once to the
      // combined Arnoldi vector rather than storing restart full-sized Z's.
      combination = 0.;
      for (int j = 0; j < count; ++j) {
        Vector v(basis.GetData()+j*n,n);
        combination.Add(coefficients[j],v);
      }
      Precondition(combination,correction);
      out.x += correction;
      a.Mult(out.x,action);
      residual = rhs; residual -= action;
      require(std::isfinite(Norm(residual)), "Nonfinite GMRES residual");
    }
    out.count = iterations;
    require(Norm(residual) <= tolerance*rhs_norm, "GMRES iteration limit reached");
  }
};
LinearSolver::LinearSolver(int n, int block, const Options &o, const Vector &metric)
    : impl(std::make_unique<Impl>(n,block,o,metric)) {}
LinearSolver::~LinearSolver() = default;
void LinearSolver::Solve(const ElementOperator &op, const Vector &rhs,
                         const Options &o, LinearResult &out, bool refresh) {
  auto &w = *impl;
  out.x.SetSize(w.n);
  out.count = 0;
  out.residual = out.assembly = out.setup = out.solve = 0.;
  const double rhs_norm = w.Norm(rhs);
  if (rhs_norm == 0.) {
    out.x = 0.;
    return;
  }

  double start = clock_s();
  if (refresh && (!o.matrix_free || o.solver != "gmres" || w.ilu)) {
    w.assembled = op.Assemble();
  }
  out.assembly = clock_s() - start;

  start = clock_s();
  if (refresh) {
    if (o.solver == "gmres") {
      if (w.ilu) {
        w.ilu->SetOperator(*w.assembled);
      } else {
        op.DiagonalBlocks(w.block, w.diagonal);
        for (auto &lu : w.diagonal_lu) {
          lu->Factor();
        }
      }
    } else if (o.solver == "dense") {
      w.assembled->ToDenseMatrix(w.dense_matrix);
      w.dense_lu->Factor();
    } else {
#ifdef MFEM_USE_SUITESPARSE
      w.direct->SetOperator(*w.assembled);
#endif
    }
  }
  out.setup = clock_s() - start;

  start = clock_s();
  if (o.solver == "gmres") {
    const mfem::Operator &a = o.matrix_free
        ? static_cast<const mfem::Operator &>(op)
        : static_cast<const mfem::Operator &>(*w.assembled);
    w.GMRES(a, rhs, o.linear_tol, o.max_krylov, out, refresh);
  } else if (o.solver == "dense") {
    w.dense_lu->Mult(rhs, out.x);
  } else {
#ifdef MFEM_USE_SUITESPARSE
    w.direct->Mult(rhs, out.x);
#endif
  }
  out.solve = clock_s() - start;

  // Always validate with the requested unpreconditioned operator and physical
  // dual residual metric, including when a global CSR was used by GMRES.
  op.Mult(out.x, w.action);
  w.action -= rhs;
  out.residual = w.Norm(w.action) / rhs_norm;
  require(std::isfinite(out.residual) &&
              out.residual <= std::max(10*o.linear_tol, 1e-12),
          "True inverse-mass weighted linear residual too large");
}
LinearResult Solve(const ElementOperator &op, const Vector &rhs, int block,
                   const Options &o, const Vector &metric) {
  // Convenient one-shot entry point for independent verification/reference code.
  // The time integrator instead reuses one LinearSolver for the whole run.
  LinearSolver solver(rhs.Size(),block,o,metric);
  LinearResult result;
  solver.Solve(op,rhs,o,result);
  return result;
}
} // namespace euler
