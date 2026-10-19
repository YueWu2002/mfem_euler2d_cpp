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
  Vector metric, residual, action, preconditioned, combination, correction,
      weighted_action, projection;
  DenseMatrix basis, h;
  Vector cosine, sine, rotated_rhs, coefficients;
  std::vector<DenseMatrix> diagonal;
  std::vector<std::vector<int>> diagonal_pivots;
  std::unique_ptr<FixedBlockILU> ilu;
  std::unique_ptr<mfem::SparseMatrix> assembled;
  DenseMatrix dense_matrix;
  std::unique_ptr<mfem::DenseMatrixInverse> dense_lu;
#ifdef MFEM_USE_MKL_PARDISO
  std::unique_ptr<mfem::PardisoSolver> direct;
#endif

  Impl(int size, int b, const Options &o, const Vector &weights)
      : n(size), block(b), restart(std::min(size,o.restart)), metric(weights),
        residual(n), action(n), preconditioned(n), combination(n), correction(n),
        weighted_action(n), projection(restart+1),
        basis(n,restart+1), h(restart+1,restart), cosine(restart), sine(restart),
        rotated_rhs(restart+1), coefficients(restart) {
    require(block > 0 && n % block == 0 && metric.Size() == n,
            "Invalid linear workspace dimensions");
    for (int i = 0; i < n; ++i) {
      require(std::isfinite(metric[i]) && metric[i] > 0, "Invalid dual residual metric");
    }
    if (o.solver == "gmres" && o.precond == "block-jacobi") {
      diagonal.resize(n/block);
      diagonal_pivots.resize(diagonal.size(), std::vector<int>(block));
      for (auto &a : diagonal) {
        a.SetSize(block);

      }
    }
    if (o.solver == "dense") {
      require(n <= 5000, "Dense direct solve limited to 5000 unknowns");
      dense_matrix.SetSize(n);
      dense_lu = std::make_unique<mfem::DenseMatrixInverse>(&dense_matrix);
    }
#ifdef MFEM_USE_MKL_PARDISO
    // Construct on numerical refresh below: MFEM's PardisoSolver owns the
    // symbolic/numeric factors for exactly one SparseMatrix.
#else
    require(o.solver != "direct", "--solver direct requires MFEM with MKL PARDISO");
#endif
  }
  double Dot(const Vector &a, const Vector &b) const {
    double sum = 0.;
    const double *ap = a.GetData(), *bp = b.GetData(), *mp = metric.GetData();
#pragma omp simd reduction(+:sum)
    for (int i = 0; i < n; ++i) { sum += mp[i]*ap[i]*bp[i]; }
    return sum;
  }
  double Norm(const Vector &a) const { return std::sqrt(Dot(a,a)); }
  void Precondition(const Vector &b, Vector &x) const {
    if (ilu) { ilu->Mult(b,x); return; }
    // Non-owning block views: no allocation, gather or scatter in block Jacobi.
    for (int e = 0; e < n/block; ++e) {
      std::copy_n(b.GetData()+e*block,block,x.GetData()+e*block);
      mfem::LUFactors lu(const_cast<double*>(diagonal[e].GetData()),
                         const_cast<int*>(diagonal_pivots[e].data()));
      lu.Solve(block,1,x.GetData()+e*block);
    }
  }
  double GMRES(const mfem::Operator &a, const Vector &rhs, double tolerance,
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
    if (rhs_norm == 0.) { return 0.; }
    int iterations = 0;
    double residual_norm = Norm(residual);
    while (iterations < limit) {
      const double beta = residual_norm;
      if (beta <= tolerance*rhs_norm) { break; }
      Vector first(basis.GetData(),n);
      first = residual; first *= 1/beta;
      h = 0.; rotated_rhs = 0.; rotated_rhs[0] = beta;
      int count = 0;
      for (int j = 0; j < restart && iterations < limit; ++j) {
        Vector v(basis.GetData()+j*n,n);
        Precondition(v,preconditioned);
        a.Mult(preconditioned,action);
        // Twice classical Gram-Schmidt (CGS2) retains reorthogonalization while
        // replacing j separate full-vector updates with two dense block actions.
        // This substantially reduces writes to the large Arnoldi work vector.
        DenseMatrix active_basis(basis.GetData(), n, j+1);
        for (int pass = 0; pass < 2; ++pass) {
          const double *ap = action.GetData(), *mp = metric.GetData();
          double *wp = weighted_action.GetData();
#pragma omp simd
          for (int k = 0; k < n; ++k) { wp[k] = mp[k]*ap[k]; }
          active_basis.MultTranspose(weighted_action.GetData(), projection.GetData());
          for (int i = 0; i <= j; ++i) {
            h(i,j) += projection[i];
          }
          active_basis.Mult(projection.GetData(), combination.GetData());
          action -= combination;
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
      DenseMatrix active_basis(basis.GetData(), n, count);
      active_basis.Mult(coefficients.GetData(), combination.GetData());
      Precondition(combination,correction);
      out.x += correction;
      a.Mult(out.x,action);
      residual = rhs; residual -= action;
      residual_norm = Norm(residual);
      require(std::isfinite(residual_norm), "Nonfinite GMRES residual");
    }
    out.count = iterations;
    require(residual_norm <= tolerance*rhs_norm, "GMRES iteration limit reached");
    return residual_norm;
  }
};
LinearSolver::LinearSolver(int n, int block, const Options &o, const Vector &metric)
    : impl(std::make_unique<Impl>(n,block,o,metric)) {}
LinearSolver::~LinearSolver() = default;
void LinearSolver::Prepare(const std::vector<std::vector<int>> &graph,const Options &o) {
  if(!impl->ilu && o.solver=="gmres" && o.precond=="block-ilu0") {
    impl->ilu=std::make_unique<FixedBlockILU>(graph,impl->block,o.ilu_ordering);
  }
}
void LinearSolver::Solve(const ElementOperator &op, const Vector &rhs,
                         const Options &o, LinearResult &out, bool refresh) {
  auto &w = *impl;
  out.x.SetSize(w.n);
  out.equation_residual.SetSize(w.n);
  out.count = 0;
  out.residual = out.assembly = out.setup = out.solve = 0.;
  const double rhs_norm = w.Norm(rhs);
  if (rhs_norm == 0.) {
    out.x = 0.;
    out.equation_residual = 0.;
    return;
  }

  // Production prepares once before Run's timer; isolated reference calls lazily prepare.
  if(o.solver=="gmres" && o.precond=="block-ilu0" && !w.ilu) {
    Prepare(op.BlockGraph(w.block),o);
  }
  double start = clock_s();
  if (refresh && (!o.matrix_free || o.solver != "gmres")) {
    w.assembled = op.Assemble();
  }
  if(refresh && w.ilu) { op.FillBlocks(*w.ilu); }
  out.assembly = clock_s() - start;

  start = clock_s();
  if (refresh) {
    if (o.solver == "gmres") {
      if (w.ilu) {
        w.ilu->Factor();
      } else {
        op.DiagonalBlocks(w.block, w.diagonal);
        for (size_t e=0; e<w.diagonal.size(); ++e) {
          mfem::LUFactors lu(w.diagonal[e].GetData(),w.diagonal_pivots[e].data());
          require(lu.Factor(w.block), "Singular block Jacobi factor");
        }
      }
    } else if (o.solver == "dense") {
      w.assembled->ToDenseMatrix(w.dense_matrix);
      w.dense_lu->Factor();
    } else {
#ifdef MFEM_USE_MKL_PARDISO
      w.direct = std::make_unique<mfem::PardisoSolver>();
      w.direct->SetMatrixType(mfem::PardisoSolver::REAL_NONSYMMETRIC);
      w.direct->SetOperator(*w.assembled);
#endif
    }
  }
  out.setup = clock_s() - start;

  start = clock_s();
  double gmres_residual_norm = -1.;
  if (o.solver == "gmres") {
    const mfem::Operator &a = o.matrix_free
        ? static_cast<const mfem::Operator &>(op)
        : static_cast<const mfem::Operator &>(*w.assembled);
    gmres_residual_norm = w.GMRES(a, rhs, o.linear_tol, o.max_krylov, out, refresh);
  } else if (o.solver == "dense") {
    w.dense_lu->Mult(rhs, out.x);
  } else {
#ifdef MFEM_USE_MKL_PARDISO
    w.direct->Mult(rhs, out.x);
#endif
  }
  out.solve = clock_s() - start;

  // Matrix-free GMRES has just formed the true residual rhs-op*x. Reuse it for
  // validation and for the caller's recovered DG residual check. When GMRES
  // used an assembled CSR (or a direct solver), retain the independent action
  // of the requested elemental operator.
  if (o.solver == "gmres" && o.matrix_free) {
    out.equation_residual = w.residual;
    out.equation_residual *= -1.;
    out.residual = gmres_residual_norm / rhs_norm;
  } else {
    op.Mult(out.x, out.equation_residual);
    out.equation_residual -= rhs;
    out.residual = w.Norm(out.equation_residual) / rhs_norm;
  }
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
