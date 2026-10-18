#pragma once
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <mfem-performance.hpp>
#include <stdexcept>
#include <string>
#include <vector>
namespace euler {
inline constexpr double pi = 3.141592653589793238462643383279502884;
// Native fixed-size tensors contain their data inline, with no dynamic storage.
using State = mfem::TVector<4>;
using Matrix4 = mfem::TMatrix<4, 4>;
using mfem::DenseMatrix;
using mfem::Vector;
inline void require(bool ok, const char *msg) {
  if (!ok) {
    throw std::runtime_error(msg);
  }
}
inline double clock_s() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}
double pressure(const State &u);
bool positive(const State &u);
State flux(const State &u, int axis);
void jacobian(const State &u, int axis, Matrix4 &j);
std::array<double, 2> speeds(const State &l, const State &r, int axis, double safety);
State exact(double x, double y, double time, double alpha);
// Body layout: (cell, GL node, conservative component).
// Trace layout: (periodic face, GL node, conservative component).
struct Options {
  int nx = 100, ny = 100, k = 1, volume_order = -1, max_newton = 25;
  int restart = 20, max_krylov = 1500, newton_steps = -1, steps = -1;
  double alpha = 45, dt = .001, end = .2, cfl = -1, safety = 1.1, tol = 1e-10,
         linear_tol = 1e-11;
  std::string assembly = "partial";
  std::string form = "both", solver = "gmres", precond = "block-jacobi",
              output = "results/run";
  bool matrix_free = true, verify = false, save_solution = true, hdg_first = false;
  bool check_update = false, timing_only = false;
  std::string ilu_ordering = "none";
};
struct VolumeTest {
  int node;
  double dx, dy;
};
struct Face {
  int l = -1, r = -1, ls = -1, rs = -1, axis = 0;
  double length = 0.;
  std::vector<double> a, b;
};
// Borrowed factors: a condensed operator must not outlive its Hybrid and X.
struct FactorEntry {
  std::vector<int> ids;
  const DenseMatrix *left = nullptr, *right = nullptr;
  double sign = -1;
};
struct Entry {
  std::vector<int> rows, cols;
  DenseMatrix block;
};
// Elemental operator: no global sparse matrix is needed for Mult().
class ElementOperator : public mfem::Operator {
  mutable Vector local_scratch, output_scratch, factor_scratch;
public:
  std::vector<Entry> entries;
  std::vector<FactorEntry> factors;
  explicit ElementOperator(int n) : Operator(n) {}
  void Add(const std::vector<int> &r, const std::vector<int> &c, const DenseMatrix &b);
  void Mult(const Vector &x, Vector &y) const override;
  virtual std::unique_ptr<mfem::SparseMatrix> Assemble() const;
  virtual DenseMatrix DiagonalBlock(int start, int size) const;
  virtual std::vector<DenseMatrix> DiagonalBlocks(int size) const;
  virtual void DiagonalBlocks(int size, std::vector<DenseMatrix> &blocks) const;
  virtual void Refresh(const Vector &) { throw std::runtime_error("Operator is not refreshable"); }
  void PrepareScratch(int gather, int output, int intermediate);
};
// Local Newton blocks only: no globally assembled augmented HDG operator.
struct Hybrid {
  std::vector<DenseMatrix> A, B, C, D;
  std::vector<std::vector<int>> ids;
  int freeze_generation = -1;
};
// X_e = A_e^{-1} B_e; y_e = A_e^{-1} r_e.
// X stores solves with local LU, never an explicitly formed inverse of A_e.
struct Condensed {
  std::vector<DenseMatrix> X;
  std::vector<Vector> y;
  Vector rhs, cell_scratch, trace_scratch, local_scratch, eps_scratch, scaled_scratch;
  DenseMatrix product_scratch;
  std::vector<std::unique_ptr<mfem::DenseMatrixInverse>> lu;
  // Destroy the borrowed-factor operator before its X storage.
  std::unique_ptr<ElementOperator> op;
  // Destroy recovery before the borrowed X and Hybrid storage.
  std::unique_ptr<mfem::Operator> recovery;
};
class Problem {
public:
  Options opt;
  double alpha, lx, ly, hx, hy, area, dt = 0;
  // n: cells; p: scalar body nodes; m=4*p: cell block size.
  // t=4*(k+1): face block size; nf: periodic faces; nu/nl: total DOFs.
  int n, p, m, t, nf, nu, nl;
  mfem::Mesh mesh;
  mfem::L2_FECollection fec, trace_fec;
  mfem::FiniteElementSpace fes;
  const mfem::FiniteElement *fe;
  DenseMatrix mass, inv_mass, V, Gx, Gy, psi;
  std::array<DenseMatrix, 4> fv;
  std::vector<double> w, ew, ex, origx, origy;
  std::vector<int> volume_node;
  std::array<std::vector<std::vector<std::pair<int, double>>>, 4> face_basis;
  std::vector<std::vector<VolumeTest>> volume_tests;
  bool collocated = false;
  Vector cell_metric, trace_metric;
  std::vector<Face> faces;
  std::vector<std::array<int, 4>> cell_faces;
  int freeze_generation = 0;
  Vector old;
  mutable Matrix4 jac_x{}, jac_y{}, face_left{}, face_right{};
  explicit Problem(const Options &o);
  std::vector<int> Uids(int e) const;
  std::vector<int> Lids(int f) const;
  State Eval(const Vector &u, int e, const DenseMatrix &basis, int q) const;
  bool Admissible(const Vector &u) const;
  Vector Initial() const;
  void Freeze(const Vector &u, double step);
  double WaveRate() const;
  double CellNorm(const Vector &r) const;
  Vector Trace(const Vector &u) const;
  Vector Skeleton(const Vector &u, const Vector &lam) const;
  Vector Residual(const Vector &u, const Vector *lam = nullptr) const;
  void Residual(const Vector &u, Vector &r, const Vector *lam = nullptr) const;
  std::vector<DenseMatrix> VolumeJacobian(const Vector &u) const;
  void VolumeJacobian(const Vector &u, std::vector<DenseMatrix> &a,
                      const std::vector<Matrix4> *cache_x = nullptr,
                      const std::vector<Matrix4> *cache_y = nullptr) const;
  std::unique_ptr<ElementOperator> DG(const Vector &u) const;
  std::unique_ptr<ElementOperator> PartialDG(const Vector &u, bool refresh = true) const;
  Hybrid HDG(const Vector &u) const;
  void PrepareHybrid(Hybrid &h, Condensed &s) const;
  void HDG(const Vector &u, Hybrid &h) const;
  Condensed Condense(const Hybrid &h, const Vector &r, const Vector *g = nullptr) const;
  void Condense(const Hybrid &h, const Vector &r, Condensed &s,
                const Vector *g = nullptr) const;
  Vector Recover(const Condensed &s, const Vector &dl) const;
  void Recover(const Condensed &s, const Vector &dl, Vector &du) const;
  std::array<double, 4> Errors(const Vector &u, double time) const;
  std::array<double, 4> Integral(const Vector &u) const;
  double Difference(const Vector &u, const Vector &v) const;
  void Save(const Vector &u, const std::string &prefix) const;

  void MassSolve(DenseMatrix &a) const;
  void MassSolve(Vector &r) const;
};
struct LinearResult {
  Vector x;
  int count = 0;
  double residual = 0, assembly = 0, setup = 0, solve = 0;
};
LinearResult Solve(const ElementOperator &op, const Vector &rhs, int block,
                   const Options &o, const Vector &metric);
// Persistent solver storage. Construct before the time loop. Numerical refresh,
// factorization, assembly and all applications remain inside Newton timings.
class LinearSolver {
  struct Impl;
  std::unique_ptr<Impl> impl;
public:
  LinearSolver(int n, int block, const Options &o, const Vector &metric);
  ~LinearSolver();
  void Solve(const ElementOperator &op, const Vector &rhs,
             const Options &o, LinearResult &result, bool refresh = true);
};
struct RunResult {
  Vector u;
  std::array<double, 4> errors{};
  double seconds = 0, drift = 0, newton_seconds = 0, final_time = 0;
  int updates = 0, time_steps = 0;
};
RunResult Run(Problem &p, const std::string &form, const std::string &prefix);
void Verify(const Options &o);
Options ParseOptions(int argc, char **argv);
} // namespace euler
