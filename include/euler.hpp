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
class EulerGas {
protected:
  const double gamma, gamma_minus_one;
public:
  explicit EulerGas(double value) : gamma(value), gamma_minus_one(value-1.) {
    require(std::isfinite(value) && value > 1., "Euler gamma must be finite and > 1");
  }
  double Gamma() const { return gamma; }
  double pressure(const State &u) const;
  bool positive(const State &u) const;
  State flux(const State &u, int axis) const;
  void jacobian(const State &u, int axis, Matrix4 &j) const;
  std::array<double, 2> speeds(const State &l, const State &r, int axis, double safety) const;
  State exact(double x, double y, double time, double alpha) const;
};
// Body layout: (cell, GL node, conservative component).
// Trace layout: (periodic face, GL node, conservative component).
struct Options {
  int nx = 100, ny = 100, k = 1, volume_order = -1, max_newton = 30;
  int restart = 20, max_krylov = 1500, newton_steps = -1, steps = -1;
  double alpha = 45, dt = .001, end = .2, cfl = -1, safety = 1.1, tol = 1e-12,
         linear_tol = 1e-12;
  std::string assembly = "partial";
  std::string form = "both", solver = "gmres", precond = "block-jacobi",
              integrator = "ssp-dirk3", output = "code_test/results/run";
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
  // Frozen-bound HLL weights shared by residual and Jacobian kernels.
  std::vector<double> inverse_sum, left_weight, right_weight, penalty;
};
// Qk face restriction: one normal-direction contraction per tangential GL node.
// Shared interpolation tables; C stores 4x4 blocks, B scalar identity factors.
class FaceCoupling {
public:
  const std::array<std::vector<std::vector<std::pair<int,double>>>,4> *basis = nullptr;
  const DenseMatrix *inverse_mass = nullptr;
  int nodes = 0, cell_size = 0;
  bool lift = false;
  std::vector<Matrix4> flux;
  // Lift B is scalar times I_4; restriction C retains full flux blocks.
  std::vector<double> scalar_flux;
  int Height() const { return lift ? cell_size : 16*nodes; }
  int Width() const { return lift ? 16*nodes : cell_size; }
  void Mult(const Vector &x, Vector &y) const;
  void Mult(const DenseMatrix &x, DenseMatrix &y) const;
  void Dense(DenseMatrix &a) const;
  // Only one output face and one input X panel; no other Schur blocks formed.
  void Panel(const DenseMatrix &x, int side, int first, int count,
             DenseMatrix &out) const;
};
// Fixed topology, block ordering and elimination schedule. Numeric blocks are
// filled directly and overwritten by ILU(0), including diagonal pivoted LU.
class FixedBlockILU {
  struct Update { int target, upper; };
  std::vector<std::vector<Update>> updates;
  std::vector<int> row, col, diag, permutation, inverse;
  std::vector<DenseMatrix> data;
  std::vector<std::vector<int>> pivots;
  mutable Vector work;
public:
  int block;
  FixedBlockILU(const std::vector<std::vector<int>> &graph, int size,
                const std::string &ordering);
  void Zero(bool off_diagonal_only = false);
  DenseMatrix &At(int original_row, int original_col);
  void Factor();
  void Mult(const Vector &b, Vector &x) const;
};
// Borrowed factors: a condensed operator must not outlive its Hybrid and X.
struct FactorEntry {
  std::vector<int> ids;
  const FaceCoupling *left = nullptr;
  const DenseMatrix *right = nullptr;
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
  Vector diagonal;
  explicit ElementOperator(int n) : Operator(n) {}
  void Add(const std::vector<int> &r, const std::vector<int> &c, const DenseMatrix &b);
  void Mult(const Vector &x, Vector &y) const override;
  virtual std::unique_ptr<mfem::SparseMatrix> Assemble() const;
  virtual DenseMatrix DiagonalBlock(int start, int size) const;
  virtual std::vector<DenseMatrix> DiagonalBlocks(int size) const;
  virtual void DiagonalBlocks(int size, std::vector<DenseMatrix> &blocks) const;
  virtual void Refresh(const Vector &) { throw std::runtime_error("Operator is not refreshable"); }
  virtual std::vector<std::vector<int>> BlockGraph(int block) const;
  virtual void FillBlocks(FixedBlockILU &target) const;
  void PrepareScratch(int gather, int output, int intermediate);
};
// Local Newton blocks only: no globally assembled augmented HDG operator.
struct Hybrid {
  std::vector<DenseMatrix> A;
  Vector D;
  std::vector<FaceCoupling> B, C;
  std::vector<std::vector<int>> ids;
  int freeze_generation = -1;
  double cached_dt = -1.;
};
// X_e = A_e^{-1} B_e; y_e = A_e^{-1} r_e.
// X stores solves with local LU, never an explicitly formed inverse of A_e.
struct Condensed {
  std::vector<DenseMatrix> X;
  std::vector<Vector> y;
  Vector rhs, cell_scratch, trace_scratch, local_scratch, scaled_scratch;
  DenseMatrix product_scratch;
  std::vector<std::vector<int>> pivots;
  DenseMatrix cell_matrix;
  std::vector<int> cell_pivots;
  // Destroy the borrowed-factor operator before its X storage.
  std::unique_ptr<ElementOperator> op;
  // Destroy recovery before the borrowed X and Hybrid storage.
  std::unique_ptr<mfem::Operator> recovery;
};
class Problem {
public:
  Options opt;
  const EulerGas gas;
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
  // Lexicographic (y*q+x) to native volume node/quadrature ordering.
  std::vector<int> tensor_node, tensor_quad;
  DenseMatrix tensor_dx, tensor_dy;
  std::array<std::vector<std::vector<std::pair<int, double>>>, 4> face_basis;
  std::vector<std::vector<VolumeTest>> volume_tests;
  bool collocated = false;
  Vector cell_metric, trace_metric;
  std::vector<Face> faces;
  std::vector<std::array<int, 4>> cell_faces;
  int freeze_generation = 0;
  Vector old;
  mutable Matrix4 jac_x{}, jac_y{}, face_left{}, face_right{};
  explicit Problem(const Options &o, const EulerGas &gas);
  std::vector<int> Uids(int e) const;
  std::vector<int> Lids(int f) const;
  State Eval(const Vector &u, int e, const DenseMatrix &basis, int q) const;
  bool Admissible(const Vector &u) const;
  Vector Initial() const;
  void Freeze(const Vector &u, double step, bool check_admissible = true);
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
  void VolumeJacobianCell(const Vector &u, int e, DenseMatrix &a,
                      const std::vector<Matrix4> *cache_x = nullptr,
                      const std::vector<Matrix4> *cache_y = nullptr) const;
  std::unique_ptr<ElementOperator> DG(const Vector &u) const;
  std::unique_ptr<ElementOperator> PartialDG(const Vector &u, bool refresh = true) const;
  Hybrid HDG(const Vector &u) const;
  void PrepareHybrid(Hybrid &h, Condensed &s) const;
  void HDG(const Vector &u, Hybrid &h) const;
  void HybridCell(const Vector &u, int e, Hybrid &h, DenseMatrix &a) const;
  void HybridDiagonal(Hybrid &h) const;
  void StreamCondense(const Vector &u, const Vector &r, Hybrid &h, Condensed &s) const;
  std::vector<std::vector<int>> SolverGraph(bool hdg) const;
  Condensed Condense(const Hybrid &h, const Vector &r, const Vector *g = nullptr) const;
  void Condense(Hybrid &h, const Vector &r, Condensed &s,
                const Vector *g = nullptr) const;
  Vector Recover(const Condensed &s, const Vector &dl) const;
  void Recover(const Condensed &s, const Vector &dl, Vector &du) const;
  std::array<double, 4> Errors(const Vector &u, double time) const;
  std::array<double, 4> Integral(const Vector &u) const;
  std::array<double, 4> Differences(const Vector &u, const Vector &v) const;
  void Save(const Vector &u, const std::string &prefix) const;

  void MassSolve(DenseMatrix &a) const;
  void MassSolve(Vector &r) const;
};
struct LinearResult {
  Vector x, equation_residual;
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
  void Prepare(const std::vector<std::vector<int>> &graph, const Options &o);
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
void Verify(const Options &o, const EulerGas &gas);
Options ParseOptions(int argc, char **argv);
} // namespace euler
