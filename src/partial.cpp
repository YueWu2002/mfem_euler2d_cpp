#include "euler.hpp"
#include <algorithm>
#include <limits>
namespace euler {
namespace {
// Quadrature-level partial assembly: cache flux Jacobians instead of the dense
// cell and face matrices. Shape/derivative tables are shared across rectangles.
class PartialDGOperator final : public ElementOperator {
  const Problem &p;
  Vector state;
  std::vector<Matrix4> jx, jy, jl, jr;
  // Per-cell flux scratch, reused by explicit one-dimensional contractions.
  mutable std::vector<State> flux_x, flux_y;
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
    require(p.collocated, "Tensor DG action requires GL collocation");
    flux_x.resize(p.p); flux_y.resize(p.p);
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
        p.gas.jacobian(v, 0, jx[e*nq+q]); p.gas.jacobian(v, 1, jy[e*nq+q]);
      }
    }
    int id = 0;
    for (const auto &f : p.faces) {
      for (int q = 0; q <= p.opt.k; q++, ++id) {
        auto &l = jl[id]; auto &r = jr[id];
        p.gas.jacobian(p.Eval(u, f.l, p.fv[f.ls], q), f.axis, l);
        p.gas.jacobian(p.Eval(u, f.r, p.fv[f.rs], q), f.axis, r);
        const double left_weight = f.left_weight[q];
        const double right_weight = f.right_weight[q];
        const double penalty = f.penalty[q];
        for (int d = 0; d < 4; d++) {
#pragma omp simd
          for (int c = 0; c < 4; c++) {
            l(c, d) = left_weight * l(c, d) + (c == d ? penalty : 0);
            r(c, d) = right_weight * r(c, d) - (c == d ? penalty : 0);
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
      const int n1=p.opt.k+1;
      for (int point=0; point<p.p; ++point) {
        const int q=p.tensor_quad[point];
        const double *v=x.GetData()+e*p.m+4*p.tensor_node[point];
        auto &fx=flux_x[point]; auto &fy=flux_y[point];
        for (int c=0; c<4; ++c) { fx[c]=0.; fy[c]=0.; }
        const int id=e*nq+q;
        // Matrix4 is column-major; traverse one contiguous column at a time and
        // vectorize the four conservative output components.
        for (int d = 0; d < 4; d++) {
          const double value = v[d];
#pragma omp simd
          for (int c = 0; c < 4; c++) {
            fx[c] += jx[id](c, d) * value;
            fy[c] += jy[id](c, d) * value;
          }
        }
        const double weight=p.dt*p.w[q];
#pragma omp simd
        for (int c=0; c<4; ++c) { fx[c]*=weight; fy[c]*=weight; }
      }
      // Weak derivatives D_x^T F_x and D_y^T F_y along tensor lines.
      // q^2 point fluxes plus two q^3 contractions; no q^2 by q^2 action.
      for (int j=0; j<n1; ++j) {
        for (int i=0; i<n1; ++i) {
          State sum{};
          for (int a=0; a<n1; ++a) {
            const double dx=p.tensor_dx(a,i), dy=p.tensor_dy(a,j);
#pragma omp simd
            for (int c=0; c<4; ++c) {
              sum[c]+=dx*flux_x[j*n1+a][c]+dy*flux_y[a*n1+i][c];
            }
          }
          double *cell=y.GetData()+e*p.m+4*p.tensor_node[j*n1+i];
#pragma omp simd
          for (int c=0; c<4; ++c) { cell[c]=-sum[c]; }
        }
      }
    }
    int id = 0;
    for (const auto &f : p.faces) {
      for (int q = 0; q <= p.opt.k; q++, id++) {
        auto l = p.Eval(x, f.l, p.fv[f.ls], q), r = p.Eval(x, f.r, p.fv[f.rs], q);
        State z{};
        for (int d = 0; d < 4; d++) {
          const double left = l[d], right = r[d];
#pragma omp simd
          for (int c = 0; c < 4; c++) {
            z[c] += jl[id](c, d) * left + jr[id](c, d) * right;
          }
        }
        double scale = p.dt * f.length / p.area * p.ew[q];
        for (const auto &entry : p.face_basis[f.ls][q]) {
          double *cell = y.GetData() + f.l*p.m + 4*entry.first;
          const double weight = scale*entry.second;
#pragma omp simd
          for (int c = 0; c < 4; ++c) { cell[c] += weight*z[c]; }
        }
        for (const auto &entry : p.face_basis[f.rs][q]) {
          double *cell = y.GetData() + f.r*p.m + 4*entry.first;
          const double weight = scale*entry.second;
#pragma omp simd
          for (int c = 0; c < 4; ++c) { cell[c] -= weight*z[c]; }
        }
      }
    }
    p.MassSolve(y);
    y += x;
  }
  std::vector<std::vector<int>> BlockGraph(int size) const override {
    require(size==p.m,"DG BSR requires cell blocks");
    return p.SolverGraph(false);
  }
  void FillBlocks(FixedBlockILU &target) const override {
    require(target.block==p.m,"DG BSR requires cell blocks");
    // VolumeJacobianCell overwrites each diagonal; do not clear it twice.
    target.Zero(true);
    for(int e=0;e<p.n;++e) {
      p.VolumeJacobianCell(state,e,target.At(e,e),&jx,&jy);
    }
    for(int fidx=0;fidx<p.nf;++fidx) {
      const auto &f=p.faces[fidx];
      auto &ll=target.At(f.l,f.l); auto &lr=target.At(f.l,f.r);
      auto &rl=target.At(f.r,f.l); auto &rr=target.At(f.r,f.r);
      for(int q=0;q<=p.opt.k;++q) {
        const int id=fidx*(p.opt.k+1)+q;
        const double scale=p.dt*f.length/p.area*p.ew[q];
        for(int side=0;side<2;++side) {
          const int test_side=side?f.rs:f.ls;
          for(int neighbor=0;neighbor<2;++neighbor) {
            const int trial_side=neighbor?f.rs:f.ls;
            auto &a=side?(neighbor?rr:rl):(neighbor?lr:ll);
            const auto &flux=neighbor?jr[id]:jl[id];
            for(const auto &i:p.face_basis[test_side][q]) {
              for(const auto &j:p.face_basis[trial_side][q]) {
                const double weight=(side?-scale:scale)*p.inv_mass(i.first,i.first)*i.second*j.second;
                for(int d=0;d<4;++d) {
#pragma omp simd
                  for(int c=0;c<4;++c) { a(4*i.first+c,4*j.first+d)+=weight*flux(c,d); }
                }
              }
            }
          }
        }
      }
    }
  }
  std::unique_ptr<mfem::SparseMatrix> Assemble() const override {
    // Two streaming passes: count exact row lengths, then fill final CSR.
    // Only one cell's tensor contributions are live; no global dense face blocks.
    auto offsets = std::make_unique<int[]>(height+1);
    offsets[0]=0;
    std::unique_ptr<int[]> columns;
    std::unique_ptr<double[]> values;
    using Item=std::pair<int,double>;
    std::vector<std::vector<Item>> rows(p.m);
    for (int pass=0; pass<2; ++pass) {
      for (int e=0; e<p.n; ++e) {
        for (auto &row:rows) { row.clear(); }
        for (int i=0; i<p.m; ++i) { rows[i].emplace_back(e*p.m+i,1.); }
        for (int q=0; q<p.V.Height(); ++q) {
          const int j=p.volume_node[q];
          for (const auto &test:p.volume_tests[q]) {
            const double scale=-p.dt*p.w[q]*p.inv_mass(test.node,test.node);
            for (int d=0; d<4; ++d) {
              for (int c=0; c<4; ++c) {
                const double value=scale*(test.dx*jx[e*p.V.Height()+q](c,d)
                                        +test.dy*jy[e*p.V.Height()+q](c,d));
                if (value!=0.) { rows[4*test.node+c].emplace_back(e*p.m+4*j+d,value); }
              }
            }
          }
        }
        for (int side=0; side<4; ++side) {
          const int face=p.cell_faces[e][side];
          const auto &f=p.faces[face];
          const bool left=f.l==e && f.ls==side;
          const double scale=(left?1.:-1.)*p.dt*f.length/p.area;
          for (int q=0; q<=p.opt.k; ++q) {
            const int id=face*(p.opt.k+1)+q;
            for (const auto &i:p.face_basis[side][q]) {
              for (int neighbor=0; neighbor<2; ++neighbor) {
                const int cell=neighbor?f.r:f.l;
                const int trial_side=neighbor?f.rs:f.ls;
                const auto &flux=neighbor?jr[id]:jl[id];
                for (const auto &j:p.face_basis[trial_side][q]) {
                  const double weight=scale*p.ew[q]*p.inv_mass(i.first,i.first)
                                      *i.second*j.second;
                  for (int d=0; d<4; ++d) {
                    for (int c=0; c<4; ++c) {
                      const double value=weight*flux(c,d);
                      if (value!=0.) { rows[4*i.first+c].emplace_back(cell*p.m+4*j.first+d,value); }
                    }
                  }
                }
              }
            }
          }
        }
        for (int i=0; i<p.m; ++i) {
          auto &row=rows[i];
          std::stable_sort(row.begin(),row.end(),
              [](const Item &a,const Item &b){return a.first<b.first;});
          size_t count=0;
          for (size_t read=0; read<row.size();) {
            const int column=row[read].first;
            double value=0.;
            do { value+=row[read++].second; }
            while(read<row.size() && row[read].first==column);
            row[count++]={column,value};
          }
          const int global=e*p.m+i;
          if (pass==0) {
            require(count<=size_t(std::numeric_limits<int>::max()-offsets[global]),
                    "CSR exceeds MFEM int index capacity");
            offsets[global+1]=offsets[global]+int(count);
          } else {
            for (size_t j=0; j<count; ++j) {
              columns[offsets[global]+j]=row[j].first;
              values[offsets[global]+j]=row[j].second;
            }
          }
        }
      }
      if(pass==0) {
        columns=std::make_unique<int[]>(offsets[height]);
        values=std::make_unique<double[]>(offsets[height]);
      }
    }
    return std::make_unique<mfem::SparseMatrix>(offsets.release(),columns.release(),
                                               values.release(),height,height);
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
            const double weight = scale*p.inv_mass(i.first,i.first)*i.second*j.second;
            for (int d = 0; d < 4; ++d) {
#pragma omp simd
              for (int c = 0; c < 4; ++c) {
                diag[f.l](4*i.first+c,4*j.first+d) +=
                    weight*jl[id](c,d);
              }
            }
          }
        }
        for (const auto &i : p.face_basis[f.rs][q]) {
          for (const auto &j : p.face_basis[f.rs][q]) {
            const double weight = scale*p.inv_mass(i.first,i.first)*i.second*j.second;
            for (int d = 0; d < 4; ++d) {
#pragma omp simd
              for (int c = 0; c < 4; ++c) {
                diag[f.r](4*i.first+c,4*j.first+d) -=
                    weight*jr[id](c,d);
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
