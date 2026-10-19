#include "euler.hpp"
#include <algorithm>
#include <limits>

// MFEM's existing MDF implementation is defined in linalg/solvers.cpp.
// Apply it once to the unweighted topology graph, not to Newton matrix values.
namespace mfem { void MinimumDiscardedFillOrdering(SparseMatrix &, Array<int> &); }

namespace euler {
FixedBlockILU::FixedBlockILU(const std::vector<std::vector<int>> &graph,
                           int size, const std::string &ordering) : block(size) {
  const int n=graph.size();
  require(block>0, "Invalid BSR block size");
  permutation.resize(n); inverse.resize(n); diag.resize(n);
  for (int i=0; i<n; ++i) { permutation[i]=i; }
  if (ordering=="mdf") {
    mfem::SparseMatrix topology(n);
    for (int i=0; i<n; ++i) {
      topology.Add(i,i,1.);
      for (int j:graph[i]) { if(j!=i) { topology.Add(i,j,1.); } }
    }
    topology.Finalize();
    mfem::Array<int> order;
    mfem::MinimumDiscardedFillOrdering(topology,order);
    for (int i=0; i<n; ++i) { permutation[i]=order[i]; }
  }
  for (int i=0; i<n; ++i) { inverse[permutation[i]]=i; }
  row.push_back(0);
  for (int i=0; i<n; ++i) {
    std::vector<int> neighbors{i};
    for (int j:graph[permutation[i]]) { neighbors.push_back(inverse[j]); }
    std::sort(neighbors.begin(),neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(),neighbors.end()),neighbors.end());
    require(neighbors.size()<=size_t(std::numeric_limits<int>::max())-col.size(),
            "BSR block index overflow");
    for (int j:neighbors) {
      if(j==i) { diag[i]=col.size(); }
      col.push_back(j);
    }
    row.push_back(col.size());
  }
  data.resize(col.size());
  for (auto &a:data) { a.SetSize(block); }
  pivots.resize(n,std::vector<int>(block));
  work.SetSize(n*block);
  updates.resize(col.size());
  // Cache all retained ILU(0) products. No symbolic searches during Factor().
  for (int i=0; i<n; ++i) {
    for (int ik=row[i]; ik<diag[i]; ++ik) {
      const int k=col[ik];
      int ij=ik+1, kj=diag[k]+1;
      while(ij<row[i+1] && kj<row[k+1]) {
        if(col[ij]<col[kj]) { ++ij; }
        else if(col[kj]<col[ij]) { ++kj; }
        else { updates[ik].push_back({ij,kj}); ++ij; ++kj; }
      }
    }
  }
}
void FixedBlockILU::Zero(bool off_diagonal_only) {
  for(int i=0;i<int(diag.size());++i) {
    for(int j=row[i];j<row[i+1];++j) {
      if(!off_diagonal_only || j!=diag[i]) { data[j]=0.; }
    }
  }
}
DenseMatrix &FixedBlockILU::At(int original_row,int original_col) {
  const int i=inverse[original_row], j=inverse[original_col];
  auto it=std::lower_bound(col.begin()+row[i],col.begin()+row[i+1],j);
  require(it!=col.begin()+row[i+1] && *it==j,"Entry outside fixed BSR graph");
  return data[it-col.begin()];
}
void FixedBlockILU::Factor() {
  for (int i=0; i<int(diag.size()); ++i) {
    for (int ik=row[i]; ik<diag[i]; ++ik) {
      const int k=col[ik];
      mfem::LUFactors lu(data[diag[k]].GetData(),pivots[k].data());
      lu.RightSolve(block,block,data[ik].GetData());
      for (const auto &u:updates[ik]) {
        mfem::AddMult_a(-1.,data[ik],data[u.upper],data[u.target]);
      }
    }
    mfem::LUFactors lu(data[diag[i]].GetData(),pivots[i].data());
    require(lu.Factor(block),"Singular BSR ILU diagonal");
  }
}
void FixedBlockILU::Mult(const Vector &b,Vector &x) const {
  const int n=diag.size();
  x.SetSize(n*block);
  for (int i=0; i<n; ++i) {
    Vector wi(work.GetData()+i*block,block);
    std::copy_n(b.GetData()+permutation[i]*block,block,wi.GetData());
    for (int ij=row[i]; ij<diag[i]; ++ij) {
      Vector wj(work.GetData()+col[ij]*block,block);
      data[ij].AddMult_a(-1.,wj,wi);
    }
  }
  for (int i=n-1; i>=0; --i) {
    Vector wi(work.GetData()+i*block,block);
    for (int ij=diag[i]+1; ij<row[i+1]; ++ij) {
      Vector wj(work.GetData()+col[ij]*block,block);
      data[ij].AddMult_a(-1.,wj,wi);
    }
    mfem::LUFactors lu(const_cast<double*>(data[diag[i]].GetData()),
                       const_cast<int*>(pivots[i].data()));
    lu.Solve(block,1,wi.GetData());
    std::copy_n(wi.GetData(),block,x.GetData()+permutation[i]*block);
  }
}
std::vector<std::vector<int>> Problem::SolverGraph(bool hdg) const {
  std::vector<std::vector<int>> graph(hdg?nf:n);
  if(hdg) {
    for(const auto &local:cell_faces) {
      for(int i:local) { for(int j:local) { graph[i].push_back(j); } }
    }
  } else {
    for(int e=0;e<n;++e) { graph[e].push_back(e); }
    for(const auto &f:faces) { graph[f.l].push_back(f.r); graph[f.r].push_back(f.l); }
  }
  for(auto &r:graph) {
    std::sort(r.begin(),r.end()); r.erase(std::unique(r.begin(),r.end()),r.end());
  }
  return graph;
}
} // namespace euler
