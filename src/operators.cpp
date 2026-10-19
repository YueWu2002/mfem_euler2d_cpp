#include "euler.hpp"

#include <algorithm>
#include <limits>

namespace euler {
void FaceCoupling::Panel(const DenseMatrix &x,int side,int first,int count,
                         DenseMatrix &out) const {
  require(!lift,"Only restriction couplings have face output panels");
  out.SetSize(4*nodes,count);
  for(int j=0;j<count;++j) {
    for(int q=0;q<nodes;++q) {
      double value[4]={}, result[4]={};
      for(const auto &entry:(*basis)[side][q]) {
#pragma omp simd
        for(int d=0;d<4;++d) { value[d]+=entry.second*x(4*entry.first+d,first+j); }
      }
      const auto &a=flux[side*nodes+q];
      for(int d=0;d<4;++d) {
#pragma omp simd
        for(int c=0;c<4;++c) { result[c]+=a(c,d)*value[d]; }
      }
      for(int c=0;c<4;++c) { out(4*q+c,j)=result[c]; }
    }
  }
}
std::vector<std::vector<int>> ElementOperator::BlockGraph(int block) const {
  std::vector<std::vector<int>> graph(height/block);
  for(int i=0;i<height/block;++i) { graph[i].push_back(i); }
  for(const auto &e:entries) {
    std::vector<int> rs,cs;
    for(int r:e.rows) { rs.push_back(r/block); }
    for(int c:e.cols) { cs.push_back(c/block); }
    std::sort(rs.begin(),rs.end()); rs.erase(std::unique(rs.begin(),rs.end()),rs.end());
    std::sort(cs.begin(),cs.end()); cs.erase(std::unique(cs.begin(),cs.end()),cs.end());
    for(int r:rs) { graph[r].insert(graph[r].end(),cs.begin(),cs.end()); }
  }
  for(const auto &f:factors) {
    for(int i=0;i<int(f.ids.size());i+=block) {
      for(int j=0;j<int(f.ids.size());j+=block) { graph[f.ids[i]/block].push_back(f.ids[j]/block); }
    }
  }
  for(auto &r:graph) {
    std::sort(r.begin(),r.end()); r.erase(std::unique(r.begin(),r.end()),r.end());
  }
  return graph;
}
void ElementOperator::FillBlocks(FixedBlockILU &target) const {
  const int b=target.block;
  target.Zero();
  if(diagonal.Size()) {
    for(int r=0;r<height;r+=b) {
      auto &a=target.At(r/b,r/b);
      for(int i=0;i<b;++i) { a(i,i)=diagonal[r+i]; }
    }
  }
  for(const auto &e:entries) {
    for(int j=0;j<int(e.cols.size());++j) {
      for(int i=0;i<int(e.rows.size());++i) {
        target.At(e.rows[i]/b,e.cols[j]/b)(e.rows[i]%b,e.cols[j]%b)+=e.block(i,j);
      }
    }
  }
  DenseMatrix panel;
  for(const auto &f:factors) {
    require(b==4*f.left->nodes,"Schur BSR expects face blocks");
    for(int side=0;side<4;++side) {
      for(int other=0;other<4;++other) {
        f.left->Panel(*f.right,side,other*b,b,panel);
        target.At(f.ids[side*b]/b,f.ids[other*b]/b).Add(f.sign,panel);
      }
    }
  }
}
void FaceCoupling::Mult(const Vector &x, Vector &y) const {
  y.SetSize(Height()); y = 0.;
  for (int side=0; side<4; ++side) {
    for (int q=0; q<nodes; ++q) {
      const auto &restriction = (*basis)[side][q];
      const auto &a = flux[side*nodes+q];
      double value[4] = {}, result[4] = {};
      const int offset = 4*(side*nodes+q);
      if (lift) {
        for (int c=0; c<4; ++c) { value[c]=x[offset+c]; }
      } else {
        for (const auto &entry : restriction) {
#pragma omp simd
          for (int c=0; c<4; ++c) { value[c]+=entry.second*x[4*entry.first+c]; }
        }
      }
      for (int d=0; d<4; ++d) {
#pragma omp simd
        for (int c=0; c<4; ++c) { result[c]+=a(c,d)*value[d]; }
      }
      if (lift) {
        for (const auto &entry : restriction) {
          const double weight=entry.second*(*inverse_mass)(entry.first,entry.first);
#pragma omp simd
          for (int c=0; c<4; ++c) { y[4*entry.first+c]+=weight*result[c]; }
        }
      } else {
        for (int c=0; c<4; ++c) { y[offset+c]=result[c]; }
      }
    }
  }
}
void FaceCoupling::Mult(const DenseMatrix &x, DenseMatrix &y) const {
  y.SetSize(Height(),x.Width());
  for (int j=0; j<x.Width(); ++j) {
    Vector input(const_cast<double*>(x.GetData())+j*x.Height(),x.Height());
    Vector output(y.GetData()+j*y.Height(),y.Height());
    Mult(input,output);
  }
}
void FaceCoupling::Dense(DenseMatrix &a) const {
  a.SetSize(Height(),Width()); a=0.;
  for (int side=0; side<4; ++side) {
    for (int q=0; q<nodes; ++q) {
      const int offset=4*(side*nodes+q);
      for (const auto &entry : (*basis)[side][q]) {
        double weight=entry.second;
        if (lift) { weight*=(*inverse_mass)(entry.first,entry.first); }
        for (int d=0; d<4; ++d) {
          for (int c=0; c<4; ++c) {
            if (lift) { a(4*entry.first+c,offset+d)=weight*flux[side*nodes+q](c,d); }
            else { a(offset+c,4*entry.first+d)=weight*flux[side*nodes+q](c,d); }
          }
        }
      }
    }
  }
}
void ElementOperator::Add(const std::vector<int> &r, const std::vector<int> &c,
                          const DenseMatrix &b) {
  entries.push_back({r, c, b});
}
void ElementOperator::Mult(const Vector &x, Vector &y) const {
  y.SetSize(height);
  y = 0.;
  if (diagonal.Size()) {
#pragma omp simd
    for(int i=0;i<height;++i) { y[i]=diagonal[i]*x[i]; }
  }
  Vector &local = local_scratch, &z = output_scratch, &mid = factor_scratch;
  for (const auto &e : entries) {
    local.SetSize(e.cols.size());
    z.SetSize(e.rows.size());
    for (int j = 0; j < local.Size(); j++) {
      local[j] = x[e.cols[j]];
    }
    e.block.Mult(local, z);
    for (int i = 0; i < z.Size(); i++) {
      y[e.rows[i]] += z[i];
    }
  }
  for (const auto &f : factors) {
    local.SetSize(f.ids.size());
    mid.SetSize(f.right->Height());
    z.SetSize(f.ids.size());
    for (int i = 0; i < local.Size(); i++) {
      local[i] = x[f.ids[i]];
    }
    f.right->Mult(local, mid);
    f.left->Mult(mid, z);
    for (int i = 0; i < local.Size(); i++) {
      y[f.ids[i]] += f.sign * z[i];
    }
  }
}
std::unique_ptr<mfem::SparseMatrix> ElementOperator::Assemble() const {
  // Batch contribution reduction avoids MFEM LIL insertion's repeated row
  // scans. The finalized storage and matrix action remain native MFEM CSR.
  using Item = std::pair<int, double>;
  std::vector<std::vector<Item>> rows(height);
  for (auto &row : rows) {
    row.reserve(64);
  }
  if(diagonal.Size()) { for(int i=0;i<height;++i) { rows[i].emplace_back(i,diagonal[i]); } }
  auto append = [&](const std::vector<int> &ri, const std::vector<int> &ci,
                    const DenseMatrix &b, double sign) {
    for (int i = 0; i < b.Height(); i++) {
      for (int j = 0; j < b.Width(); j++) {
        if (b(i, j) != 0) {
          rows[ri[i]].emplace_back(ci[j], sign * b(i, j));
        }
      }
    }
  };
  for (const auto &e : entries) {
    append(e.rows, e.cols, e.block, 1);
  }
  for (const auto &f : factors) {
    DenseMatrix b(f.ids.size());
    f.left->Mult(*f.right, b);
    append(f.ids, f.ids, b, f.sign);
  }
  auto indices = std::make_unique<int[]>(height + 1);
  indices[0] = 0;
  for (int i = 0; i < height; i++) {
    auto &row = rows[i];
    std::stable_sort(row.begin(), row.end(),
                     [](const Item &a, const Item &b) { return a.first < b.first; });
    size_t write = 0, read = 0;
    while (read < row.size()) {
      int col = row[read].first;
      double value = 0;
      while (read < row.size() && row[read].first == col) {
        value += row[read++].second;
      }
      row[write++] = {col, value};
    }
    row.resize(write);
    require(write <= size_t(std::numeric_limits<int>::max() - indices[i]),
            "CSR exceeds MFEM int index capacity");
    indices[i + 1] = indices[i] + int(write);
  }
  int nnz = indices[height];
  auto columns = std::make_unique<int[]>(nnz);
  auto data = std::make_unique<double[]>(nnz);
  for (int i = 0; i < height; i++) {
    for (size_t j = 0; j < rows[i].size(); j++) {
      int offset = indices[i] + int(j);
      columns[offset] = rows[i][j].first;
      data[offset] = rows[i][j].second;
    }
  }
  return std::make_unique<mfem::SparseMatrix>(indices.release(), columns.release(),
                                              data.release(), height, height);
}
DenseMatrix ElementOperator::DiagonalBlock(int start, int size) const {
  DenseMatrix a(size);
  a = 0.;
  if(diagonal.Size()) { for(int i=0;i<size;++i) { a(i,i)=diagonal[start+i]; } }
  for (const auto &e : entries) {
    for (int i = 0; i < e.block.Height(); i++) {
      if (e.rows[i] >= start && e.rows[i] < start + size) {
        for (int j = 0; j < e.block.Width(); j++) {
          if (e.cols[j] >= start && e.cols[j] < start + size) {
            a(e.rows[i] - start, e.cols[j] - start) += e.block(i, j);
          }
        }
      }
    }
  }
  for (const auto &f : factors) {
    DenseMatrix product;
    f.left->Mult(*f.right,product);
    for (int i=0; i<int(f.ids.size()); ++i) {
      for (int j=0; j<int(f.ids.size()); ++j) {
        if (f.ids[i]>=start && f.ids[i]<start+size &&
            f.ids[j]>=start && f.ids[j]<start+size) {
          a(f.ids[i]-start,f.ids[j]-start)+=f.sign*product(i,j);
        }
      }
    }
  }
  return a;
}
std::vector<DenseMatrix> ElementOperator::DiagonalBlocks(int size) const {
  require(height % size == 0, "Invalid diagonal block size");
  std::vector<DenseMatrix> result;
  DiagonalBlocks(size, result);
  return result;
}
void ElementOperator::DiagonalBlocks(int size, std::vector<DenseMatrix> &result) const {
  require(height % size == 0, "Invalid diagonal block size");
  result.resize(height / size);
  for (auto &a : result) {
    a.SetSize(size);
    a = 0.;
  }
  if(diagonal.Size()) {
    for(int i=0;i<height;++i) { result[i/size](i%size,i%size)=diagonal[i]; }
  }
  for (const auto &e : entries) {
    for (int i = 0; i < e.block.Height(); i++) {
      for (int j = 0; j < e.block.Width(); j++) {
        if (e.rows[i] / size == e.cols[j] / size) {
          result[e.rows[i] / size](e.rows[i] % size, e.cols[j] % size) += e.block(i, j);
        }
      }
    }
  }
  DenseMatrix panel;
  for (const auto &f : factors) {
    require(size==4*f.left->nodes, "Schur Jacobi requires complete face blocks");
    for(int side=0;side<4;++side) {
      f.left->Panel(*f.right,side,side*size,size,panel);
      auto &a=result[f.ids[side*size]/size];
      a.Add(f.sign,panel);
    }
  }
}
void ElementOperator::PrepareScratch(int gather, int output, int intermediate) {
  local_scratch.SetSize(gather);
  output_scratch.SetSize(output);
  factor_scratch.SetSize(intermediate);
}
} // namespace euler
