#include "euler.hpp"

#include <algorithm>
#include <limits>

namespace euler {
void ElementOperator::Add(const std::vector<int> &r, const std::vector<int> &c,
                          const DenseMatrix &b) {
  entries.push_back({r, c, b});
}
void ElementOperator::Mult(const Vector &x, Vector &y) const {
  y.SetSize(height);
  y = 0.;
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
    mfem::Mult(*f.left, *f.right, b);
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
    for (int i = 0; i < int(f.ids.size()); i++) {
      if (f.ids[i] >= start && f.ids[i] < start + size) {
        for (int j = 0; j < int(f.ids.size()); j++) {
          if (f.ids[j] >= start && f.ids[j] < start + size) {
            double z = 0;
            for (int h = 0; h < f.left->Width(); h++) {
              z += (*f.left)(i, h) * (*f.right)(h, j);
            }
            a(f.ids[i] - start, f.ids[j] - start) += f.sign * z;
          }
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
  for (const auto &e : entries) {
    for (int i = 0; i < e.block.Height(); i++) {
      for (int j = 0; j < e.block.Width(); j++) {
        if (e.rows[i] / size == e.cols[j] / size) {
          result[e.rows[i] / size](e.rows[i] % size, e.cols[j] % size) += e.block(i, j);
        }
      }
    }
  }
  for (const auto &f : factors) {
    // A cell restriction lists its four complete face blocks consecutively.
    // Multiply only matching face panels. Column-major i-inner traversal reads
    // C's columns contiguously instead of taking a strided dot per entry.
    bool panels = f.ids.size() % size == 0;
    for (int first = 0; panels && first < int(f.ids.size()); first += size) {
      for (int i = 0; i < size; ++i) {
        panels = panels && f.ids[first] % size == 0 && f.ids[first+i] == f.ids[first]+i;
      }
    }
    if (panels) {
      for (int first = 0; first < int(f.ids.size()); first += size) {
        auto &block = result[f.ids[first]/size];
        for (int j = 0; j < size; ++j) {
          for (int h = 0; h < f.left->Width(); ++h) {
            const double x = f.sign*(*f.right)(h,first+j);
            for (int i = 0; i < size; ++i) { block(i,j) += (*f.left)(first+i,h)*x; }
          }
        }
      }
    } else {
      // General gather/scatter fallback for independent algebra checks.
      for (int i = 0; i < int(f.ids.size()); ++i) {
        for (int j = 0; j < int(f.ids.size()); ++j) {
          if (f.ids[i]/size != f.ids[j]/size) { continue; }
          double value = 0.;
          for (int h = 0; h < f.left->Width(); ++h) { value += (*f.left)(i,h)*(*f.right)(h,j); }
          result[f.ids[i]/size](f.ids[i]%size,f.ids[j]%size) += f.sign*value;
        }
      }
    }
  }
}
void ElementOperator::PrepareScratch(int gather, int output, int intermediate) {
  local_scratch.SetSize(gather);
  output_scratch.SetSize(output);
  factor_scratch.SetSize(intermediate);
}
} // namespace euler
