# Q8 Block-Jacobi 存储估算

适用算例：二维 Euler DG，`256 x 256` 网格，`Q8`，SSP-DIRK(4,3)，
GMRES restart 20，双精度，以及每个单元一个完整稠密 Block-Jacobi 块。
这是一份容量估算；它不替代 Slurm 的最终 `MaxRSS` 测量。

## 基本规模

每单元有 `4 x (8 + 1)^2 = 324` 个自由度，共有 `256^2 = 65,536` 个
单元，因此全局自由度数为

```
N = 65,536 x 324 = 21,233,664.
```

一个全局双精度向量占 `8N = 162 MiB`。GMRES(20) 的 Arnoldi 基有 21 个
向量，约占 `3.32 GiB`。

## 主导项与总量

每个单元的完整 Jacobi 块为 `324 x 324` 个 `double`。全部块的存储量是

```
65,536 x 324^2 x 8 bytes = 51.26 GiB.
```

| 组成 | 原位 LU | Jacobian 与 LU 分开保存 |
|---|---:|---:|
| 单元 Block-Jacobi 矩阵或 LU 因子 | 51.26 GiB | 102.52 GiB |
| LU pivot | 0.08 GiB | 0.08 GiB |
| GMRES(20) Arnoldi 基 | 3.32 GiB | 3.32 GiB |
| 其他 GMRES、DIRK 与 Newton 工作向量 | 2.69–3.95 GiB | 2.69–3.95 GiB |
| **合计（未计额外局部工作区）** | **57.35–58.61 GiB** | **108.61–109.87 GiB** |

因此，内存瓶颈是完整单元块的 `O(N_elem (k+1)^4)` 稠密存储，而不是
GMRES restart 或四阶段 DIRK 本身。

## 与当前实现和实测的关系

当前 partial-DG 路径不显式保存全局 Jacobian；Block-Jacobi 的单元对角块
在 LU 分解时原位覆盖。因而上表的“原位 LU”列对应当前实现的主要存储模型。

已观察到 `k=8`、Block-Jacobi 算例的进程峰值 RSS 为 **58.99 GiB**。这与
57–59 GiB 的原位-LU估算相符，但只能说明数量级一致，不能单独证明不存在
额外临时副本。应以完成后的 `sacct`/`MaxRSS` 为准；申请内存仍需为分配器、
局部工作区与运行波动留出余量。

## 边界条件

- 若额外保存全局稀疏 Jacobian、原始块和 LU 的两份副本，或同时保留多个阶段
  的预条件器，必须在上表基础上增加相应存储。
- 该估算不适用于 point-block 或 line-block Jacobi；它们会改变预条件器的
  存储与收敛性质。

## Speed-test 运行记录

每个作业分别执行 DG 与 HDG。结果收集时，从对应 case 目录的
`cfl5_dg_resources.json`、`cfl5_hdg_resources.json` 和 `resources.csv` 填写
峰值 RSS 与墙钟时间；完成后再用 `sacct` 核对作业级 `MaxRSS`。`进程最大值`
取 DG、HDG 两项峰值中的较大者。

同一份 JSON 还记录每个求解器进程的用户 CPU 时间和系统 CPU 时间。CPU 效率
按 `(user_cpu_s + system_cpu_s) / wall_s` 计算；接近 `100%` 表示单线程计算
过程基本持续占满一个逻辑 CPU。

### 节点 CPU 信息

| node | CPU 型号 | 逻辑 CPU 数 | sockets | cores/socket | threads/core | NUMA nodes | 采集日期 |
|---|---|---:|---:|---:|---:|---:|---|
| node1601 | Intel(R) Xeon(R) Gold 6242 CPU @ 2.80GHz | 32 | 2 | 16 | 1 | 2 | 2026-10-19 |

| job ID | k | preconditioner | node | 申请内存 (GiB) | DG peak RSS (GiB) | HDG peak RSS (GiB) | 进程最大值 (GiB) | DG CPU/wall (s) | HDG CPU/wall (s) | DG/HDG CPU效率 | 状态/退出码 | 备注 |
|---:|---:|---|---|---:|---:|---:|---:|---|---|---|---|---|
| 7583428 | 8 | block-jacobi | node1601 | 364 | — | — | 58.99 (运行中采样) | — | — | — | RUNNING | 最终值待作业结束后填写 |
| 7583429 | 8 | block-ilu0 | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583430 | 7 | block-jacobi | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583431 | 7 | block-ilu0 | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583432 | 6 | block-jacobi | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583433 | 6 | block-ilu0 | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583434 | 5 | block-jacobi | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583435 | 5 | block-ilu0 | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583436 | 4 | block-jacobi | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583437 | 4 | block-ilu0 | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583438 | 3 | block-jacobi | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583439 | 3 | block-ilu0 | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583440 | 2 | block-jacobi | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583441 | 2 | block-ilu0 | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583442 | 1 | block-jacobi | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583443 | 1 | block-ilu0 | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583444 | 0 | block-jacobi | — | 120 | — | — | — | — | — | — | PENDING | |
| 7583445 | 0 | block-ilu0 | — | 120 | — | — | — | — | — | — | PENDING | |
