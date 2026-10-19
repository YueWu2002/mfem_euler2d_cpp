---
name: debug
description: 调查和修复本仓库 MFEM Euler DG/HDG 求解器、Newton 更新一致性和构建问题，按项目限定的两项代码测试验证。
---

先阅读根目录 `AGENTS.md`。本技能的文件路径均相对项目根目录。

## 定位与复现

查看工作区差异、实际失败命令、日志及二进制来源。阅读 `DESIGN.md` 的模块表，只深入与故障相关的实现：

- 参数、输出与最终解比较：`src/options.cpp`、`src/main.cpp`。
- CFL、阶段状态、Newton 停止与更新检查：`src/stepper.cpp`。
- 初始投影、物理状态和残差：`src/problem.cpp`、`src/physics.cpp`。
- DG 作用与 HDG 恢复：`src/partial.cpp`、`src/operators.cpp`、`src/condensation.cpp`。
- GMRES 与预条件：`src/linear.cpp`、`src/block_ilu.cpp`。

需要数学判断时查 `MATHEMATICS.md`，不要通过放宽容差掩盖不一致。若只请求诊断，先给出原因和证据；实现修复需处于用户授权范围。

## 构建与最小验证

先核对 `intel-build.sh`、`CMakeLists.txt` 和现有构建的编译器/MFEM 配置。需要编译时可使用 `BUILD_DIR=build-debug-session bash intel-build.sh`，实际目录应独立且不被作业引用。它是优化的 Intel 构建入口，目录名含 debug 不代表开启了调试编译选项。

运行现有 `code_test/test_solver.py`，明确指定本次二进制与新的输出目录。保持常用 BLAS/OpenMP 线程数为 1。套件范围严格为：

1. CFL 选步的 smoke test。
2. backward Euler、一个时间步、Newton 按收敛准则停止，比较 DG 与 HDG 最终解，并用 `--check-update` 对比 HDG 恢复更新与独立 DG 更新。

核对 BE 的 `integrator=be`、`time_steps=1`、`fixed_newton_steps=-1`、成功退出和非空更新记录。不能使用 `--newton-steps` 代替 Newton 收敛。`update_check_relative` 是每次线性更新的比较；最终 `*_comparison.csv` 是独立非线性演化后的解差，两者含义不同。

目前驱动使用 8×8 网格。更粗网格可能在初始投影就不满足正性；遇到该错误应先检查初始化条件。迭代 CSV 未必记录终止检查的残差，不要直接将最后一次更新前残差称为最终收敛残差。

若额外的定向诊断确实有助于定位故障，可显式执行并解释目的；不要因此扩大长期维护的 code_test 套件。报告哪些检查针对当前源码，哪些仅针对现有二进制。
