---
name: speed
description: 调查本仓库 DG/HDG 速度和内存实验、Slurm 运行状态与计时统计，保护在途作业及其依赖。
---

先阅读根目录 `AGENTS.md`。阅读当前 `speed_test/submit.sh`、`job.sh`、`run_case.py`、`measure.py`，再按需查 `collect.py`、`benchmark.py`、`DESIGN.md` 和 `MEMORY_ESTIMATE.md`。路径均相对项目根目录。

## 确认在途任务

以提交清单中的作业 ID、Slurm 队列/作业详情/记账结果和对应日志交叉确认。登录节点 `ps` 只能看到本机；调度命令不存在、超时或无输出都不能独立证明远程任务已结束。无法确认时按仍在使用处理，并报告状态未知。

保护运行与排队作业随后可能读取的文件：作业脚本、Python 驱动及导入模块、`speed_test/build/` 内的二进制、提交清单、`speed_test/slurm/` 和 case 输出目录。不要重建该二进制、移动路径或修改驱动以影响后续 case。源码调试另用独立构建和输出目录。

只有任务明确要求时才提交、重提或取消作业；不能为了查询状态再次执行 `submit.sh`。先阅读其当前实现，注意提交清单可能被覆盖。不要固定记忆节点名、内存申请、作业 ID 或度数顺序。

## 解释结果

`mean_newton_s` 衡量计时边界内每次接受 Newton 更新的平均耗时。它与进程 user/system CPU 时间、time-loop wall time、进程总 wall time 是不同指标；即使表名包含 CPU，也应按实际字段解释。保持 DG/HDG 工作量、线程数、容差与方法相同。

速度任务可以使用固定更新次数作为受控工作量，这不能证明 Newton 已收敛。`--check-update` 的额外求解会影响缓存和耗时，因此不混入正式速度实验。比较前检查期望更新数、成功状态及两种方法的完整结果。

每个进程的峰值 RSS 由 `measure.py` 单独测量，包含 setup。失败进程的峰值只代表该失败运行，不能当作完成运行的需求。区分矩阵 payload 估算与完整峰值 RSS。无结果或正在写入的文件不能填成零值。

`collect.py` 会写汇总文件，带 `--update-results` 时还写 `RESULTS.md`；只读状态查询不要调用它。需要用户要求的汇总时，先确认 case 文件完整，再执行并报告未完成项。性能判断引用实际原始文件和机器/线程条件，不以 DOF 数或理论估算代替测量。

## 结果落盘与后处理

每个 speed case 写入独立目录
`speed_test/results/k<k>_<preconditioner>/`，避免不同 degree 或预条件器互相覆盖。
该目录中的 `result.csv` 提供 DG/HDG 的正式性能行，`resources.csv` 提供两种 form
的资源汇总；`cfl5_<dg|hdg>_resources.json` 记录单个求解器进程的
`peak_rss_kib`、`peak_rss_gib`、user/system CPU 时间、wall time、状态和命令。
求解器的 summary、iteration CSV 与文本日志也留在同一 case 目录。

Slurm 标准输出和错误日志位于 `speed_test/slurm/%x-%j.out` 与 `.err`。当前
`run_case.py` 捕获子进程输出，完成一个 form 后才把它写入 case 日志；因此运行中
`.out` 只有 prolog 或没有阶段信息并不表示求解器没有进展。

当用户明确表示全部任务已结束并要求收集结果时：

1. 用 Slurm 终态和每个 case 的 `result.csv`、`resources.csv` 交叉确认完成、退出码与
   DG/HDG 双行完整性；保留失败或缺失项，不填零。
2. 在项目根目录执行 `python speed_test/collect.py`。它写入
   `speed_test/results/CPU time per Newton step.csv` 与同名 `.md`。
3. 只有用户要求同步根目录总结果页时，执行
   `python speed_test/collect.py --update-results`；该操作会更新 `RESULTS.md` 的
   speed-test 标记区。
4. 将实际值回填根目录 `MEMORY_ESTIMATE.md` 的 Speed-test 运行记录：节点、DG/HDG
   峰值 RSS、两者最大值、DG/HDG CPU 时间与 wall time、单线程 CPU 效率、作业状态及
   最终 `sacct` 的 `MaxRSS`。CPU 效率为
   `(user_cpu_s + system_cpu_s) / wall_s`；不要把作业申请内存当作实际峰值。

若用户要查看运行中的资源，优先读取 job step 的 `sstat`；作业结束后再使用 `sacct`
取得最终 `MaxRSS`、`TotalCPU` 和退出码。节点硬件信息只需首次记录一次，放入
`MEMORY_ESTIMATE.md` 的节点 CPU 信息表。
