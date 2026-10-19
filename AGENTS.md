# For agents

## 项目约束

- 项目比较串行 MFEM C++17 Euler DG 与代数等价的凝聚 HDG。阅读 `MATHEMATICS.md` 理解离散化，阅读 `DESIGN.md` 理解模块与计时边界，阅读 `README.md` 查运行入口。
- 保持两种方法的时间积分、HLL 边界、求积和求解容差可比。优化不能悄悄改动数值方法或降低收敛标准。DG/HDG 自由度减少本身不是加速证据。
- 正确性、精度、性能分别使用 `code_test/`、`accuracy_test/`、`speed_test/`；不重新创建顶层 `results/`。运行脚本放在对应测试目录，通用静态检查保留在 `scripts/`。
- `code_test` 仅维护 smoke 测试和 backward Euler 单时间步、Newton 收敛后的 DG/HDG 更新一致性测试。不要恢复旧的代数、CLI、BSR 等回归矩阵。源码里的诊断能力不等于应加入日常测试套件；确有需要可针对故障调用。
- 精度研究采用同一 CFL 下的时空同步加密，保持终止时间与其他参数一致。不使用 RK4、固定超短 dt、仅空间或仅时间加密实验。不要仅凭 CSV 的 `dt` 列判定实际步长：CFL 会覆盖它，应检查 `mean_dt`、实际 CFL 和迭代记录。
- 运行中的速度作业需要保护。修改前检查提交入口、Slurm 队列、日志和依赖；本机 `ps` 无进程不能证明 Slurm 计算节点无作业。保护在途及待启动作业依赖的可执行文件、Python 模块和输出路径，不能在其构建目录重编译或清理。
- 不自动重新提交或取消 Slurm 作业，不覆盖历史实验输出。为新诊断选用独立输出目录。删除资料前按内容和参数确认范围，不以名字含 `smoke`、`pilot` 或 `T005` 为唯一依据。

## 按任务加载技能

本项目 `.agents/` 当前为只读目录，因此项目技能保存在 `skills/`，通过这里的明确引用加载；这不代表它们已安装为全局技能。

- 调查、修复求解器或构建问题：阅读 [debug](skills/debug/SKILL.md)。
- 调查收敛阶、终止时间误差或 DG/HDG 解差：阅读 [accuracy](skills/accuracy/SKILL.md)。
- 调查速度、内存、Slurm 作业或性能统计：阅读 [speed](skills/speed/SKILL.md)。
- 实施计算速度或内存优化：阅读 [optimize](skills/optimize/SKILL.md)。

## 常用验证与交付

从项目根目录运行命令。日常代码检查入口为 `python3 code_test/test_solver.py --exe <明确的可执行文件> --output code_test/results/<本次运行名>`。该驱动只有两项测试，不支持旧的 `--quick` / `--full-time`。配置和构建不执行测试；当前 CTest 仍含独立的 `algebra` 项，不要将其误认为上述两项测试。

Intel 构建入口是 `intel-build.sh`，支持 `BUILD_DIR` 与 `MFEM_BUILD` 覆盖。先检查脚本：它使用模块环境和 `cmake --fresh`。调试时选择独立构建目录，确保结果来自本次源码；运行旧二进制通过不能证明新源码通过。

记录实际命令、可执行文件来源、输出路径、成功或失败的检查及限制。区分已完成的实验、部分输出、历史数据与推断。不要把某一次迭代数或误差值设为未来所有版本必须满足的常数。保护已有工作区修改，以中文简洁说明结果。

## AGENTS.md 和 skills 修改规范

这两部分内容一般保持不动，作为 agent 的行为规范和指导方针。当 agent 认为很必要更新这两部分内容时（例如这里的规定或指导与未来的开发路径产生分歧时），应当向用户提出申请；在得到批准后进行更新。
