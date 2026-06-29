# FSSH_by_C — Fewest-Switches Surface Hopping 的 C 语言实现

## 项目简介

本项目是 **表面跳跃法（Surface Hopping）** 中最经典算法 **FSSH（Fewest-Switches Surface Hopping）** 的纯 C 语言实现，用于模拟两态非绝热动力学过程中电子-核耦合的量子经典动力学。项目从 Python 原型重写为 C，并引入了 **多线程线程池** 以实现大规模批量轨迹的并行计算，大幅提升算力吞吐。

核心物理场景：粒子从左侧入射，经过两态交叉区域（avoided crossing），最终可能 **透射** 或 **反射**，并停留在 **高能态** 或 **低能态**。程序统计四种结果概率（TL / RL / TH / RH）随入射动量 k 的变化。

## 物理模型

当前激活的是 **模型一**（Tully 经典双态模型），在 `H_quantities.c` 中定义，另有模型二、模型三以注释形式保留，可切换：

| 模型 | 描述 | 参数 |
|------|------|------|
| 模型一（当前） | Tully Model I — 简单避免交叉 | A=0.01, B=1.6, C=0.005, D=1.0 |
| 模型二 | 双势阱型避免交叉 | A=0.1, B=0.28, E0=0.05, C=0.015, D=0.06 |
| 模型三 | 非对称耦合模型 | A=6e-4, B=0.1, C=0.9 |

透热哈密顿矩阵：

```
V11 = ±A·(1 - exp(∓Bx)),  V22 = -V11,  V12 = V21 = C·exp(-D·x²)
```

## 算法流程

每条 FSSH 轨迹的单步演化分为三个阶段（对应 Python 原型的 stepII / stepIII / stepIV）：

1. **Step II — 核动力学 + 电子振幅演化**
   - Velocity-Verlet 更新核坐标 `x` 和动量 `k`
   - 刷新当前坐标下的 H_di、U、H_ad、d12、H_eff
   - RKG 四阶方法积分电子振幅向量 `C(t)`

2. **Step III — 跃迁判定**
   - 计算跃迁概率 g12（低→高）或 g21（高→低）
   - 生成高精度随机数 ξ，若 g > ξ 则发生面跃迁

3. **Step IV — 跃迁修正**
   - 若发生跃迁：计算跃迁后动量 `k_aft`，动能不足则撤销跃迁
   - 更新核受力 F 和体系总能量

**终止条件**：核坐标超出 ±X_Max 边界，或步数超过 Step_Max。

## 项目结构

```
FSSH_c/
├── Makefile                    # 构建脚本（跨平台 Windows/Linux/macOS）
├── input_data                  # 输入文件：每行一个动量 k 值
├── output/
│   ├── result_data.txt         # 输出结果文件
│   └── main.exe                # 编译产物
├── .vscode/                    # VS Code 配置
│   ├── c_cpp_properties.json   # IntelliSense（w64devkit gcc）
│   ├── launch.json             # 调试配置
│   ├── tasks.json              # 构建/运行/清理任务
│   └── settings.json           # 编译器警告设置
├── include/                    # 头文件
│   ├── fssh_config.h           # 全局超参数结构体与接口
│   ├── fssh_traj.h             # 单条轨迹上下文 TrajContext 定义
│   ├── fssh_batch.h            # 批量运行与统计结果结构体
│   ├── fssh_thread_pool.h      # 线程池与并行批量接口
│   ├── H_quantities.h          # 物理量计算函数声明
│   └── matrix_math.h           # 2×2 复矩阵运算定义
└── src/                        # 源文件
    ├── main.c                  # 程序入口
    ├── fssh_config.c           # 全局配置默认值与 setter
    ├── fssh_traj.c             # 单条轨迹初始化、三步演化、终止判定
    ├── fssh_batch.c            # 批量模拟（串行/并行）、区间遍历、文件读取
    ├── fssh_thread_pool.c      # pthread 线程池实现
    ├── H_quantities.c          # 哈密顿量、受力、RKG积分、跃迁概率等核心物理
    └── matrix_math.c           # 矩阵乘法、共轭转置、对角化
```

## 模块说明

### matrix_math — 复矩阵数学基础
定义 `Mat2c`（2×2 复矩阵）和 `Vec2c`（2 维复向量）结构体，实现矩阵乘法、矩阵乘向量、共轭转置（dagger）、以及 2×2 厄米矩阵的解析对角化。

### H_quantities — 物理量计算
项目核心物理模块，包含：
- `H_di(x)`：构造透热哈密顿矩阵
- `dH_di_dx(x)`：中心差分计算 dH/dx
- `d12(x)`：非绝热耦合常数 = [U†·dH/dx·U]₁₂ / (E₂-E₁)
- `F(x, C)`：绝热态势能面上的核受力（中心差分）
- `x_change_by_t` / `k_change_by_t`：Velocity-Verlet 核动力学更新
- `H_eff`：有效哈密顿量 H_eff = H_ad - iħ·v·d12
- `C_change_by_t_RKG`：RKG 四阶方法积分电子振幅 ODE
- `g12` / `g21`：FSSH 跃迁概率
- `k_aft`：跃迁后动量修正（能量守恒），动能不足返回 -10086 拒绝跃迁

### fssh_traj — 单条轨迹演化
封装一条 FSSH 轨迹的完整生命周期：
- `traj_init`：初始化电子振幅、核坐标、动量、受力、总能量
- `traj_step_II/III/IV`：三阶段单步演化
- `traj_run_full`：循环演化至终止
- `traj_is_stop`：边界/步数终止判定
- `traj_get_result`：返回 [低能/高能, 透射/反射]

### fssh_config — 全局配置
全局唯一实例 `g_fssh_cfg`，默认参数：

| 参数 | 默认值 | 含义 |
|------|--------|------|
| Step_Max | 400000 | 单轨迹最大步数 |
| X_Max | 5.0 | 坐标逃逸边界 ±5.0 |
| N_Single_K_Traj | 4000 | 每个动量下模拟轨迹数 |
| M_Quality | 2000.0 | 粒子质量 |
| Delta_T | 0.1 | 时间步长 |

### fssh_batch — 批量模拟
- `batch_run_single_k(k)`：单动量串行批量模拟
- `batch_run_k_range(k1, k2, gap, path)`：动量区间遍历
- `batch_run_from_file(in_path, out_path)`：从文件读取动量列表

统计四类事件概率：**TL**（低能透射）、**RL**（低能反射）、**TH**（高能透射）、**RH**（高能反射）、**other**（超步数未出边界）。

### fssh_thread_pool — 多线程并行
基于 pthread 实现的任务队列线程池：
- 自动检测 CPU 核心数创建工作线程
- 每个线程独立随机种子（解决 `rand()` 线程不安全问题）
- 互斥锁保护任务队列和结果数组
- `batch_run_single_k_parallel(k)`：单动量下多核并行批量模拟

## 编译与运行

### 依赖
- GCC（支持 C99 + pthread）
- Windows 推荐 [w64devkit](https://github.com/skeeto/w64devkit)，Linux/macOS 使用系统 gcc

### 编译

```bash
make          # 编译，生成 output/main（Linux/macOS）或 output/main.exe（Windows）
make clean    # 清理编译产物
make run      # 编译并运行
```

### 运行

入口函数（`src/main.c`）当前配置：

```c
int main(){
    fssh_cfg_set_ntraj(4000);                              // 每个动量 4000 条轨迹
    batch_run_from_file("input_data", "output/result_data.txt");  // 从 input_data 读取动量
}
```

`input_data` 文件格式：每行一个动量值（浮点数），例如当前内容为 `13`。

### 输出格式

`output/result_data.txt`：

```
k	TL	RL	TH	RH	other
13.000000	0.738000	0.000000	0.262000	0.000000	0.000000
```

各列含义：
- **k**：入射动量
- **TL**：低能态透射概率
- **RL**：低能态反射概率
- **TH**：高能态透射概率
- **RH**：高能态反射概率
- **other**：超出最大步数未逃逸的比例

## VS Code 集成

项目附带完整的 VS Code 配置：
- `Ctrl+Shift+B`：构建（Windows 调用 `mingw32-make`，Linux/macOS 调用 `make`）
- F5 调试：使用 gdb，自动先执行 build 任务
- IntelliSense 指向 w64devkit 的 gcc

## 技术细节

- **随机数精度**：使用两次 `rand()` 组合，将精度从 15 位提升到 30 位（最小间隔 ≈ 9.3e-10）
- **数值微分**：受力和非绝热耦合均采用中心差分，步长 1e-7 ~ 1e-8
- **RKG 四阶方法**：使用 √2 相关系数的四阶 Runge-Kutta 变体积分电子振幅 ODE
- **能量守恒**：跃迁时检查动能是否足以支付势能差，不足则拒绝跃迁（返回 -10086 标记）
- **线程安全**：每个工作线程独立 srand 种子，结果写入加互斥锁
