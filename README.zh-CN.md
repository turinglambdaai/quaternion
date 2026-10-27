# quaternion

[English](README.md) · **中文**

[![CI](https://github.com/turinglambdaai/quaternion/actions/workflows/ci.yml/badge.svg)](https://github.com/turinglambdaai/quaternion/actions/workflows/ci.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

纯 Racket 实现的 Float64 四元数代数，用于三维旋转。面向图形、机器人与航空航天场景：所有计算都运行在 IEEE-754 双精度浮点数（flonum）上，除 `base` 外无其他依赖。

本包涵盖哈密顿代数、轴角构造、内旋 ZYX 欧拉角、3x3 旋转矩阵与最短路径插值（`slerp`/`nlerp`）——全部以直白的函数契约收口，在数学上要求单位四元数的地方强制前置条件。

## 安装

从本地检出的仓库安装：

```console
$ raco pkg install --auto --link /path/to/quaternion
```

或直接从 GitHub 安装：

```console
$ raco pkg install --auto --name quaternion https://github.com/turinglambdaai/quaternion
```

## 示例

```racket
#lang racket/base
(require racket/flonum quaternion)

;; 90-degree rotation about Z
(define q90z (q-from-axis-angle (flvector 0.0 0.0 1.0) 1.5707963267948966))

(q-rotate (flvector 1.0 0.0 0.0) q90z)
;; => #fl(2.220446049250313e-16 1.0 0.0) — i.e. (0, 1, 0), up to flonum noise
```

## 约定

- **乘积**：哈密顿积（`q*`）。旋转是**右手**坐标系下的**主动旋转（active rotation）**：`q-rotate` 计算 `q v q*`。
- **单位**：`q-rotate`、`q-to-axis-angle`、`q-to-euler-zyx`、`q-to-matrix`、`q-slerp` 与 `q-nlerp` 要求**单位四元数**，并以 `unit-quaternion/c` 契约强制（范数与 `1.0` 相差在 `1e-8` 以内）；它们还会防御性地对输入重新归一化。纯算术函数（`q+`、`q-`、`q*`、`q-conjugate`、`q-norm`、`q-normalize`）接受任意四元数。对单位四元数而言，共轭即逆。
- **欧拉角**：只用一种约定——**内旋 Z-Y-X（偏航-俯仰-滚转）**，`R = Rz(yaw) Ry(pitch) Rx(roll)`。`q-to-euler-zyx` 返回的 yaw/roll 位于 `(-pi, pi]`，pitch 位于 `[-pi/2, pi/2]`。pitch 到达 ±90° 时该约定进入万向节锁；分解结果将 roll 报告为 `0.0` 并精确还原确定的 yaw，因此经 `q-from-euler-zyx` 的往返仍能复现原旋转。
- **矩阵**：单个**行主序** `flvector`，长度为九（`[m00 m01 m02 m10 m11 m12 m20 m21 m22]`），按列向量方式应用为 `v' = M v`。`q-from-matrix` 假定输入是*正规*（proper）旋转矩阵（正交且行列式为 +1）；只检查形状。
- **轴角**：弧度制；轴向量会替你归一化（零轴是错误）。`q-to-axis-angle` 把角度规范化到 `[0, pi]`。
- **命名**：共轭只有 `q-conjugate` 一个名字（无 `q-conj` 别名）；`q-unit` 是 `q-normalize` 有文档记载的别名；`q-identity` 是乘法单位元。
- **比较**：`q=?` 按分量带容差比较（默认 `1e-12`）。注意 `q` 与 `-q` 表示同一个旋转——当双覆盖（double cover）有影响时，应把旋转当作旋转来比较（例如通过 `q-rotate`）。

## 开发

[权威 Scribble API 手册](quaternion/scribblings/quaternion.scrbl)定义函数契约、
轴角/欧拉角约定、矩阵和插值行为。注意 `flvector` 属于
`racket/flonum`，并非 `racket/base`；复制示例时必须引入该模块。
[可执行入门示例](examples/quickstart.rkt)会以浮点容差验证 90° 旋转，
并验证非单位四元数的契约错误。

安装 Racket 后，在仓库根目录执行：

```sh
raco pkg install --auto --link "$(pwd)"
raco test -p quaternion
racket examples/quickstart.rkt
raco setup -l quaternion --check-pkg-deps
```

CI 在 Linux、macOS、Windows 上运行上述步骤。反馈时请在
[GitHub Issues](https://github.com/turinglambdaai/quaternion/issues) 附上
脱敏最小样例与 Racket 版本。

## 许可证

MIT——见 [LICENSE](LICENSE)。
