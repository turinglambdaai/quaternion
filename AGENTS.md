# AGENTS.md

Float64 四元数代数，面向 3D 旋转的纯 Racket 实现（MIT）。

## 结构

Multi-collection 包：collection 即 `quaternion/` 子目录，scribblings 声明在 `quaternion/info.rkt`（collection 级）。

## 常用命令

```bash
raco pkg install --auto --link "$(pwd)"
raco test -p quaternion
```

## 验证与发版

- CI（`ci.yml`）：即上面两条命令，三平台跑。
- 发版 = tag `v<版本>` + CHANGELOG 同步；catalog 侧为手动提交（pkgd）。
