# quaternion

**English** · [中文](README.zh-CN.md)

[![CI](https://github.com/turlamai/quaternion/actions/workflows/ci.yml/badge.svg)](https://github.com/turlamai/quaternion/actions/workflows/ci.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Float64 quaternion algebra for 3D rotations, in pure Racket. Built for
graphics, robotics, and aerospace work: everything runs on IEEE-754 double
precision floats (flonums) with no dependencies beyond `base`.

The package covers Hamilton algebra, axis-angle construction, intrinsic
ZYX Euler angles, 3x3 rotation matrices, and shortest-path interpolation
(`slerp`/`nlerp`) — all behind straightforward function contracts that
enforce unit-quaternion preconditions where the math demands them.

## Install

From a local checkout:

```console
$ raco pkg install --auto --link /path/to/quaternion
```

Or straight from GitHub:

```console
$ raco pkg install --auto --name quaternion https://github.com/turlamai/quaternion
```

## Example

```racket
#lang racket/base
(require racket/flonum quaternion)

;; 90-degree rotation about Z
(define q90z (q-from-axis-angle (flvector 0.0 0.0 1.0) 1.5707963267948966))

(q-rotate (flvector 1.0 0.0 0.0) q90z)
;; => #fl(2.220446049250313e-16 1.0 0.0) — i.e. (0, 1, 0), up to flonum noise
```

## Conventions

- **Product**: Hamilton product (`q*`). Rotations are **active rotations**
  in a **right-handed** coordinate system: `q-rotate` computes `q v q*`.
- **Units**: `q-rotate`, `q-to-axis-angle`, `q-to-euler-zyx`, `q-to-matrix`,
  `q-slerp`, and `q-nlerp` require **unit quaternions** and enforce it with
  the `unit-quaternion/c` contract (norm within `1e-8` of `1.0`); they
  re-normalize their input defensively. Pure arithmetic (`q+`, `q-`, `q*`,
  `q-conjugate`, `q-norm`, `q-normalize`) accepts any quaternion. For unit
  quaternions, the conjugate is the inverse.
- **Euler**: exactly one convention — **intrinsic Z-Y-X (yaw-pitch-roll)**,
  `R = Rz(yaw) Ry(pitch) Rx(roll)`. `q-to-euler-zyx` returns yaw/roll in
  `(-pi, pi]` and pitch in `[-pi/2, pi/2]`. At ±90° pitch the convention
  hits gimbal lock; the decomposition reports roll `0.0` and recovers the
  determined yaw exactly, so the round trip through `q-from-euler-zyx`
  still reproduces the original rotation.
- **Matrices**: a single **row-major** `flvector` of length nine
  (`[m00 m01 m02 m10 m11 m12 m20 m21 m22]`), applied as `v' = M v` with
  column vectors. `q-from-matrix` assumes a *proper* rotation matrix
  (orthonormal, determinant +1); only the shape is checked.
- **Axis-angle**: radians; the axis is normalized for you (zero axis is an
  error). `q-to-axis-angle` canonicalizes the angle to `[0, pi]`.
- **Naming**: `q-conjugate` is the one conjugate name (no `q-conj` alias);
  `q-unit` is a documented alias for `q-normalize`; `q-identity` is the
  multiplicative identity.
- **Comparison**: `q=?` compares component-wise with a tolerance (default
  `1e-12`). Remember `q` and `-q` are the same rotation — compare
  rotations as rotations (e.g. via `q-rotate`) when the double cover
  matters.

## Development

The [Scribble API manual](quaternion/scribblings/quaternion.scrbl) is the
normative reference for contracts, axis-angle/Euler conventions, matrices and
interpolation. `flvector` comes from `racket/flonum` (not
`racket/base`), so the import above is necessary when copying the example.
The [executable quickstart](examples/quickstart.rkt) checks the 90-degree
rotation within a floating-point tolerance and the non-unit contract error.

From the repository root after installing Racket:

```sh
raco pkg install --auto --link "$(pwd)"
raco test -p quaternion
racket examples/quickstart.rkt
raco setup --check-pkg-deps -l quaternion
```

CI runs all steps on Linux, macOS and Windows. Report a minimal reproduction
and Racket version at [GitHub Issues](https://github.com/turlamai/quaternion/issues).

## License

MIT — see [LICENSE](LICENSE).
