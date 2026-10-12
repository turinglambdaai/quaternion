#lang racket/base

;; Conversions between quaternions and Euler angles / rotation matrices.
;;
;; Euler convention (the only one implemented): intrinsic Z-Y-X, a.k.a.
;; yaw-pitch-roll,   R = Rz(yaw) * Ry(pitch) * Rx(roll)
;; with yaw about Z, pitch about Y, roll about X (right-handed, radians).

(require racket/flonum
         "core.rkt")

;; Double precision closest to pi. Kept local so the package depends only on
;; `base` (racket/math lives in math-lib).
(define pi 3.141592653589793)
(define half-pi (* 0.5 pi))

;; Intrinsic Z-Y-X (yaw-pitch-roll) to quaternion:
;; q = qz(yaw) * qy(pitch) * qx(roll)
(define (q-from-euler-zyx yaw pitch roll)
  (define hy (fl* 0.5 (real->double-flonum yaw)))
  (define hp (fl* 0.5 (real->double-flonum pitch)))
  (define hr (fl* 0.5 (real->double-flonum roll)))
  (q* (q (flcos hy) 0.0 0.0 (flsin hy))
      (q* (q (flcos hp) 0.0 (flsin hp) 0.0)
          (q (flcos hr) (flsin hr) 0.0 0.0))))

;; Quaternion to intrinsic Z-Y-X (yaw-pitch-roll) Euler angles.
;;
;; Returns (values yaw pitch roll) with yaw and roll in (-pi, pi] and pitch
;; in [-pi/2, pi/2]. At pitch = +/- pi/2 the convention hits gimbal lock (yaw
;; and roll are not independent); this reports roll = 0.0 and recovers the
;; one determined combination exactly, so that round-tripping through
;; q-from-euler-zyx reproduces the original rotation. The lock branch
;; triggers when |cos(pitch)| < 1e-6.
(define (q-to-euler-zyx a)
  (define u (q-normalize a))
  (define w (q-w u))
  (define x (q-x u))
  (define y (q-y u))
  (define z (q-z u))
  (define sin-pitch (flmin 1.0 (flmax -1.0 (fl* 2.0 (fl- (fl* w y) (fl* x z))))))
  (define cos-pitch (flsqrt (flmax 0.0 (fl- 1.0 (fl* sin-pitch sin-pitch)))))
  (cond
    [(fl< cos-pitch 1e-6)
     ;; Gimbal lock: only yaw - roll (pitch = +pi/2) or yaw + roll
     ;; (pitch = -pi/2) is determined. By convention roll := 0.0.
     (values (atan (fl* sin-pitch (fl* 2.0 (fl- (fl* y z) (fl* w x))))
                   (fl- 1.0 (fl* 2.0 (fl+ (fl* x x) (fl* z z)))))
             (if (fl>= sin-pitch 0.0) half-pi (- half-pi))
             0.0)]
    [else
     (values (atan (fl* 2.0 (fl+ (fl* w z) (fl* x y)))
                   (fl- 1.0 (fl* 2.0 (fl+ (fl* y y) (fl* z z)))))
             (flasin sin-pitch)
             (atan (fl* 2.0 (fl+ (fl* w x) (fl* y z)))
                   (fl- 1.0 (fl* 2.0 (fl+ (fl* x x) (fl* y y))))))]))

;; Unit quaternion to 3x3 rotation matrix as a row-major flvector of length
;; 9: [m00 m01 m02 m10 m11 m12 m20 m21 m22], applying as v' = M v (column
;; vectors) — the same rotation as q-rotate.
(define (q-to-matrix a)
  (define u (q-normalize a))
  (define w (q-w u))
  (define x (q-x u))
  (define y (q-y u))
  (define z (q-z u))
  (define xx (fl* x x))
  (define yy (fl* y y))
  (define zz (fl* z z))
  (define xy (fl* x y))
  (define xz (fl* x z))
  (define yz (fl* y z))
  (define wx (fl* w x))
  (define wy (fl* w y))
  (define wz (fl* w z))
  (flvector (fl- 1.0 (fl* 2.0 (fl+ yy zz)))
            (fl* 2.0 (fl- xy wz))
            (fl* 2.0 (fl+ xz wy))
            (fl* 2.0 (fl+ xy wz))
            (fl- 1.0 (fl* 2.0 (fl+ xx zz)))
            (fl* 2.0 (fl- yz wx))
            (fl* 2.0 (fl- xz wy))
            (fl* 2.0 (fl+ yz wx))
            (fl- 1.0 (fl* 2.0 (fl+ xx yy)))))

;; 3x3 rotation matrix (row-major flvector of length 9) to quaternion via
;; Shepperd's method, branching on the largest component for numerical
;; stability. The public contract checks the shape only: m is assumed to be
;; a proper rotation matrix (orthonormal, determinant +1); results are
;; otherwise undefined.
(define (q-from-matrix m)
  (define m00 (flvector-ref m 0))
  (define m01 (flvector-ref m 1))
  (define m02 (flvector-ref m 2))
  (define m10 (flvector-ref m 3))
  (define m11 (flvector-ref m 4))
  (define m12 (flvector-ref m 5))
  (define m20 (flvector-ref m 6))
  (define m21 (flvector-ref m 7))
  (define m22 (flvector-ref m 8))
  (define tr (fl+ m00 m11 m22))
  (cond
    [(fl> tr 0.0)
     (define s (fl* 2.0 (flsqrt (flmax 0.0 (fl+ 1.0 tr)))))
     (q-normalize (q (fl* 0.25 s)
                     (fl/ (fl- m21 m12) s)
                     (fl/ (fl- m02 m20) s)
                     (fl/ (fl- m10 m01) s)))]
    [(and (fl> m00 m11) (fl> m00 m22))
     (define s (fl* 2.0 (flsqrt (flmax 0.0 (fl- (fl+ 1.0 m00) m11 m22)))))
     (q-normalize (q (fl/ (fl- m21 m12) s)
                     (fl* 0.25 s)
                     (fl/ (fl+ m01 m10) s)
                     (fl/ (fl+ m02 m20) s)))]
    [(fl> m11 m22)
     (define s (fl* 2.0 (flsqrt (flmax 0.0 (fl- (fl+ 1.0 m11) m00 m22)))))
     (q-normalize (q (fl/ (fl- m02 m20) s)
                     (fl/ (fl+ m01 m10) s)
                     (fl* 0.25 s)
                     (fl/ (fl+ m12 m21) s)))]
    [else
     (define s (fl* 2.0 (flsqrt (flmax 0.0 (fl- (fl+ 1.0 m22) m00 m11)))))
     (q-normalize (q (fl/ (fl- m10 m01) s)
                     (fl/ (fl+ m02 m20) s)
                     (fl/ (fl+ m12 m21) s)
                     (fl* 0.25 s)))]))

(provide (all-defined-out))
