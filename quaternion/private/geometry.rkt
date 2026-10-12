#lang racket/base

;; Axis-angle construction and active rotation of 3-vectors.

(require racket/flonum
         "core.rkt")

;; Builds a unit quaternion for a rotation of `angle` radians about `axis`.
;; The axis is normalized here, so it need not be a unit vector; a
;; zero-length axis is an error.
(define (q-from-axis-angle axis angle)
  (define ax (flvector-ref axis 0))
  (define ay (flvector-ref axis 1))
  (define az (flvector-ref axis 2))
  (define n (flsqrt (fl+ (fl* ax ax) (fl* ay ay) (fl* az az))))
  (when (fl= n 0.0)
    (raise-argument-error 'q-from-axis-angle
                          "flvector of length 3 with nonzero norm"
                          axis))
  (define half (fl* 0.5 (real->double-flonum angle)))
  (define s (flsin half))
  (define c (flcos half))
  (q c
     (fl* s (fl/ ax n))
     (fl* s (fl/ ay n))
     (fl* s (fl/ az n))))

;; Decomposes a unit quaternion into (values axis angle) with the angle in
;; [0, pi]. The double cover is resolved by picking the representative with a
;; non-negative scalar part. For a (numerically) zero rotation the axis is
;; (flvector 1.0 0.0 0.0) and the angle 0.0.
(define (q-to-axis-angle a)
  (define u (q-normalize a))
  (define rep (if (fl< (q-w u) 0.0) (q-negate u) u))
  (define w (q-w rep))
  (if (fl>= w 1.0)
      (values (flvector 1.0 0.0 0.0) 0.0)
      (let* ([angle (fl* 2.0 (flacos w))]
             [s (flsqrt (flmax 0.0 (fl- 1.0 (fl* w w))))])
        (if (fl< s 1e-12)
            (values (flvector 1.0 0.0 0.0) angle)
            (values (flvector (fl/ (q-x rep) s)
                              (fl/ (q-y rep) s)
                              (fl/ (q-z rep) s))
                    angle)))))

;; Rotates the 3-vector v by the unit quaternion a, i.e. computes a v a*
;; using the explicit expansion (no intermediate quaternions).
(define (q-rotate v a)
  (define u (q-normalize a))
  (define w (q-w u))
  (define x (q-x u))
  (define y (q-y u))
  (define z (q-z u))
  (define vx (flvector-ref v 0))
  (define vy (flvector-ref v 1))
  (define vz (flvector-ref v 2))
  ;; t = 2 * cross(u_vector, v);  v' = v + w * t + cross(u_vector, t)
  (define tx (fl* 2.0 (fl- (fl* y vz) (fl* z vy))))
  (define ty (fl* 2.0 (fl- (fl* z vx) (fl* x vz))))
  (define tz (fl* 2.0 (fl- (fl* x vy) (fl* y vx))))
  (flvector (fl+ vx (fl* w tx) (fl- (fl* y tz) (fl* z ty)))
            (fl+ vy (fl* w ty) (fl- (fl* z tx) (fl* x tz)))
            (fl+ vz (fl* w tz) (fl- (fl* x ty) (fl* y tx)))))

(provide (all-defined-out))
