#lang racket/base

;; Interpolation between unit quaternions.

(require racket/flonum
         "core.rkt")

;; Dot product; internal helper.
(define (q-dot* a b)
  (fl+ (fl* (q-w a) (q-w b))
       (fl* (q-x a) (q-x b))
       (fl* (q-y a) (q-y b))
       (fl* (q-z a) (q-z b))))

;; Blend a and b component-wise with weight t and renormalize. Internal
;; helper; t must already be a flonum in [0.0, 1.0].
(define (nlerp* a b t)
  (q-normalize
   (q (fl+ (q-w a) (fl* t (fl- (q-w b) (q-w a))))
      (fl+ (q-x a) (fl* t (fl- (q-x b) (q-x a))))
      (fl+ (q-y a) (fl* t (fl- (q-y b) (q-y a))))
      (fl+ (q-z a) (fl* t (fl- (q-z b) (q-z a)))))))

;; Spherical linear interpolation with constant angular velocity. Takes the
;; shortest path: when dot(a, b) < 0, b is negated first (the double cover).
;; When the endpoints are very close (dot > 0.9995) the trig formulation is
;; numerically ill-conditioned and normalized linear interpolation is used
;; instead. t = 0.0 and t = 1.0 return a and b unchanged. The result is a
;; unit quaternion.
(define (q-slerp a b t)
  (define tt (real->double-flonum t))
  (cond
    [(fl= tt 0.0) a]
    [(fl= tt 1.0) b]
    [else
     (define af (q-normalize a))
     (define b0 (q-normalize b))
     (define d (q-dot* af b0))
     (define bf (if (fl< d 0.0) (q-negate b0) b0))
     (define cos-theta (flmin 1.0 (flabs d)))
     (if (fl> cos-theta 0.9995)
         (nlerp* af bf tt)
         (let* ([theta (flacos cos-theta)]
                [sin-theta (flsin theta)]
                [wa (fl/ (flsin (fl* (fl- 1.0 tt) theta)) sin-theta)]
                [wb (fl/ (flsin (fl* tt theta)) sin-theta)])
           (q-normalize
            (q (fl+ (fl* wa (q-w af)) (fl* wb (q-w bf)))
               (fl+ (fl* wa (q-x af)) (fl* wb (q-x bf)))
               (fl+ (fl* wa (q-y af)) (fl* wb (q-y bf)))
               (fl+ (fl* wa (q-z af)) (fl* wb (q-z bf)))))))]))

;; Normalized linear interpolation: the cheapest interpolation between unit
;; quaternions (no trigonometry), at the cost of non-constant angular
;; velocity. Takes the shortest path like q-slerp. t = 0.0 and t = 1.0
;; return a and b unchanged.
(define (q-nlerp a b t)
  (define tt (real->double-flonum t))
  (cond
    [(fl= tt 0.0) a]
    [(fl= tt 1.0) b]
    [else
     (define af (q-normalize a))
     (define b0 (q-normalize b))
     (define bf (if (fl< (q-dot* af b0) 0.0) (q-negate b0) b0))
     (nlerp* af bf tt)]))

(provide (all-defined-out))
