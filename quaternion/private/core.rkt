#lang racket/base

;; Core quaternion type, arithmetic, and shared contracts. Everything here is
;; internal; the public API is curated by quaternion/main.rkt.

(require racket/contract
         racket/flonum)

;; ---- struct ----

;; Raw transparent record. The type is named `quat` because Racket's `struct`
;; form binds its name, and the public constructor (below) needs that name.
;; A custom writer keeps the printed form aligned with the public API.
(struct quat (w x y z) #:transparent
  #:constructor-name make-quat
  #:methods gen:custom-write
  [(define (write-proc v port mode)
     (fprintf port "(q ~a ~a ~a ~a)"
              (quat-w v) (quat-x v) (quat-y v) (quat-z v)))])

;; Public constructor: accepts any real numbers, normalizes every component
;; to a double-precision flonum.
(define (q w x y z)
  (make-quat (real->double-flonum w)
             (real->double-flonum x)
             (real->double-flonum y)
             (real->double-flonum z)))

(define q? quat?)

(define (q-w a) (quat-w a))
(define (q-x a) (quat-x a))
(define (q-y a) (quat-y a))
(define (q-z a) (quat-z a))

;; Multiplicative identity.
(define q-identity (q 1.0 0.0 0.0 0.0))

;; ---- arithmetic ----

(define (q+ a b)
  (q (+ (q-w a) (q-w b))
     (+ (q-x a) (q-x b))
     (+ (q-y a) (q-y b))
     (+ (q-z a) (q-z b))))

(define (q- a b)
  (q (- (q-w a) (q-w b))
     (- (q-x a) (q-x b))
     (- (q-y a) (q-y b))
     (- (q-z a) (q-z b))))

;; Hamilton product.
(define (q* a b)
  (q (- (* (q-w a) (q-w b))
        (* (q-x a) (q-x b))
        (* (q-y a) (q-y b))
        (* (q-z a) (q-z b)))
     (- (+ (* (q-w a) (q-x b))
           (* (q-x a) (q-w b))
           (* (q-y a) (q-z b)))
        (* (q-z a) (q-y b)))
     (- (+ (* (q-w a) (q-y b))
           (* (q-y a) (q-w b))
           (* (q-z a) (q-x b)))
        (* (q-x a) (q-z b)))
     (- (+ (* (q-w a) (q-z b))
           (* (q-z a) (q-w b))
           (* (q-x a) (q-y b)))
        (* (q-y a) (q-x b)))))

(define (q-conjugate a)
  (q (q-w a) (- (q-x a)) (- (q-y a)) (- (q-z a))))

;; Negation. Internal helper (not part of the public API): -q represents the
;; same rotation as q (the double cover).
(define (q-negate a)
  (q (- (q-w a)) (- (q-x a)) (- (q-y a)) (- (q-z a))))

(define (q-norm a)
  (flsqrt (fl+ (fl* (q-w a) (q-w a))
               (fl* (q-x a) (q-x a))
               (fl* (q-y a) (q-y a))
               (fl* (q-z a) (q-z a)))))

(define (q-normalize a)
  (define n (q-norm a))
  (if (fl= n 0.0)
      (raise-argument-error 'q-normalize "quaternion with nonzero norm" a)
      (q (fl/ (q-w a) n)
         (fl/ (q-x a) n)
         (fl/ (q-y a) n)
         (fl/ (q-z a) n))))

;; Documented alias for q-normalize.
(define q-unit q-normalize)

(define (q=? a b [tolerance 1e-12])
  (and (<= (abs (- (q-w a) (q-w b))) tolerance)
       (<= (abs (- (q-x a) (q-x b))) tolerance)
       (<= (abs (- (q-y a) (q-y b))) tolerance)
       (<= (abs (- (q-z a) (q-z b))) tolerance)))

;; ---- contracts ----

;; How far from 1.0 the norm of a "unit" quaternion may drift while still
;; passing the unit-quaternion contract.
(define unit-tolerance 1e-8)

(define (unit-q? v)
  (and (q? v)
       (<= (abs (- (q-norm v) 1.0)) unit-tolerance)))

(define unit-quaternion/c
  (flat-named-contract
   "unit quaternion (norm within 1e-8 of 1.0)"
   unit-q?))

(define (flvector3? v)
  (and (flvector? v) (= 3 (flvector-length v))))

(define (flvector9? v)
  (and (flvector? v) (= 9 (flvector-length v))))

(define flvector3/c (flat-named-contract "flvector of length 3" flvector3?))
(define flvector9/c (flat-named-contract "flvector of length 9" flvector9?))

(provide (all-defined-out))
