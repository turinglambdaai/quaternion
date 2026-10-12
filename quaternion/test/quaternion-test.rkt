#lang racket/base

;; Test suite for the quaternion package. Run with `raco test -p quaternion`.

(require rackunit
         racket/contract
         racket/flonum
         quaternion)

(define pi 3.141592653589793)
(define half-pi (* 0.5 pi))

;; ---- helpers ----

(define (check-q= actual expected [tol 1e-12])
  (check-true (q=? actual expected tol)
              (format "expected ~a, got ~a" expected actual)))

(define (flv->list v)
  (for/list ([i (in-range (flvector-length v))])
    (flvector-ref v i)))

(define (check-flvector= actual expected [tol 1e-12])
  (define a (flv->list actual))
  (define e (flv->list expected))
  (check-equal? (length a) (length e))
  (for ([av (in-list a)] [ev (in-list e)] [i (in-naturals)])
    (check-true (<= (abs (- av ev)) tol)
                (format "component ~a: expected ~a, got ~a" i ev av))))

(define (dot a b)
  (+ (* (q-w a) (q-w b))
     (* (q-x a) (q-x b))
     (* (q-y a) (q-y b))
     (* (q-z a) (q-z b))))

;; Double-cover-tolerant equality for unit quaternions: two representations
;; describe the same rotation exactly when |dot| = 1.
(define (same-rotation? a b [tol 1e-12])
  (<= (abs (- (abs (dot a b)) 1.0)) tol))

(define (neg a)
  (q (- (q-w a)) (- (q-x a)) (- (q-y a)) (- (q-z a))))

(define (matrix-apply m v)
  (define (row i)
    (+ (* (flvector-ref m (+ (* 3 i) 0)) (flvector-ref v 0))
       (* (flvector-ref m (+ (* 3 i) 1)) (flvector-ref v 1))
       (* (flvector-ref m (+ (* 3 i) 2)) (flvector-ref v 2))))
  (flvector (row 0) (row 1) (row 2)))

(define q-ident (q 1.0 0.0 0.0 0.0))
(define q-basis-i (q 0.0 1.0 0.0 0.0))
(define q-basis-j (q 0.0 0.0 1.0 0.0))
(define q-basis-k (q 0.0 0.0 0.0 1.0))
;; Arbitrary unit quaternion with all components nonzero.
(define q-abc (q-normalize (q 1.0 2.0 3.0 4.0)))

;; ---- 1. identities ----

(test-case "constructor normalizes to flonums"
  (define p (q 1 0 0 0))
  (check-true (q? p))
  (check-true (flonum? (q-w p)))
  (check-true (flonum? (q-x p)))
  (check-true (flonum? (q-y p)))
  (check-true (flonum? (q-z p)))
  (check-equal? p (q 1.0 0.0 0.0 0.0))
  (check-equal? (q-w (q 1/3 0 0 0)) 0.3333333333333333))

(test-case "identity element"
  (check-equal? q-identity (q 1.0 0.0 0.0 0.0))
  (check-q= (q* q-ident q-ident) q-ident)
  (check-q= (q* q-ident q-abc) q-abc)
  (check-q= (q* q-abc q-ident) q-abc))

(test-case "Hamilton basis products"
  (check-q= (q* q-basis-i q-basis-i) (q -1.0 0.0 0.0 0.0))
  (check-q= (q* q-basis-j q-basis-j) (q -1.0 0.0 0.0 0.0))
  (check-q= (q* q-basis-k q-basis-k) (q -1.0 0.0 0.0 0.0))
  (check-q= (q* q-basis-i q-basis-j) q-basis-k)
  (check-q= (q* q-basis-j q-basis-k) q-basis-i)
  (check-q= (q* q-basis-k q-basis-i) q-basis-j)
  (check-q= (q* q-basis-j q-basis-i) (q 0.0 0.0 0.0 -1.0))
  (check-q= (q* q-basis-k q-basis-j) (q 0.0 -1.0 0.0 0.0))
  (check-q= (q* q-basis-i q-basis-k) (q 0.0 0.0 -1.0 0.0)))

(test-case "q times its conjugate"
  (check-q= (q* q-abc (q-conjugate q-abc)) q-ident)
  ;; norm(q * conj(q)) = norm(q)^2, also for non-unit quaternions
  (define p (q 1.0 2.0 3.0 4.0))
  (check-true (<= (abs (- (q-norm (q* p (q-conjugate p)))
                          (* (q-norm p) (q-norm p))))
                  1e-9)))

(test-case "arithmetic basics"
  (check-q= (q+ (q 1.0 2.0 3.0 4.0) (q 0.5 -1.0 0.25 0.0))
            (q 1.5 1.0 3.25 4.0))
  (check-q= (q- (q 1.0 2.0 3.0 4.0) (q 1.0 2.0 3.0 4.0))
            (q 0.0 0.0 0.0 0.0)))

(test-case "q-normalize and the q-unit alias"
  (check-q= (q-normalize (q 2.0 0.0 0.0 0.0)) q-ident)
  (check-q= (q-unit (q 2.0 0.0 0.0 0.0)) q-ident)
  (check-q= (q-unit (q 0.0 3.0 4.0 0.0)) (q 0.0 0.6 0.8 0.0))
  (check-true (<= (abs (- (q-norm (q-normalize q-abc)) 1.0)) 1e-12)))

(test-case "double cover: q and -q rotate identically"
  (define v (flvector 0.3 -1.7 2.2))
  (check-flvector= (q-rotate v q-abc) (q-rotate v (neg q-abc)))
  (define p (q-from-axis-angle (flvector 1.0 1.0 0.0) 0.8))
  (check-flvector= (q-rotate v p) (q-rotate v (neg p)))
  (check-true (same-rotation? q-abc (neg q-abc))))

;; ---- 2. known rotations ----

(test-case "90 degrees about Z maps (1,0,0) to (0,1,0)"
  (define q90z (q-from-axis-angle (flvector 0.0 0.0 1.0) half-pi))
  (check-flvector= (q-rotate (flvector 1.0 0.0 0.0) q90z)
                   (flvector 0.0 1.0 0.0))
  (check-flvector= (q-rotate (flvector 0.0 1.0 0.0) q90z)
                   (flvector -1.0 0.0 0.0))
  (check-flvector= (q-rotate (flvector 1.0 2.0 3.0) q90z)
                   (flvector -2.0 1.0 3.0)))

(test-case "180 degrees about X maps (0,1,0) to (0,-1,0)"
  (define q180x (q-from-axis-angle (flvector 1.0 0.0 0.0) pi))
  (check-flvector= (q-rotate (flvector 0.0 1.0 0.0) q180x)
                   (flvector 0.0 -1.0 0.0))
  (check-flvector= (q-rotate (flvector 0.0 0.0 1.0) q180x)
                   (flvector 0.0 0.0 -1.0))
  (check-flvector= (q-rotate (flvector 1.0 0.0 0.0) q180x)
                   (flvector 1.0 0.0 0.0)))

;; ---- 3. round trips ----

(test-case "axis-angle round trip"
  (for* ([axis (in-list (list (flvector 1.0 0.0 0.0)
                              (flvector 0.0 1.0 0.0)
                              (flvector 0.0 0.0 1.0)
                              (flvector 1.0 -2.0 3.0)))]
         [angle (in-list (list 0.3 1.0 2.5))])
    (define p (q-from-axis-angle axis angle))
    (define-values (axis2 angle2) (q-to-axis-angle p))
    (check-true (<= (abs (- angle2 angle)) 1e-12)
                (format "angle: expected ~a, got ~a" angle angle2))
    (check-q= (q-from-axis-angle axis2 angle2) p)))

(test-case "to-axis-angle edge cases"
  ;; identity: zero angle, documented arbitrary axis
  (define-values (axis0 angle0) (q-to-axis-angle q-ident))
  (check-true (<= angle0 1e-12))
  (check-flvector= axis0 (flvector 1.0 0.0 0.0))
  ;; 180 degrees: axis recovered exactly
  (define-values (axisp anglep)
    (q-to-axis-angle (q-from-axis-angle (flvector 0.0 1.0 0.0) pi)))
  (check-true (<= (abs (- anglep pi)) 1e-12))
  (check-flvector= axisp (flvector 0.0 1.0 0.0))
  ;; q and -q decompose to the same axis-angle pair
  (define p (q-from-axis-angle (flvector 1.0 1.0 0.0) 0.8))
  (define-values (ax1 an1) (q-to-axis-angle p))
  (define-values (ax2 an2) (q-to-axis-angle (neg p)))
  (check-flvector= ax1 ax2)
  (check-true (<= (abs (- an1 an2)) 1e-12)))

(test-case "quaternion -> matrix -> quaternion"
  (for ([p (in-list (list (q-from-axis-angle (flvector 1.0 0.0 0.0) 0.3)
                          (q-from-axis-angle (flvector 0.0 1.0 0.0) 1.0)
                          (q-from-axis-angle (flvector 0.0 0.0 1.0) 2.5)
                          (q-from-axis-angle (flvector 1.0 -2.0 3.0) 0.7)
                          (q-from-euler-zyx 0.5 0.3 0.2)))])
    (define m (q-to-matrix p))
    (define p2 (q-from-matrix m))
    (check-true (same-rotation? p p2)
                (format "not the same rotation: ~a vs ~a" p p2))
    (check-flvector= (q-to-matrix p2) m)))

(test-case "matrix rotation agrees with q-rotate"
  (define p (q-from-axis-angle (flvector 1.0 -2.0 3.0) 0.7))
  (define v (flvector 0.4 -1.0 2.0))
  (check-flvector= (matrix-apply (q-to-matrix p) v) (q-rotate v p)))

(test-case "identity matrix"
  (check-flvector= (q-to-matrix q-ident)
                   (flvector 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0)))

;; ---- 4. Euler angles (intrinsic Z-Y-X, yaw-pitch-roll) ----

(test-case "euler ZYX matches Rz*Ry*Rx composition"
  (define y 0.7)
  (define p 0.4)
  (define r -0.2)
  (check-q= (q-from-euler-zyx y p r)
            (q* (q-from-axis-angle (flvector 0.0 0.0 1.0) y)
                (q* (q-from-axis-angle (flvector 0.0 1.0 0.0) p)
                    (q-from-axis-angle (flvector 1.0 0.0 0.0) r)))))

(test-case "euler ZYX round trip (no gimbal lock)"
  (for ([e (in-list (list (list 0.0 0.0 0.0)
                          (list 0.5 0.3 0.2)
                          (list -2.5 0.9 1.1)
                          (list 1.0 -0.75 3.0)
                          (list 3.1 0.4 -3.1)))])
    (define p (q-from-euler-zyx (list-ref e 0) (list-ref e 1) (list-ref e 2)))
    (define-values (yaw pitch roll) (q-to-euler-zyx p))
    (check-true (<= (abs (- yaw (list-ref e 0))) 1e-12)
                (format "yaw: expected ~a, got ~a" (list-ref e 0) yaw))
    (check-true (<= (abs (- pitch (list-ref e 1))) 1e-12)
                (format "pitch: expected ~a, got ~a" (list-ref e 1) pitch))
    (check-true (<= (abs (- roll (list-ref e 2))) 1e-12)
                (format "roll: expected ~a, got ~a" (list-ref e 2) roll))
    (check-q= (q-from-euler-zyx yaw pitch roll) p)))

(test-case "gimbal lock at +90 degrees pitch: exact recovery with roll := 0"
  (define p1 (q-from-euler-zyx 0.9 half-pi 0.0))
  (define-values (yaw1 pitch1 roll1) (q-to-euler-zyx p1))
  (check-true (<= (abs (- pitch1 half-pi)) 1e-12))
  (check-true (<= roll1 1e-12) (format "roll should be 0.0, got ~a" roll1))
  (check-flvector= (q-to-matrix (q-from-euler-zyx yaw1 pitch1 roll1))
                   (q-to-matrix p1)))

(test-case "gimbal lock at -90 degrees pitch: exact recovery with roll := 0"
  (define p1 (q-from-euler-zyx 0.9 (- half-pi) 0.4))
  (define-values (yaw1 pitch1 roll1) (q-to-euler-zyx p1))
  (check-true (<= (abs (+ pitch1 half-pi)) 1e-12))
  (check-true (<= roll1 1e-12))
  (check-flvector= (q-to-matrix (q-from-euler-zyx yaw1 pitch1 roll1))
                   (q-to-matrix p1)))

;; ---- 5. interpolation ----

(define rot-x-90 (q-from-axis-angle (flvector 1.0 0.0 0.0) half-pi))

(test-case "slerp and nlerp endpoints"
  (define a (q-from-axis-angle (flvector 0.0 0.0 1.0) 0.4))
  (define b (q-from-axis-angle (flvector 1.0 1.0 0.0) 1.2))
  (check-q= (q-slerp a b 0.0) a)
  (check-q= (q-slerp a b 1.0) b)
  (check-q= (q-nlerp a b 0.0) a)
  (check-q= (q-nlerp a b 1.0) b))

(test-case "slerp midpoint of perpendicular half-quaternions has sqrt(1/2) angles"
  ;; Between the identity and a 180-degree rotation about X the quaternion
  ;; angle is pi/2; halfway lies the 90-degree rotation about X, whose scalar
  ;; part is cos(pi/4) = sqrt(1/2).
  (check-q= (q-slerp q-ident (q 0.0 1.0 0.0 0.0) 0.5) rot-x-90)
  (check-q= (q-nlerp q-ident (q 0.0 1.0 0.0 0.0) 0.5) rot-x-90))

(test-case "slerp results are unit quaternions"
  (define a (q-from-axis-angle (flvector 1.0 2.0 3.0) 0.7))
  (define b (q-from-axis-angle (flvector 0.0 1.0 1.0) 2.0))
  (for ([t (in-list (list 0.0 0.1 0.25 0.5 0.75 0.9 1.0))])
    (define m (q-slerp a b t))
    (check-true (<= (abs (- (q-norm m) 1.0)) 1e-12)
                (format "t=~a not unit: norm ~a" t (q-norm m)))))

(test-case "slerp and nlerp take the shortest path (double cover)"
  (define a (q-from-axis-angle (flvector 0.0 0.0 1.0) 0.2))
  (define b (q-from-axis-angle (flvector 0.0 0.0 1.0) 1.0))
  (check-q= (q-slerp a (neg b) 0.5) (q-slerp a b 0.5))
  (check-q= (q-nlerp a (neg b) 0.5) (q-nlerp a b 0.5))
  ;; halfway between rotations of 0.2 and 1.0 about the same axis is 0.6
  (check-q= (q-slerp a b 0.5) (q-from-axis-angle (flvector 0.0 0.0 1.0) 0.6)))

(test-case "nlerp stays normalized"
  (define a (q-from-axis-angle (flvector 1.0 2.0 3.0) 0.7))
  (define b (q-from-axis-angle (flvector 0.0 1.0 1.0) 2.0))
  (for ([t (in-list (list 0.0 0.2 0.5 0.8 1.0))])
    (check-true (<= (abs (- (q-norm (q-nlerp a b t)) 1.0)) 1e-12))))

;; ---- 6. contracts and error messages ----

(test-case "non-unit quaternion inputs are rejected"
  (check-exn #rx"unit quaternion"
             (lambda () (q-rotate (flvector 1.0 0.0 0.0) (q 2.0 0.0 0.0 0.0))))
  (check-exn #rx"unit quaternion"
             (lambda () (q-slerp (q 2.0 0.0 0.0 0.0) q-ident 0.5)))
  (check-exn #rx"unit quaternion"
             (lambda () (q-slerp q-ident (q 1.0 1.0 1.0 1.0) 0.5)))
  (check-exn #rx"unit quaternion"
             (lambda () (q-nlerp (q 0.0 0.0 0.0 0.0) q-ident 0.5)))
  (check-exn #rx"unit quaternion"
             (lambda () (q-to-axis-angle (q 0.0 0.0 0.0 0.0))))
  (check-exn #rx"unit quaternion"
             (lambda () (q-to-euler-zyx (q 2.0 0.0 0.0 0.0))))
  (check-exn #rx"unit quaternion"
             (lambda () (q-to-matrix (q 2.0 0.0 0.0 0.0)))))

(test-case "malformed arguments raise contract errors"
  (check-exn exn:fail:contract?
             (lambda () (q-rotate (flvector 1.0 0.0) q-ident)))
  (check-exn exn:fail:contract?
             (lambda () (q-from-axis-angle (flvector 0.0 0.0 0.0) 1.0)))
  (check-exn exn:fail:contract?
             (lambda () (q-normalize (q 0.0 0.0 0.0 0.0))))
  (check-exn exn:fail:contract?
             (lambda () (q-slerp q-ident q-ident 1.5)))
  (check-exn exn:fail:contract?
             (lambda () (q=? q-ident q-ident -1.0))))

(test-case "q=? tolerance"
  (check-true (q=? q-ident (q 1.0 0.0 0.0 0.0)))
  (check-false (q=? (q 1.0 0.0 0.0 0.0) (q 1.0 1e-9 0.0 0.0)))
  (check-true (q=? (q 1.0 0.0 0.0 0.0) (q 1.0 1e-9 0.0 0.0) 1e-8))
  (check-true (q=? q-abc (q-normalize (q 1.0 2.0 3.0 4.0)))))

(test-case "exported contracts"
  (check-pred contract? unit-quaternion/c)
  (check-pred contract? flvector3/c)
  (check-pred contract? flvector9/c))
