#lang racket/base
;; Runnable form of the README's 90-degree Z rotation.
(require racket/flonum quaternion)

(define q90z
  (q-from-axis-angle (flvector 0.0 0.0 1.0)
                     1.5707963267948966))
(define out (q-rotate (flvector 1.0 0.0 0.0) q90z))
(unless (and (<= (abs (flvector-ref out 0)) 1e-12)
             (<= (abs (- (flvector-ref out 1) 1.0)) 1e-12)
             (<= (abs (flvector-ref out 2)) 1e-12))
  (error 'quaternion-quickstart "90-degree rotation result differs"))

;; The public contract must reject a non-unit quaternion.
(define rejected?
  (with-handlers ([exn:fail:contract? (lambda (_) #t)])
    (q-rotate (flvector 1.0 0.0 0.0) (q 2.0 0.0 0.0 0.0))
    #f))
(unless rejected?
  (error 'quaternion-quickstart "unit-quaternion contract not enforced"))
(displayln "Quaternion quickstart OK")
