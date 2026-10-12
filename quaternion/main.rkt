#lang racket/base

;; quaternion: float64 quaternion algebra and 3D rotation utilities.
;;
;; This module is the public API. Every export is guarded by a contract;
;; functions whose mathematics require unit quaternions enforce that
;; precondition via `unit-quaternion/c`. See
;; quaternion/scribblings/quaternion.scrbl for the full documentation.

(require racket/contract
         "private/core.rkt"
         "private/geometry.rkt"
         "private/convert.rkt"
         "private/interp.rkt")

(provide (contract-out
          ;; ---- construction and components ----
          [q (-> real? real? real? real? q?)]
          [q? (-> any/c boolean?)]
          [q-w (-> q? flonum?)]
          [q-x (-> q? flonum?)]
          [q-y (-> q? flonum?)]
          [q-z (-> q? flonum?)]
          [q-identity q?]

          ;; ---- arithmetic (any quaternion) ----
          [q+ (-> q? q? q?)]
          [q- (-> q? q? q?)]
          [q* (-> q? q? q?)]
          [q-conjugate (-> q? q?)]
          [q-norm (-> q? flonum?)]
          [q-normalize (-> q? q?)]
          [q-unit (-> q? q?)]
          [q=? (->* (q? q?) ((>=/c 0.0)) boolean?)]

          ;; ---- rotation (unit quaternion required) ----
          [q-rotate (-> flvector3/c unit-quaternion/c flvector3/c)]
          [q-from-axis-angle (-> flvector3/c real? q?)]
          [q-to-axis-angle (-> unit-quaternion/c (values flvector3/c (>=/c 0.0)))]

          ;; ---- Euler angles: intrinsic Z-Y-X (yaw-pitch-roll) ----
          [q-from-euler-zyx (-> real? real? real? q?)]
          [q-to-euler-zyx (-> unit-quaternion/c (values real? real? real?))]

          ;; ---- rotation matrices: row-major flvector of length 9 ----
          [q-from-matrix (-> flvector9/c q?)]
          [q-to-matrix (-> unit-quaternion/c flvector9/c)]

          ;; ---- interpolation (unit quaternions required) ----
          [q-slerp (-> unit-quaternion/c unit-quaternion/c (between/c 0.0 1.0) q?)]
          [q-nlerp (-> unit-quaternion/c unit-quaternion/c (between/c 0.0 1.0) q?)]

          ;; ---- exported contracts ----
          [unit-quaternion/c contract?]
          [flvector3/c contract?]
          [flvector9/c contract?]))
