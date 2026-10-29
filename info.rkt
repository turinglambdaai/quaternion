#lang info

(define collection 'multi)
(define version "0.1")
(define deps '("base" "rackunit-lib"))
(define build-deps '("scribble-lib" "racket-doc" "rackunit-lib"))
(define pkg-desc "Float64 quaternion algebra for 3D rotations: Hamilton product, axis-angle, Euler ZYX, rotation matrices, slerp.")
(define pkg-authors '("turlamai"))
(define license 'MIT)
(define test-omit-paths '("quaternion/scribblings"))
