#lang scribble/manual

@(require scribble/eval
          (for-label quaternion
                     racket/base
                     racket/contract
                     racket/flonum))

@(define the-eval (make-base-eval '(require quaternion racket/flonum)))

@title{quaternion: Float64 Quaternions for 3D Rotation}
@author{turinglambdaai}

@defmodule[quaternion]

This library provides double-precision (flonum) quaternion algebra for 3D
rotation tasks common in graphics, robotics, and aerospace: Hamilton
products, axis-angle construction, intrinsic Z-Y-X Euler angles, 3x3
rotation matrices, and shortest-path interpolation.

Everything is implemented in pure Racket on top of IEEE-754 double
precision floats; the package depends on nothing beyond
@racketmodname[racket/base].

@section[#:tag "conventions"]{Conventions}

@itemlist[
 @item{A quaternion is @racket[(q w x y z)]: @racket[w] is the scalar part
  and @racket[(x y z)] the vector part. All components are flonums; the
  @racket[q] constructor converts arbitrary real numbers.}
 @item{Multiplication follows Hamilton's convention. Rotations are active
  rotations in a right-handed coordinate system: rotating vector @tt{v} by
  unit quaternion @tt{a} computes @tt{a v a*}, where @tt{a*} is the
  conjugate.}
 @item{Angles are in radians.}
 @item{Functions whose mathematics require unit quaternions take a
  @racket[unit-quaternion/c] contract: the norm must be within @tt{1e-8} of
  @racket[1.0]. Those functions additionally re-normalize their input
  defensively. Pure arithmetic (@racket[q+], @racket[q-], @racket[q*],
  @racket[q-conjugate], @racket[q-norm], @racket[q-normalize]) accepts any
  quaternion.}
 @item{Euler angles use the intrinsic Z-Y-X convention (yaw-pitch-roll);
  see @secref["euler"].}
 @item{Rotation matrices are a single row-major @racket[flvector] of length
  nine: @tt{[m00 m01 m02 m10 m11 m12 m20 m21 m22]}; see
  @secref["matrices"].}
]

@section[#:tag "construction"]{Construction and Components}

@defproc[(q [w real?] [x real?] [y real?] [z real?]) q?]{
  Constructs a quaternion from four real numbers. Every component is
  converted to a double-precision flonum.

@examples[#:eval the-eval
  (q 1 0 0 0)
  (q 1/2 0.25 2 3)]
}

@defthing[q-identity q?]{
  The multiplicative identity @racket[(q 1.0 0.0 0.0 0.0)].
}

@defproc[(q? [v any/c]) boolean?]{
  Returns @racket[#t] if @racket[v] is a quaternion.
}

@defproc*[
  (((q-w [a q?]) flonum?)
   ((q-x [a q?]) flonum?)
   ((q-y [a q?]) flonum?)
   ((q-z [a q?]) flonum?))
]{
  Component accessors: @racket[q-w] returns the scalar part; the others
  return the vector part.
}

@section{Arithmetic}

@defproc[(q+ [a q?] [b q?]) q?]{Component-wise sum.}

@defproc[(q- [a q?] [b q?]) q?]{Component-wise difference @tt{a - b}.}

@defproc[(q* [a q?] [b q?]) q?]{
  Hamilton product. For unit quaternions it composes rotations:
  @racket[(q* a b)] applies @racket[b] first, then @racket[a].

@examples[#:eval the-eval
  (define i (q 0.0 1.0 0.0 0.0))
  (define j (q 0.0 0.0 1.0 0.0))
  (q* i i)
  (q* i j)]
}

@defproc[(q-conjugate [a q?]) q?]{
  Conjugate @tt{(w -x -y -z)}. For a unit quaternion the conjugate is also
  the multiplicative inverse. This library standardizes on the name
  @racket[q-conjugate]; there is no @tt{q-conj} alias.
}

@defproc[(q-norm [a q?]) flonum?]{Euclidean length of the four components.}

@defproc[(q-normalize [a q?]) q?]{
  Scales @racket[a] to unit length. Raises @racket[exn:fail:contract] when
  given the zero quaternion.
}

@defproc[(q-unit [a q?]) q?]{
  Documented alias for @racket[q-normalize].
}

@defproc[(q=? [a q?] [b q?] [tolerance (>=/c 0.0) 1e-12]) boolean?]{
  Component-wise equality: @racket[#t] when the absolute differences of all
  four components are at most @racket[tolerance]. Note that a quaternion
  and its negation represent the same rotation even though they differ
  component-wise; compare rotations as rotations via @racket[q-rotate] or
  the matrix conversions when the double cover matters.
}

@section[#:tag "rotation"]{Axis-Angle and Rotation}

@defproc[(q-rotate [v flvector3/c] [a unit-quaternion/c]) flvector3/c]{
  Rotates the 3-vector @racket[v] by the unit quaternion @racket[a]
  (equivalent to @tt{a v a*}). Returns a fresh @racket[flvector].

@examples[#:eval the-eval
  (q-rotate (flvector 1.0 0.0 0.0)
            (q-from-axis-angle (flvector 0.0 0.0 1.0) 1.5707963267948966))]
}

@defproc[(q-from-axis-angle [axis flvector3/c] [angle real?]) q?]{
  Builds the unit quaternion for a rotation of @racket[angle] radians about
  @racket[axis]. The axis is normalized by this function, so it need not be
  a unit vector; a zero-length axis raises @racket[exn:fail:contract].
}

@defproc[(q-to-axis-angle [a unit-quaternion/c])
         (values flvector3/c (>=/c 0.0))]{
  Decomposes the unit quaternion @racket[a] into a unit axis and an angle
  in @tt{[0, pi]}. Because a quaternion and its negation encode the same
  rotation, the decomposition is canonicalized to the representative with a
  non-negative scalar part. For a (numerically) zero rotation the angle is
  @racket[0.0] and the axis is @racket[(flvector 1.0 0.0 0.0)].
}

@section[#:tag "euler"]{Euler Angles}

This library implements exactly one Euler convention: @bold{intrinsic
Z-Y-X}, also known as @bold{yaw-pitch-roll},

@verbatim|{
  R = Rz(yaw) * Ry(pitch) * Rx(roll)
}|

with yaw about Z, pitch about Y, and roll about X (right-handed axes,
radians).

@bold{Gimbal lock.} At @tt{pitch = +pi/2} or @tt{pitch = -pi/2} yaw and
roll are not independent (the convention degenerates). Rather than reject
such input, @racket[q-to-euler-zyx] reports @racket[0.0] for the roll and
recovers the single determined combination of yaw and roll exactly, so that
@racket[(q-from-euler-zyx (q-to-euler-zyx a))] reproduces the original
rotation even in the locked case. The lock branch triggers when
@tt{|cos(pitch)| < 1e-6}.

@defproc[(q-from-euler-zyx [yaw real?] [pitch real?] [roll real?]) q?]{
  Builds the unit quaternion for @tt{R = Rz(yaw) Ry(pitch) Rx(roll)}.

@examples[#:eval the-eval
  (define q180x (q-from-euler-zyx 0.0 0.0 3.141592653589793))
  (q-rotate (flvector 0.0 1.0 0.0) q180x)]
}

@defproc[(q-to-euler-zyx [a unit-quaternion/c]) (values real? real? real?)]{
  Decomposes the unit quaternion into three angles
  @racket[(values yaw pitch roll)] with yaw and roll in @tt{(-pi, pi]} and
  pitch in @tt{[-pi/2, pi/2]}. See above for the gimbal-lock behavior.
}

@section[#:tag "matrices"]{Rotation Matrices}

Matrices are represented as a single row-major @racket[flvector] of length
nine:

@verbatim|{
  (flvector m00 m01 m02
            m10 m11 m12
            m20 m21 m22)
}|

applied as @tt{v' = M v} with column vectors — the same rotation that
@racket[q-rotate] performs.

@defproc[(q-to-matrix [a unit-quaternion/c]) flvector9/c]{
  Converts a unit quaternion to its 3x3 rotation matrix.
}

@defproc[(q-from-matrix [m flvector9/c]) q?]{
  Converts a rotation matrix to a quaternion using Shepperd's method
  (branching on the largest component for numerical stability). The result
  is normalized.

  The contract checks only the shape of @racket[m]; it is assumed to be a
  @emph{proper} rotation matrix (orthonormal with determinant @tt{+1}).
  Results are undefined for scaled, reflected, or otherwise malformed
  matrices.
}

@section{Interpolation}

@defproc[(q-slerp [a unit-quaternion/c] [b unit-quaternion/c]
                  [t (between/c 0.0 1.0)]) q?]{
  Spherical linear interpolation between unit quaternions with constant
  angular velocity along the great circle between @racket[a] and
  @racket[b].

  The double cover is handled the way users expect: when @tt{dot(a, b) < 0}
  the endpoint @racket[b] is negated first, so interpolation always takes
  the shortest path. When the endpoints are very close
  (@tt{dot(a, b) > 0.9995}) the trigonometric formulation is numerically
  ill-conditioned and normalized linear interpolation is used instead.
  @racket[t] = @racket[0.0] and @racket[t] = @racket[1.0] return @racket[a]
  and @racket[b] unchanged. The result is a unit quaternion.

@examples[#:eval the-eval
  (q-slerp q-identity (q 0.0 1.0 0.0 0.0) 0.0)
  (q-slerp q-identity (q 0.0 1.0 0.0 0.0) 1.0)]
}

@defproc[(q-nlerp [a unit-quaternion/c] [b unit-quaternion/c]
                  [t (between/c 0.0 1.0)]) q?]{
  Normalized linear interpolation: blends component-wise and renormalizes.
  Cheaper than @racket[q-slerp] (no trigonometry) at the cost of
  non-constant angular velocity; reasonable for small angles or when frame
  smoothness matters more than exact speed. Takes the shortest path like
  @racket[q-slerp]; endpoints are returned unchanged.
}

@section{Contracts}

@defthing[unit-quaternion/c contract?]{
  Flat contract recognizing quaternions whose norm is within @tt{1e-8} of
  @racket[1.0]. Use it to guard your own functions that expect unit
  quaternions.
}

@defthing[flvector3/c contract?]{
  Flat contract for an @racket[flvector] of length three (a 3-vector).
}

@defthing[flvector9/c contract?]{
  Flat contract for an @racket[flvector] of length nine (a 3x3 matrix in
  row-major order).
}

@(close-eval the-eval)
