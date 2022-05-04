;; Prologue
;; :PROPERTIES:
;; :header-args+: :tangle ../src/fdg/prologue.clj :comments org
;; :END:


(ns fdg.prologue
  (:refer-clojure :exclude [+ - * / = compare zero? ref partial
                            numerator denominator])
  (:require [sicmutils.env
             :as e :refer :all
             :exclude [Lagrange-equations Gamma]]))
;; Functional Abstraction

;; But this corrected use of Leibniz notation is ugly. We had to introduce
;; extraneous symbols ($q$ and $\dot{q}$) in order to indicate the argument
;; position specifying the partial derivative. Nothing would change here if we
;; replaced $q$ and $\dot{q}$ by $a$ and $b$.[fn:3] We can simplify the notation by
;; admitting that the partial derivatives of the Lagrangian are themselves new
;; functions, and by specifying the particular partial derivative by the position
;; of the argument that is varied

;; $$\frac{d}{d t}\left(\left(\partial_{2} L\right)\left(t, w(t), \frac{d}{d t}
;; w(t)\right)\right)-\left(\partial_{1} L\right)\left(t, w(t), \frac{d}{d t}
;; w(t)\right)=0,$$

;; where $\partial_{i}L$ is the function which is the partial derivative of the
;; function L with respect to the ith argument.[fn:4]

;; Two different notions of derivative appear in this expression. The functions
;; $\partial_2 L$ $\partial_1 L$, constructed from the Lagrangian $L$, have the
;; same arguments as $L$.

;; The derivative $d/dt$ is an expression derivative. It applies to an expression
;; that involves the variable $t$ and it gives the rate of change of the value of
;; the expression as the value of the variable $t$ is varied.

;; These are both useful interpretations of the idea of a derivative. But functions
;; give us more power. There are many equivalent ways to write expressions that
;; compute the same value. For example $1/(1/r_1 + 1/r_2)=(r_1r_2)/(r_1 + r_2)$.
;; These expressions compute the same function of the two variables $r_1$ and
;; $r_2$. The first expression fails if $r_1 = 0$ but the second one gives the
;; right value of the function. If we abstract the function, say as $\Pi(r_1,
;; r_2)$, we can ignore the details of how it is computed. The ideas become clearer
;; because they do not depend on the detailed shape of the expressions.

;; So let’s get rid of the expression derivative $d/dt$ and replace it with an
;; appropriate functional derivative. If $f$ is a function then we will write $Df$
;; as the new function that is the derivative of $f$:[fn:5]

;; $$(D f)(t)=\left.\frac{d}{d x} f(x)\right|_{x=t}.$$

;; To do this for the Lagrange equation we need to construct a function to take the
;; derivative of.

;; Given a configuration-space path $w$, there is a standard way to make the
;; state-space path. We can abstract this method as a mathematical function
;; $\Gamma$:

;; $$\Gamma[w](t)=\left(t, w(t), \frac{d}{d t} w(t)\right).$$

;; Using $\Gamma$ we can write:

;; $$\frac{d}{dt}\left(\left(\partial_{2} L\right) \left(\Gamma[w](t)\right)
;; \right) - \left(\partial_{1} L\right) \left(\Gamma[w](t)\right)=0.$$

;; If we now define composition of functions $(f \circ g)(x) = f(g(x))$, we can
;; express the Lagrange equations entirely in terms of functions:

;; $$D\left(\left(\partial_{2} L\right) \circ \left(\Gamma[w]\right)\right)
;; \\ -\left(\partial_{1} L\right) \circ \left(\Gamma[w]\right)=0.$$

;; The functions $\partial_1 L$ and $\partial_2 L$ are partial derivatives of the
;; function $L$. Composition with $\Gamma[w]$ evaluates these partials with
;; coordinates and velocites appropriate for the path $w$, making functions of
;; time. Applying $D$ takes the time derivative. The Lagrange equation states that
;; the difference of the resulting functions of time must be zero. This statement
;; of the Lagrange equation is complete, unambiguous, and functional. It is not
;; encumbered with the particular choices made in expressing the Lagrangian. For
;; example, it doesn’t matter if the time is named $t$ or $\tau$, and it has an
;; explicit place for the path to be tested.

;; This expression is equivalent to a computer program:[fn:6]


(declare Gamma)


;; #+RESULTS:
;; : #'fdg.prologue/Gamma


(defn Lagrange-equations [Lagrangian]
  (fn [w]
    (- (D (compose ((partial 2) Lagrangian) (Gamma w)))
       (compose ((partial 1) Lagrangian) (Gamma w)))))


;; #+RESULTS:
;; : #'fdg.prologue/Lagrange-equations

;; In the Lagrange equations procedure the parameter =Lagrangian= is a procedure
;; that implements the Lagrangian. The derivatives of the Lagrangian, for example
;; =((partial 2) Lagrangian)=, are also procedures. The state-space path procedure
;; =(Gamma w)= is constructed from the configuration-space path procedure =w= by
;; the procedure =Gamma=:


(defn Gamma [w]
  (fn [t]
    (up t (w t) ((D w) t))))


;; #+RESULTS:
;; : #'fdg.prologue/Gamma

;; where =up= is a constructor for a data structure that represents a state of the
;; dynamical system (time, coordinates, velocities).

;; The result of applying the =Lagrange-equations= procedure to a procedure
;; =Lagrangian= that implements a Lagrangian function is a procedure that takes a
;; configuration-space path procedure =w= and returns a procedure that gives the
;; residual of the Lagrange equations for that path at a time.

;; For example, consider the harmonic oscillator, with Lagrangian

;; $$L(t, q, v) = \frac{1}{2}mv^2 - \frac{1}{2}kq^2,$$

;; for mass $m$ and spring constant $k$. this lagrangian is implemented by


(defn L-harmonic [m k]
  (fn [[_ q v]]
      (- (* 1/2 m (square v))
         (* 1/2 k (square q)))))


;; #+RESULTS:
;; : #'fdg.prologue/L-harmonic

;; We know that the motion of a harmonic oscillator is a sinusoid with a given
;; amplitude $a$, frequency $\omega$, and phase $\varphi$:

;; $$x(t) = a \cos(\omega t + \varphi).$$

;; Suppose we have forgotten how the constants in the solution relate to the
;; physical parameters of the oscillator. Let’s plug in the proposed solution and
;; look at the residual:


(defn proposed-solution [t]
  (* 'a (cos (+ (* 'omega t) 'phi))))


;; #+RESULTS:
;; : #'fdg.prologue/proposed-solution


(tex$$
 (((Lagrange-equations (L-harmonic 'm 'k))
   proposed-solution)
  't))


;; #+RESULTS:
;; :results:
;; $$- a\,m\,{\omega}^{2}\,\cos\left(\omega\,t + \phi\right) + a\,k\,\cos\left(\omega\,t + \phi\right)$$
;; :end:

;; The residual here shows that for nonzero amplitude, the only solutions allowed
;; are ones where $(k - m\omega^2) = 0$ or $\omega = \sqrt{k/m}$.

;; But, suppose we had no idea what the solution looks like. We could propose a
;; literal function for the path:

(tex$$
 (((Lagrange-equations (L-harmonic 'm 'k))
   (literal-function 'x))
  't))
