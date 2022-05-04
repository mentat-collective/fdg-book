;; Chapter 1: Introduction
;; :PROPERTIES:
;; :header-args+: :tangle ../src/fdg/ch1.clj :comments org
;; :END:


(ns fdg.ch1
  (:refer-clojure :exclude [+ - * / = compare zero? ref partial
                            numerator denominator])
  (:require [sicmutils.env :as e :refer :all :exclude [F->C]]))

(define-coordinates t e/R1-rect)

;; The Metric

;; Let's now take a step into the geometry. A surface has a metric which tells us
;; how to measure sizes and angles at every point on the surface. (Metrics are
;; introduced in Chapter 9.)

;; The metric is a symmetric function of two vector fields that gives a number for
;; every point on the manifold. (Vector fields are introduced in Chapter 3).
;; Metrics may be used to compute the length of a vector field at each point, or
;; alternatively to compute the inner product of two vector fields at each point.
;; For example, the metric for the sphere of radius $R$ is

;; \begin{equation}
;; \mathsf{g}(\mathsf{u}, \mathsf{v})=R^{2} \mathsf{d} \theta(\mathsf{u})
;; \mathsf{d} \theta(\mathsf{v})+R^{2}(\sin \theta)^{2} \mathsf{d}
;; \phi(\mathsf{u}) \mathsf{d} \phi(\mathsf{v}),
;; \end{equation}

;; where $\mathsf{u}$ and $\mathsf{v}$ are vector fields, and $\mathsf{d}\theta$
;; and $\mathsf{d}\phi$ are one-form fields that extract the named components of
;; the vector-field argument. (One-form fields are introduced in Chapter 3.) We can
;; think of $\mathsf{d}\theta(\mathsf{u})$ as a function of a point that gives the
;; size of the vector field $\mathsf{u}$ in the $\theta$ direction at the point.
;; Notice that $\mathsf{g}(\mathsf{u}, \mathsf{u})$ is a weighted sum of the
;; squares of the components of $\mathsf{u}$. In fact, if we identify

;; \begin{align*}
;; &\mathsf{d} \theta(\mathsf{v})=\dot{\theta} \\
;; &\mathsf{d} \phi(\mathsf{v})=\dot{\phi},
;; \end{align*}

;; then the coefficients in the metric are the same as the coefficients in the
;; value of the Lagrangian, equation (1.1), apart from a factor of $m/2$.

;; We can generalize this result and write a Lagrangian for free motion of a
;; particle of mass $m$ on a manifold with metric $\mathsf{g}$:

;; \begin{equation}
;; L_{2}(x, v)=\sum_{i j} \frac{1}{2} m g_{i j}(x) v^{i} v^{j}
;; \end{equation}

;; This is written using indexed variables to indicate components of the geometric
;; objects expressed with respect to an unspecified coordinate system. The metric
;; coefficients $g_{ij}$ are, in general, a function of the position coordinates
;; $x$, because the properties of the space may vary from place to place.

;; We can capture this geometric statement as a program:


(defn L2 [mass metric]
  (fn [place velocity]
    (* 1/2 mass ((metric velocity velocity) place))))



;; #+RESULTS:
;; : #'chapter1/L2

;; This program gives the Lagrangian in a coordinate-independent, geometric way. It
;; is entirely in terms of geometric objects, such as a place on the configuration
;; manifold, the velocity at that place, and the metric that describes the local
;; shape of the manifold. But to compute we need a coordinate system. We express
;; the dynamical state in terms of coordinates and velocity components in the
;; coordinate system. For each coordinate system there is a natural vector basis
;; and the geometric velocity vectors can be constructed by contracting the basis
;; with the components of the velocity. Thus, we can form a coordinate
;; representation of the Lagrangian.


(defn Lc [mass metric coordsys]
  (let [e (coordinate-system->vector-basis coordsys)]
    (fn [[_ x v]]
      ((L2 mass metric) ((point coordsys) x) (* e v)))))



;; #+RESULTS:
;; : #'chapter1/Lc

;; The manifold point $\mathsf{m}$ represented by the coordinates $x$ is given by
;; =(define m ((point coordsys) x))=. The coordinates of $\mathsf{m}$ in a
;; different coordinate system are given by =((chart coordsys2) m)=. The manifold
;; point $\mathsf{m}$ is a geometric object that is the same point independent of
;; how it is specified. Similarly, the velocity vector $\mathsf{e}v$ is a geometric
;; object, even though it is specified using components $v$ with respect to the
;; basis $\mathsf{e}$. Both $v$ and $\mathsf{e}$ have as many components as the
;; dimension of the space so their product is interpreted as a contraction.

;; Let's make a general metric on a 2-dimensional real manifold:[fn:4]


(def the-metric (literal-metric 'g R2-rect))



;; #+RESULTS:
;; : #'chapter1/the-metric

;; The metric is expressed in rectangular coordinates, so the coordinate system is
;; =R2-rect=.[fn:5] The component functions will be labeled as subscripted ~g~s.

;; We can now make the Lagrangian for the system:


(def L (Lc 'm the-metric R2-rect))



;; #+RESULTS:
;; : #'chapter1/L

;; And we can apply our Lagrangian to an arbitrary state:


(simplify
 (L (up 't (up 'x 'y) (up 'vx 'vy))))

;; Euler-Lagrange Residuals

;; The Euler-Lagrange equations are satisfied on realizable paths. Let $\gamma$ be
;; a path on the manifold of configurations. (A path is a map from the
;; 1-dimensional real line to the configuration manifold. We introduce maps between
;; manifolds in Chapter 6.) Consider an arbitrary path:[fn:6]


(def gamma (literal-manifold-map 'q R1-rect R2-rect))



;; #+RESULTS:
;; : #'chapter1/gamma

;; The values of $\gamma$ are points on the manifold, not a coordinate
;; representation of the points. We may evaluate =gamma= only on points of the
;; real-line manifold; =gamma= produces points on the $\mathbb{R}^2$ manifold. So
;; to go from the literal real-number coordinate ='t= to a point on the real line
;; we use =((point R1-rect) 't)= and to go from a point =m= in $\mathbb{R}^2$ to
;; its coordinate representation we use =((chart R2-rect) m)=. (The procedures
;; point and chart are introduced in Chapter 2.) Thus


((chart R2-rect) (gamma ((point R1-rect) 't)))



;; #+RESULTS[e358f13540c3718dda00f1f02ee0bcfe3dcbeafd]:
;; : (up (q↑0 t) (q↑1 t))


(def coordinate-path
  (compose (chart R2-rect) gamma (point R1-rect)))



;; #+RESULTS:
;; : #'chapter1/coordinate-path


(coordinate-path 't)



;; #+RESULTS[dfcb38a21d28a56f5ab7d03d663824ba30315508]:
;; : (up (q↑0 t) (q↑1 t))

;; Now we can compute the residuals of the Euler-Lagrange equations, but we get a
;; large messy expression that we will not show.[fn:7] However, we will save it to
;; compare with the residuals of the geodesic equations.


(def Lagrange-residuals
  (((Lagrange-equations L) coordinate-path) 't))

;; Geodesic Equations

;; Now we get deeper into the geometry. The traditional way to write the geodesic
;; equations is
;; \begin{equation}
;; \nabla_{\mathsf{v}} \mathsf{v}=0
;; \end{equation}
;; where $\nabla$ is a covariant derivative operator. Roughly, $\nabla_{\mathsf{v}}
;; \mathsf{w}$ is a directional derivative. It gives a measure of the variation of
;; the vector field $\mathsf{w}$ as you walk along the manifold in the direction of
;; $\mathsf{v}$. (We will explain this in depth in Chapter 7.) $\nabla_{\mathsf{v}}
;; \mathsf{v}=0$ is intended to convey that the velocity vector is
;; parallel-transported by itself. When you walked East on the Equator you had to
;; hold the stick so that it was parallel to the Equator. But the stick is
;; constrained to the surface of the Earth, so moving it along the Equator required
;; turning it in three dimensions. The $\nabla$ thus must incorporate the
;; 3-dimensional shape of the Earth to provide a notion of "parallel" appropriate
;; for the denizens of the surface of the Earth. This information will appear as
;; the "Christoffel coefficients" in the coordinate representation of the geodesic
;; equations.

;; The trouble with the traditional way to write the geodesic equations (1.4) is
;; that the arguments to the covariant derivative are vector fields and the
;; velocity along the path is not a vector field. A more precise way of stating
;; this relation is:
;; \begin{equation}
;; \nabla^\gamma_{\partial/\partial\mathsf{t}} d\gamma\left(\partial/\partial \mathsf{t}\right) = 0.
;; \end{equation}
;; (We know that this may be unfamiliar notation, but we will explain it in
;; Chapter 7.)

;; In coordinates, the geodesic equations are expressed
;; \begin{equation}
;; D^{2} q^{i}(t)+\sum_{j k} \Gamma_{j k}^{i}(\gamma(t)) D q^{j}(t) D q^{k}(t)=0,
;; \end{equation}
;; where $q(t)$ is the coordinate path corresponding to the manifold path $\gamma$,
;; and $\Gamma^i_{jk}\left(\mathsf{m}\right)$ are Christoffel coefficients. The
;; $\Gamma^i_{jk}\left(\mathsf{m}\right)$ describe the "shape" of the manifold
;; close to the manifold point $\mathsf{m}$. They can be derived from the metric
;; $g$.

;; We can get and save the geodesic equation residuals by:

;; #+name: Cartan

(def Cartan
  (Christoffel->Cartan
   (metric->Christoffel-2
    the-metric
    (coordinate-system->basis R2-rect))))



;; #+RESULTS: Cartan
;; : #'chapter1/Cartan


(def geodesic-equation-residuals
  (((((covariant-derivative Cartan gamma) d:dt)
     ((differential gamma) d:dt))
    (chart R2-rect))
   ((point R1-rect) 't)))



;; #+RESULTS:
;; : #'chapter1/geodesic-equation-residuals

;; where =d/dt= is a vector field on the real line[fn:8] and =Cartan= is a way of
;; encapsulating the geometry, as specified by the Christoffel coefficients. The
;; Christoffel coefficients are computed from the metric:


(def Cartan
  (Christoffel->Cartan
   (metric->Christoffel-2
    the-metric
    (coordinate-system->basis R2-rect))))



;; #+RESULTS:
;; : #'chapter1/Cartan

;; The two messy residual results that we did not show are related by the metric.
;; If we change the representation of the geodesic equations by "lowering" them
;; using the mass and the metric, we see that the residuals are equal:


(def metric-components
  (metric->components
   the-metric
   (coordinate-system->basis R2-rect)))



;; #+RESULTS:
;; : #'chapter1/metric-components


(simplify
 (- Lagrange-residuals
    (* (* 'm (metric-components (gamma ((point R1-rect) 't))))
       geodesic-equation-residuals)))
