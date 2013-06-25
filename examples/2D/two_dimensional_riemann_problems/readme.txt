/**
@page TWOD-RP 2D Example: Two dimensional Riemann problems

How to simulate 2D Riemann Problems

@section TWOD-RiemannP Two dimensional Riemann problems

The examples reported here concern with 2D Riemann problems. The initial conditions are the same as the ones reported by Kurganov and Tadmor [\ref kurganov "1"].
A 2D quadrant [0,1]x[0,1] is subdivided into 4 sub-quadrants:

@code
    ^ Y
    |
  1 -----------------------------------
    |                |                |
    |                |                |
    |                |                |
    |                |                |
    |      (2)       |      (1)       |
    |                |                |
    |                |                |
    |                |                |
0.5 -----------------------------------
    |                |                |
    |                |                |
    |                |                |
    |                |                |
    |      (3)       |      (4)       |
    |                |                |
    |                |                |
    |                |                |  X
    ------------------------------------->
    0                0.5              1
@endcode

The initial conditions are imposed as following:

@code
             |(p,r,u,v)(1) if x > 0.5 and y > 0.5
             |(p,r,u,v)(2) if x < 0.5 and y > 0.5
(p,r,u,v)(0)=|(p,r,u,v)(3) if x < 0.5 and y < 0.5
             |(p,r,u,v)(4) if x > 0.5 and y < 0.5
@endcode

@section CONDITIONS Initial condition of the problems considered

The Kurganov and Tadmor configurations considered are:

@subsection kt-c03 Configuration 3:
@code
         (p,r,u,v)(1) = (1.500,1.5000,0.000,0.000)
         (p,r,u,v)(2) = (0.300,0.5323,1.206,0.000)
         (p,r,u,v)(3) = (0.029,0.1380,1.206,1.206)
         (p,r,u,v)(4) = (0.300,0.5323,0.000,1.206)
@endcode
An example of solution is reported in figure \ref fig-kt-c03 "1".
\anchor fig-kt-c03 \image html kt-c03-scaled-r.png "Figure 1: density profiles for t=0.3"

@subsection kt-c04 Configuration 4:
@code
         (p,r,u,v)(1) = (1.100,1.1000,0.0000,0.0000)
         (p,r,u,v)(2) = (0.350,0.5065,0.8939,0.0000)
         (p,r,u,v)(3) = (1.100,1.1000,0.8939,0.8939)
         (p,r,u,v)(4) = (0.350,0.5065,0.0000,0.8939)
@endcode
An example of solution is reported in figure \ref fig-kt-c04 "2".
\anchor fig-kt-c04 \image html kt-c04-scaled-r.png "Figure 2: density profiles for t=0.3"

@subsection kt-c05 Configuration 5:
@code
         (p,r,u,v)(1) = (1.0,1.0,-0.75,-0.5)
         (p,r,u,v)(2) = (1.0,2.0,-0.75, 0.5)
         (p,r,u,v)(3) = (1.0,1.0, 0.75, 0.5)
         (p,r,u,v)(4) = (1.0,3.0, 0.75,-0.5)
@endcode
An example of solution is reported in figure \ref fig-kt-c05 "3".
\anchor fig-kt-c05 \image html kt-c05-scaled-r.png "Figure 3: density profiles for t=0.3"

@subsection kt-c06 Configuration 6:
@code
         (p,r,u,v)(1) = (1.0,1.0, 0.75,-0.5)
         (p,r,u,v)(2) = (1.0,2.0, 0.75, 0.5)
         (p,r,u,v)(3) = (1.0,1.0,-0.75, 0.5)
         (p,r,u,v)(4) = (1.0,3.0,-0.75,-0.5)
@endcode
An example of solution is reported in figure \ref fig-kt-c06 "4".
\anchor fig-kt-c06 \image html kt-c06-scaled-r.png "Figure 4: density profiles for t=0.3"

@subsection kt-c12 Configuration 12:
@code
         (p,r,u,v)(1) = (0.4,0.5313,0.0000,0.0000)
         (p,r,u,v)(2) = (1.0,1.0000,0.7276,0.0000)
         (p,r,u,v)(3) = (1.0,0.8000,0.0000,0.0000)
         (p,r,u,v)(4) = (1.0,1.0000,0.0000,0.7276)
@endcode
An example of solution is reported in figure \ref fig-kt-c12 "5".
\anchor fig-kt-c12 \image html kt-c12-scaled-r.png "Figure 5: density profiles for t=0.3"

@subsection kt-c17 Configuration 17:
@code
         (p,r,u,v)(1) = (1.0,1.0000,0.0,-0.4000)
         (p,r,u,v)(2) = (1.0,2.0000,0.0,-0.3000)
         (p,r,u,v)(3) = (0.4,1.0625,0.0, 0.2145)
         (p,r,u,v)(4) = (0.4,0.5197,0.0,-1.1259)
@endcode
An example of solution is reported in figure \ref fig-kt-c17 "6".
\anchor fig-kt-c17 \image html kt-c17-scaled-r.png "Figure 6: density profiles for t=0.3"

@section rEFERENCES References

\anchor kurganov [1] Kurganov, A. and Tadmor, E., <em>Solution of Two Dimensional Riemann Problems for Gas Dynamics without Riemann Problem Solvers</em>, Numerical Methods for Partial Differential Equations, vol. 18, n.5, 584--608, 2002.

*/
