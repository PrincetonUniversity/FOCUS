
!title (bnormal) ! Calculate total bnormal and its derivatives.

!latex \briefly{Calculate the total bnormal of all coils on plasma surface and the derivatives with respect to coil geometry and currents, including the first and second dirivatives. 
!latex          Calling \emph{bnormal(0), bnormal(1), bnormal(2)} calculates the $0-order$, $1^{st}-order$ and $2^{nd}-order$ derivatives respectively.}

!latex \calledby{\link{costfun}}
!latex \calls{\link{bfield}}

!latex \tableofcontents

!latex \newcommand{\cmt}{\cos(mt)}
!latex \newcommand{\smt}{\sin(mt)}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!latex \subsection{Bnormal on plasma surface (0-order)}
!latex \bi
!latex \item[1.] The ``energy'' of normal magnetic field is defined as the ``quadratic-flux'' on a given ``plasma boundary'',
!latex           \be \label{eq:bnormal}
!latex               bnorm & = &  \int_{\cal S} \frac{1}{2}( B_n)^2 \; ds \nonumber \\
!latex                     & = &  \Delta \t \, \Delta \z \sum_{j,k} \sqrt g_{j,k}
!latex                            \frac{1}{2} \left( B_{n,j,k} \right)^2 
!latex           \ee
!latex           where \be B_{n,j,k} \equiv n^x_{j,k} \sum\nolimits_i B^{x,i}_{j,k} 
!latex                                    + n^y_{j,k} \sum\nolimits_i B^{y,i}_{j,k} 
!latex                                    + n^z_{j,k} \sum\nolimits_i B^{z,i}_{j,k}, \ee
!latex           where $B^{x,i}_{j,k}$, $B^{y,i}_{j,k}$ and $B^{z,i}_{j,k}$ are the Cartesian components of the magnetic field,
!latex           which depend explicity on the geometry of and the current in the $i$-th coil.
!latex           The normal vector to the plasma boundary at the angular location $(\t_{j,k},\z_{j,k})$ 
!latex           is \mbox{${\bf n}_{j,k} \equiv n^x_{j,k} \, {\bf i} + n^y_{j,k} \, {\bf j} + n^z_{j,k} \, {\bf k}$}.
!latex           (This is pre-computed in \link{surface}.)
!latex \item[2.] The resolution of the discretized surface integral is given by \inputvar{Nteta} and \inputvar{Nzeta} (see \link{global} and \link{surface}).
!latex \item[3.] The magnetic field at ${\bf x} \equiv x {\bf i} + y {\bf j} + z {\bf k}$ is given by the Biot-Savart integral,
!latex           \be \label{eq:BiotSavart}
!latex              {\bf B} \equiv I \int_{\cal C} \frac{d{\bf l}\times {\bf r}}{r^3}
!latex           \ee
!latex           where ${\bf r}\equiv {\bf x}-\bar{\bf x}$, where $\bar {\bf x}$ is a point on the plasma boundary
!latex \item[4.] In component form, \Eqn{BiotSavart} is
!latex           \be B^x & \equiv & I \int_0^{2\pi} \frac{\dot {\bar y} (z-\bar z) - \dot {\bar z} (y-\bar y) }{r^3} dt, \\
!latex               B^y & \equiv & I \int_0^{2\pi} \frac{\dot {\bar z} (x-\bar x) - \dot {\bar x} (z-\bar z) }{r^3} dt, \\
!latex               B^z & \equiv & I \int_0^{2\pi} \frac{\dot {\bar x} (y-\bar y) - \dot {\bar y} (x-\bar x) }{r^3} dt,
!latex           \ee
!latex           where $\dot {\bar x} \equiv d\bar x/dt$, etc. and the Fourier representation of the curves is given in \link{iccoil}.
!latex \ei
!latex \subsection{First derivatives}
!latex \bi
!latex  \item[1.] The derivatives of $bnorm$ with respect to the coil parameters, e.g. $\alpha^i$, take the form:
!latex           \be \frac{\partial bnorm}{\partial \alpha^i} & = & \Delta \t \, \Delta \z \sum_{j,k} \sqrt g_{j,k}
!latex                        \left( 
!latex                        n^x_{j,k} \sum\nolimits_i B^{x,i}_{j,k} + n^y_{j,k} \sum\nolimits_i B^{y,i}_{j,k} + n^z_{j,k} \sum\nolimits_i B^{z,i}_{j,k}
!latex                                          \right)
!latex                        \left(
!latex                        n^x_{j,k} \frac{\partial B^{x,i}_{j,k}}{\partial \alpha^i} 
!latex                      + n^y_{j,k} \frac{\partial B^{y,i}_{j,k}}{\partial \alpha^i} 
!latex                      + n^z_{j,k} \frac{\partial B^{z,i}_{j,k}}{\partial \alpha^i} 
!latex                        \right)
!latex           \ee
!latex \item[2.] The integrals over $t$ are provided by \nag{}{D01EAF}. Only the integrands are required:
!latex \be \begin{array}{ccccccccccccccccccccccccccccccccccccccccccc}
!latex                       B^x                    & \equiv & \ds \frac{ \Delta y \; \dot z - \Delta z \; \dot y}{r^3}                               \\ \\
!latex    \ds \frac{\partial B^x}{\partial x_{c,m}} & \equiv &                                                          & - B^x 3 \Delta x \cmt / r^2 \\ \\
!latex    \ds \frac{\partial B^x}{\partial y_{c,m}} & \equiv & \ds \frac{  \cmt    \; \dot z + \Delta z \; \smt m}{r^3} & - B^x 3 \Delta y \cmt / r^2 \\ \\
!latex    \ds \frac{\partial B^x}{\partial z_{c,m}} & \equiv & \ds \frac{-\Delta y \; \smt m - \cmt     \; \dot y}{r^3} & - B^x 3 \Delta z \cmt / r^2 \\ \\
!latex    \ds \frac{\partial B^x}{\partial x_{s,m}} & \equiv &                                                          & - B^x 3 \Delta x \smt / r^2 \\ \\
!latex    \ds \frac{\partial B^x}{\partial y_{s,m}} & \equiv & \ds \frac{  \smt    \; \dot z - \Delta z \; \cmt m}{r^3} & - B^x 3 \Delta y \smt / r^2 \\ \\
!latex    \ds \frac{\partial B^x}{\partial z_{s,m}} & \equiv & \ds \frac{ \Delta y \; \cmt m - \smt     \; \dot y}{r^3} & - B^x 3 \Delta z \smt / r^2

!latex \end{array} \ee
!latex \ei

!latex \subsection{The second derivatives}
!latex \bi
!latex \item[1.] The integrands in $B^x$, $B^y$ and $B^z$ can be expressed in a concise way:
!latex \be B^i \equiv g^i r^{-3} \ee
!latex where,
!latex \be g^i \equiv \varepsilon^{ijk} \; \Delta l^j \; \dot l^k \ee
!latex Here $\varepsilon^{ijk}$ is Levi-Civita symbol and for simplification I will omit this symbol in later $i,j,k$ cases.
!latex \item[2.] Therefore, the derivatives of $B^x$, $B^y$ and $B^z$ integrands (both the first and second derivatives) can be expressed as:
!latex \be 
!latex \begin{array}{ccc}
!latex \ds \frac{\partial B^i}{\partial x_l} & \equiv & \ds \frac{\partial g^i}{\partial x_l} r^{-3} - \frac{3}{r^{4}}\frac{\partial r}{\partial x_l}g^i\\ \\
!latex \ds \frac{\partial^2 B^i}{\partial x_l \partial x_m} & \equiv & \ds   \frac{\partial^2 g^i}{\partial x_l \partial x_m} r^{-3} - \frac{3}{r^{4}}\frac{\partial r}{\partial x_m}\frac{\partial g^i}{\partial x_l} 
!latex                                                          \ds + \frac{12}{r^{5}}\frac{\partial r}{\partial x_l}\frac{\partial r}{\partial x_m}g^i - \frac{3}{r^{-4}}\frac{\partial^2 r}{\partial x_l \partial x_m}g^i
!latex                                                          \ds - \frac{3}{r^{4}}\frac{\partial r}{\partial x_l}\frac{\partial g^i}{\partial x_m}
!latex \end{array}
!latex \ee
!latex \item[3.] In that case, $\frac{\partial^2 B^i}{\partial x_l \partial x_m}$ is just related to the derivatives of $g^i$ and $r$.So we can also write out all the derivatives of $g^i$ and $r$.
!latex \be
!latex \begin{array}{ccc}
!latex \ds \frac{\partial g^i}{\partial x_l} & \equiv & \ds \frac{\partial \Delta l^j}{\partial x_l} \dot l^k + \Delta l^j \frac{\partial \dot l^k}{\partial x_l}
!latex                                          - \frac{\partial \Delta l^k}{\partial x_l} \dot l^j - \Delta l^k \frac{\partial \dot l^j}{\partial x_l}\\ \\
!latex
!latex \ds \frac{\partial^2 g^i}{\partial x_l \partial x_m} & \equiv & \ds \frac{\partial^2 \Delta l^j}{\partial x_l \partial x_m} \dot l^k + \frac{\partial \Delta l^j}{\partial x_l} \frac{\partial \dot l^k}{\partial x_m}
!latex                                                         + \frac{\partial \Delta l^j}{\partial x_m} \frac{\partial \dot l^k}{\partial x_l} + \Delta l^j \frac{\partial^2 \dot l^k}{\partial x_l \partial x_m}\\ \\
!latex                                                  &&   \ds -    \frac{\partial^2 \Delta l^k}{\partial x_l \partial x_m} \dot l^j - \frac{\partial \Delta l^k}{\partial x_l} \frac{\partial \dot l^j}{\partial x_m}
!latex                                                         - \frac{\partial \Delta l^k}{\partial x_m} \frac{\partial \dot l^j}{\partial x_l} - \Delta l^k \frac{\partial^2 \dot l^j}{\partial x_l \partial x_m}\\ \\
!latex                \ds \frac{\partial r}{\partial x^i_l} & \equiv & \ds \frac{\Delta l^i}{r} \frac{\partial \Delta l^i}{\partial x_l}\\ \\
!latex \ds \frac{\partial^2 r}{\partial x^i_l \partial x^j_m} & \equiv & \ds \delta^i_j \frac{1}{r} \frac{\partial \Delta l^i}{\partial x_l}  \frac{\partial \Delta l^j}{\partial x_m}
!latex                                                           - \frac{\Delta l^i}{r^2} \frac{\partial \Delta r}{\partial x_m} \frac{1}{r} \frac{\partial \Delta l^i}{\partial x_l}\\ \\
!latex \end{array}
!latex \ee
!latex \item[4.] The derivatives of $bnorm$ in \Eqn{bnormal} can be written as,
!latex \be
!latex \ds \frac{\partial{bnorm}}{\partial{I^i}}   & \equiv & \int_{\cal S} \sum_j I^j(B_x^j n_x + B_y^j n_y + B_z^j n_z) \, (B_x^i n_x + B_y^i n_y + B_z^i n_z) \, ds\\
!latex \ds \frac{\partial{bnorm}}{\partial{x^i_n}} & \equiv & \int_{\cal S} \sum_j I^j(B_x^j n_x + B_y^j n_y + B_z^j n_z) \, (\frac{\partial{B_x^i}}{\partial{x_n^i}} n_x +
!latex                                                                  \frac{\partial{B_y^i}}{\partial{x_n^i}} n_y + \frac{\partial{B_z^i}}{\partial{x_n^i}} n_z) I^i \, ds\\

!latex \ds \frac{\partial^2{bnorm}}{\partial{I^i}\partial{I^j}} & \equiv & \int_{\cal S} (B_x^j n_x + B_y^j n_y + B_z^j n_z) (B_x^i n_x + B_y^i n_y + B_z^i n_z) \, ds\\
!latex \ds \frac{\partial^2{bnorm}}{\partial{I^i}\partial{x^j_n}} & \equiv & \int_{\cal S} (B_x^i n_x + B_y^i n_y + B_z^i n_z) (\frac{\partial{B_x^j}}{\partial{x_n^j}} n_x +
!latex                                                                  \frac{\partial{B_y^j}}{\partial{x_n^j}} n_y + \frac{\partial{B_z^j}}{\partial{x_n^j}} n_z)I^i \, ds\\
!latex \ds \frac{\partial^2{bnorm}}{\partial{x^i_m}\partial{x^j_n}} & \equiv & \int_{\cal S} (\frac{\partial{B_x^j}}{\partial{x_n^j}} n_x + \frac{\partial{B_y^j}}{\partial{x_n^j}} n_y 
!latex     + \frac{\partial{B_z^j}}{\partial{x_n^j}} n_z)(\frac{\partial{B_x^i}}{\partial{x_n^i}} n_x + \frac{\partial{B_y^i}}{\partial{x_n^i}} n_y + \frac{\partial{B_z^i}}{\partial{x_n^i}} n_z)I^i I^j \, ds\\
!latex \ee
!latex \ei

!latex \subsection{Normalization}
!latex \bi
!latex \item[1.] It's recommended to normalize all the cost functions, even the weights may need to be normalized. While dealing with bnormal function, it's normalized as,
!latex \be
!latex  Bnorm \equiv \int_s \frac{1}{2} \frac{(\vec{B} \cdot \vec{n})^2}{|B|^2} ds
!latex \ee
!latex \item[2.] For simplification, we can denote $Bn$ for $\vec{B} \cdot \vec{n}$ and $Bm$ for $|B|^2$. And their derivatives can be written as,
!latex \be \begin{array}{cccccccc}
!latex \ds \frac{\partial{Bn}}{\partial{I^i  }} & = &  B_x^i n_x \,& +  & \, B_y^i n_y \, & +  & \, B_z^i n_z \\
!latex \ds \frac{\partial{Bn}}{\partial{x_m^i}} & = &    ( \frac{\partial{B_x^i}}{\partial{x_m^i}} n_x \, & + & \, \frac{\partial{B_y^i}}{\partial{x_m^i}} n_y \, & + & \, \frac{\partial{B_z^i}}{\partial{x_m^i}} n_z ) \, I^i \\
!latex \ds \frac{\partial{Bm}}{\partial{I^i  }} & = &  2 ( B_x^i B_x \, & + & \, B_y^i B_y \, & + & \, B_z^i B_z )  \\
!latex \ds \frac{\partial{Bm}}{\partial{x_m^i}} & = &  2 ( \frac{\partial{B_x^i}}{\partial{x_m^i}} B_x \, & + & \, \frac{\partial{B_y^i}}{\partial{x_m^i}} B_y \, & + & \, \frac{\partial{B_z^i}}{\partial{x_m^i}} B_z ) \, I^i \\
!latex \end{array} \ee
!latex Here, the superscript $i$ is denoted as the $i^th$ coil's current ($I$) or geometric variables ($x$). And $B_x$ (or $B_y$ and $B_z$) means the total magnetic field at the surface point, while $B_x^i$ means 
!latex the magnetic field generated by the $i^th$ coil without timing current and Biot-Savart constant $\frac{\mu}{4\pi}$.\\
!latex \item[3.] Similarily, we can also write down the second derivatives for $Bn$ and $Bm$.
!latex \be \begin{array}{ccccccc}
!latex \ds \frac{\partial^2{Bn}}{\partial{I^i  }\partial{I^j  }} & = & 0 & \\
!latex \ds \frac{\partial^2{Bn}}{\partial{I^i  }\partial{x_n^j}} & = & \left \{ \begin{array}{ll} 0 & \mathrm{, \; if \; i \neq j;}\\                                                               
!latex      \frac{\partial{B_x^i}}{\partial{x_m^i}} n_x \, + \, \frac{\partial{B_y^i}}{\partial{x_m^i}} n_y \, + \, \frac{\partial{B_z^i}}{\partial{x_m^i}} n_z & \mathrm{, \; if \; i = j} \end{array} \right . \\
!latex \ds \frac{\partial^2{Bn}}{\partial{x_m^i}\partial{x_n^j}} & = & \left \{ \begin{array}{ll} 0 & \mathrm{, \; if \; i \neq j;}\\
!latex       ( \frac{\partial^2{B_x^i}}{\partial{x_m^i}\partial{x_n^i}} n_x + \frac{\partial^2{B_y^i}}{\partial{x_m^i}\partial{x_n^i}} n_y + \frac{\partial^2{B_z^i}}{\partial{x_m^i}\partial{x_n^i}} n_z  )
!latex      & \mathrm{, \; if \; i = j} \end{array} \right .
!latex \end{array} \ee
!latex \\
!latex \be \begin{array}{cccccccccccccccccccccccc}
!latex \ds \frac{\partial^2{Bm}}{\partial{I^i  }\partial{I^j  }} & = & 2 ( B_x^i B_x^j \,+ \, B_y^i B_y^j \,+ \, B_z^i B_z^j ) \\
!latex \ds \frac{\partial^2{Bm}}{\partial{I^i  }\partial{x_m^j}} & = & \left \{ 
!latex      \begin{array}{lllll} 2 ( B_x^i \frac{\partial{B_x^j}}{\partial{x_m^j}} \,+ \, B_y^i \frac{\partial{B_y^j}}{\partial{x_m^j}}\,+ \, B_z^i\frac{\partial{B_z^j}}{\partial{x_m^j}} ) I^j 
!latex                                                                                                                                                &  \mathrm{, \; if \; i \neq j;}\\                  
!latex      2 ( B_x^i \frac{\partial{B_x^j}}{\partial{x_m^j}} \,+ \, B_y^i \frac{\partial{B_y^j}}{\partial{x_m^j}}\,+ \, B_z^i\frac{\partial{B_z^j}}{\partial{x_m^j}} ) I^j \, + \,
!latex      2 ( \frac{\partial{B_x^i}}{\partial{x_m^i}} B_x \,+ \, \frac{\partial{B_y^i}}{\partial{x_m^i}} B_y \,+ \, \frac{\partial{B_z^i}}{\partial{x_m^i}} B_z )
!latex                                                                                                                                                & \mathrm{, \; if \; i = j} \end{array} \right . \\
!latex \ds \frac{\partial^2{Bm}}{\partial{x_m^i}\partial{x_n^j}} & = & \left \{ 
!latex      \begin{array}{ll} 2 ( \frac{\partial{B_x^i}}{\partial{x_n^i}}  \frac{\partial{B_x^j}}{\partial{x_m^j}} \,+ \, \frac{\partial{B_y^i}}{\partial{x_n^i}} \frac{\partial{B_y^j}}{\partial{x_m^j}} \,
!latex                       + \, \frac{\partial{B_z^i}}{\partial{x_n^i}} \frac{\partial{B_z^j}}{\partial{x_m^j}}) I^i I^j & \mathrm{, \; if \; i \neq j;}\\
!latex                         2 ( \frac{\partial{B_x^i}}{\partial{x_n^i}}  \frac{\partial{B_x^j}}{\partial{x_m^j}} \,+ \, \frac{\partial{B_y^i}}{\partial{x_n^i}} \frac{\partial{B_y^j}}{\partial{x_m^j}} \,
!latex                        + \, \frac{\partial{B_z^i}}{\partial{x_n^i}} \frac{\partial{B_z^j}}{\partial{x_m^j}}) I^i I^j \, + \, 
!latex                        2( \frac{\partial^2{B_x^i}}{\partial{x_m^i}\partial{x_n^i}} B_x + \frac{\partial^2{B_y^i}}{\partial{x_m^i}\partial{x_n^i}} B_y + \frac{\partial^2{B_z^i}}{\partial{x_m^i}\partial{x_n^i}} B_z  ) I^i
!latex      & \mathrm{, \; if \; i = j} \end{array} \right .
!latex \end{array} \ee
!latex \\
!latex \item[4.] So the derivatives of bnormal cost function can be represented as,
!latex \be \begin{array}{cccccccccccccccccccccccc}
!latex \ds \frac{\partial{bnorm}}{\partial{x}} & = & \frac{Bn}{Bm} \frac{\partial{Bn}}{\partial{x}} & - & \frac{1}{2} \frac{Bn^2}{Bm^2} \frac{\partial{Bm}}{\partial{x}} \\
!latex \ds \frac{\partial^2{bnorm}}{\partial{x_1}\partial{x_2}} & = & \frac{1}{Bm} \frac{\partial{Bn}}{\partial{x_1}} \frac{\partial{Bn}}{\partial{x_2}} & - & 
!latex          \frac{Bn}{Bm^2} \frac{\partial{Bn}}{\partial{x_1}} \frac{\partial{Bm}}{\partial{x_2}} & + & \frac{Bn}{Bm} \frac{\partial^2{Bn}}{\partial{x_1}\partial{x_2}} \\
!latex       && \frac{Bn^2}{Bm^3} \frac{\partial{Bm}}{\partial{x_1}} \frac{\partial{Bm}}{\partial{x_2}} & - & \frac{Bn}{Bm^2} \frac{\partial{Bm}}{\partial{x_1}} \frac{\partial{Bn}}{\partial{x_2}}
!latex          & - & \frac{1}{2} \frac{Bn^2}{Bm^2} \frac{\partial^2{Bm}}{\partial{x_1}\partial{x_2}} \\
!latex \end{array} \ee
!latex Here,$x, x_1, x_2$ can be both currents and geometric variables.
!latex \item[5.] \textbf{Note: } \emph{This normalization method works for 0-order cost function. But derivatives on current are actually divided by corresponding currents additionally.}
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$subroutine bnormal(nderiv)                              ! The subroutine takes a long time; need to be optimized later; 06/09/2016 !bnormal2 is faster without debug flag; 06/22/2016
!!$  use kmodule, only: zero, half, one, pi2, sqrtmachprec,  &
!!$       coil, surf, NFcoil, cmt, smt, NDcoil, Ncoils, Nteta, Nzeta, discretefactor, &
!!$       bnorm, t1B, t2B, bn, bm, &
!!$       ncpu, myid, ounit
!!$  implicit none
!!$  include "mpif.h"
!!$
!!$  INTEGER           :: nderiv
!!$
!!$
!!$  INTEGER           :: astat, ierr
!!$  INTEGER           :: icoil, mm, NN, Cdof, iteta, jzeta, kseg, ii, jj, kk, ll, c1, c2, n1, n2   ! icoil, iteta, jzeta, kseg are local
!!$  REAL              :: r, rm2, rm3, rm4, bx, by, bz, lm, c12, lbnorm
!!$  REAL, allocatable :: dlx(:), dly(:), dlz(:), ltx(:), lty(:), ltz(:), rp( :, :), l1B( :, :), l2B( :, :, :, :), b1n( :, :), b2n( :, :, :, :), b1m( :, :), b2m( :, :, :, :)
!!$  INTEGER           :: array2size, array4size
!!$  REAL              :: bsconstant = 1.0E-7
!!$
!!$
!!$
!!$  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$  NN   = NFcoil
!!$  Cdof  = 6 * NN + 6
!!$  bnorm = zero
!!$
!!$  SALLOCATE( dlx, (0:Cdof),               zero )
!!$  SALLOCATE( dly, (0:Cdof),               zero )
!!$  SALLOCATE( dlz, (0:Cdof),               zero )
!!$  SALLOCATE( ltx, (0:Cdof),               zero )
!!$  SALLOCATE( lty, (0:Cdof),               zero )
!!$  SALLOCATE( ltz, (0:Cdof),               zero )
!!$  if ( .not. allocated(bn) ) allocate( bn(0:Nteta-1, 0:Nzeta-1) )
!!$  if ( .not. allocated(bm) ) allocate( bm(0:Nteta-1, 0:Nzeta-1) )
!!$  bn = zero; bm = zero
!!$
!!$  lbnorm = 0
!!$
!!$  select case ( nderiv )
!!$
!!$   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$  case(0)
!!$
!!$   do jzeta = 0, Nzeta - 1
!!$    do iteta = 0, Nteta - 1
!!$     if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
!!$     bx = zero; by = zero; bz = zero
!!$     do icoil = 1, Ncoils
!!$
!!$      if(.not. allocated(coil(icoil)%Bx)) then 
!!$       SALLOCATE( coil(icoil)%Bx, (0:Cdof, 0:Cdof), zero )
!!$       SALLOCATE( coil(icoil)%By, (0:Cdof, 0:Cdof), zero )
!!$       SALLOCATE( coil(icoil)%Bz, (0:Cdof, 0:Cdof), zero )
!!$      else
!!$       coil(icoil)%Bx = zero
!!$       coil(icoil)%By = zero
!!$       coil(icoil)%Bz = zero
!!$      endif
!!$
!!$      do kseg = 0, NDcoil - 1
!!$       dlx = zero; ltx = zero;
!!$       dly = zero; lty = zero;          
!!$       dlz = zero; ltz = zero;          
!!$
!!$       dlx(0) = coil(icoil)%xx(kseg) - surf(1)%xx(iteta,jzeta)
!!$       dly(0) = coil(icoil)%yy(kseg) - surf(1)%yy(iteta,jzeta)
!!$       dlz(0) = coil(icoil)%zz(kseg) - surf(1)%zz(iteta,jzeta)
!!$       r = sqrt(dlx(0)**2 + dly(0)**2 + dlz(0)**2)
!!$       rm3 = r**(-3)
!!$
!!$       ltx(0) = coil(icoil)%xt(kseg)
!!$       lty(0) = coil(icoil)%yt(kseg)
!!$       ltz(0) = coil(icoil)%zt(kseg)
!!$
!!$       coil(icoil)%Bx(0,0) = coil(icoil)%Bx(0,0) + ( dly(0)*ltz(0) - dlz(0)*lty(0) ) * rm3
!!$       coil(icoil)%By(0,0) = coil(icoil)%By(0,0) + ( dlz(0)*ltx(0) - dlx(0)*ltz(0) ) * rm3
!!$       coil(icoil)%Bz(0,0) = coil(icoil)%Bz(0,0) + ( dlx(0)*lty(0) - dly(0)*ltx(0) ) * rm3
!!$
!!$      enddo    ! enddo kseg
!!$!      coil(icoil)%Bx = coil(icoil)%Bx * pi2 / NDcoil * bsconstant
!!$!      coil(icoil)%By = coil(icoil)%By * pi2 / NDcoil * bsconstant
!!$!      coil(icoil)%Bz = coil(icoil)%Bz * pi2 / NDcoil * bsconstant
!!$
!!$      bx = bx + coil(icoil)%Bx( 0, 0) * coil(icoil)%I * bsconstant
!!$      by = by + coil(icoil)%By( 0, 0) * coil(icoil)%I * bsconstant
!!$      bz = bz + coil(icoil)%Bz( 0, 0) * coil(icoil)%I * bsconstant
!!$
!!$     enddo ! end do icoil
!!$     bn(iteta, jzeta) = bx * surf(1)%nx(iteta,jzeta) + by * surf(1)%ny(iteta,jzeta) + bz * surf(1)%nz(iteta,jzeta)
!!$     bm(iteta, jzeta) = bx * bx + by * by + bz * bz ! magnitude square of B
!!$     lbnorm = lbnorm +  bn(iteta, jzeta)**2 * surf(1)%ds(iteta,jzeta) / bm(iteta, jzeta)
!!$
!!$    enddo ! end do iteta
!!$   enddo ! end do jzeta
!!$
!!$   call MPI_REDUCE( lbnorm, bnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
!!$   RlBCAST( bnorm, 1, 0 )
!!$   bnorm  = bnorm * half * discretefactor
!!$
!!$
!!$   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$  case(1)
!!$
!!$   SALLOCATE( rp , (  0:Cdof, 0:Cdof), zero )
!!$   SALLOCATE( l1B, (1:Ncoils, 0:Cdof), zero )
!!$   SALLOCATE( b1n, (1:Ncoils, 0:Cdof), zero )
!!$   SALLOCATE( b1m, (1:Ncoils, 0:Cdof), zero )
!!$
!!$
!!$   do jzeta = 0, Nzeta - 1
!!$    do iteta = 0, Nteta - 1
!!$     if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
!!$     bx = zero; by = zero; bz = zero
!!$     do icoil = 1, Ncoils
!!$
!!$      if(.not. allocated(coil(icoil)%Bx)) then 
!!$       SALLOCATE( coil(icoil)%Bx, (0:Cdof, 0:Cdof), zero )
!!$       SALLOCATE( coil(icoil)%By, (0:Cdof, 0:Cdof), zero )
!!$       SALLOCATE( coil(icoil)%Bz, (0:Cdof, 0:Cdof), zero )
!!$      else
!!$       coil(icoil)%Bx = zero
!!$       coil(icoil)%By = zero
!!$       coil(icoil)%Bz = zero
!!$      endif
!!$
!!$      do kseg = 0, NDcoil - 1
!!$       dlx = zero; ltx = zero; rp = zero
!!$       dly = zero; lty = zero;          
!!$       dlz = zero; ltz = zero;          
!!$
!!$       dlx(0) = coil(icoil)%xx(kseg) - surf(1)%xx(iteta,jzeta)
!!$       dly(0) = coil(icoil)%yy(kseg) - surf(1)%yy(iteta,jzeta)
!!$       dlz(0) = coil(icoil)%zz(kseg) - surf(1)%zz(iteta,jzeta)
!!$       r = sqrt(dlx(0)**2 + dly(0)**2 + dlz(0)**2)
!!$       rm3 = r**(-3)
!!$       rm4 = r**(-4)
!!$       ltx(0) = coil(icoil)%xt(kseg)
!!$       lty(0) = coil(icoil)%yt(kseg)
!!$       ltz(0) = coil(icoil)%zt(kseg)
!!$
!!$       ! The first derivatives of delta x and x'(t) on x,y,z
!!$       do ll = 0, NN
!!$        dlx( ll        + 1 ) = cmt(kseg,ll)
!!$        dlx( ll +   NN + 2 ) = smt(kseg,ll)
!!$        dly( ll + 2*NN + 3 ) = cmt(kseg,ll)
!!$        dly( ll + 3*NN + 4 ) = smt(kseg,ll)
!!$        dlz( ll + 4*NN + 5 ) = cmt(kseg,ll)
!!$        dlz( ll + 5*NN + 6 ) = smt(kseg,ll)
!!$
!!$        ltx(ll        + 1) = -ll*smt(kseg, ll)
!!$        ltx(ll +   NN + 2) =  ll*cmt(kseg, ll)
!!$        lty(ll + 2*NN + 3) = -ll*smt(kseg, ll)
!!$        lty(ll + 3*NN + 4) =  ll*cmt(kseg, ll)
!!$        ltz(ll + 4*NN + 5) = -ll*smt(kseg, ll)
!!$        ltz(ll + 5*NN + 6) =  ll*cmt(kseg, ll)
!!$       enddo
!!$
!!$       ! The first derivatives of r 
!!$       do ll = 1, 2*NN+2
!!$        rp(ll       , 0) = dlx(0) * dlx(ll       ) / r
!!$        rp(ll+2*NN+2, 0) = dly(0) * dly(ll+2*NN+2) / r
!!$        rp(ll+4*NN+4, 0) = dlz(0) * dlz(ll+4*NN+4) / r
!!$       enddo
!!$
!!$
!!$       coil(icoil)%Bx(0,0) = coil(icoil)%Bx(0,0) + ( dly(0)*ltz(0) - dlz(0)*lty(0) ) * rm3
!!$       coil(icoil)%By(0,0) = coil(icoil)%By(0,0) + ( dlz(0)*ltx(0) - dlx(0)*ltz(0) ) * rm3
!!$       coil(icoil)%Bz(0,0) = coil(icoil)%Bz(0,0) + ( dlx(0)*lty(0) - dly(0)*ltx(0) ) * rm3
!!$
!!$       do ll = 1, Cdof
!!$        coil(icoil)%Bx(ll, 0) = coil(icoil)%Bx(ll, 0) + ( dly(ll)*ltz( 0) + dly( 0)*ltz(ll) - dlz(ll)*lty( 0) - dlz( 0)*lty(ll) ) * rm3 &
!!$             - 3.0 * rp(ll,0) * (  dly( 0)*ltz( 0) - dlz( 0)*lty( 0) ) * rm4
!!$        coil(icoil)%By(ll, 0) = coil(icoil)%By(ll, 0) + ( dlz(ll)*ltx( 0) + dlz( 0)*ltx(ll) - dlx(ll)*ltz( 0) - dlx( 0)*ltz(ll) ) * rm3 &
!!$             - 3.0 * rp(ll,0) * (  dlz( 0)*ltx( 0) - dlx( 0)*ltz( 0) ) * rm4
!!$        coil(icoil)%Bz(ll, 0) = coil(icoil)%Bz(ll, 0) + ( dlx(ll)*lty( 0) + dlx( 0)*lty(ll) - dly(ll)*ltx( 0) - dly( 0)*ltx(ll) ) * rm3 &
!!$             - 3.0 * rp(ll,0) * (  dlx( 0)*lty( 0) - dly( 0)*ltx( 0) ) * rm4    
!!$       enddo
!!$
!!$      enddo    ! enddo kseg
!!$      coil(icoil)%Bx = coil(icoil)%Bx * pi2 / NDcoil * bsconstant
!!$      coil(icoil)%By = coil(icoil)%By * pi2 / NDcoil * bsconstant
!!$      coil(icoil)%Bz = coil(icoil)%Bz * pi2 / NDcoil * bsconstant
!!$
!!$      bx = bx + coil(icoil)%Bx( 0, 0) * coil(icoil)%I
!!$      by = by + coil(icoil)%By( 0, 0) * coil(icoil)%I
!!$      bz = bz + coil(icoil)%Bz( 0, 0) * coil(icoil)%I
!!$
!!$     enddo ! end do icoil
!!$     bn(iteta, jzeta) = bx * surf(1)%nx(iteta,jzeta) + by * surf(1)%ny(iteta,jzeta) + bz * surf(1)%nz(iteta,jzeta)
!!$     bm(iteta, jzeta) = bx * bx + by * by + bz * bz ! magnitude square of B
!!$
!!$     b1n = zero; b1m = zero
!!$
!!$     do c1 = 1, Ncoils
!!$      b1n(c1, 0       )  =       coil(c1)%Bx( 0, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By( 0, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz( 0, 0)*surf(1)%nz(iteta,jzeta)
!!$      b1m(c1, 0       )  = 2 * ( coil(c1)%Bx( 0, 0)*bx                      + coil(c1)%By( 0, 0)*by                      + coil(c1)%Bz( 0, 0)*bz                      )
!!$      do n1 = 1, Cdof
!!$       b1n(c1,n1      )  =     ( coil(c1)%Bx(n1, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By(n1, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz(n1, 0)*surf(1)%nz(iteta,jzeta) ) * coil(c1)%I
!!$       b1m(c1,n1      )  = 2 * ( coil(c1)%Bx(n1, 0)*bx                      + coil(c1)%By(n1, 0)*by                      + coil(c1)%Bz(n1, 0)*bz                      ) * coil(c1)%I
!!$      enddo
!!$     enddo
!!$
!!$     ! derivatives for bnomal cost function
!!$     lbnorm = lbnorm +  bn(iteta, jzeta)**2 * surf(1)%ds(iteta,jzeta) / bm(iteta, jzeta)
!!$
!!$     do n1 = 0, Cdof
!!$      do c1 = 1, Ncoils
!!$       l1B(c1,n1) = l1B(c1,n1) + ( bn(iteta,jzeta)/bm(iteta,jzeta)*b1n(c1,n1) - half*bn(iteta,jzeta)**2/bm(iteta,jzeta)**2*b1m(c1,n1) ) * surf(1)%ds(iteta,jzeta)
!!$      enddo
!!$     enddo
!!$
! comment out; symmetric stuff; 06/23/2016
!!$     lbnorm = lbnorm +  bn(iteta, jzeta)**2 * surf(1)%ds(iteta,jzeta)
!!$
!!$     ! first derivatives of energy on DOFs
!!$     do c1 = 1, Ncoils
!!$      l1B(c1, 0) = l1B(c1, 0) + bn(iteta, jzeta) * ( coil(c1)%Bx(  0, 0) * surf(1)%nx(iteta,jzeta) &
!!$           + coil(c1)%By(  0, 0) * surf(1)%ny(iteta,jzeta) &
!!$           + coil(c1)%Bz(  0, 0) * surf(1)%nz(iteta,jzeta) ) * surf(1)%ds(iteta,jzeta)              ! first  derivatives on currents
!!$     enddo
!!$
!!$     do n1 = 1, Cdof
!!$      do c1 = 1, Ncoils
!!$       l1B(c1,n1) = l1B(c1,n1) + bn(iteta, jzeta) * ( coil(c1)%Bx( n1, 0) * surf(1)%nx(iteta,jzeta) &
!!$            + coil(c1)%By( n1, 0) * surf(1)%ny(iteta,jzeta) &
!!$            + coil(c1)%Bz( n1, 0) * surf(1)%nz(iteta,jzeta) ) * coil(c1)%I * surf(1)%ds(iteta,jzeta) ! first  derivatives on geometry
!!$      enddo ! end do n1
!!$     enddo ! end do c1
!!$
!!$
!!$    enddo ! end do iteta
!!$   enddo ! end do jzeta
!!$
!!$   call MPI_REDUCE( lbnorm, bnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
!!$   RlBCAST( bnorm, 1, 0 )
!!$   bnorm  = bnorm * half * discretefactor
!!$
!!$   if ( .not. allocated(t1B) ) allocate(t1B(1:Ncoils, 0:Cdof), stat=astat)
!!$   t1B = zero
!!$   ! SALLOCATE( t1B, (1:Ncoils, 0:Cdof), zero )
!!$
!!$   array2size = Ncoils * ( Cdof + 1 )
!!$   call MPI_REDUCE( l1B, t1B, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
!!$   RlBCAST( t1B, array2size, 0 )
!!$   t1B = t1B * discretefactor
!!$
!!$   DALLOCATE(l1B)
!!$   DALLOCATE(b1n)
!!$   DALLOCATE(b1m)
!!$   DALLOCATE( rp)
!!$
!!$
!!$   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$  case(2)
!!$
!!$   SALLOCATE( rp , (0:Cdof  , 0:Cdof), zero )
!!$   SALLOCATE( l1B, (1:Ncoils, 0:Cdof), zero )
!!$   SALLOCATE( b1n, (1:Ncoils, 0:Cdof), zero )
!!$   SALLOCATE( b1m, (1:Ncoils, 0:Cdof), zero )
!!$   SALLOCATE( l2B, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero )
!!$   SALLOCATE( b2n, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero )
!!$   SALLOCATE( b2m, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero )
!!$
!!$   do jzeta = 0, Nzeta - 1
!!$    do iteta = 0, Nteta - 1
!!$     if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
!!$     bx = zero; by = zero; bz = zero
!!$     do icoil = 1, Ncoils
!!$
!!$      if(.not. allocated(coil(icoil)%Bx)) then 
!!$       SALLOCATE( coil(icoil)%Bx, (0:Cdof, 0:Cdof), zero )
!!$       SALLOCATE( coil(icoil)%By, (0:Cdof, 0:Cdof), zero )
!!$       SALLOCATE( coil(icoil)%Bz, (0:Cdof, 0:Cdof), zero )
!!$      else
!!$       coil(icoil)%Bx = zero
!!$       coil(icoil)%By = zero
!!$       coil(icoil)%Bz = zero
!!$      endif
!!$
!!$      do kseg = 0, NDcoil - 1
!!$       dlx = zero; ltx = zero; rp = zero
!!$       dly = zero; lty = zero;          
!!$       dlz = zero; ltz = zero;          
!!$
!!$       dlx(0) = coil(icoil)%xx(kseg) - surf(1)%xx(iteta,jzeta)
!!$       dly(0) = coil(icoil)%yy(kseg) - surf(1)%yy(iteta,jzeta)
!!$       dlz(0) = coil(icoil)%zz(kseg) - surf(1)%zz(iteta,jzeta)
!!$       r = sqrt(dlx(0)**2 + dly(0)**2 + dlz(0)**2)
!!$       rm2 = r**(-2)
!!$       rm3 = r**(-3)
!!$       rm4 = r**(-4)
!!$       ltx(0) = coil(icoil)%xt(kseg)
!!$       lty(0) = coil(icoil)%yt(kseg)
!!$       ltz(0) = coil(icoil)%zt(kseg)
!!$
!!$       ! The first derivatives of delta x and x'(t) on x,y,z
!!$       do ll = 0, NN
!!$        dlx( ll        + 1 ) = cmt(kseg,ll)
!!$        dlx( ll +   NN + 2 ) = smt(kseg,ll)
!!$        dly( ll + 2*NN + 3 ) = cmt(kseg,ll)
!!$        dly( ll + 3*NN + 4 ) = smt(kseg,ll)
!!$        dlz( ll + 4*NN + 5 ) = cmt(kseg,ll)
!!$        dlz( ll + 5*NN + 6 ) = smt(kseg,ll)
!!$
!!$        ltx(ll        + 1) = -ll*smt(kseg, ll)
!!$        ltx(ll +   NN + 2) =  ll*cmt(kseg, ll)
!!$        lty(ll + 2*NN + 3) = -ll*smt(kseg, ll)
!!$        lty(ll + 3*NN + 4) =  ll*cmt(kseg, ll)
!!$        ltz(ll + 4*NN + 5) = -ll*smt(kseg, ll)
!!$        ltz(ll + 5*NN + 6) =  ll*cmt(kseg, ll)
!!$       enddo
!!$
!!$       ! The first derivatives of r 
!!$       do ll = 1, 2*NN+2
!!$        rp(ll       , 0) = dlx(0) * dlx(ll       ) / r
!!$        rp(ll+2*NN+2, 0) = dly(0) * dly(ll+2*NN+2) / r
!!$        rp(ll+4*NN+4, 0) = dlz(0) * dlz(ll+4*NN+4) / r
!!$       enddo
!!$
!!$       ! The second derivatives of r   
!!$       do mm = 1, 2*NN+2
!!$        do ll = 1, 2*NN+2
!!$         rp(ll       ,mm       ) = dlx(ll       ) * dlx(mm       ) / r - dlx(0) * dlx(ll       ) * rp(mm,       0) * rm2
!!$         rp(ll       ,mm+2*NN+2) =                                     - dlx(0) * dlx(ll       ) * rp(mm+2*NN+2,0) * rm2
!!$         rp(ll       ,mm+4*NN+4) =                                     - dlx(0) * dlx(ll       ) * rp(mm+4*NN+4,0) * rm2
!!$
!!$         rp(ll+2*NN+2,mm       ) =                                     - dly(0) * dly(ll+2*NN+2) * rp(mm       ,0) * rm2
!!$         rp(ll+2*NN+2,mm+2*NN+2) = dly(ll+2*NN+2) * dly(mm+2*NN+2) / r - dly(0) * dly(ll+2*NN+2) * rp(mm+2*NN+2,0) * rm2
!!$         rp(ll+2*NN+2,mm+4*NN+4) =                                     - dly(0) * dly(ll+2*NN+2) * rp(mm+4*NN+4,0) * rm2
!!$
!!$         rp(ll+4*NN+4,mm       ) =                                     - dlz(0) * dlz(ll+4*NN+4) * rp(mm       ,0) * rm2
!!$         rp(ll+4*NN+4,mm+2*NN+2) =                                     - dlz(0) * dlz(ll+4*NN+4) * rp(mm+2*NN+2,0) * rm2
!!$         rp(ll+4*NN+4,mm+4*NN+4) = dlz(ll+4*NN+4) * dlz(mm+4*NN+4) / r - dlz(0) * dlz(ll+4*NN+4) * rp(mm+4*NN+4,0) * rm2
!!$        enddo
!!$       enddo
!!$
!!$       coil(icoil)%Bx(0,0) = coil(icoil)%Bx(0,0) + ( dly(0)*ltz(0) - dlz(0)*lty(0) ) * rm3
!!$       coil(icoil)%By(0,0) = coil(icoil)%By(0,0) + ( dlz(0)*ltx(0) - dlx(0)*ltz(0) ) * rm3
!!$       coil(icoil)%Bz(0,0) = coil(icoil)%Bz(0,0) + ( dlx(0)*lty(0) - dly(0)*ltx(0) ) * rm3
!!$
!!$       do ll = 1, Cdof
!!$        coil(icoil)%Bx(ll, 0) = coil(icoil)%Bx(ll, 0) + ( dly(ll)*ltz( 0) + dly( 0)*ltz(ll) - dlz(ll)*lty( 0) - dlz( 0)*lty(ll) ) * rm3 &
!!$                                     - 3.0 * rp(ll,0) * ( dly( 0)*ltz( 0)                   - dlz( 0)*lty( 0)                   ) * rm4
!!$        coil(icoil)%By(ll, 0) = coil(icoil)%By(ll, 0) + ( dlz(ll)*ltx( 0) + dlz( 0)*ltx(ll) - dlx(ll)*ltz( 0) - dlx( 0)*ltz(ll) ) * rm3 &
!!$                                     - 3.0 * rp(ll,0) * ( dlz( 0)*ltx( 0)                   - dlx( 0)*ltz( 0)                   ) * rm4
!!$        coil(icoil)%Bz(ll, 0) = coil(icoil)%Bz(ll, 0) + ( dlx(ll)*lty( 0) + dlx( 0)*lty(ll) - dly(ll)*ltx( 0) - dly( 0)*ltx(ll) ) * rm3 &
!!$                                     - 3.0 * rp(ll,0) * ( dlx( 0)*lty( 0)                   - dly( 0)*ltx( 0)                   ) * rm4
!!$       enddo
!!$
!!$       do mm = 1, Cdof
!!$        do ll = 1, mm
!!$         coil(icoil)%Bx(ll,mm) = coil(icoil)%Bx(ll,mm) + ( dly(ll)*ltz(mm) + dly(mm)*ltz(ll) - dlz(ll)*lty(mm) - dlz(mm)*lty(ll) ) * rm3 &
!!$                             - 3.0 * rm4 * ( rp(mm, 0) * ( dly(ll)*ltz( 0) + dly( 0)*ltz(ll) - dlz(ll)*lty( 0) - dlz( 0)*lty(ll) )       &
!!$                                           + rp(ll,mm) * ( dly( 0)*ltz( 0)                   - dlz( 0)*lty( 0)                   )       &
!!$                                           + rp(ll, 0) * ( dly(mm)*ltz( 0) + dly( 0)*ltz(mm) - dlz(mm)*lty( 0) - dlz( 0)*lty(mm) )       &
!!$                         - 4.0 * rp(ll, 0) * rp(mm, 0) * ( dly( 0)*ltz( 0)                   - dlz( 0)*lty( 0)                   ) / r   )
!!$
!!$         coil(icoil)%By(ll,mm) = coil(icoil)%By(ll,mm) + ( dlz(ll)*ltx(mm) + dlz(mm)*ltx(ll) - dlx(ll)*ltz(mm) - dlx(mm)*ltz(ll) ) * rm3 &
!!$                             - 3.0 * rm4 * ( rp(mm, 0) * ( dlz(ll)*ltx( 0) + dlz( 0)*ltx(ll) - dlx(ll)*ltz( 0) - dlx( 0)*ltz(ll) )       &
!!$                                           + rp(ll,mm) * ( dlz( 0)*ltx( 0)                   - dlx( 0)*ltz( 0)                   )       &
!!$                                           + rp(ll, 0) * ( dlz(mm)*ltx( 0) + dlz( 0)*ltx(mm) - dlx(mm)*ltz( 0) - dlx( 0)*ltz(mm) )       &
!!$                         - 4.0 * rp(ll, 0) * rp(mm, 0) * ( dlz( 0)*ltx( 0)                   - dlx( 0)*ltz( 0)                   ) / r   )
!!$
!!$         coil(icoil)%Bz(ll,mm) = coil(icoil)%Bz(ll,mm) + ( dlx(ll)*lty(mm) + dlx(mm)*lty(ll) - dly(ll)*ltx(mm) - dly(mm)*ltx(ll) ) * rm3 &
!!$                             - 3.0 * rm4 * ( rp(mm, 0) * ( dlx(ll)*lty( 0) + dlx( 0)*lty(ll) - dly(ll)*ltx( 0) - dly( 0)*ltx(ll) )       &
!!$                                           + rp(ll,mm) * ( dlx( 0)*lty( 0)                   - dly( 0)*ltx( 0)                   )       &
!!$                                           + rp(ll, 0) * ( dlx(mm)*lty( 0) + dlx( 0)*lty(mm) - dly(mm)*ltx( 0) - dly( 0)*ltx(mm) )       &
!!$                         - 4.0 * rp(ll, 0) * rp(mm, 0) * ( dlx( 0)*lty( 0)                   - dly( 0)*ltx( 0)                   ) / r   )
!!$        enddo
!!$       enddo
!!$
!!$
!!$
!!$       ! Derivatives for Bn(1:3, 0:Cdof, 0:Cdof) ; sum over segments 
!!$
!!$      enddo    ! enddo kseg
!!$
!!$      do mm = 1, Cdof-1
!!$       do ll = mm+1, Cdof
!!$        coil(icoil)%Bx(ll,mm) = coil(icoil)%Bx(mm,ll)
!!$        coil(icoil)%By(ll,mm) = coil(icoil)%By(mm,ll)
!!$        coil(icoil)%Bz(ll,mm) = coil(icoil)%Bz(mm,ll)
!!$       enddo
!!$      enddo
!!$
!!$      coil(icoil)%Bx = coil(icoil)%Bx * pi2 / NDcoil * bsconstant
!!$      coil(icoil)%By = coil(icoil)%By * pi2 / NDcoil * bsconstant
!!$      coil(icoil)%Bz = coil(icoil)%Bz * pi2 / NDcoil * bsconstant
!!$
!!$      bx = bx + coil(icoil)%Bx( 0, 0) * coil(icoil)%I
!!$      by = by + coil(icoil)%By( 0, 0) * coil(icoil)%I
!!$      bz = bz + coil(icoil)%Bz( 0, 0) * coil(icoil)%I
!!$
!!$     enddo ! end do icoil
!!$     bn(iteta, jzeta) = bx * surf(1)%nx(iteta,jzeta) + by * surf(1)%ny(iteta,jzeta) + bz * surf(1)%nz(iteta,jzeta)
!!$     bm(iteta, jzeta) = bx * bx + by * by + bz * bz ! magnitude square of B    
!!$
!!$     ! derivatives for bn & bm terms
!!$     b1n = zero; b1m = zero
!!$     b2n = zero; b2m = zero
!!$
!!$     ! b1n( 0, 0) = bn(iteta, jzeta)
!!$     ! b1m( 0, 0) = bm(iteta, jzeta)
!!$
!!$     do c1 = 1, Ncoils
!!$      b1n(c1, 0       )  =       coil(c1)%Bx( 0, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By( 0, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz( 0, 0)*surf(1)%nz(iteta,jzeta)
!!$      b1m(c1, 0       )  = 2 * ( coil(c1)%Bx( 0, 0)*bx                      + coil(c1)%By( 0, 0)*by                      + coil(c1)%Bz( 0, 0)*bz                      )
!!$      do n1 = 1, Cdof
!!$       b1n(c1,n1      )  =     ( coil(c1)%Bx(n1, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By(n1, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz(n1, 0)*surf(1)%nz(iteta,jzeta) ) * coil(c1)%I
!!$       b1m(c1,n1      )  = 2 * ( coil(c1)%Bx(n1, 0)*bx                      + coil(c1)%By(n1, 0)*by                      + coil(c1)%Bz(n1, 0)*bz                      ) * coil(c1)%I
!!$      enddo
!!$     enddo
!!$     
!!$     ! current & current
!!$     do c2 = 1, Ncoils
!!$      do c1 = 1, Ncoils
!!$       b2m(c1, 0,c2, 0 ) = 2 * ( coil(c1)%Bx( 0, 0)*coil(c2)%Bx( 0, 0)      + coil(c1)%By( 0, 0)*coil(c2)%By( 0, 0)      + coil(c1)%Bz( 0, 0)*coil(c2)%Bz( 0, 0)      )
!!$      enddo
!!$     enddo
!!$
!!$     ! current & geometry
!!$     do c2 = 1, Ncoils
!!$      do n2 = 1, Cdof
!!$      do c1 = 1, Ncoils
!!$
!!$        b2m(c1, 0,c2,n2) = 2 * ( coil(c1)%Bx( 0, 0)*coil(c2)%Bx(n2, 0)      + coil(c1)%By( 0, 0)*coil(c2)%By(n2, 0)      + coil(c1)%Bz( 0, 0)*coil(c2)%Bz(n2, 0)      ) * coil(c2)%I
!!$        b2m(c1,n2,c2, 0) = 2 * ( coil(c1)%Bx(n2, 0)*coil(c2)%Bx( 0, 0)      + coil(c1)%By(n2, 0)*coil(c2)%By( 0, 0)      + coil(c1)%Bz(n2, 0)*coil(c2)%Bz( 0, 0)      ) * coil(c1)%I
!!$
!!$       enddo
!!$
!!$       b2n(c2, 0,c2,n2)  =       coil(c2)%Bx(n2, 0)*surf(1)%nx(iteta,jzeta) + coil(c2)%By(n2, 0)*surf(1)%ny(iteta,jzeta) + coil(c2)%Bz(n2, 0)*surf(1)%nz(iteta,jzeta)
!!$       b2n(c2,n2,c2, 0)  =       b2n(c2, 0,c2,n2)
!!$       b2m(c2, 0,c2,n2)  = 2 * ( coil(c2)%Bx(n2, 0)*bx                      + coil(c2)%By(n2, 0)*by                      + coil(c2)%Bz(n2, 0)*bz                      ) + b2m(c2, 0,c2,n2)
!!$       b2m(c2,n2,c2, 0)  =       b2m(c2, 0,c2,n2)
!!$       
!!$      enddo
!!$     enddo
!!$     
!!$      
!!$     ! geometry & geometry 
!!$     do c2 = 1, Ncoils
!!$      do n2 = 1, Cdof
!!$       do n1 = 1, Cdof
!!$        do c1 = 1, Ncoils
!!$
!!$         b2m(c1,n1,c2,n2) = 2 * ( coil(c1)%Bx(n1, 0)*coil(c2)%Bx(n2, 0)      + coil(c1)%By(n1, 0)*coil(c2)%By(n2, 0)      + coil(c1)%Bz(n1, 0)*coil(c2)%Bz(n2, 0)      ) * coil(c1)%I * coil(c2)%I
!!$
!!$        enddo
!!$        
!!$        b2n(c2,n1,c2,n2)  =       coil(c2)%Bx(n1,n2)*surf(1)%nx(iteta,jzeta) + coil(c2)%By(n1,n2)*surf(1)%ny(iteta,jzeta) + coil(c2)%Bz(n1,n2)*surf(1)%nz(iteta,jzeta)
!!$        b2m(c2,n1,c2,n2)  = 2 * ( coil(c2)%Bx(n1,n2)*bx                      + coil(c2)%By(n1,n2)*by                      + coil(c2)%Bz(n1,n2)*bz                      ) * coil(c2)%I + b2m(c2,n1,c2,n2)
!!$
!!$       enddo
!!$      enddo
!!$     enddo
!!$     
!!$        
!!$        
!!$
!!$          
!!$
!!$
!!$
!!$
!!$!!!$     do c2 = 1, Ncoils
!!$!!!$      do c1 = 1, Ncoils!c1 <= c2; symmetric matrix
!!$!!!$
!!$!!!$       !b2n(c1, 0,c2, 0) = zero
!!$!!!$       b2m(c1, 0,c2, 0 ) = 2 * ( coil(c1)%Bx( 0, 0)*coil(c2)%Bx( 0, 0)      + coil(c1)%By( 0, 0)*coil(c2)%By( 0, 0)      + coil(c1)%Bz( 0, 0)*coil(c2)%Bz( 0, 0)      )
!!$!!!$
!!$!!!$       do n2 = 1, Cdof
!!$!!!$
!!$!!!$        b2m(c1, 0,c2,n2) = 2 * ( coil(c1)%Bx( 0, 0)*coil(c2)%Bx(n2, 0)      + coil(c1)%By( 0, 0)*coil(c2)%By(n2, 0)      + coil(c1)%Bz( 0, 0)*coil(c2)%Bz(n2, 0)      ) * coil(c2)%I
!!$!!!$       b2m(c1,n2,c2, 0) = 2 * ( coil(c1)%Bx(n2, 0)*coil(c2)%Bx( 0, 0)      + coil(c1)%By(n2, 0)*coil(c2)%By( 0, 0)      + coil(c1)%Bz(n2, 0)*coil(c2)%Bz( 0, 0)      ) * coil(c1)%I
!!$!!!$
!!$!!!$        do n1 = 1, Cdof
!!$!!!$
!!$!!!$         b2m(c1,n1,c2,n2)= 2 * ( coil(c1)%Bx(n1, 0)*coil(c2)%Bx(n2, 0)      + coil(c1)%By(n1, 0)*coil(c2)%By(n2, 0)      + coil(c1)%Bz(n1, 0)*coil(c2)%Bz(n2, 0)      ) * coil(c1)%I * coil(c2)%I
!!$!!!$
!!$!!!$        enddo ! enddo n1
!!$!!!$       enddo  !enddo n2
!!$!!!$      enddo !enddo c1
!!$!!!$
!!$!!!$      do n2 = 1, Cdof
!!$!!!$
!!$!!!$       b2n(c2, 0,c2,n2)  =       coil(c2)%Bx(n2, 0)*surf(1)%nx(iteta,jzeta) + coil(c2)%By(n2, 0)*surf(1)%ny(iteta,jzeta) + coil(c2)%Bz(n2, 0)*surf(1)%nz(iteta,jzeta)
!!$!!!$       b2m(c2, 0,c2,n2)  = 2 * ( coil(c2)%Bx(n2, 0)*bx                      + coil(c2)%By(n2, 0)*by                      + coil(c2)%Bz(n2, 0)*bz                      ) * coil(c2)%I + b2m(c2, 0,c2,n2)
!!$!!!$
!!$!!!$       do n1 = 1, Cdof
!!$!!!$
!!$!!!$        b2n(c2,n1,c2,n2)  =       coil(c2)%Bx(n1,n2)*surf(1)%nx(iteta,jzeta) + coil(c2)%By(n1,n2)*surf(1)%ny(iteta,jzeta) + coil(c2)%Bz(n1,n2)*surf(1)%nz(iteta,jzeta)
!!$!!!$        b2m(c2,n1,c2,n2)  = 2 * ( coil(c2)%Bx(n1,n2)*bx                      + coil(c2)%By(n1,n2)*by                      + coil(c2)%Bz(n1,n2)*bz                      ) * coil(c2)%I + b2m(c2,n1,c2,n2)
!!$!!!$
!!$!!!$       enddo !enddo n1
!!$!!!$      enddo  !enddo n2
!!$!!!$
!!$!!!$     enddo   !enddo c2
!!$
!!$     ! derivatives for bnomal cost function
!!$     lbnorm = lbnorm +  bn(iteta, jzeta)**2 * surf(1)%ds(iteta,jzeta) / bm(iteta, jzeta)
!!$
!!$     do n1 = 0, Cdof
!!$      do c1 = 1, Ncoils
!!$       l1B(c1,n1) = l1B(c1,n1) + ( bn(iteta,jzeta)/bm(iteta,jzeta)*b1n(c1,n1) - half*bn(iteta,jzeta)**2/bm(iteta,jzeta)**2*b1m(c1,n1) ) * surf(1)%ds(iteta,jzeta)
!!$      enddo
!!$     enddo
!!$
!!$     do n2 = 0, Cdof
!!$      do c2 = 1, Ncoils
!!$       do n1 = 0, Cdof
!!$        do c1 = 1, Ncoils    
!!$         l2B(c1,n1,c2,n2) = l2B(c1,n1,c2,n2) + ( one/bm(iteta,jzeta)*   b1n(c1,n1)*b1n(c2,n2) - bn(iteta,jzeta)/bm(iteta,jzeta)**2*b1n(c1,n1)*b1m(c2,n2) +      bn(iteta,jzeta)   /bm(iteta,jzeta)   *b2n(c1,n1,c2,n2) &
!!$                               +  bn(iteta,jzeta)**2/bm(iteta,jzeta)**3*b1m(c1,n1)*b1m(c2,n2) - bn(iteta,jzeta)/bm(iteta,jzeta)**2*b1m(c1,n1)*b1n(c2,n2) - half*bn(iteta,jzeta)**2/bm(iteta,jzeta)**2*b2m(c1,n1,c2,n2) ) &
!!$                               * surf(1)%ds(iteta,jzeta)
!!$        enddo
!!$       enddo
!!$      enddo
!!$     enddo
!!$
!!$
!!$
!!$    enddo ! end do iteta
!!$   enddo ! end do jzeta
!!$
!!$
!!$   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$   call MPI_REDUCE( lbnorm, bnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
!!$   RlBCAST( bnorm, 1, 0)
!!$   bnorm  = bnorm * half * discretefactor
!!$
!!$   if ( .not. allocated(t1B) ) allocate(t1B(1:Ncoils, 0:Cdof), stat=astat)
!!$   if ( .not. allocated(t2B) ) allocate(t2B(1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), stat=astat)
!!$   t1B = zero; t2B = zero
!!$   ! SALLOCATE( t1B, (1:Ncoils, 0:Cdof), zero )
!!$   ! SALLOCATE( t2B, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero )
!!$
!!$   array2size = Ncoils * ( Cdof + 1 )
!!$   call MPI_REDUCE( l1B, t1B, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
!!$   RlBCAST( t1B, array2size, 0 )
!!$
!!$!comment out; symmetric stuffs; 06/23/2016
!!$   do c2 = 1, Ncoils-1
!!$    do c1 = c2+1, Ncoils
!!$     l2B(c1,0,c2,0) = l2B(c2,0,c1,0)
!!$    enddo
!!$   enddo
!!$
!!$   do c2 = 1, Ncoils
!!$    do n1 = 1, Cdof
!!$     do c1 = 1, Ncoils
!!$      l2B(c1,n1,c2,0) = l2B(c2,0,c1,n1)
!!$     enddo
!!$    enddo
!!$   enddo
!!$
!!$   do n2 = 0, Cdof
!!$    do c2 = 1, Ncoils-1
!!$     do n1 = 0, Cdof
!!$      do c1 = c2+1, Ncoils
!!$       l2B(c1,n1,c2,n2) = l2B(c2,n2,c1,n1)
!!$      enddo
!!$     enddo ! end n1
!!$    enddo    ! end c2     
!!$   enddo       ! end n2
!!$
!!$   array4size = array2size * array2size
!!$   call MPI_REDUCE( l2B, t2B, array4size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
!!$   RlBCAST( t2B, array4size, 0 )
!!$
!!$
!!$
!!$   t1B = t1B * discretefactor
!!$   t2B = t2B * discretefactor
!!$
!!$   DALLOCATE( rp)
!!$   DALLOCATE(l1B)
!!$   DALLOCATE(b1n)
!!$   DALLOCATE(b1m)
!!$   DALLOCATE(l2B)
!!$   DALLOCATE(b2n)
!!$   DALLOCATE(b2m)
!!$
!!$  end select
!!$  DALLOCATE( dlx )
!!$  DALLOCATE( dly )
!!$  DALLOCATE( dlz )
!!$  DALLOCATE( ltx )
!!$  DALLOCATE( lty )
!!$  DALLOCATE( ltz )
!!$  return
!!$end subroutine bnormal









!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bnormal( nderiv )

  use kmodule, only: zero, half, one, pi2, sqrtmachprec, bsconstant, &
       coil, surf, NFcoil, cmt, smt, NDcoil, Ncoils, Nteta, Cdof, Nzeta, discretefactor, &
       bnorm, t1B, t2B, bn, bm, tbn, SaveBx, SaveBy, SaveBz, &
       ncpu, myid, ounit

  implicit none
  include "mpif.h"

  INTEGER           :: nderiv


  INTEGER           :: astat, ierr
  INTEGER           :: icoil, NN,  iteta, jzeta, c1, c2, n1, n2   ! icoil, iteta, jzeta, kseg are local
  REAL              :: r, rm2, rm3, rm4, bx, by, bz, lm, c12, lbnorm, start, finish
  REAL, allocatable :: l1B( :, :), l2B( :, :, :, :), b1n( :, :), b2n( :, :, :, :), b1m( :, :), b2m( :, :, :, :)
  INTEGER           :: array2size, array4size
  REAL              :: lBx(0:Nteta, 0:Nzeta), lBy(0:Nteta, 0:Nzeta), lBz(0:Nteta, 0:Nzeta)

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if ( .not. allocated(SaveBx) ) then
     SALLOCATE(SaveBx, (0:Nteta, 0:Nzeta), zero)
     SALLOCATE(SaveBy, (0:Nteta, 0:Nzeta), zero)
     SALLOCATE(SaveBz, (0:Nteta, 0:Nzeta), zero)
  endif

  lBx = zero; lBy = zero; lBz = zero

  !if ( .not. allocated(bn) ) allocate( bn(0:Nteta-1, 0:Nzeta-1) )
  if ( .not. allocated(bm) ) allocate( bm(0:Nteta-1, 0:Nzeta-1) )
  !if ( .not. allocated(tbn)) allocate(tbn(0:Nteta-1, 0:Nzeta-1) )
  if ( .not. allocated(bn) ) allocate( bn(0:Nteta, 0:Nzeta) )  ! for saving the Bn data; 2017/02/10
  if ( .not. allocated(tbn)) allocate(tbn(0:Nteta, 0:Nzeta) )

  if ( .not. allocated(surf(1)%bnt) ) then
     allocate( surf(1)%bnt(0:Nteta, 0:Nzeta) ) ! Bnorm target from plasma currents or other coils
     surf(1)%bnt = zero
  endif
  
  NN = NFcoil; bnorm = zero; bn = zero; bm = zero; lbnorm = zero; tbn = zero

  do icoil = 1, Ncoils
   if(.not. allocated(coil(icoil)%Bx)) then 
    SALLOCATE( coil(icoil)%Bx, (0:Cdof, 0:Cdof), zero )
    SALLOCATE( coil(icoil)%By, (0:Cdof, 0:Cdof), zero )
    SALLOCATE( coil(icoil)%Bz, (0:Cdof, 0:Cdof), zero )
   else
    coil(icoil)%Bx = zero
    coil(icoil)%By = zero
    coil(icoil)%Bz = zero
   endif
  enddo

  select case ( nderiv )

   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  case(0)

   do jzeta = 0, Nzeta - 1
    do iteta = 0, Nteta - 1

     if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
     bx = zero; by = zero; bz = zero

     do icoil = 1, Ncoils

      call bfield0(icoil, iteta, jzeta, coil(icoil)%Bx(0,0), coil(icoil)%By(0,0), coil(icoil)%Bz(0,0))

      bx = bx + coil(icoil)%Bx( 0, 0) * coil(icoil)%I * pi2/coil(icoil)%D * bsconstant    
      by = by + coil(icoil)%By( 0, 0) * coil(icoil)%I * pi2/coil(icoil)%D * bsconstant
      bz = bz + coil(icoil)%Bz( 0, 0) * coil(icoil)%I * pi2/coil(icoil)%D * bsconstant

     enddo ! end do icoil

     lBx(iteta, jzeta) = bx; lBy(iteta, jzeta) = by; lBz(iteta, jzeta) = bz;
     bn(iteta, jzeta) = bx * surf(1)%nx(iteta,jzeta) + by * surf(1)%ny(iteta,jzeta) + bz * surf(1)%nz(iteta,jzeta)
     bm(iteta, jzeta) = bx * bx + by * by + bz * bz ! magnitude square of B
     !bm(iteta, jzeta) = one  ! no normalization;
#ifdef BNORM
     lbnorm = lbnorm +  (bn(iteta, jzeta) - surf(1)%bnt(iteta, jzeta))**2 / bm(iteta, jzeta) * surf(1)%ds(iteta,jzeta)  !change to minus on 2017/02/07
#else
     lbnorm = lbnorm +  bn(iteta, jzeta)**2 / bm(iteta, jzeta) * surf(1)%ds(iteta,jzeta)
#endif

    enddo ! end do iteta
   enddo ! end do jzeta

   call MPI_REDUCE( lbnorm, bnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( bnorm, 1, 0 )
   bnorm  = bnorm * half * discretefactor

   call MPI_REDUCE( bn, tbn, (1+Nteta)*(1+Nzeta), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( tbn,(1+Nteta)*(1+Nzeta), 0)

   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  case(1)

   SALLOCATE( l1B, (1:Ncoils, 0:Cdof), zero )
   SALLOCATE( b1n, (1:Ncoils, 0:Cdof), zero )
   SALLOCATE( b1m, (1:Ncoils, 0:Cdof), zero )

   do jzeta = 0, Nzeta - 1
    do iteta = 0, Nteta - 1

     if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
     bx = zero; by = zero; bz = zero

     do icoil = 1, Ncoils

!!$      if (icoil .eq. 1 .and. iteta .eq. 0 .and. jzeta .eq. 0) then
!!$      call CPU_TIME( start  )
!!$      write(ounit,*) coil(icoil)%xx(0), coil(icoil)%yy(0), coil(icoil)%zz(0)
!!$      write(ounit,*) coil(icoil)%xt(0), coil(icoil)%yt(0), coil(icoil)%zt(0)
!!$      call bfield1(icoil, iteta, jzeta, coil(icoil)%Bx(0:Cdof,0), coil(icoil)%By(0:Cdof,0), coil(icoil)%Bz(0:Cdof,0))
!!$      write(ounit,*) "Bx/x = ", coil(icoil)%Bx(0:Cdof,0)
!!$      call CPU_TIME( finish )
!!$      write(ounit,*) "Time for one single call is ", finish-start
!!$      endif

      call bfield1(icoil, iteta, jzeta, coil(icoil)%Bx(0:Cdof,0), coil(icoil)%By(0:Cdof,0), coil(icoil)%Bz(0:Cdof,0))

      coil(icoil)%Bx = coil(icoil)%Bx * pi2 / coil(icoil)%D * bsconstant
      coil(icoil)%By = coil(icoil)%By * pi2 / coil(icoil)%D * bsconstant
      coil(icoil)%Bz = coil(icoil)%Bz * pi2 / coil(icoil)%D * bsconstant

      bx = bx + coil(icoil)%Bx( 0, 0) * coil(icoil)%I
      by = by + coil(icoil)%By( 0, 0) * coil(icoil)%I
      bz = bz + coil(icoil)%Bz( 0, 0) * coil(icoil)%I

     enddo ! end do icoil

     lBx(iteta, jzeta) = bx; lBy(iteta, jzeta) = by; lBz(iteta, jzeta) = bz; 
     bn(iteta, jzeta) = bx * surf(1)%nx(iteta,jzeta) + by * surf(1)%ny(iteta,jzeta) + bz * surf(1)%nz(iteta,jzeta)
     bm(iteta, jzeta) = bx * bx + by * by + bz * bz ! magnitude square of B

     b1n = zero; b1m = zero

     do c1 = 1, Ncoils
      b1n(c1, 0       )  =       coil(c1)%Bx( 0, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By( 0, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz( 0, 0)*surf(1)%nz(iteta,jzeta)
      b1m(c1, 0       )  = 2 * ( coil(c1)%Bx( 0, 0)*bx                      + coil(c1)%By( 0, 0)*by                      + coil(c1)%Bz( 0, 0)*bz                      )
      do n1 = 1, Cdof
       b1n(c1,n1      )  =     ( coil(c1)%Bx(n1, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By(n1, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz(n1, 0)*surf(1)%nz(iteta,jzeta) ) * coil(c1)%I
       b1m(c1,n1      )  = 2 * ( coil(c1)%Bx(n1, 0)*bx                      + coil(c1)%By(n1, 0)*by                      + coil(c1)%Bz(n1, 0)*bz                      ) * coil(c1)%I
      enddo
     enddo

     ! derivatives for bnomal cost function
#ifdef BNORM
     lbnorm = lbnorm +  (bn(iteta, jzeta) - surf(1)%bnt(iteta, jzeta))**2 / bm(iteta, jzeta) * surf(1)%ds(iteta,jzeta)
#else
     lbnorm = lbnorm +  bn(iteta, jzeta)**2 / bm(iteta, jzeta) * surf(1)%ds(iteta,jzeta)
#endif

     do n1 = 0, Cdof
      do c1 = 1, Ncoils
#ifdef BNORM
        l1B(c1,n1) = l1B(c1,n1) + ( (bn(iteta,jzeta)-surf(1)%bnt(iteta,jzeta))/bm(iteta,jzeta)*b1n(c1,n1) &
                                - half*(bn(iteta,jzeta)-surf(1)%bnt(iteta,jzeta))**2/bm(iteta,jzeta)**2*b1m(c1,n1) ) * surf(1)%ds(iteta,jzeta)
#else
        l1B(c1,n1) = l1B(c1,n1) + ( bn(iteta,jzeta)/bm(iteta,jzeta)*b1n(c1,n1) - half*bn(iteta,jzeta)**2/bm(iteta,jzeta)**2*b1m(c1,n1) ) * surf(1)%ds(iteta,jzeta)
#endif
      enddo
     enddo

!!$     ! first derivatives of energy on DOFs
!!$     do c1 = 1, Ncoils
!!$      l1B(c1, 0) = l1B(c1, 0) + bn(iteta, jzeta) * ( coil(c1)%Bx(  0, 0) * surf(1)%nx(iteta,jzeta) &
!!$           + coil(c1)%By(  0, 0) * surf(1)%ny(iteta,jzeta) &
!!$           + coil(c1)%Bz(  0, 0) * surf(1)%nz(iteta,jzeta) ) * surf(1)%ds(iteta,jzeta)              ! first  derivatives on currents
!!$     enddo
!!$
!!$     do n1 = 1, Cdof
!!$      do c1 = 1, Ncoils
!!$       l1B(c1,n1) = l1B(c1,n1) + bn(iteta, jzeta) * ( coil(c1)%Bx( n1, 0) * surf(1)%nx(iteta,jzeta) &
!!$            + coil(c1)%By( n1, 0) * surf(1)%ny(iteta,jzeta) &
!!$            + coil(c1)%Bz( n1, 0) * surf(1)%nz(iteta,jzeta) ) * coil(c1)%I * surf(1)%ds(iteta,jzeta) ! first  derivatives on geometry
!!$      enddo ! end do n1
!!$     enddo ! end do c1


    enddo ! end do iteta
   enddo ! end do jzeta

   call MPI_REDUCE( lbnorm, bnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( bnorm, 1, 0 )
   bnorm  = bnorm * half * discretefactor

   call MPI_REDUCE( bn, tbn, (1+Nteta)*(1+Nzeta), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( tbn,(1+Nteta)*(1+Nzeta), 0)

   if ( .not. allocated(t1B) ) allocate(t1B(1:Ncoils, 0:Cdof), stat=astat)
   t1B = zero

   array2size = Ncoils * ( Cdof + 1 )
   call MPI_REDUCE( l1B, t1B, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( t1B, array2size, 0 )

   t1B = t1B * discretefactor

   DALLOCATE(l1B)
   DALLOCATE(b1n)
   DALLOCATE(b1m)

   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case(2)

   SALLOCATE( l1B, (1:Ncoils, 0:Cdof), zero )
   SALLOCATE( b1n, (1:Ncoils, 0:Cdof), zero )
   SALLOCATE( b1m, (1:Ncoils, 0:Cdof), zero )
   SALLOCATE( l2B, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero )
   SALLOCATE( b2n, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero )
   SALLOCATE( b2m, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero )

   do jzeta = 0, Nzeta - 1
    do iteta = 0, Nteta - 1

     if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
     bx = zero; by = zero; bz = zero

     do icoil = 1, Ncoils

      call bfield2(icoil, iteta, jzeta, coil(icoil)%Bx(0:Cdof,0:Cdof), coil(icoil)%By(0:Cdof,0:Cdof), coil(icoil)%Bz(0:Cdof,0:Cdof))

      coil(icoil)%Bx = coil(icoil)%Bx * pi2/coil(icoil)%D * bsconstant
      coil(icoil)%By = coil(icoil)%By * pi2/coil(icoil)%D * bsconstant
      coil(icoil)%Bz = coil(icoil)%Bz * pi2/coil(icoil)%D * bsconstant

      bx = bx + coil(icoil)%Bx( 0, 0) * coil(icoil)%I
      by = by + coil(icoil)%By( 0, 0) * coil(icoil)%I
      bz = bz + coil(icoil)%Bz( 0, 0) * coil(icoil)%I

     enddo ! end do icoil

     lBx(iteta, jzeta) = bx; lBy(iteta, jzeta) = by; lBz(iteta, jzeta) = bz; 
     bn(iteta, jzeta) = bx * surf(1)%nx(iteta,jzeta) + by * surf(1)%ny(iteta,jzeta) + bz * surf(1)%nz(iteta,jzeta)
     bm(iteta, jzeta) = bx * bx + by * by + bz * bz ! magnitude square of B    

     ! derivatives for bn & bm terms
     b1n = zero; b1m = zero
     b2n = zero; b2m = zero

     ! b1n( 0, 0) = bn(iteta, jzeta)
     ! b1m( 0, 0) = bm(iteta, jzeta)

     do c1 = 1, Ncoils
      b1n(c1, 0       )  =       coil(c1)%Bx( 0, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By( 0, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz( 0, 0)*surf(1)%nz(iteta,jzeta)
      b1m(c1, 0       )  = 2 * ( coil(c1)%Bx( 0, 0)*bx                      + coil(c1)%By( 0, 0)*by                      + coil(c1)%Bz( 0, 0)*bz                      )
      do n1 = 1, Cdof
       b1n(c1,n1      )  =     ( coil(c1)%Bx(n1, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By(n1, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz(n1, 0)*surf(1)%nz(iteta,jzeta) ) * coil(c1)%I
       b1m(c1,n1      )  = 2 * ( coil(c1)%Bx(n1, 0)*bx                      + coil(c1)%By(n1, 0)*by                      + coil(c1)%Bz(n1, 0)*bz                      ) * coil(c1)%I
      enddo
     enddo

     
     ! current & current
     do c2 = 1, Ncoils
      do c1 = 1, Ncoils
       b2m(c1, 0,c2, 0 ) = 2 * ( coil(c1)%Bx( 0, 0)*coil(c2)%Bx( 0, 0)      + coil(c1)%By( 0, 0)*coil(c2)%By( 0, 0)      + coil(c1)%Bz( 0, 0)*coil(c2)%Bz( 0, 0)      )
      enddo
     enddo

     ! current & geometry
     do c2 = 1, Ncoils
      do n2 = 1, Cdof
      do c1 = 1, Ncoils

        b2m(c1, 0,c2,n2) = 2 * ( coil(c1)%Bx( 0, 0)*coil(c2)%Bx(n2, 0)      + coil(c1)%By( 0, 0)*coil(c2)%By(n2, 0)      + coil(c1)%Bz( 0, 0)*coil(c2)%Bz(n2, 0)      ) * coil(c2)%I
        b2m(c1,n2,c2, 0) = 2 * ( coil(c1)%Bx(n2, 0)*coil(c2)%Bx( 0, 0)      + coil(c1)%By(n2, 0)*coil(c2)%By( 0, 0)      + coil(c1)%Bz(n2, 0)*coil(c2)%Bz( 0, 0)      ) * coil(c1)%I

       enddo

       b2n(c2, 0,c2,n2)  =       coil(c2)%Bx(n2, 0)*surf(1)%nx(iteta,jzeta) + coil(c2)%By(n2, 0)*surf(1)%ny(iteta,jzeta) + coil(c2)%Bz(n2, 0)*surf(1)%nz(iteta,jzeta)
       b2n(c2,n2,c2, 0)  =       b2n(c2, 0,c2,n2)
       b2m(c2, 0,c2,n2)  = 2 * ( coil(c2)%Bx(n2, 0)*bx                      + coil(c2)%By(n2, 0)*by                      + coil(c2)%Bz(n2, 0)*bz                      ) + b2m(c2, 0,c2,n2)
       b2m(c2,n2,c2, 0)  =       b2m(c2, 0,c2,n2)
       
      enddo
     enddo
     
      
     ! geometry & geometry 
     do c2 = 1, Ncoils
      do n2 = 1, Cdof
       do n1 = 1, Cdof
        do c1 = 1, Ncoils

         b2m(c1,n1,c2,n2) = 2 * ( coil(c1)%Bx(n1, 0)*coil(c2)%Bx(n2, 0)      + coil(c1)%By(n1, 0)*coil(c2)%By(n2, 0)      + coil(c1)%Bz(n1, 0)*coil(c2)%Bz(n2, 0)      ) * coil(c1)%I * coil(c2)%I

        enddo
        
        b2n(c2,n1,c2,n2)  =     ( coil(c2)%Bx(n1,n2)*surf(1)%nx(iteta,jzeta) + coil(c2)%By(n1,n2)*surf(1)%ny(iteta,jzeta) + coil(c2)%Bz(n1,n2)*surf(1)%nz(iteta,jzeta) ) * coil(c2)%I
        b2m(c2,n1,c2,n2)  = 2 * ( coil(c2)%Bx(n1,n2)*bx                      + coil(c2)%By(n1,n2)*by                      + coil(c2)%Bz(n1,n2)*bz                      ) * coil(c2)%I + b2m(c2,n1,c2,n2)

       enddo
      enddo
     enddo

!!$     do c2 = 1, Ncoils
!!$      do c1 = 1, Ncoils!c1 <= c2; symmetric matrix
!!$
!!$       !b2n(c1, 0,c2, 0) = zero
!!$       b2m(c1, 0,c2, 0 ) = 2 * ( coil(c1)%Bx( 0, 0)*coil(c2)%Bx( 0, 0)      + coil(c1)%By( 0, 0)*coil(c2)%By( 0, 0)      + coil(c1)%Bz( 0, 0)*coil(c2)%Bz( 0, 0)      )
!!$
!!$       do n2 = 1, Cdof
!!$
!!$        b2m(c1, 0,c2,n2) = 2 * ( coil(c1)%Bx( 0, 0)*coil(c2)%Bx(n2, 0)      + coil(c1)%By( 0, 0)*coil(c2)%By(n2, 0)      + coil(c1)%Bz( 0, 0)*coil(c2)%Bz(n2, 0)      ) * coil(c2)%I
!!$!       b2m(c1,n2,c2, 0) = 2 * ( coil(c1)%Bx(n2, 0)*coil(c2)%Bx( 0, 0)      + coil(c1)%By(n2, 0)*coil(c2)%By( 0, 0)      + coil(c1)%Bz(n2, 0)*coil(c2)%Bz( 0, 0)      ) * coil(c1)%I
!!$
!!$        do n1 = 1, Cdof
!!$
!!$         b2m(c1,n1,c2,n2)= 2 * ( coil(c1)%Bx(n1, 0)*coil(c2)%Bx(n2, 0)      + coil(c1)%By(n1, 0)*coil(c2)%By(n2, 0)      + coil(c1)%Bz(n1, 0)*coil(c2)%Bz(n2, 0)      ) * coil(c1)%I * coil(c2)%I
!!$
!!$        enddo ! enddo n1
!!$       enddo  !enddo n2
!!$      enddo !enddo c1
!!$
!!$      do n2 = 1, Cdof
!!$
!!$       b2n(c2, 0,c2,n2)  =       coil(c2)%Bx(n2, 0)*surf(1)%nx(iteta,jzeta) + coil(c2)%By(n2, 0)*surf(1)%ny(iteta,jzeta) + coil(c2)%Bz(n2, 0)*surf(1)%nz(iteta,jzeta)
!!$       b2m(c2, 0,c2,n2)  = 2 * ( coil(c2)%Bx(n2, 0)*bx                      + coil(c2)%By(n2, 0)*by                      + coil(c2)%Bz(n2, 0)*bz                      ) * coil(c2)%I + b2m(c2, 0,c2,n2)
!!$
!!$       do n1 = 1, Cdof
!!$
!!$        b2n(c2,n1,c2,n2) =        coil(c2)%Bx(n1,n2)*surf(1)%nx(iteta,jzeta) + coil(c2)%By(n1,n2)*surf(1)%ny(iteta,jzeta) + coil(c2)%Bz(n1,n2)*surf(1)%nz(iteta,jzeta)
!!$        b2m(c2,n1,c2,n2)  = 2 * ( coil(c2)%Bx(n1,n2)*bx                      + coil(c2)%By(n1,n2)*by                      + coil(c2)%Bz(n1,n2)*bz                      ) * coil(c2)%I + b2m(c2,n1,c2,n2)
!!$
!!$       enddo !enddo n1
!!$      enddo  !enddo n2
!!$
!!$     enddo   !enddo c2

     ! derivatives for bnomal cost function
#ifdef BNORM
     lbnorm = lbnorm +  (bn(iteta, jzeta) - surf(1)%bnt(iteta, jzeta))**2 / bm(iteta, jzeta) * surf(1)%ds(iteta,jzeta)
#else
     lbnorm = lbnorm +  bn(iteta, jzeta)**2 / bm(iteta, jzeta) * surf(1)%ds(iteta,jzeta)
#endif

     do n1 = 0, Cdof
        do c1 = 1, Ncoils
#ifdef BNORM
            l1B(c1,n1) = l1B(c1,n1) + ( (bn(iteta,jzeta)-surf(1)%bnt(iteta,jzeta))/bm(iteta,jzeta)*b1n(c1,n1) &
                                    - half*(bn(iteta,jzeta)-surf(1)%bnt(iteta,jzeta))**2/bm(iteta,jzeta)**2*b1m(c1,n1) ) * surf(1)%ds(iteta,jzeta)
#else
            l1B(c1,n1) = l1B(c1,n1) + ( bn(iteta,jzeta)/bm(iteta,jzeta)*b1n(c1,n1) - half*bn(iteta,jzeta)**2/bm(iteta,jzeta)**2*b1m(c1,n1) ) * surf(1)%ds(iteta,jzeta)
#endif
      enddo
     enddo

     do n2 = 0, Cdof
      do c2 = 1, Ncoils
       do n1 = 0, Cdof
        do c1 = 1, Ncoils
#ifdef BNORM
         l2B(c1,n1,c2,n2) = l2B(c1,n1,c2,n2) + ( one/bm(iteta,jzeta)* b1n(c1,n1)*b1n(c2,n2) &
                                             - (bn(iteta,jzeta)-surf(1)%bnt(iteta,jzeta))/bm(iteta,jzeta)**2*b1n(c1,n1)*b1m(c2,n2) &
                                             + (bn(iteta,jzeta)-surf(1)%bnt(iteta,jzeta))/bm(iteta,jzeta)*b2n(c1,n1,c2,n2) &
                                             + (bn(iteta,jzeta)-surf(1)%bnt(iteta,jzeta))**2/bm(iteta,jzeta)**3*b1m(c1,n1)*b1m(c2,n2) &
                                             - (bn(iteta,jzeta)-surf(1)%bnt(iteta,jzeta))/bm(iteta,jzeta)**2*b1m(c1,n1)*b1n(c2,n2) &
                                             - half*(bn(iteta,jzeta)-surf(1)%bnt(iteta,jzeta))**2/bm(iteta,jzeta)**2*b2m(c1,n1,c2,n2) ) &
                                             * surf(1)%ds(iteta,jzeta)
#else
           l2B(c1,n1,c2,n2) = l2B(c1,n1,c2,n2) + ( one/bm(iteta,jzeta)*b1n(c1,n1)*b1n(c2,n2) &
                                               - bn(iteta,jzeta)/bm(iteta,jzeta)**2*b1n(c1,n1)*b1m(c2,n2) &
                                               + bn(iteta,jzeta)/bm(iteta,jzeta)   *b2n(c1,n1,c2,n2) &
                                               + bn(iteta,jzeta)**2/bm(iteta,jzeta)**3*b1m(c1,n1)*b1m(c2,n2) &
                                               - bn(iteta,jzeta)/bm(iteta,jzeta)**2*b1m(c1,n1)*b1n(c2,n2) &
                                               - half*bn(iteta,jzeta)**2/bm(iteta,jzeta)**2*b2m(c1,n1,c2,n2) ) &
                                               * surf(1)%ds(iteta,jzeta)
#endif
        enddo
       enddo
      enddo
     enddo


!!$     ! first  derivatives on currents
!!$     do c1 = 1, Ncoils
!!$      l1B(c1, 0) = l1B(c1, 0) + ( bn(iteta, jzeta) / bm(iteta,jzeta) * &
!!$      ( coil(c1)%Bx( 0, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By( 0, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz( 0, 0)*surf(1)%nz(iteta,jzeta) ) - bn**2/bm**2 * &
!!$      ( coil(c1)%Bx( 0, 0)*bx                      + coil(c1)%By( 0, 0)*by                      + coil(c1)%Bz( 0, 0)*bz                     ))              * surf(1)%ds(iteta,jzeta)
!!$     enddo
!!$
!!$     ! first  derivatives on geometry
!!$     do n1 = 1, Cdof
!!$      do c1 = 1, Ncoils
!!$        l1B(c1,n1) = l1B(c1,n1) + ( bn(iteta, jzeta) / bm(iteta,jzeta) * &
!!$      ( coil(c1)%Bx(n1, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By(n1, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz(n1, 0)*surf(1)%nz(iteta,jzeta) ) - bn**2/bm**2 * &
!!$      ( coil(c1)%Bx(n1, 0)*bx                      + coil(c1)%By(n1, 0)*by                      + coil(c1)%Bz(n1, 0)*bz                     )) * coil(c1)%I * surf(1)%ds(iteta,jzeta)
!!$      enddo
!!$     enddo
!!$
!!$     ! second derivatives on current & current
!!$     do c2 = 1, Ncoils
!!$      do c1 = 1, c2
!!$       l2B(c1,0,c2,0) = l2B(c1,0,c2,0) + 1/bm(iteta, jzeta)*( coil(c1)%Bx( 0, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By( 0, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz( 0, 0)*surf(1)%nz(iteta,jzeta) ) &
!!$                                                           *( coil(c2)%Bx( 0, 0)*surf(1)%nx(iteta,jzeta) + coil(c2)%By( 0, 0)*surf(1)%ny(iteta,jzeta) + coil(c2)%Bz( 0, 0)*surf(1)%nz(iteta,jzeta) ) &
!!$                   - 2*bn(iteta, jzeta)/bm(iteta, jzeta)**2*( coil(c1)%Bx( 0, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By( 0, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz( 0, 0)*surf(1)%nz(iteta,jzeta) ) &
!!$                                                           *( coil(c2)%Bx( 0, 0)*bx                      + coil(c2)%By( 0, 0)*by                      + coil(c2)%Bz( 0, 0)*bz                      )
!!$      enddo
!!$     enddo
!!$     
!!$     ! second derivatives on current & geometry
!!$     do n2 = 1, Cdof
!!$      do c2 = 1, Ncoils
!!$       do c1 = 1, Ncoils
!!$
!!$        l2B(c1, 0,c2,n2) = l2B(c1, 0,c2,n2) + 1/bm(iteta, jzeta)*( coil(c1)%Bx( 0, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By( 0, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz( 0, 0)*surf(1)%nz(iteta,jzeta) ) &
!!$                                                     *coil(c2)%I*( coil(c2)%Bx(n2, 0)*surf(1)%nx(iteta,jzeta) + coil(c2)%By(n2, 0)*surf(1)%ny(iteta,jzeta) + coil(c2)%Bz(n2, 0)*surf(1)%nz(iteta,jzeta) ) &
!!$                        - 2*bn(iteta, jzeta)/bm(iteta, jzeta)**2*( coil(c1)%Bx( 0, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By( 0, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz( 0, 0)*surf(1)%nz(iteta,jzeta) ) &
!!$                                                     *coil(c2)%I*( coil(c2)%Bx(n2, 0)*bx                      + coil(c2)%By(n2, 0)*by                      + coil(c2)%Bz(n2, 0)*bz                      )
!!$
!!$       enddo ! end c1
!!$
!!$       l2B(c2, 0,c2,n2) = l2B(c2, 0,c2,n2) + bn(iteta, jzeta)/bm(iteta, jzeta)*&
!!$                                                                 ( coil(c2)%Bx(n2, 0)*surf(1)%nx(iteta,jzeta) + coil(c2)%By(n2, 0)*surf(1)%ny(iteta,jzeta) + coil(c2)%Bz(n2, 0)*surf(1)%nz(iteta,jzeta) )
!!$
!!$      enddo ! end c2
!!$     enddo
!!$
!!$     ! second derivatives on geometry & geometry
!!$     do n2 = 1, Cdof
!!$      do c2 = 1, Ncoils
!!$       do n1 = 1, Cdof
!!$        do c1 = 1, c2
!!$
!!$         l2B(c1,n1,c2,n2) = l2B(c1,n1,c2,n2) + 1/bm(iteta, jzeta)*( coil(c1)%Bx(n1, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By(n1, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz(n1, 0)*surf(1)%nz(iteta,jzeta) ) &
!!$                                          *coil(c1)%I *coil(c2)%I*( coil(c2)%Bx(n2, 0)*surf(1)%nx(iteta,jzeta) + coil(c2)%By(n2, 0)*surf(1)%ny(iteta,jzeta) + coil(c2)%Bz(n2, 0)*surf(1)%nz(iteta,jzeta) ) &
!!$                         - 2*bn(iteta, jzeta)/bm(iteta, jzeta)**2*( coil(c1)%Bx( 0, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By( 0, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz( 0, 0)*surf(1)%nz(iteta,jzeta) ) &
!!$                                                      *coil(c2)%I*( coil(c2)%Bx(n2, 0)*bx                      + coil(c2)%By(n2, 0)*by                      + coil(c2)%Bz(n2, 0)*bz                      )
!!$ 
!!$        enddo ! end do c1
!!$        l2B(c2,n1,c2,n2) = l2B(c2,n1,c2,n2) +   bn(iteta, jzeta) * ( coil(c2)%Bx( n1,n2) * surf(1)%nx(iteta,jzeta) &
!!$             + coil(c2)%By( n1,n2) * surf(1)%ny(iteta,jzeta) &
!!$             + coil(c2)%Bz( n1,n2) * surf(1)%nz(iteta,jzeta) )  * coil(c2)%I * surf(1)%ds(iteta,jzeta)
!!$       enddo ! end do n1
!!$      enddo ! end do c2
!!$     enddo ! end do n2



!!$  lbnorm = lbnorm +  bn(iteta, jzeta)**2 * surf(1)%ds(iteta,jzeta)
!!$  do c1 = 1, Ncoils
!!$  l1B(c1, 0) = l1B(c1, 0) + bn(iteta, jzeta) * ( coil(c1)%Bx(  0, 0) * surf(1)%nx(iteta,jzeta) &
!!$                                               + coil(c1)%By(  0, 0) * surf(1)%ny(iteta,jzeta) &
!!$                                               + coil(c1)%Bz(  0, 0) * surf(1)%nz(iteta,jzeta) ) * surf(1)%ds(iteta,jzeta)              ! first  derivatives on currents
!!$  enddo
!!$
!!$  do n1 = 1, Cdof
!!$  do c1 = 1, Ncoils
!!$  l1B(c1,n1) = l1B(c1,n1) + bn(iteta, jzeta) * ( coil(c1)%Bx( n1, 0) * surf(1)%nx(iteta,jzeta) &
!!$                                               + coil(c1)%By( n1, 0) * surf(1)%ny(iteta,jzeta) &
!!$                                               + coil(c1)%Bz( n1, 0) * surf(1)%nz(iteta,jzeta) ) * coil(c1)%I * surf(1)%ds(iteta,jzeta) ! first  derivatives on geometry
!!$  enddo
!!$  enddo
!!$
!!$     do c2 = 1, Ncoils
!!$      do c1 = 1, c2
!!$       l2B(c1,0,c2,0) = l2B(c1,0,c2,0) + ( coil(c2)%Bx(  0, 0) * surf(1)%nx(iteta,jzeta) &
!!$            + coil(c2)%By(  0, 0) * surf(1)%ny(iteta,jzeta) &
!!$            + coil(c2)%Bz(  0, 0) * surf(1)%nz(iteta,jzeta) ) &
!!$            * ( coil(c1)%Bx(  0, 0) * surf(1)%nx(iteta,jzeta) &
!!$            + coil(c1)%By(  0, 0) * surf(1)%ny(iteta,jzeta) &
!!$            + coil(c1)%Bz(  0, 0) * surf(1)%nz(iteta,jzeta) ) * surf(1)%ds(iteta,jzeta)                      ! second derivatives on c & c
!!$      enddo
!!$     enddo
!!$
!!$     do n2 = 1, Cdof
!!$      do c2 = 1, Ncoils
!!$       do c1 = 1, Ncoils
!!$        l2B(c1, 0,c2,n2) = l2B(c1, 0,c2,n2) + (                    ( coil(c2)%Bx( n2, 0) * surf(1)%nx(iteta,jzeta) &
!!$             + coil(c2)%By( n2, 0) * surf(1)%ny(iteta,jzeta) &
!!$             + coil(c2)%Bz( n2, 0) * surf(1)%nz(iteta,jzeta) )  * coil(c2)%I &
!!$             * ( coil(c1)%Bx(  0, 0) * surf(1)%nx(iteta,jzeta) &
!!$             + coil(c1)%By(  0, 0) * surf(1)%ny(iteta,jzeta) &
!!$             + coil(c1)%Bz(  0, 0) * surf(1)%nz(iteta,jzeta) ) ) * surf(1)%ds(iteta,jzeta) ! second derivatives on c & g 
!!$       enddo ! end c1
!!$       l2B(c2, 0,c2,n2) = l2B(c2, 0,c2,n2) + ( bn(iteta, jzeta) * ( coil(c2)%Bx( n2, 0) * surf(1)%nx(iteta,jzeta) &
!!$            + coil(c2)%By( n2, 0) * surf(1)%ny(iteta,jzeta) &
!!$            + coil(c2)%Bz( n2, 0) * surf(1)%nz(iteta,jzeta) ) ) * surf(1)%ds(iteta,jzeta) ! second derivatives on c & g 
!!$      enddo ! end c2
!!$     enddo
!!$
!!$
!!$
!!$
!!$     do n2 = 1, Cdof
!!$      do c2 = 1, Ncoils
!!$       do n1 = 1, Cdof
!!$        do c1 = 1, c2
!!$         l2B(c1,n1,c2,n2) = l2B(c1,n1,c2,n2) + (                    ( coil(c2)%Bx( n2, 0) * surf(1)%nx(iteta,jzeta) &
!!$              + coil(c2)%By( n2, 0) * surf(1)%ny(iteta,jzeta) &
!!$              + coil(c2)%Bz( n2, 0) * surf(1)%nz(iteta,jzeta) )  * coil(c2)%I &
!!$              * ( coil(c1)%Bx( n1, 0) * surf(1)%nx(iteta,jzeta) &
!!$              + coil(c1)%By( n1, 0) * surf(1)%ny(iteta,jzeta) &
!!$              + coil(c1)%Bz( n1, 0) * surf(1)%nz(iteta,jzeta) ) ) * coil(c1)%I * surf(1)%ds(iteta,jzeta) ! second derivatives on g & g
!!$
!!$        enddo ! end do c1
!!$        l2B(c2,n1,c2,n2) = l2B(c2,n1,c2,n2) +   bn(iteta, jzeta) * ( coil(c2)%Bx( n1,n2) * surf(1)%nx(iteta,jzeta) &
!!$             + coil(c2)%By( n1,n2) * surf(1)%ny(iteta,jzeta) &
!!$             + coil(c2)%Bz( n1,n2) * surf(1)%nz(iteta,jzeta) )  * coil(c2)%I * surf(1)%ds(iteta,jzeta)
!!$       enddo ! end do n1
!!$      enddo ! end do c2
!!$     enddo ! end do n2


    enddo ! end do iteta
   enddo ! end do jzeta


   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   call MPI_REDUCE( lbnorm, bnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( bnorm, 1, 0)
   bnorm  = bnorm * half * discretefactor

   call MPI_REDUCE( bn, tbn, (1+Nteta)*(1+Nzeta), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( tbn,(1+Nteta)*(1+Nzeta), 0)

   if ( .not. allocated(t1B) ) allocate(t1B(1:Ncoils, 0:Cdof), stat=astat)
   if ( .not. allocated(t2B) ) allocate(t2B(1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), stat=astat)
   t1B = zero; t2B = zero

   array2size = Ncoils * ( Cdof + 1 )
   call MPI_REDUCE( l1B, t1B, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( t1B, array2size, 0 )

!!$   ! symmetric terms
!!$   do c2 = 1, Ncoils-1
!!$    do c1 = c2+1, Ncoils
!!$     l2B(c1,0,c2,0) = l2B(c2,0,c1,0)
!!$    enddo
!!$   enddo
!!$
!!$   do c2 = 1, Ncoils
!!$    do n1 = 1, Cdof
!!$     do c1 = 1, Ncoils
!!$      l2B(c1,n1,c2,0) = l2B(c2,0,c1,n1)
!!$     enddo
!!$    enddo
!!$   enddo
!!$
!!$   do n2 = 1, Cdof
!!$    do c2 = 1, Ncoils-1
!!$     do n1 = 1, Cdof
!!$      do c1 = c2+1, Ncoils
!!$       l2B(c1,n1,c2,n2) = l2B(c2,n2,c1,n1)
!!$      enddo
!!$     enddo ! end n1
!!$    enddo    ! end c2     
!!$   enddo       ! end n2

   array4size = array2size * array2size
   call MPI_REDUCE( l2B, t2B, array4size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( t2B, array4size, 0 )

   t1B = t1B * discretefactor
   t2B = t2B * discretefactor

   DALLOCATE(l1B)
   DALLOCATE(b1n)
   DALLOCATE(b1m)
   DALLOCATE(l2B)
   DALLOCATE(b2n)
   DALLOCATE(b2m)

  end select

  call MPI_REDUCE( lBx, SaveBx, (1+Nteta)*(1+Nzeta), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  RlBCAST( SaveBx,(1+Nteta)*(1+Nzeta), 0)
  call MPI_REDUCE( lBy, SaveBy, (1+Nteta)*(1+Nzeta), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  RlBCAST( SaveBy,(1+Nteta)*(1+Nzeta), 0)
  call MPI_REDUCE( lBz, SaveBz, (1+Nteta)*(1+Nzeta), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  RlBCAST( SaveBz,(1+Nteta)*(1+Nzeta), 0)

  do icoil = 1, Ncoils
   if( allocated(coil(icoil)%Bx)) then
    DALLOCATE(coil(icoil)%Bx)
    DALLOCATE(coil(icoil)%By)
    DALLOCATE(coil(icoil)%Bz)
   endif
  enddo

  tbn(Nteta,0:Nzeta-1) = tbn(0, 0:Nzeta-1)  !unassigned arrays; 2017/02/10
  tbn(0:Nteta-1,Nzeta) = tbn(0:Nteta-1, 0)
  tbn(Nteta, Nzeta) = tbn(0, 0)

  SaveBx(Nteta,0:Nzeta-1) = SaveBx(0, 0:Nzeta-1)  !unassigned arrays; 2017/05/04
  SaveBx(0:Nteta-1,Nzeta) = SaveBx(0:Nteta-1, 0)
  SaveBx(Nteta, Nzeta) = SaveBx(0, 0)
  SaveBy(Nteta,0:Nzeta-1) = SaveBy(0, 0:Nzeta-1)  !unassigned arrays; 2017/05/04
  SaveBy(0:Nteta-1,Nzeta) = SaveBy(0:Nteta-1, 0)
  SaveBy(Nteta, Nzeta) = SaveBy(0, 0)
  SaveBz(Nteta,0:Nzeta-1) = SaveBz(0, 0:Nzeta-1)  !unassigned arrays; 2017/05/04
  SaveBz(0:Nteta-1,Nzeta) = SaveBz(0:Nteta-1, 0)
  SaveBz(Nteta, Nzeta) = SaveBz(0, 0)
! bn(iteta,jzeta) & bm(iteta,jzeta) are still allocated.
  return
end subroutine bnormal

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bnormal2( nderiv )

  use kmodule, only: zero, half, one, pi2, sqrtmachprec, bsconstant, &
       coil, surf, NFcoil, cmt, smt, NDcoil, Ncoils, Nteta, Cdof, Nzeta, discretefactor, &
       bnorm, t1B, t2B, bn, bm, tbn, SaveBx, SaveBy, SaveBz, &
       ncpu, myid, ounit

  implicit none
  include "mpif.h"

  INTEGER           :: nderiv


  INTEGER           :: astat, ierr
  INTEGER           :: icoil, NN,  iteta, jzeta, c1, c2, n1, n2   ! icoil, iteta, jzeta, kseg are local
  REAL              :: r, rm2, rm3, rm4, bx, by, bz, lm, c12, lbnorm, start, finish
  REAL, allocatable :: l1B( :, :), l2B( :, :, :, :), b1n( :, :), b2n( :, :, :, :), b1m( :, :), b2m( :, :, :, :)
  INTEGER           :: array2size, array4size
  REAL              :: lBx(0:Nteta, 0:Nzeta), lBy(0:Nteta, 0:Nzeta), lBz(0:Nteta, 0:Nzeta)

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if ( .not. allocated(SaveBx) ) then
     SALLOCATE(SaveBx, (0:Nteta, 0:Nzeta), zero)
     SALLOCATE(SaveBy, (0:Nteta, 0:Nzeta), zero)
     SALLOCATE(SaveBz, (0:Nteta, 0:Nzeta), zero)
  endif

  lBx = zero; lBy = zero; lBz = zero
  
  !if ( .not. allocated(bn) ) allocate( bn(0:Nteta-1, 0:Nzeta-1) )
  !if ( .not. allocated(bm) ) allocate( bm(0:Nteta-1, 0:Nzeta-1) )
  !if ( .not. allocated(tbn)) allocate(tbn(0:Nteta-1, 0:Nzeta-1) )
  if ( .not. allocated(bn) ) allocate( bn(0:Nteta, 0:Nzeta) )  ! for saving the Bn data; 2017/02/10
  if ( .not. allocated(tbn)) allocate(tbn(0:Nteta, 0:Nzeta) )

  if ( .not. allocated(surf(1)%bnt) ) then
     allocate( surf(1)%bnt(0:Nteta, 0:Nzeta) ) ! Bnorm target from plasma currents or other coils
     surf(1)%bnt = zero
  endif
  
  NN = NFcoil; bnorm = zero; bn = zero; lbnorm = zero; tbn = zero

  do icoil = 1, Ncoils
   if(.not. allocated(coil(icoil)%Bx)) then 
    SALLOCATE( coil(icoil)%Bx, (0:Cdof, 0:Cdof), zero )
    SALLOCATE( coil(icoil)%By, (0:Cdof, 0:Cdof), zero )
    SALLOCATE( coil(icoil)%Bz, (0:Cdof, 0:Cdof), zero )
   else
    coil(icoil)%Bx = zero
    coil(icoil)%By = zero
    coil(icoil)%Bz = zero
   endif
  enddo

  select case ( nderiv )

   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  case(0)

   do jzeta = 0, Nzeta - 1
    do iteta = 0, Nteta - 1

     if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
     bx = zero; by = zero; bz = zero

     do icoil = 1, Ncoils

      call bfield0(icoil, iteta, jzeta, coil(icoil)%Bx(0,0), coil(icoil)%By(0,0), coil(icoil)%Bz(0,0))

      bx = bx + coil(icoil)%Bx( 0, 0) * coil(icoil)%I * pi2/coil(icoil)%D * bsconstant    
      by = by + coil(icoil)%By( 0, 0) * coil(icoil)%I * pi2/coil(icoil)%D * bsconstant
      bz = bz + coil(icoil)%Bz( 0, 0) * coil(icoil)%I * pi2/coil(icoil)%D * bsconstant

     enddo ! end do icoil

     !SaveBx(iteta, jzeta) = bx; SaveBy(iteta, jzeta) = by; SaveBz(iteta, jzeta) = bz
     lBx(iteta, jzeta) = bx; lBy(iteta, jzeta) = by; lBz(iteta, jzeta) = bz;
     bn(iteta, jzeta) = bx * surf(1)%nx(iteta,jzeta) + by * surf(1)%ny(iteta,jzeta) + bz * surf(1)%nz(iteta,jzeta)
#ifdef BNORM
     lbnorm = lbnorm +  (bn(iteta, jzeta) - surf(1)%bnt(iteta, jzeta))**2 * surf(1)%ds(iteta,jzeta)  !change to minus on 2017/02/07
#else
     lbnorm = lbnorm +  bn(iteta, jzeta)**2 * surf(1)%ds(iteta,jzeta)
#endif

    enddo ! end do iteta
   enddo ! end do jzeta

   call MPI_REDUCE( lbnorm, bnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( bnorm, 1, 0 )
   bnorm  = bnorm * half * discretefactor

   call MPI_REDUCE( bn, tbn, (1+Nteta)*(1+Nzeta), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( tbn,(1+Nteta)*(1+Nzeta), 0)

   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  case(1)

   SALLOCATE( l1B, (1:Ncoils, 0:Cdof), zero )
   SALLOCATE( b1n, (1:Ncoils, 0:Cdof), zero )
   !SALLOCATE( b1m, (1:Ncoils, 0:Cdof), zero )

   do jzeta = 0, Nzeta - 1
    do iteta = 0, Nteta - 1

     if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
     bx = zero; by = zero; bz = zero

     do icoil = 1, Ncoils

      call bfield1(icoil, iteta, jzeta, coil(icoil)%Bx(0:Cdof,0), coil(icoil)%By(0:Cdof,0), coil(icoil)%Bz(0:Cdof,0))

      coil(icoil)%Bx = coil(icoil)%Bx * pi2 / coil(icoil)%D * bsconstant
      coil(icoil)%By = coil(icoil)%By * pi2 / coil(icoil)%D * bsconstant
      coil(icoil)%Bz = coil(icoil)%Bz * pi2 / coil(icoil)%D * bsconstant

      bx = bx + coil(icoil)%Bx( 0, 0) * coil(icoil)%I
      by = by + coil(icoil)%By( 0, 0) * coil(icoil)%I
      bz = bz + coil(icoil)%Bz( 0, 0) * coil(icoil)%I

     enddo ! end do icoil

     lBx(iteta, jzeta) = bx; lBy(iteta, jzeta) = by; lBz(iteta, jzeta) = bz;
     bn(iteta, jzeta) = bx * surf(1)%nx(iteta,jzeta) + by * surf(1)%ny(iteta,jzeta) + bz * surf(1)%nz(iteta,jzeta)

     b1n = zero!; b1m = zero

     do c1 = 1, Ncoils
      b1n(c1, 0       )  =       coil(c1)%Bx( 0, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By( 0, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz( 0, 0)*surf(1)%nz(iteta,jzeta)
      do n1 = 1, Cdof
       b1n(c1,n1      )  =     ( coil(c1)%Bx(n1, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By(n1, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz(n1, 0)*surf(1)%nz(iteta,jzeta) ) * coil(c1)%I
      enddo
     enddo

     ! derivatives for bnomal cost function
#ifdef BNORM
     lbnorm = lbnorm +  (bn(iteta, jzeta) - surf(1)%bnt(iteta, jzeta))**2 * surf(1)%ds(iteta,jzeta)
#else
     lbnorm = lbnorm +  bn(iteta, jzeta)**2 * surf(1)%ds(iteta,jzeta)
#endif

     do n1 = 0, Cdof
      do c1 = 1, Ncoils
#ifdef BNORM
        l1B(c1,n1) = l1B(c1,n1) + ( (bn(iteta,jzeta)-surf(1)%bnt(iteta,jzeta))*b1n(c1,n1) ) * surf(1)%ds(iteta,jzeta)                                
#else
        l1B(c1,n1) = l1B(c1,n1) + ( bn(iteta,jzeta)*b1n(c1,n1) ) * surf(1)%ds(iteta,jzeta)
#endif
      enddo
     enddo

!!$     ! first derivatives of energy on DOFs
!!$     do c1 = 1, Ncoils
!!$      l1B(c1, 0) = l1B(c1, 0) + bn(iteta, jzeta) * ( coil(c1)%Bx(  0, 0) * surf(1)%nx(iteta,jzeta) &
!!$           + coil(c1)%By(  0, 0) * surf(1)%ny(iteta,jzeta) &
!!$           + coil(c1)%Bz(  0, 0) * surf(1)%nz(iteta,jzeta) ) * surf(1)%ds(iteta,jzeta)              ! first  derivatives on currents
!!$     enddo
!!$
!!$     do n1 = 1, Cdof
!!$      do c1 = 1, Ncoils
!!$       l1B(c1,n1) = l1B(c1,n1) + bn(iteta, jzeta) * ( coil(c1)%Bx( n1, 0) * surf(1)%nx(iteta,jzeta) &
!!$            + coil(c1)%By( n1, 0) * surf(1)%ny(iteta,jzeta) &
!!$            + coil(c1)%Bz( n1, 0) * surf(1)%nz(iteta,jzeta) ) * coil(c1)%I * surf(1)%ds(iteta,jzeta) ! first  derivatives on geometry
!!$      enddo ! end do n1
!!$     enddo ! end do c1


    enddo ! end do iteta
   enddo ! end do jzeta

   call MPI_REDUCE( lbnorm, bnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( bnorm, 1, 0 )
   bnorm  = bnorm * half * discretefactor

   call MPI_REDUCE( bn, tbn, (1+Nteta)*(1+Nzeta), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( tbn,(1+Nteta)*(1+Nzeta), 0)

   if ( .not. allocated(t1B) ) allocate(t1B(1:Ncoils, 0:Cdof), stat=astat)
   t1B = zero

   array2size = Ncoils * ( Cdof + 1 )
   call MPI_REDUCE( l1B, t1B, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( t1B, array2size, 0 )

   t1B = t1B * discretefactor

   DALLOCATE(l1B)
   DALLOCATE(b1n)
   !DALLOCATE(b1m)

   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case(2)

   SALLOCATE( l1B, (1:Ncoils, 0:Cdof), zero )
   SALLOCATE( b1n, (1:Ncoils, 0:Cdof), zero )
   !SALLOCATE( b1m, (1:Ncoils, 0:Cdof), zero )
   SALLOCATE( l2B, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero )
   SALLOCATE( b2n, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero )
   !SALLOCATE( b2m, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero )

   do jzeta = 0, Nzeta - 1
    do iteta = 0, Nteta - 1

     if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
     bx = zero; by = zero; bz = zero

     do icoil = 1, Ncoils

      call bfield2(icoil, iteta, jzeta, coil(icoil)%Bx(0:Cdof,0:Cdof), coil(icoil)%By(0:Cdof,0:Cdof), coil(icoil)%Bz(0:Cdof,0:Cdof))

      coil(icoil)%Bx = coil(icoil)%Bx * pi2/coil(icoil)%D * bsconstant
      coil(icoil)%By = coil(icoil)%By * pi2/coil(icoil)%D * bsconstant
      coil(icoil)%Bz = coil(icoil)%Bz * pi2/coil(icoil)%D * bsconstant

      bx = bx + coil(icoil)%Bx( 0, 0) * coil(icoil)%I
      by = by + coil(icoil)%By( 0, 0) * coil(icoil)%I
      bz = bz + coil(icoil)%Bz( 0, 0) * coil(icoil)%I

     enddo ! end do icoil

     lBx(iteta, jzeta) = bx; lBy(iteta, jzeta) = by; lBz(iteta, jzeta) = bz;
     bn(iteta, jzeta) = bx * surf(1)%nx(iteta,jzeta) + by * surf(1)%ny(iteta,jzeta) + bz * surf(1)%nz(iteta,jzeta)
     !bm(iteta, jzeta) = bx * bx + by * by + bz * bz ! magnitude square of B    

     ! derivatives for bn & bm terms
     b1n = zero; !b1m = zero
     b2n = zero; !b2m = zero

     ! b1n( 0, 0) = bn(iteta, jzeta)
     ! b1m( 0, 0) = bm(iteta, jzeta)

     do c1 = 1, Ncoils
      b1n(c1, 0       )  =       coil(c1)%Bx( 0, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By( 0, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz( 0, 0)*surf(1)%nz(iteta,jzeta)
      do n1 = 1, Cdof
       b1n(c1,n1      )  =     ( coil(c1)%Bx(n1, 0)*surf(1)%nx(iteta,jzeta) + coil(c1)%By(n1, 0)*surf(1)%ny(iteta,jzeta) + coil(c1)%Bz(n1, 0)*surf(1)%nz(iteta,jzeta) ) * coil(c1)%I 
      enddo
     enddo

     ! current & geometry
     do c2 = 1, Ncoils
      do n2 = 1, Cdof

       b2n(c2, 0,c2,n2)  =       coil(c2)%Bx(n2, 0)*surf(1)%nx(iteta,jzeta) + coil(c2)%By(n2, 0)*surf(1)%ny(iteta,jzeta) + coil(c2)%Bz(n2, 0)*surf(1)%nz(iteta,jzeta)
       b2n(c2,n2,c2, 0)  =       b2n(c2, 0,c2,n2)
       
      enddo
     enddo
     
      
     ! geometry & geometry 
     do c2 = 1, Ncoils
      do n2 = 1, Cdof
       do n1 = 1, Cdof
        
        b2n(c2,n1,c2,n2)  =     ( coil(c2)%Bx(n1,n2)*surf(1)%nx(iteta,jzeta) + coil(c2)%By(n1,n2)*surf(1)%ny(iteta,jzeta) + coil(c2)%Bz(n1,n2)*surf(1)%nz(iteta,jzeta) ) * coil(c2)%I

       enddo
      enddo
     enddo

     ! derivatives for bnomal cost function
#ifdef BNORM
     lbnorm = lbnorm +  (bn(iteta, jzeta) - surf(1)%bnt(iteta, jzeta))**2 * surf(1)%ds(iteta,jzeta)
#else
     lbnorm = lbnorm +  bn(iteta, jzeta)**2 * surf(1)%ds(iteta,jzeta)
#endif

     do n1 = 0, Cdof
        do c1 = 1, Ncoils
#ifdef BNORM
            l1B(c1,n1) = l1B(c1,n1) + ( (bn(iteta,jzeta)-surf(1)%bnt(iteta,jzeta))*b1n(c1,n1) ) * surf(1)%ds(iteta,jzeta)
#else
            l1B(c1,n1) = l1B(c1,n1) + ( bn(iteta,jzeta)*b1n(c1,n1) ) * surf(1)%ds(iteta,jzeta)
#endif
      enddo
     enddo

     do n2 = 0, Cdof
      do c2 = 1, Ncoils
       do n1 = 0, Cdof
        do c1 = 1, Ncoils
#ifdef BNORM
          l2B(c1,n1,c2,n2) = l2B(c1,n1,c2,n2) + ( b1n(c1,n1)*b1n(c2,n2) + (bn(iteta,jzeta)-surf(1)%bnt(iteta,jzeta))*b2n(c1,n1,c2,n2) ) * surf(1)%ds(iteta,jzeta)
#else
         l2B(c1,n1,c2,n2) = l2B(c1,n1,c2,n2) + ( b1n(c1,n1)*b1n(c2,n2) + bn(iteta,jzeta) *b2n(c1,n1,c2,n2) ) * surf(1)%ds(iteta,jzeta)
#endif
        enddo
       enddo
      enddo
     enddo

    enddo ! end do iteta
   enddo ! end do jzeta


   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   call MPI_REDUCE( lbnorm, bnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( bnorm, 1, 0)
   bnorm  = bnorm * half * discretefactor

   call MPI_REDUCE( bn, tbn, (1+Nteta)*(1+Nzeta), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( tbn,(1+Nteta)*(1+Nzeta), 0)

   if ( .not. allocated(t1B) ) allocate(t1B(1:Ncoils, 0:Cdof), stat=astat)
   if ( .not. allocated(t2B) ) allocate(t2B(1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), stat=astat)
   t1B = zero; t2B = zero

   array2size = Ncoils * ( Cdof + 1 )
   call MPI_REDUCE( l1B, t1B, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( t1B, array2size, 0 )

!!$   ! symmetric terms
!!$   do c2 = 1, Ncoils-1
!!$    do c1 = c2+1, Ncoils
!!$     l2B(c1,0,c2,0) = l2B(c2,0,c1,0)
!!$    enddo
!!$   enddo
!!$
!!$   do c2 = 1, Ncoils
!!$    do n1 = 1, Cdof
!!$     do c1 = 1, Ncoils
!!$      l2B(c1,n1,c2,0) = l2B(c2,0,c1,n1)
!!$     enddo
!!$    enddo
!!$   enddo
!!$
!!$   do n2 = 1, Cdof
!!$    do c2 = 1, Ncoils-1
!!$     do n1 = 1, Cdof
!!$      do c1 = c2+1, Ncoils
!!$       l2B(c1,n1,c2,n2) = l2B(c2,n2,c1,n1)
!!$      enddo
!!$     enddo ! end n1
!!$    enddo    ! end c2     
!!$   enddo       ! end n2

   array4size = array2size * array2size
   call MPI_REDUCE( l2B, t2B, array4size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( t2B, array4size, 0 )

   t1B = t1B * discretefactor
   t2B = t2B * discretefactor

   DALLOCATE(l1B)
   DALLOCATE(b1n)
   !DALLOCATE(b1m)
   DALLOCATE(l2B)
   DALLOCATE(b2n)
   !DALLOCATE(b2m)

  end select

  call MPI_REDUCE( lBx, SaveBx, (1+Nteta)*(1+Nzeta), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  RlBCAST( SaveBx,(1+Nteta)*(1+Nzeta), 0)
  call MPI_REDUCE( lBy, SaveBy, (1+Nteta)*(1+Nzeta), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  RlBCAST( SaveBy,(1+Nteta)*(1+Nzeta), 0)
  call MPI_REDUCE( lBz, SaveBz, (1+Nteta)*(1+Nzeta), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  RlBCAST( SaveBz,(1+Nteta)*(1+Nzeta), 0)

  SaveBx(Nteta,0:Nzeta-1) = SaveBx(0, 0:Nzeta-1)  !unassigned arrays; 2017/05/04
  SaveBx(0:Nteta-1,Nzeta) = SaveBx(0:Nteta-1, 0)
  SaveBx(Nteta, Nzeta) = SaveBx(0, 0)
  SaveBy(Nteta,0:Nzeta-1) = SaveBy(0, 0:Nzeta-1)  !unassigned arrays; 2017/05/04
  SaveBy(0:Nteta-1,Nzeta) = SaveBy(0:Nteta-1, 0)
  SaveBy(Nteta, Nzeta) = SaveBy(0, 0)
  SaveBz(Nteta,0:Nzeta-1) = SaveBz(0, 0:Nzeta-1)  !unassigned arrays; 2017/05/04
  SaveBz(0:Nteta-1,Nzeta) = SaveBz(0:Nteta-1, 0)
  SaveBz(Nteta, Nzeta) = SaveBz(0, 0)

  do icoil = 1, Ncoils
   if( allocated(coil(icoil)%Bx)) then
    DALLOCATE(coil(icoil)%Bx)
    DALLOCATE(coil(icoil)%By)
    DALLOCATE(coil(icoil)%Bz)
   endif
  enddo

  tbn(Nteta,0:Nzeta-1) = tbn(0, 0:Nzeta-1)  !unassigned arrays; 2017/02/10
  tbn(0:Nteta-1,Nzeta) = tbn(0:Nteta-1, 0)
  tbn(Nteta, Nzeta) = tbn(0, 0)
! bn(iteta,jzeta) & bm(iteta,jzeta) are still allocated.
  return
end subroutine bnormal2
