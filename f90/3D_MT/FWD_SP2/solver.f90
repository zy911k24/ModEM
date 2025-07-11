! *****************************************************************************
! module containing iterative equation solvers. Uses operators and
! pre-conditioners defined in SG3DFWC1 to solve equations for divergence
! correction, induction operator. 
! Source code is completely general; only the
! module interface is  specific to implementation of operators in SG3DFWC1
! 
! modified to translate cvectors into normal arrays to be used with sparse 
! matrix operation
module solver

   use math_constants   ! math/ physics constants
   use utilities, only: isnan
   use spoptools        ! for sparse-matrix operations
#if defined(FG) 
   use mpi
#if (defined(CUDA) || defined(HIP))
   use Declaration_MPI, only: comm_nccl
#endif
#endif
#if defined(CUDA)
   use cudaFortMap      ! cuda GPU api bindings for fortran
#elif defined(HIP)
   use hipFortMap       ! hip GPU api bindings for fortran
#endif
   !use griddef	! staggered grid definitions
   !use sg_scalar
   !use sg_vector
   implicit none

  ! solverControl_t as a data type controls diagnostic tools for joint forward
  ! modeling-inversion scheme
  type  :: solverControl_t
     ! maximum number of iterations in one call to iterative solver
     integer                                               :: maxIt
     ! convergence criteria: return from solver if relative error < tol
     real (kind=prec)                                      :: tol
     ! actual number of iterations before return
     integer                                               :: niter
     ! relative error for each iteration
     real (kind=prec), pointer, dimension(:)               :: rerr
     ! logical variable indicating if algorithm "failed"
     logical                                               :: failed = .false.
  end type solverControl_t
! currently only BiCG is supported for FG and GPU versions
#if defined(FG)
  interface BICG
      module procedure BiCG
      module procedure BiCGfg
#if defined(CUDA) 
      module procedure cuBiCGfg
#endif 
  end interface
#else
  interface BICG
      module procedure BiCG
#if defined(CUDA)
      module procedure cuBiCG
#elif defined(HIP)
      module procedure hipBiCG
#endif 
  end interface
#endif

Contains


! *****************************************************************************
! Solver contains subroutines for: 
! a) PCG- a quasi-generic pre-conditioned congugate gradient, and 
! b) QMR - Quasi-Minimal Residual method (pre-conditioned, no look-ahead)
! c) BICG - bicgstab Stablilized Biconjugate Gradient method
! d) TFQMR - Transpose-free QMR (don't use it with CC-DC, as TFQMR doesn't 
!            work well with frequent interuption) 
! *****************************************************************************

subroutine PCG(b,x,PCGiter)
  ! Purpose: a quasi-generic pre-conditioned conjugate gradient
  ! routine.
  ! solves Ax = b
  ! 
  ! modified to translate cvectors/scalars into normal arrays to be used with 
  ! sparse matrix operation. That said, the implementation is still the old way
  ! used in matrix free method:
  ! the real A is already defined in the DivCgrad function for divergence 
  ! correction.
    use modeloperator3d, only: A => DivCgrad, Minv => DivCgradILU

  implicit none
  complex (kind=prec), intent(in), dimension(:)    :: b
  complex (kind=prec), intent(inout),dimension(:)  :: x
  type (solverControl_t), intent(inout)            :: PCGiter

  ! local variables
  complex(kind=prec),allocatable,dimension(:)  :: r,s,p,q
  complex(kind=prec)    :: beta,alpha,delta,deltaOld
  real(kind=prec)       :: bnorm, rnorm
  integer               :: i
  integer               :: xsize


  xsize=size(x,1)
  ! Allocation of r, z, p, q
  ! residual
  allocate(r(xsize))
  allocate(s(xsize))
  allocate(p(xsize))
  allocate(q(xsize))

  ! r = b-A*x
  call A(x,r) ! in DivCorr, x should be zero (no initial guess)
  r = b - r   ! hence r should be indentical to b
  bnorm = sqrt(dot_product(b,b))
  if (isnan(bnorm)) then
  ! this usually means an inadequate model, in which case Maxwell's fails
      write(0,*) 'Error: b in PCG contains NaNs; exiting...'
      stop
  endif
  rnorm = sqrt(dot_product(r,r))
  i = 1
  PCGiter%rerr(i) =real(rnorm/bnorm)

  loop: do while ((PCGiter%rerr(i).gt.PCGiter%tol).and.(i.lt.PCGiter%maxIt))
     ! precondiction first
     call Minv(r,s)
     delta = dot_product(r,s)
     if(i.eq.1) then
        p = s
     else
        beta = delta/deltaOld
        p = s + beta*p
     end if
     call A(p,q)
     alpha = delta/dot_product(p,q)
     x = x + p * alpha
     r = r - q * alpha
     deltaOld = delta
     i = i + 1
     rnorm = sqrt(dot_product(r,r))
     if (isnan(rnorm)) then
         write(0,*) 'Error: residual in PCG contains NaNs; exiting...'
         stop
     endif
     PCGiter%rerr(i) = real(rnorm/bnorm)
  end do loop

  PCGiter%niter = i

  ! explicitly deallocate the temperary work arrays
  ! seems not necesary in fortran
  deallocate(r)
  deallocate(s)
  deallocate(p)
  deallocate(q)

end subroutine PCG ! PCG


! *****************************************************************************
subroutine QMR(b,x,QMRiter)
  ! Purpose ... a quasi-minimal residual method routine, set up for solving
  ! A x = b using routines in  mult_A. Actual code is generic, but interface
  ! is fairly specific

  ! interface...........
  ! redefining some of the interfaces for our convenience (locally)
  ! generic routines for vector operations for edge/ face nodes
  ! in a staggered grid
  ! 
  ! modified to be used with sparse matrix operators
  ! Again this is not that generic as the A is already defined in Mult_Aii... 
  use modeloperator3d, only: A => mult_Aii, M1solve => PC_Lsolve,        &
     &                         M2solve => PC_Usolve

  implicit none
  !  b is right hand side
  complex (kind=prec), intent(in), dimension(:)    :: b
  !  solution vector is x ... on input is provided with the initial
  !  guess, on output is the most recent iterate
  complex (kind=prec), intent(inout),dimension(:)  :: x
  type (solverControl_t), intent(inout) :: QMRiter

   ! local variables
  complex(kind=prec),allocatable,dimension(:)  :: XMIN,R,VT
  complex(kind=prec),allocatable,dimension(:)  :: Y,Z,WT,V,W,YT,ZT
  complex(kind=prec),allocatable,dimension(:)  :: P,Q,PT,D,S
  complex (kind=prec)          :: ETA,PDE,EPSIL,RDE,BETA,DELTA,RHO
  complex (kind=prec)          :: PSI,RHO1,GAMM,GAMM1,THET,THET1,TM2
  real (kind=prec)             :: bnorm,rnorm,rnormin
  complex (kind=prec)          :: rhoInv,psiInv
  integer                      :: iter, xsize
  logical                      :: adjoint, ilu_adjt

  xsize=size(x,1)
  ! Allocate work arrays
  allocate(XMIN(xsize))
  allocate(R(xsize))
  allocate(VT(xsize))
  allocate(Y(xsize))
  allocate(Z(xsize))
  allocate(WT(xsize))
  allocate(V(xsize))
  allocate(W(xsize))
  allocate(YT(xsize))
  allocate(ZT(xsize))
  allocate(P(xsize))
  allocate(Q(xsize))
  allocate(PT(xsize))
  allocate(D(xsize))
  allocate(S(xsize))

  ! NOTE: this iterative solver is QMR without look-ahead
  ! patterned after the scheme given on page 24 of Barrett et al.
  ! "Templates for the solution of linear systems of equations:
  ! Building blocks for iterative methods"
  ! Note that there are a couple of small differences, due to
  ! the fact that our system is complex

  ! R=b-Ax 
  adjoint = .false.
  call A(x,adjoint,R)
  ! R= b - Ax, for inital guess x, that has been inputted to the routine
  R = b - R

  ! Norm of rhs, residual
  bnorm = CDSQRT(dot_product(b, b))
  rnorm = CDSQRT(dot_product(R, R))
  rnormin = rnorm
  XMIN = x

  ! this usually means an inadequate model, in which case Maxwell's fails
  if (isnan(abs(bnorm))) then
      write(0,*) 'Error: b in QMR contains NaNs; exiting...'
      stop
  endif

  !  iter is iteration counter
  iter = 1
  QMRiter%rerr(iter) = real(rnorm/bnorm)
  ! write(6,*) 'initial residual:', QMRiter%rerr(iter) 

  ! L
  VT = R
  ilu_adjt = .false.
  call M1solve(VT,ilu_adjt,Y)
!  Y = VT
  RHO = CDSQRT(dot_product(Y,Y))
  ! U
  WT = R
  ilu_adjt = .true.
  call M2solve(WT,ilu_adjt,Z)
!  Z = WT
  PSI  = CDSQRT(dot_product(Z,Z))
  GAMM = C_ONE
  ETA  = C_MinusONE

  ! the do loop goes on while the relative error is greater than the tolerance
  ! and the iterations are less than maxIt
  loop: do while ((QMRiter%rerr(iter).gt.QMRiter%tol).and.&
       (iter.lt.QMRiter%maxIt))
      if ((RHO.eq.C_ZERO).or.(PSI.eq.C_ZERO)) then
        QMRiter%failed = .true.
        write(0,*) 'QMR FAILED TO CONVERGE : RHO'
        write(0,*) 'QMR FAILED TO CONVERGE : PSI'
        stop
      endif

      rhoInv = (1/RHO)*cmplx(1.0, 0.0, 8)
      psiInv = (1/PSI)*cmplx(1.0, 0.0, 8)
      V = VT * rhoInv
      W = WT * psiInv
      Y = Y * rhoinv
      Z = Z * psiinv

      DELTA = dot_product(Z,Y)
      if(DELTA.eq.C_ZERO) then
        QMRiter%failed = .true.
        write(0,*) 'QMR FAILS TO CONVERGE : DELTA'
        exit
      endif

      ilu_adjt = .false.
      call M2solve(Y,ilu_adjt,YT)
!      YT = Y
      ilu_adjt = .true.
      Call M1solve(Z,ilu_adjt,ZT)
!      ZT = Z

      if (iter.eq.1) then
        P = YT
        Q = ZT
      else
      ! these calculations are only done when iter greater than 1
        PDE = PSI*DELTA/EPSIL
        RDE = RHO*CONJG(DELTA/EPSIL)
        P = YT - PDE * P
        Q = ZT - RDE * Q
      endif

      adjoint = .false.
      Call A(P, adjoint, PT)
      EPSIL = dot_product(Q,PT)
      if (EPSIL.eq.C_ZERO) then
        QMRiter%failed = .true.
        write(0,*) 'QMR FAILED TO CONVERGE : EPSIL'
        exit
      endif

      BETA = EPSIL/DELTA
      if (BETA.eq.C_ZERO) then
        QMRiter%failed = .true.
        write(0,*) 'QMR FAILED TO CONVERGE : BETA'
        exit
      endif
      VT = PT - BETA * V

      RHO1 = RHO
      ilu_adjt = .false.
      Call M1solve(VT, ilu_adjt, Y)
!      Y = VT
      RHO = CDSQRT(dot_product(Y,Y))

      adjoint = .true.
      Call A(Q, adjoint, WT)
      WT = WT - conjg(BETA) * W

      ilu_adjt = .true.
      Call M2solve(WT,ilu_adjt,Z)
!      Z = WT
      PSI = CDSQRT(dot_product(Z,Z))

      if (iter.gt.1) then
        THET1 = THET
      endif
      THET = RHO/(GAMM*CDABS(BETA))
      GAMM1 = GAMM
      GAMM = C_ONE/CDSQRT(C_ONE + THET*THET)
      if (GAMM.eq.C_ZERO) then
        write(0,*) 'QMR FAILS TO CONVERGE : GAMM'
        exit
      endif

      ETA = -ETA*RHO1*GAMM*GAMM/(BETA*GAMM1*GAMM1)
      if (iter.eq.1) then
        D = ETA * P
        S = ETA * PT
      else
        TM2 = THET1*THET1*GAMM*GAMM
        D =  ETA*P + TM2*D
        S =  ETA*PT + TM2*S
      endif

      x = x + D
      R = R - S
      rnorm = CDSQRT(dot_product(R,R))
      if (rnorm .lt. rnormin) then! update the best solution so far
         rnormin = rnorm
         XMIN = X
      end if
      iter = iter + 1

      ! Keeping track of errors
      ! QMR book-keeping between divergence correction calls
      QMRiter%rerr(iter) = real(rnorm/bnorm)
      !write(6,*) 'iter # ', iter, 'residual:', QMRiter%rerr(iter)
  end do loop

  QMRiter%niter = iter
  ! x = XMIN ! use the last instead of the best
  QMRiter%rerr(iter) = real(rnormin/bnorm)

  ! deallocate all the work arrays
  ! seems not necesary here
  deallocate(XMIN)
  deallocate(R)
  deallocate(VT)
  deallocate(Y)
  deallocate(Z)
  deallocate(WT)
  deallocate(V)
  deallocate(W)
  deallocate(YT)
  deallocate(ZT)
  deallocate(P)
  deallocate(Q)
  deallocate(PT)
  deallocate(D)
  deallocate(S)

end subroutine QMR ! qmr

! *****************************************************************************
subroutine TFQMR(b,x,KSPiter,adjt)
  ! a transpose-free version of Quasi-Minimum Residue Algorithm,
  ! set up for solving
  ! A x = b using routines in  mult_Aii.
  ! solves for the interior (edge) field
  ! see: 
  ! Freund, Roland, A transpose-free quasi-minimal residual algorithm for
  ! non-Hermitian linear systems, SIAM J. Sci. Comp., 14 (1993), 470--482.
  !
  ! to me TFQMR is something like a middle ground between the Conjugate 
  ! Gradient and the Minimal Residual methods - the good part (like BiCG) 
  ! is this does not need the y = A^T*x operation like the original QMR
  ! therefore it is a good choice for the SP2 framework - as the CCGD 
  ! matrix is no longer symmetric 
  ! 
  ! the TFQMR also converges smoothly, which seems more stable comparing 
  ! with the BiCG (which is literally like a roller-coaster)
  ! the code is therefore easier to coupe with when the mixed precision
  ! method is considered. 
  ! 
  ! the downside of TFQMR, is that it doesn't work well when interupted 
  ! frequently (needs a long Krylov subspace) 
  ! so it won't work well with CC-DC - you have been warned
  ! 
  ! modified from my matlab version of TFQMR...
  ! so the naming might sound a little different from conventional ones
  ! also added the optional adjoint to solve adjoint system A^Tx = b 
  ! 
  ! NOTE: like BICG, TFQMR performs two sub line searches within a 
  !      iteration, but here we only store the relerr for the second sub
  !      just to be compatitive with QMR
  ! 
  ! interface...........
  ! redefining some of the interfaces for our convenience (locally)
  ! generic routines for vector operations for edge/ face nodes
  ! in a staggered grid
  !
  ! NOTE: this has not been extensively tested 
  ! if you have time reading this, test it!
  use modeloperator3d, only: A => mult_Aii, M1solve => PC_Lsolve,        &
     &                                      M2solve => PC_Usolve

  implicit none
  !  b is right hand side
  complex (kind=prec), intent(in), dimension(:)    :: b
  !  solution vector is x ... on input is provided with the initial
  !  guess, on output is the iterate with smallest residual. 
  !
  complex (kind=prec), intent(inout),dimension(:)  :: x
  type (solverControl_t), intent(inout)            :: KSPiter
  logical,intent(in),optional                      :: adjt

  ! local variables
  complex (kind=prec),allocatable,dimension(:)  :: R, R0, V, W, D
  complex (kind=prec),allocatable,dimension(:)  :: AX, AY, AD, Y, YP1
  complex (kind=prec),allocatable,dimension(:)  :: PY, PY1 ! preconditioned Y
  complex (kind=prec),allocatable,dimension(:)  :: xmin, xhalf
  real    (kind=prec)                           :: rnorm, rnorm1, bnorm,rnormin
  real    (kind=prec)                           :: xnorm, dnorm, btol
  real    (kind=prec)                           :: THETA, TAU, C
  complex (kind=prec)                           :: RHO, RHO1, ALPHA, BETA
  complex (kind=prec)                           :: ETA, SIGMA
  integer                                       :: iter, xsize, imin
  integer                                       :: maxiter
  logical                                       :: adjoint, ilu_adjt, converged

  ! stagnant check
  integer                                       :: maxstagsteps, stag
  integer                                       :: restarted, maxrestarts
  integer                                       :: last_restart, interval
  logical                                       :: restart
  real (kind=prec)                              :: eps = R_tiny
 
  if (present(adjt)) then
      adjoint = adjt
      ilu_adjt = adjt
  else
      adjoint = .false.
      ilu_adjt = .false.
  endif
  xsize = size(x,1)
  ! Norm of rhs
  bnorm = SQRT(dot_product(b, b))
  if (isnan(bnorm)) then
      write(0,*) 'Error: b in TFQMR contains NaNs; exiting...'
      stop
  else if ( bnorm .eq. 0.0) then ! zero rhs -> zero solution
      write(0,*) 'Warning: b in TFQMR has all zeros, returning zero solution'
      x = b 
      KSPiter%niter=1
      KSPiter%failed=.false.
      KSPiter%rerr=0.0
      return
  endif
  ! allocate the local variables
  allocate(R(xsize))
  ! now calculate the (original) residual
  call A(x,adjoint,R)
  ! R= b - Ax, for inital guess x, that has been input to the routine
  R = b - R
  ! Norm of residual
  rnorm = SQRT(dot_product(R, R))
  btol = KSPiter%tol * bnorm
  if ( rnorm .le. btol ) then ! the first guess is already good enough
     ! returning
      KSPiter%niter=1
      KSPiter%failed=.false.
      KSPiter%rerr(1)=real(rnorm/bnorm)
      deallocate(R)
      return 
  else
      ! allocate the rest here
      allocate(xhalf(xsize))
      allocate(xmin(xsize))
      allocate(AX(xsize))
      allocate(Y(xsize))
      allocate(YP1(xsize))
      allocate(AY(xsize))
      allocate(R0(xsize))
      allocate(D(xsize))
      allocate(AD(xsize))
      allocate(PY(xsize))
      allocate(PY1(xsize))
      allocate(W(xsize))
      allocate(V(xsize))
  end if 
!================= Now start configuring the iteration ===================!
  ! and store the best residual so far
  rnormin = rnorm
  xmin = x
  imin = 1
  KSPiter%rerr(1) = real(rnormin/bnorm)
  ! write(6,'(A, ES12.6)') ' initial relative residual = ', rnorm/bnorm
  converged = .false.
  maxiter = KSPiter%maxit 
  ! parameters for stagnant check
  stag = 0
  maxstagsteps = 5
  maxrestarts = 3
  last_restart = 0
  restarted = 0
  restart = .TRUE.
  ! interval = 120 !hard coded here
!============================== looooops! ================================!
  do iter= 1, maxIter
      ! if (mod((iter-last_restart), interval).eq.1) then
      !     restart = .true.
      ! end if
      if (restart) then
          ! store the first residual
          R0 = R
          W = R0
          Y = W
          ! L
          call M1solve(Y,ilu_adjt,PY1)
          ! U
          call M2solve(PY1,ilu_adjt,PY)
          ! A*Y = AY
          call A(PY,adjoint,AY)
          ! V = AY
          V = AY
          D = C_ZERO
          AD = D
          THETA = 0
          ETA = C_ZERO
          TAU = rnorm
          RHO = dot_product(R, R)
          RHO1 = RHO
          restart = .FALSE.
      end if
      alpha = RHO / dot_product(R0, V)
      YP1 = Y - alpha * V
      ! =================== first half of the TFQMR iteration ===========!
      W = W - ALPHA * AY
      ! SIGMA 0 should be zero
      SIGMA = (THETA*THETA/ALPHA)*ETA
      D = PY + SIGMA * D
      AD = AY + SIGMA * AD
      THETA = SQRT(dot_product(W, W))/TAU
      C = 1.0D0 / SQRT(1.0D0 + THETA * THETA)
      TAU = TAU * THETA * C ! = norm(r)/norm(b)
      ETA = C*C * ALPHA
      xnorm = SQRT(dot_product(x,x))
      dnorm = SQRT(dot_product(D,D))
      if (abs(ETA)*dnorm .lt. eps*xnorm) then
          stag = stag + 1
      else 
          stag = 0
      end if
      ! update the xhalf - first half of the iteration
      xhalf = x + eta*d
      ! update the residual
      R = R - ETA*AD
      rnorm = SQRT(dot_product(R,R))
      KSPiter%rerr(iter)=real(rnorm/bnorm)
      ! stagnant check
      ! write(6,*) 'iter # ',iter,'first half relres= ', KSPiter%rerr(iter)
      if (rnorm.le.btol) then
          ! double check if the residual is really less than tol 
          call A(xhalf,adjoint,AX)
          R = b - AX
          rnorm = SQRT(dot_product(R,R))
          if (rnorm .le. btol) then
              x = xhalf
              KSPiter%rerr(iter)=real(rnorm/bnorm)
              KSPiter%failed = .false.
              KSPiter%niter = iter
              converged = .true.
              exit
          end if
      end if
      if (stag .ge. maxstagsteps) then
          stag = 0 ! bail out
          ! try restarting
          restarted = restarted + 1
          if (restarted .gt. maxrestarts) then
              ! stagnant - exiting
              converged = .false.
              KSPiter%failed = .true.
              KSPiter%niter = iter
              exit 
          else
              x = xhalf
              last_restart = iter
              restart = .TRUE.
              continue
          end if
      end if
      if (rnorm .lt. rnormin) then
          ! store the best solution so far
          rnormin = rnorm
          xmin = xhalf
          imin = iter
      end if
      ! L
      call M1solve(YP1,ilu_adjt,PY1)
      ! U
      call M2solve(PY1,ilu_adjt,PY)
      ! update AY
      call A(PY,adjoint,AY)
      Y = YP1
      ! ================== second half of the TFQMR iteration ===========!
      W = W - ALPHA * AY
      SIGMA = (THETA*THETA/ALPHA)*ETA
      D = PY + SIGMA * D
      AD = AY + SIGMA * AD
      THETA = SQRT(dot_product(W, W))/TAU
      C = 1.0D0 / SQRT(1.0D0 + THETA * THETA)
      TAU = TAU * THETA * C
      ETA = C*C * ALPHA
      ! stag check
      xnorm = SQRT(dot_product(xhalf,xhalf))
      dnorm = SQRT(dot_product(D,D))
      if (abs(ETA)*dnorm .lt. eps*xnorm) then
          stag = stag + 1
      else 
          stag = 0
      end if
      ! update the x - second half of the iteration
      x = xhalf + ETA*D
      ! update the residual
      R = R - ETA*AD
      rnorm = SQRT(dot_product(R,R))
      KSPiter%rerr(iter)=real(rnorm/bnorm)
      ! write(6,*) 'iter # ',iter,'second half relres= ', KSPiter%rerr(iter)
      if (rnorm.lt.btol) then
          ! double check if the residual is really less than tol 
          call A(x,adjoint,AX)
          R = b - AX
          rnorm = SQRT(dot_product(R,R))
          if (rnorm .le. btol) then
              KSPiter%rerr(iter)=real(rnorm/bnorm)
              KSPiter%failed = .false.
              KSPiter%niter = iter
              converged = .true.
              exit
          end if
      end if
      if (stag .ge. maxstagsteps) then
          stag = 0 ! bail out
          ! try restarting
          restarted = restarted + 1
          if (restarted .gt. maxrestarts) then
              ! stagnant - exiting
              converged = .false.
              KSPiter%failed = .true.
              KSPiter%niter = iter
              exit 
          else
              last_restart = iter
              restart = .TRUE.
              continue
          end if
      end if
      if (rnorm .lt. rnormin) then
          ! store the best solution so far
          rnormin = rnorm
          xmin = x
          imin = iter
      end if
      ! update the RHO
      RHO = dot_product(R0,W)
      BETA = RHO / RHO1
      ! store the previous RHO
      RHO1 = RHO
      ! update Y
      YP1 = W + BETA*Y
      ! L
      call M1solve(YP1,ilu_adjt,PY1)
      ! U
      call M2solve(PY1,ilu_adjt,PY)
      ! partial update of V
      V = BETA*(AY+BETA*V)
      ! update AY
      call A(PY,adjoint,AY)
      ! second part of update for V
      V = AY + V
      Y = YP1
  end do
 
  if (.not. converged) then 
      ! it should be noted that this is the way my matlab version works
      ! the TFQMR will return the 'best' (smallest residual) iteration
      ! KSPiter%niter=imin ! comment this line
      KSPiter%niter = maxiter
      ! KSPiter%rerr(KSPiter%maxit) = KSPiter%rerr(imin)  ! and this line
      ! to use the last iteration result instead of the 'best'
  end if
  deallocate(xhalf)
  deallocate(xmin)
  deallocate(AX)
  deallocate(R)
  deallocate(R0)
  deallocate(Y)
  deallocate(YP1)
  deallocate(PY)
  deallocate(PY1)
  deallocate(AY)
  deallocate(V)
  deallocate(W)
  deallocate(D)
  deallocate(AD)
end subroutine TFQMR ! tfqmr 

! *****************************************************************************
subroutine BiCG(b,x,KSPiter,adjt)
  ! Stablized version of BiConjugate Gradient, set up for solving
  ! A x = b using routines in  mult_Aii.
  ! solves for the interior (edge) field
  !
  ! modified from my matlab version of BICGstab...
  ! so the naming might sound a little different from conventional ones
  ! also added the optional adjoint to solve adjoint system A^Tx = b
  ! 
  ! NOTE: BICG actually performs two sub (or half) line searches within a 
  !       iteration, but here we only store the relerr for the second sub
  !       just to be compatitive with QMR
  !
  ! 
  ! interface...........
  ! redefining some of the interfaces for our convenience (locally)
  ! generic routines for vector operations for edge/ face nodes
  ! in a staggered grid
  !
  ! NOTE: this has not been extensively tested - I believe it feels a 
  ! little unstable (dispite the name)...
  ! if you have time reading this, test it!
  use modeloperator3d, only: A => mult_Aii, M1solve => PC_Lsolve,        &
     &                                      M2solve => PC_Usolve

  implicit none
  !  b is right hand side
  complex (kind=prec), intent(in), dimension(:)    :: b
  !  solution vector is x ... on input is provided with the initial
  !  guess, on output is the iterate with smallest residual. 
  !
  complex (kind=prec), intent(inout),dimension(:)  :: x
  type (solverControl_t), intent(inout)            :: KSPiter
  logical,intent(in),optional                      :: adjt

  ! local variables
  complex (kind=prec),allocatable,dimension(:)  :: R,RT,V,T
  complex (kind=prec),allocatable,dimension(:)  :: P,PT,PH,S,ST,SH,AX
  complex (kind=prec),allocatable,dimension(:)  :: xhalf,xmin
  real    (kind=prec)                           :: rnorm, bnorm, rnormin, btol
  complex (kind=prec)                           :: RHO, ALPHA, BETA, OMEGA
  complex (kind=prec)                           :: RTV,TT,RHO1
  integer                                       :: iter, xsize, imin
  integer                                       :: maxiter
  logical                                       :: adjoint, ilu_adjt, converged
 
  if (present(adjt)) then
      adjoint = adjt
      ilu_adjt = adjt
  else
      adjoint = .false.
      ilu_adjt = .false.
  endif
  xsize = size(x,1)
  ! Norm of rhs
  bnorm = SQRT(dot_product(b, b))
  if (isnan(bnorm)) then
  ! this usually means an inadequate model, in which case Maxwell's fails
      write(0,*) 'Error: b in BICG contains NaNs; exiting...'
      stop
  else if ( bnorm .eq. 0.0) then ! zero rhs -> zero solution
      write(0,*) 'Warning: b in BICG has all zeros, returning zero solution'
      x = b 
      KSPiter%niter=1
      KSPiter%failed=.false.
      KSPiter%rerr=0.0
      return
  endif
  ! allocate the local variables
  allocate(R(xsize))
  ! now calculate the (original) residual
  call A(x,adjoint,R)
  ! R= b - Ax, for inital guess x, that has been inputted to the routine
  R = b - R
  ! Norm of residual
  rnorm = CDSQRT(dot_product(R, R))
  btol = KSPiter%tol * bnorm
  if ( rnorm .le. btol ) then ! the first guess is already good enough
     ! returning
      KSPiter%niter=1
      KSPiter%failed=.false.
      KSPiter%rerr(1)=real(rnorm/bnorm)
      deallocate(R)
      return 
  else
      ! allocate the rest here
      ! allocate(xhalf(xsize))
      allocate(xmin(xsize))
      allocate(AX(xsize))
      allocate(RT(xsize))
      allocate(P(xsize))
      allocate(PT(xsize))
      allocate(PH(xsize))
      allocate(S(xsize))
      allocate(ST(xsize))
      allocate(SH(xsize))
      allocate(V(xsize))
      allocate(T(xsize))
  end if 
!================= Now start configuring the iteration ===================!
  ! the adjoint (shadow) residual
  rnormin = rnorm
  KSPiter%rerr(1) = real(rnormin/bnorm)
  ! write(6,*) 'initial residual',  KSPiter%rerr(1)
  converged = .false.
  maxiter = KSPiter%maxit 
  imin = 0
  RHO = C_ONE
  OMEGA = C_ONE
  RT = R
  xmin = x
  imin = 1
!============================== looooops! ================================!
  do iter= 1, maxiter
      RHO1 = RHO
      RHO = dot_product(RT,R)
      if (RHO .eq. 0.0) then
          KSPiter%failed = .true.
          exit
      end if 
      if (iter .eq. 1) then
          P = R
      else 
          BETA = (RHO/RHO1)*(ALPHA/OMEGA) 
          if (BETA .eq. 0.0) then
              KSPiter%failed = .true.
              exit
          end if
          P= R + BETA * (P - OMEGA * V);
      end if 
      ! first half of the iteration
      ! L
      call M1solve(P,ilu_adjt,PT)
      ! U
      call M2solve(PT,ilu_adjt,PH)
!      PH = P
      call A(PH,adjoint,V)
      RTV = dot_product(RT,V)
      if (RTV.eq.0.0) then
          KSPiter%failed = .true.
          exit
      end if
      ALPHA = RHO / RTV
      if (ALPHA.eq.0.0) then
          KSPiter%failed = .true.
          exit
      end if
      x = x + ALPHA * PH ! the first half 
      S = R - ALPHA*V  !residual for the 0.5 x
      ! second half of the iteration
      ! L
      call M1solve(S,ilu_adjt,ST)
      ! U
      call M2solve(ST,ilu_adjt,SH)
!      SH = S
      call A(SH,adjoint,T)
      TT = dot_product(T,T)
      if (TT.eq.0.0) then
          KSPiter%failed = .true.
          exit
      end if
      OMEGA = dot_product(T,S)/TT
      if (OMEGA.eq.0.0) then
          KSPiter%failed = .true.
          exit
      end if
      x = x + OMEGA * SH  ! the second half 
      R = S - OMEGA * T  !residual for the 1.0 x
      rnorm = CDSQRT(dot_product(R,R))
      KSPiter%rerr(iter) = real(rnorm / bnorm)
      ! write(6,*) 'iter # ',iter,' x residual: ', KSPiter%rerr(2*iter)
      if (rnorm.lt.btol) then
          KSPiter%failed = .false.
          KSPiter%niter = iter
          converged = .true.
          exit
      end if
      if (rnorm .lt. rnormin) then
          rnormin = rnorm
          xmin = x
          imin = iter
      end if
  end do
 
  if (.not. converged) then 
      ! it should be noted that this is the way my matlab version works
      ! the bicg will return the 'best' (smallest residual) iteration
      ! x = xmin;  ! comment this line 
      KSPiter%niter=maxiter
      ! KSPiter%rerr(KSPiter%maxit) = KSPiter%rerr(imin)  ! and this line
      ! to use the last iteration result instead of the 'best'
  end if

  ! deallocate(xhalf)
  deallocate(xmin)
  deallocate(AX)
  deallocate(R)
  deallocate(RT)
  deallocate(P)
  deallocate(PT)
  deallocate(PH)
  deallocate(S)
  deallocate(ST)
  deallocate(SH)
  deallocate(V)
  deallocate(T)
end subroutine BiCG ! BICG
  
#if defined(FG) 
! *****************************************************************************
subroutine BiCGfg(b,x,KSPiter,comm_local,adjt)
  ! fine-grained parallel version 
  ! Stablized version of BiConjugate Gradient, set up for solving
  ! A x = b using routines in  mult_Aii.
  ! solves for the interior (edge) field
  !
  ! modified from my matlab version of BICGstab...
  ! so the naming might sound a little different from conventional ones
  ! also added the optional adjoint to solve adjoint system A^Tx = b
  ! 
  ! NOTE: BICG actually performs two sub (or half) line searches within a 
  !       iteration, but here we only store the relerr for the second sub
  !       just to be compatitive with QMR
  !
  ! 
  ! interface...........
  ! redefining some of the interfaces for our convenience (locally)
  ! generic routines for vector operations for edge/ face nodes
  ! in a staggered grid
  !
  ! NOTE: this has not been extensively tested - I believe it feels a 
  ! little unstable (dispite the name)...
  ! if you have time reading this, test it!
  use modeloperator3d, only: A => mult_Alcl,  M1solve => PC_LsolveLcl,     &
     &                                        M2solve => PC_UsolveLcl

  implicit none
  !  b is right hand side (which is actually b_local
  complex (kind=prec), intent(in), dimension(:)    :: b
  !  solution vector is x_local ... on input is provided with the initial
  !  guess, on output is the iterate with smallest residual. 
  !
  complex (kind=prec), intent(inout),dimension(:)  :: x
  type (solverControl_t), intent(inout)            :: KSPiter
  integer,intent(in)                               :: comm_local
  logical,intent(in),optional                      :: adjt

  ! local variables
  complex (kind=prec),allocatable,dimension(:)  :: R,RT,V,T, xmin
  complex (kind=prec),allocatable,dimension(:)  :: P,PT,PH,S,ST,SH,AX
  real    (kind=prec)                           :: rnorm, bnorm, rnormin, btol
  complex (kind=prec)                           :: RHO, ALPHA, BETA, OMEGA
  complex (kind=prec)                           :: RTV,TT,RHO1
  integer                                       :: iter, imin
  integer                                       :: maxiter, i, j
  logical                                       :: adjoint, ilu_adjt, converged
  ! fine-grain parallel related
  integer                                       :: rank_local, size_local, ierr
  complex (kind=prec),allocatable, dimension(:) :: xbuff, rbuff, bbuff
  integer,allocatable,dimension(:)              :: isizes, displs
  integer                                       :: fsize, lsize
  ! buffer 
  real    (kind=prec)                           :: bnrm
  complex (kind=prec)                           :: bdot
 
  if (present(adjt)) then
      adjoint = adjt
      ilu_adjt = adjt
  else
      adjoint = .false.
      ilu_adjt = .false.
  endif
  ! now see how many workers do we have
  call MPI_COMM_RANK(comm_local,rank_local,ierr)
  call MPI_COMM_SIZE(comm_local,size_local,ierr)
  ! firstly try to figure out how the data is distributed on 
  ! each process
  allocate(isizes(size_local))
  allocate(displs(size_local))
  ! local array size
  lsize = size(x)
  ! try to get the local row sizes for each process, initialize by zero
  isizes = 0
  ! NOTE this can be used in-place, if the comm_local is a intra-communicator
  ! it should be - but I will refrain to do that as who knows our users 
  ! will configure their MPI processes...
  call MPI_ALLGATHER(lsize,  1, MPI_INTEGER, isizes, 1, MPI_INTEGER&
 &        , comm_local, ierr)
  fsize = sum(isizes)
  ! also the displacement
  displs = 0
  do i=2,size_local
      displs(i) = sum(isizes(1:i-1))
  end do
  ! write(6, *) 'local size =', lsize, 'displs =', displs(rank_local+1), &
  !&   'full size = ', fsize,  'rank =',   rank_local
  !Norm of rhs
  ! the idea is to calculate the dot product of blocal for all processes and
  ! sum the result
  bnrm = dot_product(b, b)
  ! note that ALLREDUCE is equivalent to a REDUCE and a BCAST
  call MPI_ALLREDUCE(bnrm, bnorm, 1, MPI_DOUBLE, MPI_SUM,    &
 &        comm_local, ierr)
  bnorm = sqrt(bnorm)
  if (isnan(bnorm)) then
  ! this usually means an inadequate model, in which case Maxwell's fails
      write(0,*) 'Error: b in BICG contains NaNs; exiting...'
      stop
  else if ( bnorm .eq. 0.0) then ! zero rhs -> zero solution
      write(0,*) 'Warning: b in BICG has all zeros, returning zero solution'
      x = b 
      KSPiter%niter=1
      KSPiter%failed=.false.
      KSPiter%rerr=0.0
      return
  endif
  ! firstly allocate for the buffer
  allocate(xbuff(fsize))
  allocate(R(lsize))
  ! need to gather all local results to xbuff
  call MPI_ALLGATHERV(x, lsize, MPI_DOUBLE_COMPLEX, xbuff, isizes, &
 &        displs, MPI_DOUBLE_COMPLEX, comm_local, ierr)
  ! now calculate the (original) residual
  ! R= b - Ax, for inital guess x, that has been inputted to the routine
  call A(xbuff,x,adjoint,R)
  R = b - R
  bnrm = dot_product(R, R)
  ! note that ALLREDUCE is equivalent to a REDUCE and a BCAST
  call MPI_ALLREDUCE(bnrm, rnorm, 1, MPI_DOUBLE, MPI_SUM, comm_local, ierr)
  ! Norm of residual
  rnorm = sqrt(rnorm)
  btol = KSPiter%tol * bnorm
  if ( rnorm .le. btol ) then ! the first guess is already good enough
     ! returning
      KSPiter%niter=1
      KSPiter%failed=.false.
      KSPiter%rerr(1)=real(rnorm/bnorm)
      deallocate(xbuff)
      deallocate(R)
      return 
  end if 
  ! wait for the others
  call MPI_BARRIER(comm_local, ierr)
  ! allocate the local parameters of "lsize"
  allocate(xmin(lsize))
  allocate(RT(lsize))
  allocate(P(lsize))
  allocate(PT(lsize))
  allocate(PH(lsize))
  allocate(S(lsize))
  allocate(ST(lsize))
  allocate(SH(lsize))
  allocate(V(lsize))
  allocate(T(lsize))
!================= Now start configuring the iteration ===================!
  ! the adjoint (shadow) residual
  rnormin = rnorm
  KSPiter%rerr(1) = real(rnormin/bnorm)
  ! write(6,*) 'initial residual',  KSPiter%rerr(1)
  converged = .false.
  maxiter = KSPiter%maxit 
  imin = 0
  RHO = C_ONE
  OMEGA = C_ONE
  RT = R
  xmin = x
  imin = 1
!============================== looooops! ================================!
  do iter= 1, maxiter
      RHO1 = RHO
      bdot = dot_product(RT,R)
      ! note that ALLREDUCE is equivalent to a REDUCE and a BCAST
      call MPI_ALLREDUCE(bdot, RHO, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, &
 &          comm_local, ierr)
      if (RHO .eq. 0.0) then
          KSPiter%failed = .true.
          exit
      end if 
      if (iter .eq. 1) then
          P = R
      else 
          BETA = (RHO/RHO1)*(ALPHA/OMEGA) 
          if (BETA .eq. 0.0) then
              KSPiter%failed = .true.
              exit
          end if
          ! note we only update the "local" part here
          P = R + BETA * (P - OMEGA * V);
      end if 
      ! first half of the iteration
      ! L solve PT in L * PT = P
      call M1solve(P,ilu_adjt,PT)
      ! U solve PH in L * PH = PT
      call M2solve(PT,ilu_adjt,PH)
      ! need to gather all local PH results to xbuff
      call MPI_ALLGATHERV(PH, lsize, MPI_DOUBLE_COMPLEX, xbuff, isizes, &
 &        displs, MPI_DOUBLE_COMPLEX, comm_local, ierr)
      ! calculate V = A * PH
      call A(xbuff,PH,adjoint,V)
      ! RTV = dot(RT, V)
      bdot = dot_product(RT,V)
      ! note that ALLREDUCE is equivalent to a REDUCE and a BCAST
      call MPI_ALLREDUCE(bdot, RTV, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, &
 &          comm_local, ierr)
      if (RTV.eq.0.0) then
          KSPiter%failed = .true.
          exit
      end if
      ALPHA = RHO / RTV
      if (ALPHA.eq.0.0) then
          KSPiter%failed = .true.
          exit
      end if
      ! update Xhalf - note that we need to update "each local block" of x
      x = x + ALPHA * PH ! the first half 
      S = R - ALPHA * V  !residual for the 0.5 x
      ! second half of the iteration
      ! L solve ST in L * ST = S
      call M1solve(S,ilu_adjt,ST)
      ! U solve SH in L * SH = ST
      call M2solve(ST,ilu_adjt,SH)
      ! need to gather all local SH results to xbuff
      call MPI_ALLGATHERV(SH, lsize, MPI_DOUBLE_COMPLEX, xbuff, isizes, &
 &        displs, MPI_DOUBLE_COMPLEX, comm_local, ierr)
      ! calculate T = A * SH
      call A(xbuff,SH,adjoint,T)
      bdot = dot_product(T,T)
      ! note that ALLREDUCE is equivalent to a REDUCE and a BCAST
      call MPI_ALLREDUCE(bdot, TT, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, &
 &          comm_local, ierr)
      ! wait for the others
      call MPI_BARRIER(comm_local, ierr)
      if (TT.eq.0.0) then
          KSPiter%failed = .true.
          exit
      end if
      ! omega = dot(T,S)/TT
      bdot = dot_product(T,S)
      ! note that ALLREDUCE is equivalent to a REDUCE and a BCAST
      call MPI_ALLREDUCE(bdot, OMEGA, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, &
 &          comm_local, ierr)
      OMEGA = OMEGA/TT
      if (OMEGA.eq.0.0) then
          KSPiter%failed = .true.
          exit
      end if
      x = x + OMEGA * SH  ! the second half 
      R = S - OMEGA * T   ! residual for the "0.5" iteration
      ! rnorm = sqrt(dot(R,R))
      bnrm = dot_product(R,R)
      ! note that ALLREDUCE is equivalent to a REDUCE and a BCAST
      call MPI_ALLREDUCE(bnrm, rnorm, 1, MPI_DOUBLE, MPI_SUM, &
 &          comm_local, ierr)
      rnorm = sqrt(rnorm)
      ! wait for the others
      call MPI_BARRIER(comm_local, ierr)
      KSPiter%rerr(iter) = real(rnorm / bnorm)
      ! write(6,*) 'iteration # ', iter ,' relres= ', KSPiter%rerr(iter)
      if (rnorm.lt.btol) then
          KSPiter%failed = .false.
          KSPiter%niter = iter
          converged = .true.
          exit
      end if
      if (rnorm .lt. rnormin) then
          rnormin = rnorm
          xmin = x
          imin = iter
      end if
  end do
 
  if (.not. converged) then 
      ! it should be noted that this is the way my matlab version works
      ! the bicg will return the 'best' (smallest residual) iteration
      ! x = xmin;  ! comment this line 
      KSPiter%niter=maxiter
      ! KSPiter%rerr(KSPiter%maxit) = KSPiter%rerr(imin)  ! and this line
      ! to use the last iteration result instead of the 'best'
  end if

  deallocate(xmin)
  deallocate(xbuff)
  deallocate(R)
  deallocate(RT)
  deallocate(P)
  deallocate(PT)
  deallocate(PH)
  deallocate(S)
  deallocate(ST)
  deallocate(SH)
  deallocate(V)
  deallocate(T)
end subroutine BiCGfg ! BICGp
#endif

! *****************************************************************************
#if defined(CUDA)
subroutine cuBiCG(b,x,KSPiter,device_idx,adjt)
  ! Stablized version of BiConjugate Gradient, set up for solving
  ! A x = b 
  ! solves for the interior (edge) field
  !
  ! modified (no hell no) to call the CUDA lib to calculate with GPU 
  ! I kept the naming convention from my BiCGStab (see above CPU version)
  ! for instance X -> devPtrX, T -> devPtrT to manipulate the GPU memory
  ! also added the optional adjoint to solve adjoint system A^Tx = b
  ! 
  ! NOTE: BICG actually performs two sub (or half) line searches within a 
  !       iteration, but here we only store the relerr for the second sub
  !       just to be compatitive with QMR
  ! 
  ! interface...........
  ! redefining some of the interfaces for our convenience (locally)
  ! generic routines for vector operations for edge/ face nodes
  ! in a staggered grid
  !
  ! NOTE: this has not been extensively tested - I believe it feels a 
  ! little unstable (dispite the name)...
  ! if you have time reading this, test it!

      use modeloperator3d, only:  AAii, L, U, LH, UH,vOmegaMuSig
      use spoptools
      implicit none
      !  b is right hand side
      complex (kind=prec),intent(in),target,dimension(:)    :: b
      !  solution vector is x ... on input is provided with the initial
      !  guess, on output is the iterate with smallest residual. 
      complex (kind=prec),intent(inout),target,dimension(:) :: x
      type (solverControl_t),intent(inout)                  :: KSPiter
      integer,intent(in)                                    :: device_idx
      logical,intent(in),optional                           :: adjt
    
      ! local variables
      integer                                    :: n, nnz, iter, maxIter
      complex (kind=prec),pointer,dimension(:)   :: iwus

      real    (kind=prec), target                :: rnorm, bnorm, rnorm0, btol
      real    (kind=prec), target                :: xnorm
      complex (kind=prec)                        :: RHO1, ALPHA, BETA, OMEGA
      complex (kind=prec), target                :: RTV,TT,RHO
      real    (kind=SP), target                  :: ctime
      integer                                    :: TRANS
      integer                                    :: Lfillmode,Ldiagtype
      integer                                    :: Ufillmode,Udiagtype
      integer                                    :: ierr,ierr2, idx
      logical                                    :: converged, adjoint
      integer                                    :: interval
      logical                                    :: restart

      ! size to determine corresponding buffer/memory
      integer*8                :: Arowp1_i_size,Arow_d_size
      integer*8                :: Acol_d_size,Annz_i_size,Annz_d_size
      integer*8                :: Mrow_d_size,Lnnz_i_size,Lnnz_d_size
      integer*8                :: Unnz_i_size,Unnz_d_size
      integer(c_int),target    :: zeroLoc
      ! integer(c_int)         :: mbsize
      integer(c_size_t)        :: bsize, lbsize, ubsize

      ! --------------------- pointers to *host* memory ------------------ !   
      type(c_ptr) :: cublasHandle 
      type(c_ptr) :: cusparseHandle
      type(c_ptr) :: cuStream
      type(c_ptr) :: cuEvent1
      type(c_ptr) :: cuEvent2
      type(c_ptr) :: matA
      type(c_ptr) :: matL
      type(c_ptr) :: matU
      type(c_ptr) :: vecX
      type(c_ptr) :: vecY 
      type(c_ptr) :: LsolveHandle
      type(c_ptr) :: UsolveHandle
      type(c_ptr) :: ilu_info
      type(c_ptr) :: ArowPtr
      type(c_ptr) :: AcolPtr
      type(c_ptr) :: AvalPtr
      type(c_ptr) :: iwusPtr
      type(c_ptr) :: MrowPtr
      type(c_ptr) :: McolPtr
      type(c_ptr) :: MvalPtr
      type(c_ptr) :: xPtr  
      type(c_ptr) :: bPtr
      type(c_ptr) :: bnormPtr
      type(c_ptr) :: rnormPtr
      type(c_ptr) :: rnormPtr0
      type(c_ptr) :: rhoPtr
      type(c_ptr) :: rtvPtr
      type(c_ptr) :: ttPtr
      type(c_ptr) :: zeroLocPtr
      type(c_ptr) :: timePtr
      ! -------------------- pointers to *device* memory ----------------- ! 
      type(c_ptr) :: devPtrArow
      type(c_ptr) :: devPtrAcol
      type(c_ptr) :: devPtrAval
      type(c_ptr) :: devPtriwus
      type(c_ptr) :: devPtrLrow
      type(c_ptr) :: devPtrLcol
      type(c_ptr) :: devPtrLval
      type(c_ptr) :: devPtrUrow
      type(c_ptr) :: devPtrUcol
      type(c_ptr) :: devPtrUval
      type(c_ptr) :: devPtrX
      type(c_ptr) :: devPtrX0
      type(c_ptr) :: devPtrRHS
      type(c_ptr) :: devPtrR
      type(c_ptr) :: devPtrRT
      type(c_ptr) :: devPtrP
      type(c_ptr) :: devPtrPT
      type(c_ptr) :: devPtrPH
      type(c_ptr) :: devPtrS
      type(c_ptr) :: devPtrST
      type(c_ptr) :: devPtrSH
      type(c_ptr) :: devPtrV
      type(c_ptr) :: devPtrT
      type(c_ptr) :: devPtrAX
      type(c_ptr) :: buffer
      type(c_ptr) :: bufferL
      type(c_ptr) :: bufferU
      
      ! note we don't check the input parameters compatibility 
      ierr2 = 0
      converged = .FALSE.
      zeroLoc = 0
      if (present(adjt)) then 
          ! write(6,'(A)') ' adjt = ', adjt
          adjoint = adjt
          if (adjt) then
              TRANS = CUSPARSE_OPERATION_TRANSPOSE
              Lfillmode = CUSPARSE_FILL_MODE_UPPER
              Ufillmode = CUSPARSE_FILL_MODE_LOWER
              Ldiagtype = CUSPARSE_DIAG_TYPE_NON_UNIT
              Udiagtype = CUSPARSE_DIAG_TYPE_UNIT
          else
              TRANS = CUSPARSE_OPERATION_NON_TRANSPOSE
              Lfillmode = CUSPARSE_FILL_MODE_LOWER
              Ufillmode = CUSPARSE_FILL_MODE_UPPER
              Ldiagtype = CUSPARSE_DIAG_TYPE_UNIT
              Udiagtype = CUSPARSE_DIAG_TYPE_NON_UNIT
          end if
      else
          adjoint = .FALSE.
          TRANS = CUSPARSE_OPERATION_NON_TRANSPOSE
          Lfillmode = CUSPARSE_FILL_MODE_LOWER
          Ufillmode = CUSPARSE_FILL_MODE_UPPER
          Ldiagtype = CUSPARSE_DIAG_TYPE_UNIT
          Udiagtype = CUSPARSE_DIAG_TYPE_NON_UNIT
      end if
      
      ! firstly need to translate the ModEM SP datatypes to simple
      ! CSR plain vectors
      ! note we still keep the Ae = CCGDe + iwuse 
      ! idea as in ModEM SP2, we will test if this will save some time 
      ! (probably not much)
      ! call CSR_R2Cdiag(AAii,VOmegaMuSig,Aii)
      n = AAii%nrow
      nnz = size(AAii%col)
      ArowPtr = c_loc(AAii%row)
      AcolPtr = c_loc(AAii%col)
      AvalPtr = c_loc(AAii%val)
      ! other host pointers
      xPtr  = c_loc(x)  ! x = A \ b
      bPtr  = c_loc(b)  ! b = ones(n,1)
      allocate(iwus(n))
      iwus = VOmegaMuSig*ISIGN*CMPLX(0.0,1.0,8)
      iwusPtr = c_loc(iwus) ! i V omega mu sigma
      bnormPtr = c_loc(bnorm)
      rnormPtr = c_loc(rnorm)
      rnormPtr0 = c_loc(rnorm0)
      rhoPtr = c_loc(RHO)
      rtvPtr = c_loc(rtv)
      ttPtr = c_loc(tt)
      zeroLocPtr = c_loc(zeroLoc)
      ! pointer to store the time 
      ctime = 0.0 
      timePtr = c_loc(ctime)
      ! now remove the Aii matrix structure
      ! call deall_spMatCSR(Aii)
      
      ! get the size of the vector and matrix (need to allocate on the device)
      Arowp1_i_size=sizeof(AAii%row(1:n+1))
      Arow_d_size=sizeof(b(1:n))
      Acol_d_size=sizeof(x(1:n)) ! this is useful if A is not square
      Annz_i_size=sizeof(AAii%col(1:nnz))
      Annz_d_size=sizeof(AAii%val(1:nnz))

      ! select the current cuda device 
      ! note this only works for physical devices (not working for MIG devices)
      ierr = cudaSetDevice(device_idx);
      ierr2 = ierr2 + ierr
      ! firstly define the CUDA Stream and cuda handles
      ierr = cudaStreamCreateWithFlags(cuStream, cudaStreamNonBlocking) 
      ierr2 = ierr2 + ierr
      ! initialize the cusparse lib
      ierr = cusparseCreate(cusparseHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseSetStream(cusparseHandle,cuStream) 
      ierr2 = ierr2 + ierr
      ! now initialize the cublas lib
      ierr = cublasCreate(cublasHandle)
      ierr2 = ierr2 + ierr
      ierr = cublasSetStream(cublasHandle,cuStream)
      ierr2 = ierr2 + ierr
      ! and creates two cuda events to record the time
      ierr = cudaEventCreate(cuEvent1)
      ierr2 = ierr2 + ierr
      ierr = cudaEventCreate(cuEvent2)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I4)') 'Error during cuda initialize ',ierr2
          stop
      end if 
      ! record the event before memory manipulation
      ! ierr = cudaEventRecord(cuEvent1, cuStream)
      ! ierr2 = ierr2 + ierr
      ! write(6,*) 'Allocating GPU memory'
      ierr = cudaMalloc(devPtrX,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrX0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrRHS,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrR,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrRT,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrP,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrPT,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrPH,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrS,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrST,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrSH,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrT,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrV,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrAX,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrAval,Annz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrAcol,Annz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrArow,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtriwus,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error during CUDA allocation: ',ierr2
          stop
      end if 
      ! write(6,*) 'reset GPU memory to all zeros'
      ierr = cudaMemset(devPtrX,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrX0,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrRHS,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrR,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrRT,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrP,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrPT,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrPH,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrS,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrST,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrSH,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrT,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrV,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrAX,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrAval,0,Annz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrAcol,0,Annz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrArow,0,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtriwus,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I3)') 'Error during device memory reseting : ',ierr2
          stop
      end if 
      ! transfer memory over to GPU
      ! write(6,*) 'Transferring memory to GPU'
      ! initialize A
      ierr = cudaMemcpyAsync(devPtrArow,ArowPtr,Arowp1_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrAcol,AcolPtr,Annz_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrAval,AvalPtr,Annz_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtriwus,iwusPtr,Arow_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ! finally need to create the SpMatDescr 
      ierr = cusparseCreateCsr(matA, n, n, nnz, devPtrArow, devPtrAcol, &
     &       devPtrAval, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, &
     &       CUSPARSE_INDEX_BASE_ONE, CUDA_R_64F) 
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error assembling system matrix ", ierr2
          stop
      end if
      ! initialize rhs and x
      ierr = cudaMemcpyAsync(devPtrX,xPtr,Arow_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrRHS,bPtr,Arow_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error during cuda memcpy ", ierr2
          stop
      end if
      ! now deallocate temp array
      if ( c_associated(iwusPtr) ) then
          ! this is a little tricky as iwus is associated with
          ! iwusPtr in C
          nullify(iwus)
          call cf_free(iwusPtr)
      end if
      ! firstly we need to create two dense vectors 
      ! as the user-hostile developers in Nvidia think of a new idea 
      ! to mess up the interfaces
      ierr = cusparseCreateDnVec(vecX, n, devPtrX, CUDA_C_64F)
      ierr2 = ierr2 + ierr
      ierr = cusparseCreateDnVec(vecY, n, devPtrR, CUDA_C_64F)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, *) " error creating the dense Vecs "
          stop
      end if
      ierr = cusparseSpMV_bufferSize_cmplx(cusparseHandle,              &
     &        TRANS, C_ONE, matA, vecX, C_ONE, vecY,                   &
     &        CUDA_C_64F, CUSPARSE_SPMV_CSR_ALG2, bsize)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error estimating buffer for spMV ",  ierr2
          stop
      end if
      ! write(6,*) 'spmv buffersize = ', mbsize
      ! and finally (re)allocate the buffer
      ierr = cudaMalloc(buffer, bsize)
      ierr2 = ierr2 + ierr
      ! now let's deal with L and U
      ! write(6,'(A)') ' Setup L and U preconditioners on GPU'
      ! L first
      if (adjoint) then
          nnz = size(LH%col)
          ! sizes (row size is the same as A
          Lnnz_i_size=sizeof(LH%col(1:nnz))
          Lnnz_d_size=sizeof(LH%val(1:nnz))
          MrowPtr = c_loc(LH%row)
          McolPtr = c_loc(LH%col)
          MvalPtr = c_loc(LH%val)
      else
          nnz = size(L%col)
          ! sizes (row size is the same as A
          Lnnz_i_size=sizeof(L%col(1:nnz))
          Lnnz_d_size=sizeof(L%val(1:nnz))
          MrowPtr = c_loc(L%row)
          McolPtr = c_loc(L%col)
          MvalPtr = c_loc(L%val)
      end if
      ! initialize L in GPU memory
      ierr = cudaMalloc(devPtrLval,Lnnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrLcol,Lnnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrLrow,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! set everything to zeros
      ierr = cudaMemset(devPtrLval,0,Lnnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrLcol,0,Lnnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrLrow,0,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! now copy the values to device
      ierr = cudaMemcpyAsync(devPtrLval,MvalPtr,Lnnz_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrLcol,McolPtr,Lnnz_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrLrow,MrowPtr,Arowp1_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ! finish synchronizing
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      ! finally need to create the SpMatDescr 
      ierr = cusparseCreateCsr(matL, n, n, nnz, devPtrLrow, devPtrLcol, &
     &       devPtrLval, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, &
     &       CUSPARSE_INDEX_BASE_ONE, CUDA_C_64F) 
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error assembling L matrix ", ierr2
          stop
      end if
      ! now start messing up with SpSV
      ! now estimate the buffersize needed by SpSV (Lsolve)
      ! solves y in L*y = a*x (if a=1)
      ierr = cusparseSpMatSetAttribute(matL, CUSPARSE_SPMAT_FILL_MODE, &
    &        Lfillmode, 4)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpMatSetAttribute(matL, CUSPARSE_SPMAT_DIAG_TYPE, &
    &        Ldiagtype, 4)
      ierr2 = ierr2 + ierr
      ! still need to establish a context handler for Lsolve
      ierr = cusparseSpSV_createDescr(LsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_bufferSize_dcmplx(cusparseHandle,             &
     &        CUSPARSE_OPERATION_NON_TRANSPOSE, C_ONE, matL, vecX, vecY,&
     &        CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle, lbsize)
      ierr2 = ierr2 + ierr

      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error estimating buffer for Lsolve ",  ierr2
          stop
      end if
      ierr = cudaMalloc(bufferL, lbsize)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_analysis_dcmplx(cusparseHandle,         &
     &        CUSPARSE_OPERATION_NON_TRANSPOSE, C_ONE, matL, vecX, vecY, &
     &        CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle, bufferL)
      ierr2 = ierr2 + ierr
      ! write(6,*) 'spsv Lsolve buffersize = ', lbsize
      ! then U 
      if (adjoint) then
          nnz = size(UH%col)
          ! sizes (row size is the same as A
          Unnz_i_size=sizeof(UH%col(1:nnz))
          Unnz_d_size=sizeof(UH%val(1:nnz))
          MrowPtr = c_loc(UH%row)
          McolPtr = c_loc(UH%col)
          MvalPtr = c_loc(UH%val)
      else
          nnz = size(U%col)
          ! sizes (row size is the same as A
          Unnz_i_size=sizeof(U%col(1:nnz))
          Unnz_d_size=sizeof(U%val(1:nnz))
          MrowPtr = c_loc(U%row)
          McolPtr = c_loc(U%col)
          MvalPtr = c_loc(U%val)
      end if
      ! initialize U in GPU memory
      ierr = cudaMalloc(devPtrUval,Unnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrUcol,Unnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrUrow,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! set everything to zeros
      ierr = cudaMemset(devPtrUval,0,Unnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrUcol,0,Unnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrUrow,0,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! now copy the values to device
      ierr = cudaMemcpyAsync(devPtrUval,MvalPtr,Unnz_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrUcol,McolPtr,Unnz_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrUrow,MrowPtr,Arowp1_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ! finally need to create the SpMatDescr 
      ierr = cusparseCreateCsr(matU, n, n, nnz, devPtrUrow, devPtrUcol, &
     &       devPtrUval, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, &
     &       CUSPARSE_INDEX_BASE_ONE, CUDA_C_64F) 
      ierr2 = ierr2 + ierr
      ! now estimate the buffersize needed by SpSV (Usolve)
      ! solves y in U*y = a*x (if a=1)
      ierr = cusparseSpMatSetAttribute(matU, CUSPARSE_SPMAT_FILL_MODE, &
    &        Ufillmode, 4)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpMatSetAttribute(matU, CUSPARSE_SPMAT_DIAG_TYPE, &
    &        Udiagtype, 4)
      ierr2 = ierr2 + ierr
      ! still need to establish a context handler for Usolve
      ierr = cusparseSpSV_createDescr(UsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_bufferSize_dcmplx(cusparseHandle,             &
     &        CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matU, vecX, vecY,&
     &        CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, UsolveHandle, ubsize)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error estimating buffer for Usolve ",  ierr2
          stop
      end if
      ! write(6,*) 'spsv Usolve buffersize = ', ubsize
      ierr = cudaMalloc(bufferU, ubsize)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_analysis_dcmplx(cusparseHandle,         &
     &        CUSPARSE_OPERATION_NON_TRANSPOSE, C_ONE, matU, vecX, vecY, &
     &        CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, UsolveHandle, bufferU)
      ierr2 = ierr2 + ierr
      !write(6, '(A, I8)') " allocated SpMV/SpSV buffer on GPU (MB) ",&
     !&            bsize/1024/1024
      ! see how long it takes for the intialization part
      ! ierr = cudaEventRecord(cuEvent2, cuStream)
      ! ierr2 = ierr2 + ierr
      ! ierr = cudaDeviceSynchronize()
      ! ierr = cudaEventElapsedTime(timePtr, cuEvent1, cuEvent2)
      ! ierr2 = ierr2 + ierr
      ! write(6,*) 'initial GPU memory cpy time is ', ctime, 'ms'
      ! if (ierr2 .ne. 0 ) then
      !     write(6, *) " error recording the time ", ierr2
      !     stop
      ! end if
      ! now Compute the initial residual
      ! write(6,*) 'Computing initial residual'
      ! setup the vecX and vecY
      ierr = cusparseDnVecSetValues(vecX, devPtrX)
      ierr2 = ierr2 + ierr
      ! R = Diag(iwus)*X <-- diagonal
      call kernelc_hadac(devPtrX, devPtriwus, devPtrR, n, cuStream)
      ierr = cusparseDnVecSetValues(vecY, devPtrR)
      ierr2 = ierr2 + ierr
      ! now calculate the y = Ax
      ! R = CC*X + Diag(iwus)*X
      ! thank God (or whatever deserves it) that SpMV does not 
      ! need to be analyzed
      ierr = cusparseSpMV_cmplx(cusparseHandle,TRANS,                  &
     &        C_ONE, matA, vecX, C_ONE, vecY, CUDA_C_64F,             &
     &        CUSPARSE_SPMV_CSR_ALG2, buffer)
      ierr2 = ierr2 + ierr
      ! and get the values back 
      ierr = cusparseDnVecGetValues(vecY, devPtrR)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, *) " error with SpMV operation ", ierr2
          stop
      end if
      ! bnorm = nrm2(b)
      ierr = cublasZnrm2(cublasHandle,n,devPtrRHS,1,bnormPtr)
      ierr2 = ierr2 + ierr
      if (isnan(bnorm)) then
      ! this usually means an inadequate model, in which case Maxwell's fails
          write(6,*) 'Error: b in BICG contains NaNs; exiting...'
          stop
      else if ( bnorm .eq. 0.0) then ! zero rhs -> zero solution
          write(6,*) 'Warning: b in BICG has all zeros, returning zero &
     &        solution'
          x = b 
          KSPiter%niter=1
          KSPiter%failed=.false.
          KSPiter%rerr(1)=0.0
          converged = .true.
          goto 9527 
      endif
      ! R = -Ax 
      ! ierr = cublasZscal(cublasHandle,n,C_MINUSONE,devPtrR,1)
      ! ierr2 = ierr2 + ierr
      ! R = b - Ax
      ! ierr = cublasZaxpy(cublasHandle,n,C_ONE,devPtrRHS,1,devPtrR,1)
      ! ierr2 = ierr2 + ierr
      ! TEST: R = b - Ax with an all-in-one kernel of xpby
      call kernelc_xpbyc(devPtrRHS, C_MINUSONE, devPtrR, n, cuStream)
      ! rnorm = nrm2(b - Ax)
      ierr = cublasZnrm2(cublasHandle,n,devPtrR,1,rnormPtr)
      ierr2 = ierr2 + ierr
      rnorm0 = rnorm
      if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error during residual estimation", ierr2
          stop
      end if
      ! write(6,'(A, ES12.6)') ' initial relative residual = ', rnorm/bnorm
      !================= Now start configuring the iteration ================!
      btol = KSPiter%tol*bnorm
      maxIter = KSPiter%maxIt
      KSPiter%rerr(1) = real(rnorm/bnorm)
      if ( rnorm .le. btol ) then ! the first guess is already good enough
         ! returning
          write(6, *) " The first guess is good enough, exiting..."
          iter=1
          KSPiter%failed=.false.
          converged = .true.
          goto 9527 
      end if
      ! intial the parameters for restarting
      restart = .FALSE.
      interval = 120
      ! write(6,'(A,I4)') ' maxiter = ', maxiter
      ! RHO = C_ONE
      kspLoop: do iter = 1,maxIter
        ! write(6,'(A, I4)') ' KSP iteration #', iter
        if (mod(iter, interval) .eq. 1) then ! hard coded here
            restart = .TRUE.
        end if
        if (restart) then 
            ! restart the iteration (to steepest decend) every interval times
            ! RT = R
            ierr = cublasZcopy(cublasHandle,n,devPtrR,1,devPtrRT,1)
            ierr2 = ierr2 + ierr
            ! current and previous RHO, RHO0 should be one
            RHO1 = C_ONE
            ! RHO = dot(RT,R)
            ierr = cublasZdot(cublasHandle,n,devPtrRT,1,devPtrR,1,rhoPtr)
            ierr2 = ierr2 + ierr
            ! P = R
            ierr = cublasZcopy(cublasHandle,n,devPtrR,1,devPtrP,1)
            ierr2 = ierr2 + ierr
            OMEGA = C_ONE
            restart = .FALSE.
            if (ierr2 .ne. 0 ) then
              write(6, '(A, I2)') " Error steering search direction ", ierr2
              stop
            end if
        else
            ! current and previous RHO, RHO0 should be one
            RHO1 = RHO
            ! RHO = dot(RT,R)
            ierr = cublasZdot(cublasHandle,n,devPtrRT,1,devPtrR,1,rhoPtr)
            ierr2 = ierr2 + ierr
            BETA=(RHO/RHO1)*(ALPHA/OMEGA)
            if (BETA .eq. 0.0) then !bad beta
                converged=.FALSE.
                KSPiter%failed = .true.
                exit kspLoop
            end if
            ! P = R + BETA * (P - OMEGA * V)
            ! P = P - OMEGA * V
            ! ierr = cublasZaxpy(cublasHandle,n,-(OMEGA),devPtrV,1,devPtrP,1)
            ! ierr2 = ierr2 + ierr
            ! P = BETA * P
            ! ierr = cublasZscal(cublasHandle,n,BETA,devPtrP,1)
            ! ierr2 = ierr2 + ierr
            ! P = R + P
            ! ierr = cublasZaxpy(cublasHandle,n,C_ONE,devPtrR,1,devPtrP,1)
            ! ierr2 = ierr2 + ierr
            ! TEST: P = R + BETA * P  with an all-in-one kernel 
            ! call kernelc_xpbyc(devPtrR, BETA, devPtrP, n, cuStream)
            ! TEST: update P with an all-in-one kernel 
            call kernelc_update_pc(devPtrR, devPtrV, BETA, OMEGA, devPtrP, n,&
       &            cuStream)
            if (ierr2 .ne. 0 ) then
              write(6, '(A, I2)') " Error steering search direction ", ierr2
              stop
            end if
        end if
        ! record - start of two SPSVs
        ! ierr = cudaEventRecord(cuEvent1, cuStream)
        ! ierr2 = ierr2 + ierr
        ! ============== first half of the conjugate iteration ============= !
        ! L solve --> L*PT = P
        ! write(6,'(A)') ' Lsolve '
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecX, devPtrP)
        ierr2 = ierr2 + ierr
        ierr = cusparseDnVecSetValues(vecY, devPtrPT)
        ierr2 = ierr2 + ierr
        ierr = cusparseSpSV_solve_dcmplx(cusparseHandle,            &
     &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matL, vecX, vecY, &
     &         CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecY, devPtrPT)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I2)') " Error during Lsolve ", ierr2
          stop
        end if
        ! U solve --> U*PH = PT
        ! write(6,'(A)') ' Usolve '
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecX, devPtrPT)
        ierr2 = ierr2 + ierr
        ierr = cusparseDnVecSetValues(vecY, devPtrPH)
        ierr2 = ierr2 + ierr
        ierr = cusparseSpSV_solve_dcmplx(cusparseHandle,            &
     &         CUSPARSE_OPERATION_NON_TRANSPOSE, C_ONE, matU, vecX, vecY, &
     &         CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, UsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecY, devPtrPH)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I2)') " Error during Usolve ", ierr2
          stop
        end if
        ! record - end of two SPSVs
        ! ierr = cudaEventRecord(cuEvent2, cuStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = cudaDeviceSynchronize()
        ! ierr = cudaEventElapsedTime(timePtr, cuEvent1, cuEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '2 SPSV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !    write(6,'(A,I4)') " error recording the time ", ierr2
        !     stop
        ! end if
        ! record - start of one SPMV
        ! ierr = cudaEventRecord(cuEvent1, cuStream)
        ! ierr2 = ierr2 + ierr
        ! V = A*PH
        ! write(6,'(A)') ' Axpy '
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecX, devPtrPH)
        ierr2 = ierr2 + ierr
        ! V = Diag(iwus)*PH
        call kernelc_hadac(devPtrPH, devPtriwus, devPtrV, n, cuStream)
        ierr = cusparseDnVecSetValues(vecY, devPtrV)
        ierr2 = ierr2 + ierr
        ! V = CC*PH + Diag(iwus)*PH
        ierr = cusparseSpMV_cmplx(cusparseHandle,TRANS, &
     &         C_ONE, matA, vecX, C_ONE, vecY, CUDA_C_64F,  &
     &         CUSPARSE_SPMV_CSR_ALG2, buffer)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecY, devPtrV)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
            write(6,'(A,I4)') " Error calculating A*PH", ierr2
            stop
        end if
        ! record - end of one SPMV
        ! ierr = cudaEventRecord(cuEvent2, cuStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = cudaDeviceSynchronize()
        ! ierr = cudaEventElapsedTime(timePtr, cuEvent1,cuEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '1 SPMV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !    write(6, *) " error recording the time ", ierr2
        !    stop
        ! end if
        ! now check the residual of first half of iteration
        ! rtv = dot(RT,V)
        ierr = cublasZdot(cublasHandle,n,devPtrRT,1,devPtrV,1,rtvPtr)
        ierr2 = ierr2 + ierr
        ALPHA = RHO/rtv
        ! write(6,*) 'alpha = ', ALPHA
        if (ALPHA .eq. 0.0) then !bad alpha
            KSPiter%failed = .true.
            converged=.FALSE.
            exit kspLoop
        end if
        ! x = x + ALPHA * PH
        ierr = cublasZaxpy(cublasHandle,n,ALPHA,devPtrPH,1,devPtrX,1)
        ierr2 = ierr2 + ierr
        ! S = R
        ierr = cudaMemcpyAsync(devPtrS, devPtrR, Arow_d_size, &
       &       cudaMemcpyDeviceToDevice)
        ierr2 = ierr2 + ierr
        ! S = R - ALPHA * V
        ierr = cublasZaxpy(cublasHandle,n,-(ALPHA),devPtrV,1,devPtrS,1)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error estimating the first Conj res. ", ierr2
          stop
        end if
        ! ============== second half of the conjugate iteration ============= !
        ! record - start of two SPSV
        ! ierr = cudaEventRecord(cuEvent1, cuStream)
        ! ierr2 = ierr2 + ierr
        ! L solve --> L*ST = S
        ! write(6,'(A)') ' Lsolve '
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecX, devPtrS)
        ierr2 = ierr2 + ierr
        ierr = cusparseDnVecSetValues(vecY, devPtrST)
        ierr2 = ierr2 + ierr
    !    ierr = cusparseSpSV_analysis_dcmplx(cusparseHandle,        &
    ! &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matL, vecX, vecY, &
    ! &         CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle, bufferL)
    !    ierr2 = ierr2 + ierr
        ierr = cusparseSpSV_solve_dcmplx(cusparseHandle,           &
     &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matL, vecX, vecY, &
     &         CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecY, devPtrST)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I2)') " Error during Lsolve ", ierr2
          stop
        end if
        
        ! U solve --> U*SH = ST
        ! write(6,'(A)') ' Usolve '
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecX, devPtrST)
        ierr2 = ierr2 + ierr
        ierr = cusparseDnVecSetValues(vecY, devPtrSH)
        ierr2 = ierr2 + ierr
    !    ierr = cusparseSpSV_analysis_dcmplx(cusparseHandle,         &
    ! &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matU, vecX, vecY, &
    ! &         CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, UsolveHandle, bufferU)
    !    ierr2 = ierr2 + ierr
        ierr = cusparseSpSV_solve_dcmplx(cusparseHandle,            &
     &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matU, vecX, vecY, &
     &         CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, UsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecY, devPtrSH)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I2)') " Error during Usolve ", ierr2
          stop
        end if

        ! record - end of 2 SPSVs
        ! ierr = cudaEventRecord(cuEvent2, cuStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = cudaDeviceSynchronize()
        ! ierr = cudaEventElapsedTime(timePtr, cuEvent1,cuEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '2 SPSV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !     write(6, *) " error recording the time ", ierr2
        !     stop
        ! end if
        ! record - start of 1 SPMV
        ! ierr = cudaEventRecord(cuEvent1, cuStream)
        ! ierr2 = ierr2 + ierr
        ! T = A*SH 
        ! write(6,'(A)') ' Axpy '
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecX, devPtrSH)
        ierr2 = ierr2 + ierr
        ! T = Diag(iwus)*SH
        call kernelc_hadac(devPtrSH, devPtriwus, devPtrT, n, cuStream)
        ierr = cusparseDnVecSetValues(vecY, devPtrT)
        ierr2 = ierr2 + ierr
        ! T = A*SH + Diag(iwus)*SH
        ierr = cusparseSpMV_cmplx(cusparseHandle,TRANS, &
       &        C_ONE, matA, vecX, C_ONE, vecY, CUDA_C_64F,                 & 
       &        CUSPARSE_SPMV_CSR_ALG2, buffer)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecY, devPtrT)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
            write(6,'(A,I4)') " Error calculating A*SH", ierr2
            stop
        end if

        ! record - end of 1 SPMV 
        ! ierr = cudaEventRecord(cuEvent2, cuStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = cudaDeviceSynchronize()
        ! ierr = cudaEventElapsedTime(timePtr, cuEvent1,cuEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '2 SPMV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !     write(6, *) " error recording the time ", ierr2
        !     stop
        ! end if

        ! now check the second half of iteration
        ! calculate the residual norm for the second half of iteration 
        ! tt = dot(T,T)
        ierr = cublasZdot(cublasHandle,n,devPtrT,1,devPtrT,1,ttPtr)
        ierr2 = ierr2 + ierr
        ! rtv = dot(T,S)
        ierr = cublasZdot(cublasHandle,n,devPtrT,1,devPtrS,1,rtvPtr)
        ierr2 = ierr2 + ierr
        OMEGA = rtv/tt
        ! write(6,*) 'omega = ', OMEGA
        if (OMEGA .eq. 0.0) then !bad omega
            KSPiter%failed = .true.
            converged=.false.
            exit kspLoop
        end if
        ! x = x + OMEGA * SH
        ierr = cublasZaxpy(cublasHandle,n,OMEGA,devPtrSH,1,devPtrX,1)
        ierr2 = ierr2 + ierr
        ! R = S
        ierr = cudaMemcpyAsync(devPtrR, devPtrS, Arow_d_size, &
       &       cudaMemcpyDeviceToDevice)
        ierr2 = ierr2 + ierr
        ! R = S - OMEGA * T
        ierr = cublasZaxpy(cublasHandle,n,-(OMEGA),devPtrT,1,devPtrR,1)
        ierr2 = ierr2 + ierr
        ierr = cudaDeviceSynchronize()
        ierr2 = ierr2 + ierr
        ! early second half convergence check (norm of R)
        ! rnorm = norm(R)
        ierr = cublasZnrm2(cublasHandle,n,devPtrR,1,rnormPtr)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error estimating the second Conj res. ", ierr2
          stop
        end if
        KSPiter%rerr(iter) = rnorm/bnorm
        ! write(6,'(A12,I4,A10,ES12.6)') 'iteration #', iter, ' relres= ', &
        !&           KSPiter%rerr(iter)
        ! now check the second half of iteration
        ! NOTE WE TEST AN IDEA THAT OMIT THE X, WHICH SAVES ANOTHER 
        ! SPMV FOR US
        if (rnorm.lt.btol) then
            ! double check if the residual is really less than tol
            ! R = A*X
            ! still need to use the vecX/Y types
            ierr = cusparseDnVecSetValues(vecX, devPtrX)
            ierr2 = ierr2 + ierr
            ! R = Diag(iwus)*X
            call kernelc_hadac(devPtrX, devPtriwus, devPtrR, n, cuStream)
            ierr = cusparseDnVecSetValues(vecY, devPtrR)
            ierr2 = ierr2 + ierr
            ! R = CC*X + Diag(iwus)*X
            ierr = cusparseSpMV_cmplx(cusparseHandle,TRANS, &
     &          C_ONE, matA, vecX, C_ONE, vecY, CUDA_C_64F,  &
     &          CUSPARSE_SPMV_CSR_ALG2, buffer)
            ierr2 = ierr2 + ierr
            ! and get the values back
            ierr = cusparseDnVecGetValues(vecY, devPtrR)
            ierr2 = ierr2 + ierr
            if (ierr2 .ne. 0 ) then
                write(6,'(A,I4)') " Error calculating A*X", ierr2
                stop
            end if
            ! R = -Ax
            ! ierr = cublasZscal(cublasHandle,n,C_MINUSONE,devPtrR,1)
            ! ierr2 = ierr2 + ierr
            ! R = b - Ax
            ! ierr = cublasZaxpy(cublasHandle,n,C_ONE,devPtrRHS,1,devPtrR,1)
            ! ierr2 = ierr2 + ierr
            ! TEST: R = b - Ax with an all-in-one kernel of xpby
            call kernelc_xpbyc(devPtrRHS, C_MINUSONE, devPtrR, n, cuStream)
            ! rnorm = norm(R)
            ierr = cublasZnrm2(cublasHandle,n,devPtrR,1,rnormPtr)
            ierr2 = ierr2 + ierr
            if (rnorm.lt.btol) then
                KSPiter%rerr(iter) = real(rnorm/bnorm)
                converged=.TRUE.
                KSPiter%failed = .false.
                KSPiter%niter = iter
                exit kspLoop
            end if
        end if
        if (rnorm .lt. rnorm0) then
            rnorm0 = rnorm
            ierr = cudaMemCpyAsync(devPtrX0, devPtrX, Arow_d_size, &
                cudaMemcpyDeviceToDevice)
            ierr2 = ierr2 + ierr
        end if
      end do kspLoop
 9527 continue
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error synchronizing after iterations ", ierr2
          stop
      end if
      ! write(6,'(A)') ' Copy solution from GPU to CPU'
      if (.not. converged) then ! solution not found 
          KSPiter%niter=KSPiter%maxit
          ! return the best solution so far
          ierr = cudaMemcpy(xPtr,devPtrX0,Arow_d_size,cudaMemcpyDeviceToHost)
      else ! solution found
          KSPiter%niter=iter
          ! return the last solution
          ierr = cudaMemcpy(xPtr,devPtrX,Arow_d_size,cudaMemcpyDeviceToHost)
      end if
      if (ierr .ne. 0 ) then
          write(6, '(A, I2)') " cudaMemcpy back to host error: ", ierr
          stop
      end if
      ! \activiate lightsaber
      ! clear gpu mem
      ierr = cudaFree(devPtrArow)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrAcol)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrAval)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtriwus)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrLrow)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrLcol)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrLval)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrUrow)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrUcol)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrUval)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrX)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrX0)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrRHS)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrR)
      ierr2 = ierr2 + ierr 
      ierr = cudaFree(devPtrRT)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrP)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrPT)
      ierr2 = ierr2 + ierr 
      ierr = cudaFree(devPtrPH)
      ierr2 = ierr2 + ierr 
      ierr = cudaFree(devPtrS)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrST)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrSH)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrT)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrV)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrAX)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(buffer)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(bufferL)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(bufferU)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error during cudafree: ',ierr2
          stop
      end if 
      ierr = cusparseSpSV_destroyDescr(LsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_destroyDescr(UsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroyDnVec(vecX)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroyDnVec(vecY)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroySpMat(matL)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroySpMat(matU)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroySpMat(matA)
      ierr2 = ierr2 + ierr
      ierr = cublasDestroy(cublasHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroy(cusparseHandle)
      ierr2 = ierr2 + ierr
      ierr = cudaStreamDestroy(cuStream)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error during cuda handle destruction: ',ierr2
          stop
      end if 
      ierr = cf_resetFlag(device_idx)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error setting device flags: ',ierr2
          stop
      end if 
      return
end subroutine cuBiCG ! cuBiCG

! *****************************************************************************
subroutine cuBiCGmix(b,x,KSPiter,device_idx,adjt)
  ! Stablized version of BiConjugate Gradient, set up for solving
  ! A x = b 
  ! solves for the interior (edge) field
  !
  ! modified (no hell no) to call the CUDA lib to calculate with GPU 
  ! I kept the naming convention from my BiCGStab (see above CPU version)
  ! for instance X -> devPtrX, T -> devPtrT to manipulate the GPU memory
  ! also added the optional adjoint to solve adjoint system A^Tx = b
  ! 
  ! NOTE: BICG actually performs two sub (or half) line searches within a 
  !       iteration, but here we only store the relerr for the second sub
  !       just to be compatitive with QMR
  ! 
  ! interface...........
  ! redefining some of the interfaces for our convenience (locally)
  ! generic routines for vector operations for edge/ face nodes
  ! in a staggered grid
  !
  ! NOTE: this has not been extensively tested - I believe it feels a 
  ! little unstable (dispite the name)...
  ! if you have time reading this, test it!

      use modeloperator3d, only:  AAii, L, U, LH, UH,vOmegaMuSig
      use spoptools
      implicit none
      !  b is right hand side
      complex (kind=prec),intent(in),target,dimension(:)    :: b
      !  solution vector is x ... on input is provided with the initial
      !  guess, on output is the iterate with smallest residual. 
      complex (kind=prec),intent(inout),target,dimension(:) :: x
      type (solverControl_t),intent(inout)                  :: KSPiter
      integer,intent(in)                                    :: device_idx
      logical,intent(in),optional                           :: adjt
    
      ! local variables
      integer                                    :: n, nnz, iter, maxIter
      integer                                    :: nnzL, nnzU
      complex (kind=prec),pointer,dimension(:)   :: iwus

      real    (kind=prec), target                :: rnorm, bnorm, rnorm0, btol
      real    (kind=prec), target                :: xnorm
      complex (kind=prec)                        :: RHO1, ALPHA, BETA, OMEGA
      complex (kind=prec), target                :: RTV,TT,RHO
      real    (c_float),target                   :: ctime 
      integer                                    :: TRANS
      integer (c_int)                            :: Lfillmode,Ldiagtype
      integer (c_int)                            :: Ufillmode,Udiagtype
      integer                                    :: ierr,ierr2, idx
      logical                                    :: converged, adjoint
      integer                                    :: interval
      logical                                    :: restart

      ! size to determine corresponding buffer/memory
      integer*8                :: Arowp1_i_size,Arow_d_size
      integer*8                :: Acol_d_size,Annz_i_size,Annz_d_size
      integer*8                :: Mrow_d_size,Lnnz_i_size,Lnnz_d_size
      integer*8                :: Unnz_i_size,Unnz_d_size
      integer(c_int),target    :: zeroLoc
      ! integer(c_int)           :: mbsize
      integer(c_size_t)        :: bsize, lbsize, ubsize

      ! ------------------- pointers to *host* memory ------------------- !   
      type(c_ptr) :: cublasHandle 
      type(c_ptr) :: cusparseHandle
      type(c_ptr) :: cuStream
      type(c_ptr) :: cuEvent1
      type(c_ptr) :: cuEvent2
      type(c_ptr) :: timePtr
      type(c_ptr) :: matA
      type(c_ptr) :: matL
      type(c_ptr) :: matU
      type(c_ptr) :: vecX
      type(c_ptr) :: vecY 
      type(c_ptr) :: vecXs
      type(c_ptr) :: vecYs
      type(c_ptr) :: LsolveHandle
      type(c_ptr) :: UsolveHandle
      type(c_ptr) :: ilu_info
      type(c_ptr) :: ArowPtr
      type(c_ptr) :: AcolPtr
      type(c_ptr) :: AvalPtr
      type(c_ptr) :: iwusPtr
      type(c_ptr) :: MrowPtr
      type(c_ptr) :: McolPtr
      type(c_ptr) :: MvalPtr
      type(c_ptr) :: xPtr  
      type(c_ptr) :: bPtr
      type(c_ptr) :: bnormPtr
      type(c_ptr) :: rnormPtr
      type(c_ptr) :: rnormPtr0
      type(c_ptr) :: xnormPtr ! for stagnant check
      type(c_ptr) :: rhoPtr
      type(c_ptr) :: rtvPtr
      type(c_ptr) :: ttPtr
      type(c_ptr) :: zeroLocPtr
      ! -------------------- pointers to *device* memory ----------------- ! 
      type(c_ptr) :: devPtrArow
      type(c_ptr) :: devPtrMval
      type(c_ptr) :: devPtrAcol
      type(c_ptr) :: devPtrAval
      type(c_ptr) :: devPtriwus
      type(c_ptr) :: devPtrLrow
      type(c_ptr) :: devPtrLcol
      type(c_ptr) :: devPtrLval
      type(c_ptr) :: devPtrUrow
      type(c_ptr) :: devPtrUcol
      type(c_ptr) :: devPtrUval
      type(c_ptr) :: devPtrX
      type(c_ptr) :: devPtrX0
      type(c_ptr) :: devPtrRHS
      type(c_ptr) :: devPtrR
      type(c_ptr) :: devPtrRT
      type(c_ptr) :: devPtrP
      type(c_ptr) :: devPtrPT
      type(c_ptr) :: devPtrPH
      type(c_ptr) :: devPtrS
      type(c_ptr) :: devPtrST
      type(c_ptr) :: devPtrSH
      type(c_ptr) :: devPtrV
      type(c_ptr) :: devPtrT
      type(c_ptr) :: devPtrAX
      type(c_ptr) :: devPtr32 ! single temp devptr
      type(c_ptr) :: buffer
      type(c_ptr) :: bufferL
      type(c_ptr) :: bufferU
      
      ! note we don't check the input parameters compatibility 
      ierr2 = 0
      converged = .FALSE.
      zeroLoc = 0
      if (present(adjt)) then 
          ! write(6,'(A)') ' adjt = ', adjt
          adjoint = adjt
          if (adjt) then
              TRANS = CUSPARSE_OPERATION_TRANSPOSE
              Lfillmode = CUSPARSE_FILL_MODE_UPPER
              Ufillmode = CUSPARSE_FILL_MODE_LOWER
              Ldiagtype = CUSPARSE_DIAG_TYPE_NON_UNIT
              Udiagtype = CUSPARSE_DIAG_TYPE_UNIT
          else
              TRANS = CUSPARSE_OPERATION_NON_TRANSPOSE
              Lfillmode = CUSPARSE_FILL_MODE_LOWER
              Ufillmode = CUSPARSE_FILL_MODE_UPPER
              Ldiagtype = CUSPARSE_DIAG_TYPE_UNIT
              Udiagtype = CUSPARSE_DIAG_TYPE_NON_UNIT
          end if
      else
          adjoint = .FALSE.
          TRANS = CUSPARSE_OPERATION_NON_TRANSPOSE
          Lfillmode = CUSPARSE_FILL_MODE_LOWER
          Ufillmode = CUSPARSE_FILL_MODE_UPPER
          Ldiagtype = CUSPARSE_DIAG_TYPE_UNIT
          Udiagtype = CUSPARSE_DIAG_TYPE_NON_UNIT
      end if
      
      ! firstly need to translate the ModEM SP datatypes to simple
      ! CSR plain vectors
      ! note we still keep the Ae = CCGDe + iwuse 
      ! idea as in ModEM SP2, we will test if this will save some time 
      ! (probably not much)
      n = AAii%nrow
      nnz = size(AAii%col)
      ! note here A is real now
      ArowPtr = c_loc(AAii%row)
      AcolPtr = c_loc(AAii%col)
      AvalPtr = c_loc(AAii%val)
      ! other host pointers
      xPtr  = c_loc(x)  ! x = A \ b
      bPtr  = c_loc(b)  ! b = ones(n,1)
      allocate(iwus(n))
      iwus = VOmegaMuSig*ISIGN*CMPLX(0.0,1.0,8)
      iwusPtr = c_loc(iwus) ! i V omega mu sigma
      bnormPtr = c_loc(bnorm)
      rnormPtr = c_loc(rnorm)
      rnormPtr0 = c_loc(rnorm0)
      xnormPtr = c_loc(xnorm)
      rhoPtr = c_loc(RHO)
      rtvPtr = c_loc(rtv)
      ttPtr = c_loc(tt)
      zeroLocPtr = c_loc(zeroLoc)
      ! pointer to store the time 
      ctime = 0.0
      timePtr = c_loc(ctime)
      
      ! get the size of the vector and matrix (need to allocate on the device)
      Arow_d_size=sizeof(b(1:n))
      Acol_d_size=sizeof(x(1:n)) ! this is useful if A is not square
      Arowp1_i_size=sizeof(AAii%row(1:n+1))
      Annz_i_size=sizeof(AAii%col(1:nnz))
      Annz_d_size=sizeof(AAii%val(1:nnz))
      Mrow_d_size=Arow_d_size/2

      ! select the current cuda device 
      ! note this only works for physical devices (not working for MIG devices)
      ierr = cudaSetDevice(device_idx);
      ierr2 = ierr2 + ierr
      ! firstly define the CUDA Stream and cuda handles
      ierr = cudaStreamCreateWithFlags(cuStream, cudaStreamNonBlocking) 
      ierr2 = ierr2 + ierr
      ! initialize the cusparse lib
      ierr = cusparseCreate(cusparseHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseSetStream(cusparseHandle,cuStream) 
      ierr2 = ierr2 + ierr
      ! now initialize the cublas lib
      ierr = cublasCreate(cublasHandle)
      ierr2 = ierr2 + ierr
      ierr = cublasSetStream(cublasHandle,cuStream)
      ierr2 = ierr2 + ierr
      ! and creates two cuda events to record the time
      ierr = cudaEventCreate(cuEvent1)
      ierr2 = ierr2 + ierr
      ierr = cudaEventCreate(cuEvent2)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error during cuda initialize ',ierr2
          stop
      end if 
      ! record, start allocating the device mem
      ! ierr = cudaEventRecord(cuEvent1, cuStream)
      ! ierr2 = ierr2 + ierr
      ! write(6,*) 'Allocating GPU memory'
      ierr = cudaMalloc(devPtrX,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrX0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrRHS,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrR,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrRT,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrP,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrPT,Mrow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrPH,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrS,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrST,Mrow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrSH,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrT,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrV,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrAX,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrAval,Annz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrAcol,Annz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrArow,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrMval,Annz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtriwus,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtr32,Mrow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error during CUDA allocation: ',ierr2
          stop
      end if 
      ! write(6,*) 'reset GPU memory to all zeros'
      ierr = cudaMemset(devPtrX,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrX0,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrRHS,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrR,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrRT,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrP,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrPT,0,Mrow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrPH,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrS,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrST,0,Mrow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrSH,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrT,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrV,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrAX,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrAval,0,Annz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrAcol,0,Annz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrArow,0,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrMval,0,Annz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtriwus,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtr32,0,Mrow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I3)') 'Error during device memory reseting : ',ierr2
          stop
      end if 
      ! now we need to transfer memory over to GPU
      ! initialize rhs and x
      ierr = cudaMemcpyAsync(devPtrX,xPtr,Arow_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrX0,devPtrX,Arow_d_size, &
     &        cudaMemcpyDeviceToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrRHS,bPtr,Arow_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error during cuda memcpy ", ierr2
          stop
      end if
      ! initialize A
      ! note we should not do this before L/U as we use devPtrAval when
      ! converting L/U to single
      ierr = cudaMemcpyAsync(devPtrArow,ArowPtr,Arowp1_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrAcol,AcolPtr,Annz_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrAval,AvalPtr,Annz_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtriwus,iwusPtr,Arow_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ! finally need to create the SpMatDescr (now as real)
      ierr = cusparseCreateCsr(matA, n, n, nnz, devPtrArow, devPtrAcol, &
     &       devPtrAval, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, &
     &       CUSPARSE_INDEX_BASE_ONE, CUDA_R_64F) 
      ierr2 = ierr2 + ierr
      ! finish synchronizing
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error assembling system matrix ", ierr2
          stop
      end if
      ! now deallocate temp array
      if ( c_associated(iwusPtr) ) then
          ! this is a little tricky as iwus is associated with
          ! iwusPtr in C
          nullify(iwus)
          call cf_free(iwusPtr)
      end if
      ! now start messing up with SpMV
      ierr = cusparseCreateDnVec(vecX, n, devPtrX, CUDA_C_64F)
      ierr2 = ierr2 + ierr
      ierr = cusparseCreateDnVec(vecY, n, devPtrR, CUDA_C_64F)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, *) " error creating the dense Vecs "
          stop
      end if
      ! now estimate the buffersize needed by SpMV
      ! y = a*A*x + b*y --> y = Ax
      ierr = cusparseSpMV_bufferSize_cmplx(cusparseHandle,              &
     &        TRANS, C_ONE, matA, vecX, C_ONE, vecY,                   &
     &        CUDA_C_64F, CUSPARSE_SPMV_CSR_ALG2, bsize)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error estimating buffer for spMV ",  ierr2
          stop
      end if
      ! write(6,*) 'spmv buffersize = ', mbsize
      ! and finally (re)allocate the buffer
      ierr = cudaMalloc(buffer, bsize)
      ierr2 = ierr2 + ierr
      ! now we need to deal with L and U
      ! which are already decomposed by CPU
      ! write(6,'(A)') ' Setup L and U preconditioners on GPU'
      ! L first
      if (adjoint) then
          nnzL = size(LH%col)
          ! sizes (row size is the same as A
          Lnnz_i_size=sizeof(LH%col(1:nnzL))
          MrowPtr = c_loc(LH%row)
          McolPtr = c_loc(LH%col)
          MvalPtr = c_loc(LH%val)
          Lnnz_d_size=sizeof(LH%val)/2
      else
          nnzL = size(L%col)
          ! sizes (row size is the same as A
          Lnnz_i_size=sizeof(L%col(1:nnzL))
          MrowPtr = c_loc(L%row)
          McolPtr = c_loc(L%col)
          MvalPtr = c_loc(L%val)
          Lnnz_d_size=sizeof(L%val)/2
      end if
      ! initialize L in GPU memory
      ierr = cudaMalloc(devPtrLval,Lnnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrLcol,Lnnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrLrow,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! set everything to zeros
      ierr = cudaMemset(devPtrLval,0,Lnnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrLcol,0,Lnnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrLrow,0,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! now copy the values to device
      ierr = cudaMemcpyAsync(devPtrMval,MvalPtr,Lnnz_d_size*2, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      call kernelc_d2s(devPtrMval, devPtrLval, nnzL)
      ierr = cudaMemcpyAsync(devPtrLcol,McolPtr,Lnnz_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrLrow,MrowPtr,Arowp1_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ! finish synchronizing
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      ! finally need to create the SpMatDescr 
      ierr = cusparseCreateCsr(matL,n,n,nnzL,devPtrLrow,devPtrLcol,&
     &       devPtrLval, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, &
     &       CUSPARSE_INDEX_BASE_ONE, CUDA_C_32F) 
      ierr2 = ierr2 + ierr
      ierr = cusparseSpMatSetAttribute(matL, CUSPARSE_SPMAT_FILL_MODE, &
    &        Lfillmode, 4)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpMatSetAttribute(matL, CUSPARSE_SPMAT_DIAG_TYPE, &
    &        Ldiagtype, 4)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error assembling L matrix ", ierr2
          stop
      end if
      ! now start messing up with SpSV
      ! firstly we need to create two dense vectors 
      ! as the user-hostile developers in Nvidia think of a new idea 
      ! to mess up the interfaces
      ! both should be zero at current stage
      ierr = cusparseCreateDnVec(vecXs, n, devPtrPT, CUDA_C_32F)
      ierr2 = ierr2 + ierr
      ierr = cusparseCreateDnVec(vecYs, n, devPtrST, CUDA_C_32F)
      ierr2 = ierr2 + ierr
      ! now estimate the buffersize needed by SpSV (Lsolve)
      ! solves y in L*y = a*x (if a=1)
      ! still need to establish a context handler for Lsolve
      ierr = cusparseSpSV_createDescr(LsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_bufferSize_cmplx(cusparseHandle,             &
     &        CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONEs,matL,vecXs,vecYs,&
     &        CUDA_C_32F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle, lbsize)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error estimating buffer for Lsolve ",  ierr2
          stop
      end if
      ierr = cudaMalloc(bufferL, lbsize)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_analysis_cmplx(cusparseHandle,         &
    &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONEs,matL,vecXs,vecYs,&
    &         CUDA_C_32F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle, bufferL)
      ierr2 = ierr2 + ierr
      ! write(6,*) 'spsv Lsolve buffersize = ', lbsize
      ! then U 
      if (adjoint) then
          nnzU = size(UH%col)
          ! sizes (row size is the same as A
          Unnz_i_size=sizeof(UH%col(1:nnzU))
          MrowPtr = c_loc(UH%row)
          McolPtr = c_loc(UH%col)
          MvalPtr = c_loc(UH%val)
          Unnz_d_size=sizeof(UH%val)/2
      else
          nnzU = size(U%col)
          ! sizes (row size is the same as A
          Unnz_i_size=sizeof(U%col(1:nnzU))
          MrowPtr = c_loc(U%row)
          McolPtr = c_loc(U%col)
          MvalPtr = c_loc(U%val)
          Unnz_d_size=sizeof(U%val)/2
      end if
      ! initialize U in GPU memory
      ierr = cudaMalloc(devPtrUval,Unnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrUcol,Unnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrUrow,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! set everything to zeros
      ierr = cudaMemset(devPtrUval,0,Unnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrUcol,0,Unnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrUrow,0,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! now copy the values to device
      ierr = cudaMemcpyAsync(devPtrMval,MvalPtr,Unnz_d_size*2, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      call kernelc_d2s(devPtrMval, devPtrUval, nnzU)
      ierr = cudaMemcpyAsync(devPtrUcol,McolPtr,Unnz_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrUrow,MrowPtr,Arowp1_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ! now deallocate the temp val
      ierr = cudaFree(devPtrMval)
      ierr2 = ierr2 + ierr
      ! finish synchronizing
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      ! finally need to create the SpMatDescr 
      ierr = cusparseCreateCsr(matU,n,n,nnzU,devPtrUrow,devPtrUcol,&
     &       devPtrUval, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, &
     &       CUSPARSE_INDEX_BASE_ONE, CUDA_C_32F) 
      ierr2 = ierr2 + ierr
      ierr = cusparseSpMatSetAttribute(matU, CUSPARSE_SPMAT_FILL_MODE, &
    &        Ufillmode, 4)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpMatSetAttribute(matU, CUSPARSE_SPMAT_DIAG_TYPE, &
    &        Udiagtype, 4)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error assembling U matrix ", ierr2
          stop
      end if
      ! now estimate the buffersize needed by SpSV (Usolve)
      ! solves y in U*y = a*x (if a=1)
      ! still need to establish a context handler for Usolve
      ierr = cusparseSpSV_createDescr(UsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_bufferSize_cmplx(cusparseHandle,             &
     &        CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONEs,matU,vecXs,vecYs,&
     &        CUDA_C_32F,CUSPARSE_SPSV_ALG_DEFAULT, UsolveHandle, ubsize)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error estimating buffer for Usolve ",  ierr2
          stop
      end if
      ! write(6,*) 'spsv Usolve buffersize = ', ubsize
      ierr = cudaMalloc(bufferU, ubsize)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_analysis_cmplx(cusparseHandle,         &
    &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONEs,matU,vecXs,vecYs,&
    &         CUDA_C_32F,CUSPARSE_SPSV_ALG_DEFAULT, UsolveHandle, bufferU)
      ierr2 = ierr2 + ierr
      !write(6, '(A, I8)') " allocated SpMV/SpSV buffer on GPU (MB) ",&
     !&            bsize/1024/1024
      ! see how long it takes for the intialization part
      ! ierr = cudaEventRecord(cuEvent2, cuStream)
      ! ierr2 = ierr2 + ierr
      ! ierr = cudaDeviceSynchronize()
      ! ierr = cudaEventElapsedTime(timePtr, cuEvent1,cuEvent2)
      ! ierr2 = ierr2 + ierr
      ! write(6,*) 'initial GPU memory cpy time is ', ctime, 'ms'
      ! if (ierr2 .ne. 0 ) then
      !     write(6, *) " error recording the time ", ierr2
      !     stop
      ! end if
      ! now Compute the initial residual
      ! write(6,*) 'Computing initial residual'
      ! need to deal with the dense vecs
      ierr = cusparseDnVecSetValues(vecX, devPtrX)
      ierr2 = ierr2 + ierr
      ! R = Diag(iwus)*X
      call kernelc_hadac(devPtrX, devPtriwus, devPtrR, n, cuStream)
      ierr = cusparseDnVecSetValues(vecY, devPtrR)
      ierr2 = ierr2 + ierr
      ! now calculate the y = Ax
      ! thank God (or whatever deserves it) that SpMV does not
      ! need to be analyzed
      ! R = CC*X + Diag(iwus)*X
      ierr = cusparseSpMV_cmplx(cusparseHandle,TRANS,                  &
     &        C_ONE, matA, vecX, C_ONE, vecY, CUDA_C_64F,             &
     &        CUSPARSE_SPMV_CSR_ALG2, buffer)
      ierr2 = ierr2 + ierr
      ! and get the values back
      ierr = cusparseDnVecGetValues(vecY, devPtrR)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, *) " error with SpMV operation ", ierr2
          stop
      end if
      ! bnorm = nrm2(b)
      ierr = cublasZnrm2(cublasHandle,n,devPtrRHS,1,bnormPtr)
      ierr2 = ierr2 + ierr
      if (isnan(bnorm)) then
      ! this usually means an inadequate model, in which case Maxwell's fails
          write(0,*) 'Error: b in BICG contains NaNs; exiting...'
          stop
      else if ( bnorm .eq. 0.0) then ! zero rhs -> zero solution
          write(0,*) 'Warning: b in BICG has all zeros, returning zero solution'
          x = b 
          KSPiter%niter=1
          KSPiter%failed=.false.
          KSPiter%rerr(1)=0.0
          converged = .true.
          goto 9527 
      endif
      ! R = -Ax 
      ! ierr = cublasZscal(cublasHandle,n,C_MINUSONE,devPtrR,1)
      ! ierr2 = ierr2 + ierr
      ! R = b - Ax
      ! ierr = cublasZaxpy(cublasHandle,n,C_ONE,devPtrRHS,1,devPtrR,1)
      ! ierr2 = ierr2 + ierr
      ! TEST: R = b - Ax with an all-in-one kernel of xpby
      call kernelc_xpbyc(devPtrRHS, C_MINUSONE, devPtrR, n, cuStream)
      ! rnorm = nrm2(b - Ax)
      ierr = cublasZnrm2(cublasHandle,n,devPtrR,1,rnormPtr)
      ierr2 = ierr2 + ierr
      rnorm0 = rnorm
      if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error during residual estimation", ierr2
          stop
      end if
      ! write(6,'(A, ES10.4)') ' initial relative residual = ', rnorm/bnorm
      !================= Now start configuring the iteration ================!
      btol = KSPiter%tol*bnorm
      maxIter = KSPiter%maxIt
      KSPiter%rerr(1) = real(rnorm/bnorm)
      if ( rnorm .le. btol ) then ! the first guess is already good enough
         ! returning
          write(6, *) " The first guess is good enough, exiting..."
          iter=1
          KSPiter%failed=.false.
          converged = .true.
          goto 9527 
      end if
      RHO = C_ONE
      ! parameters for stagnant check
      restart = .TRUE.
      interval = 120 ! hard coded here
      ! write(6,'(A,I4)') ' maxiter = ', maxiter
      kspLoop: do iter = 1,maxIter
        ! write(6,'(A, I4)') ' KSP iteration #', iter
        if (mod(iter, interval) .eq. 1) then ! hard coded here
            restart = .TRUE.
        end if
        if (restart) then
            ! restart the iteration (to steepest decend) every interval times
            ! RT = R 
            ierr = cublasZcopy(cublasHandle,n,devPtrR,1,devPtrRT,1)
            ierr2 = ierr2 + ierr
            ! current and previous RHO, RHO0 should be one
            RHO1 = C_ONE
            ! RHO = dot(RT,R)
            ierr = cublasZdot(cublasHandle,n,devPtrRT,1,devPtrR,1,rhoPtr)
            ierr2 = ierr2 + ierr
            ! P = R
            ierr = cublasZcopy(cublasHandle,n,devPtrR,1,devPtrP,1)
            ierr2 = ierr2 + ierr
            OMEGA = C_ONE
            restart = .FALSE.
            if (ierr2 .ne. 0 ) then
              write(6, '(A, I2)') " Error steering search direction ", ierr2
              stop
            end if
        else
            ! current and previous RHO, RHO0 should be one
            RHO1=RHO
            ! RHO = dot(RT,R)
            ierr = cublasZdot(cublasHandle,n,devPtrRT,1,devPtrR,1,rhoPtr)
            ierr2 = ierr2 + ierr
            BETA=(RHO/RHO1)*(ALPHA/OMEGA)
            if (BETA .eq. 0.0) then !bad beta
                converged=.FALSE.
                KSPiter%failed = .true.
                exit kspLoop
            end if
            ! P = R + BETA * (P - OMEGA * V)
            ! ierr = cublasZaxpy(cublasHandle,n,-(OMEGA),devPtrV,1,devPtrP,1)
            ! ierr2 = ierr2 + ierr
            ! ierr = cublasZscal(cublasHandle,n,BETA,devPtrP,1)
            ! ierr2 = ierr2 + ierr
            ! ierr = cublasZaxpy(cublasHandle,n,C_ONE,devPtrR,1,devPtrP,1)
            ! ierr2 = ierr2 + ierr
            ! TEST: P = R + BETA * P with an all-in-one kernel of xpby
            ! call kernelc_xpbyc(devPtrR, BETA, devPtrP, n, cuStream)
            ! TEST: update P with an all-in-one kernel 
            call kernelc_update_pc(devPtrR, devPtrV, BETA, OMEGA, devPtrP, n, &
      &         cuStream)
            if (ierr2 .ne. 0 ) then
              write(6, '(A, I2)') " Error steering search direction ", ierr2
              stop
            end if
        end if
        ! ============== first half of the conjugate iteration ============= !
        ! L solve --> L*PT = P
        ! write(6,'(A)') ' Lsolve '
        ! record -start of SPSVs
        ! ierr = cudaEventRecord(cuEvent1, cuStream)
        ! ierr2 = ierr2 + ierr
        ! firstly need to convert P to FP32
        call kernelc_d2s(devPtrP, devPtr32, n)
        ierr = cusparseDnVecSetValues(vecXs, devPtr32)
        ierr2 = ierr2 + ierr
        ierr = cusparseDnVecSetValues(vecYs, devPtrPT)
            restart = .FALSE.
        ierr2 = ierr2 + ierr
        ierr = cusparseSpSV_solve_cmplx(cusparseHandle,            &
     &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONEs,matL,vecXs,vecYs,&
     &         CUDA_C_32F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecYs, devPtrPT)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I2)') " Error during Lsolve ", ierr2
          stop
        end if
        ! U solve --> U*PH = PT
        ! write(6,'(A)') ' Usolve '
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecXs, devPtrPT)
        ierr2 = ierr2 + ierr
        ierr = cusparseDnVecSetValues(vecYs, devPtr32)
        ierr2 = ierr2 + ierr
        ierr = cusparseSpSV_solve_cmplx(cusparseHandle,            &
     &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONEs,matU,vecXs,vecYs,&
     &         CUDA_C_32F,CUSPARSE_SPSV_ALG_DEFAULT, UsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecYs, devPtr32)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I2)') " Error during Usolve ", ierr2
          stop
        end if
        ! still need to convert PH to FP64
        call kernelc_s2d(devPtr32, devPtrPH, n)
        ! record - end of SPSVs
        ! ierr = cudaEventRecord(cuEvent2, cuStream)
        ! ierr = cudaDeviceSynchronize()
        ! ierr2 = ierr2 + ierr
        ! ierr = cudaEventElapsedTime(timePtr, cuEvent1,cuEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '2 SPSV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !     write(6,'(A,I4)') " error recording the time ", ierr2
        !     stop
        ! end if
        ! record -start of SPMV
        ! ierr = cudaEventRecord(cuEvent1, cuStream)
        ! ierr2 = ierr2 + ierr
        ! write(6,'(A)') ' Axpy '
        ! V = A*PH
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecX, devPtrPH)
        ierr2 = ierr2 + ierr
        ! V = Diag(iwus)*PH
        call kernelc_hadac(devPtrPH, devPtriwus, devPtrV, n, cuStream)
        ierr = cusparseDnVecSetValues(vecY, devPtrV)
        ierr2 = ierr2 + ierr
        ! V = CC*PH + Diag(iwus)*PH
        ierr = cusparseSpMV_cmplx(cusparseHandle,TRANS, &
     &         C_ONE, matA, vecX, C_ONE, vecY, CUDA_C_64F,  &
     &         CUSPARSE_SPMV_CSR_ALG2, buffer)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecY, devPtrV)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
            write(6,'(A,I4)') " Error calculating A*PH", ierr2
            stop
        end if
        ! record - end of SPMV
        ! ierr = cudaEventRecord(cuEvent2, cuStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = cudaDeviceSynchronize()
        ! ierr = cudaEventElapsedTime(timePtr, cuEvent1,cuEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '1 SPMV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !     write(6,'(A,I4)') " error recording the time ", ierr2
        !     stop
        ! end if
        ! record -start of SPMV
        ! ierr = cudaEventRecord(cuEvent1, cuStream)
        ! now check the residual of first half of iteration
        ! calculate the residual norm for the first half of iteration 
        ierr = cublasZdot(cublasHandle,n,devPtrRT,1,devPtrV,1,rtvPtr)
        ierr2 = ierr2 + ierr
        ALPHA = RHO/rtv
        ! write(6,*) 'alpha = ', ALPHA
        if (ALPHA .eq. 0.0) then !bad alpha
            converged=.FALSE.
            KSPiter%failed = .true.
            exit kspLoop
        end if
        ! x = x + ALPHA * PH
        ierr = cublasZaxpy(cublasHandle,n,ALPHA,devPtrPH,1,devPtrX,1)
        ierr2 = ierr2 + ierr
        ! S = R
        ierr = cudaMemCpyAsync(devPtrS, devPtrR, Arow_d_size, &
     &       cudaMemcpyDeviceToDevice)
        ierr2 = ierr2 + ierr
        ! S = R - ALPHA*V
        ierr = cublasZaxpy(cublasHandle,n,-(ALPHA),devPtrV,1,devPtrS,1)
        ierr2 = ierr2 + ierr
        ! ============= second half of the conjugate iteration ============= !
        ! record -start of 2 SPSVs
        ! ierr = cudaEventRecord(cuEvent1, cuStream)
        ! ierr2 = ierr2 + ierr
        ! write(6,'(A)') ' Lsolve '
        ! L solve --> L*ST = S
        ! still need to convert S to FP32
        call kernelc_d2s(devPtrS, devPtr32, n)
        ierr = cusparseDnVecSetValues(vecXs, devPtr32)
        ierr2 = ierr2 + ierr
        ierr = cusparseDnVecSetValues(vecYs, devPtrST)
        ierr2 = ierr2 + ierr
    !    ierr = cusparseSpSV_analysis_cmplx(cusparseHandle,        &
    ! &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONEs,matL,vecXs,vecYs,&
    ! &         CUDA_C_32F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle, bufferL)
    !    ierr2 = ierr2 + ierr
        ierr = cusparseSpSV_solve_cmplx(cusparseHandle,           &
     &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONEs,matL,vecXs,vecYs,&
     &         CUDA_C_32F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecYs, devPtrST)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I2)') " Error during Lsolve ", ierr2
          stop
        end if
        
        ! U solve --> U*SH = ST
        ! write(6,'(A)') ' Usolve '
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecXs, devPtrST)
        ierr2 = ierr2 + ierr
        ierr = cusparseDnVecSetValues(vecYs, devPtr32)
        ierr2 = ierr2 + ierr
        ierr = cusparseSpSV_solve_cmplx(cusparseHandle,            &
     &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONEs,matU,vecXs,vecYs,&
     &         CUDA_C_32F,CUSPARSE_SPSV_ALG_DEFAULT, UsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecYs, devPtr32)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I2)') " Error during Usolve ", ierr2
          stop
        end if
        ! still need to convert SH to double
        call kernelc_s2d(devPtr32, devPtrSH, n)
        ! record - end of 2 SPSVs
        ! ierr = cudaEventRecord(cuEvent2, cuStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = cudaDeviceSynchronize()
        ! ierr = cudaEventElapsedTime(timePtr, cuEvent1,cuEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '2 SPSV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !     write(6, *) " error recording the time ", ierr2
        !     stop
        ! end if
        ! record - start of 2 SPMVs
        ! ierr = cudaEventRecord(cuEvent1, cuStream)
        ! ierr2 = ierr2 + ierr
        ! write(6,'(A)') ' SPMV '
        ! T = A*SH 
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecX, devPtrSH)
        ierr2 = ierr2 + ierr
        ! T = Diag(iwus)*SH
        call kernelc_hadac(devPtrSH, devPtriwus, devPtrT, n, cuStream)
        ierr = cusparseDnVecSetValues(vecY, devPtrT)
        ierr2 = ierr2 + ierr
        ! T = A*SH + Diag(iwus)*SH
        ierr = cusparseSpMV_cmplx(cusparseHandle,TRANS, &
       &        C_ONE, matA, vecX, C_ONE, vecY, CUDA_C_64F,                 & 
       &        CUSPARSE_SPMV_CSR_ALG2, buffer)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecY, devPtrT)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
            write(6,'(A,I4)') " Error calculating Axpy", ierr2
            stop
        end if

        ! record - end of 1 SPMV 
        ! ierr = cudaEventRecord(cuEvent2, cuStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = cudaDeviceSynchronize()
        ! ierr = cudaEventElapsedTime(timePtr, cuEvent1, cuEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '2 SPMV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !     write(6, *) " error recording the time ", ierr2
        !     stop
        ! end if
        ! now check the second half of iteration
        ! calculate the residual norm for the second half of iteration 

        ! tt = dot(T,T)
        ierr = cublasZdot(cublasHandle,n,devPtrT,1,devPtrT,1,ttPtr)
        ierr2 = ierr2 + ierr
        ! rtv = dot(T,S)
        ierr = cublasZdot(cublasHandle,n,devPtrT,1,devPtrS,1,rtvPtr)
        ierr2 = ierr2 + ierr
        OMEGA = rtv/tt
        ! write(6,*) 'omega = ', OMEGA
        if (OMEGA .eq. 0.0) then !bad omega
            KSPiter%failed = .true.
            converged=.false.
            exit kspLoop
        end if
        ! x = x + OMEGA * SH
        ierr = cublasZaxpy(cublasHandle,n,OMEGA,devPtrSH,1,devPtrX,1)
        ierr2 = ierr2 + ierr
        ! R = S
        ierr = cudaMemCpyAsync(devPtrR, devPtrS, Arow_d_size, &
     &       cudaMemcpyDeviceToDevice)
        ierr2 = ierr2 + ierr
        ! R = S - OMEGA * T
        ierr = cublasZaxpy(cublasHandle,n,-(OMEGA),devPtrT,1,devPtrR,1)
        ierr2 = ierr2 + ierr

        ! early second half convergence check (norm of R)
        ! rnorm = norm(r)
        ierr = cublasZnrm2(cublasHandle,n,devPtrR,1,rnormPtr)
        ierr2 = ierr2 + ierr
        ierr = cudaDeviceSynchronize()
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error estimating the second Conj res. ", ierr2
          stop
        end if

        KSPiter%rerr(iter) = rnorm/bnorm
        ! write(6,'(A12,I4,A10,ES12.6)') 'iteration #', iter, ' relres= ', &
        !&           KSPiter%rerr(iter)
        ! now check the second half of iteration
        ! NOTE WE TEST AN IDEA THAT OMIT THE X, WHICH SAVES ANOTHER 
        ! SPMV FOR US
        if (rnorm.lt.btol) then
            ! double check if the residual is really less than tol
            ! R = A*X
            ! still need to use the vecX/Y types
            ierr = cusparseDnVecSetValues(vecX, devPtrX)
            ierr2 = ierr2 + ierr
            ! R = Diag(iwus)*X
            call kernelc_hadac(devPtrX, devPtriwus, devPtrR, n, cuStream)
            ierr = cusparseDnVecSetValues(vecY, devPtrR)
            ierr2 = ierr2 + ierr
            ! R = CC*X + Diag(iwus)*X
            ierr = cusparseSpMV_cmplx(cusparseHandle,TRANS, &
     &          C_ONE, matA, vecX, C_ONE, vecY, CUDA_C_64F,  &
     &          CUSPARSE_SPMV_CSR_ALG2, buffer)
            ierr2 = ierr2 + ierr
            ! and get the values back
            ierr = cusparseDnVecGetValues(vecY, devPtrR)
            ierr2 = ierr2 + ierr
            if (ierr2 .ne. 0 ) then
                write(6,'(A,I4)') " Error calculating A*X", ierr2
                stop
            end if
            ! R = -Ax
            ! ierr = cublasZscal(cublasHandle,n,C_MINUSONE,devPtrR,1)
            ! ierr2 = ierr2 + ierr
            ! R = b - Ax
            ! ierr = cublasZaxpy(cublasHandle,n,C_ONE,devPtrRHS,1,devPtrR,1)
            ! ierr2 = ierr2 + ierr
            ! TEST: R = b - Ax with an all-in-one kernel of xpby
            call kernelc_xpbyc(devPtrRHS, C_MINUSONE, devPtrR, n, cuStream)
            ! rnorm = norm(R)
            ierr = cublasZnrm2(cublasHandle,n,devPtrR,1,rnormPtr)
            ierr2 = ierr2 + ierr
            if (rnorm.lt.btol) then
                KSPiter%rerr(iter) = real(rnorm/bnorm)
                converged=.TRUE.
                KSPiter%failed = .false.
                KSPiter%niter = iter
                exit kspLoop
            end if
        end if
        if (rnorm .lt. rnorm0) then
            rnorm0 = rnorm
            ierr = cudaMemCpyAsync(devPtrX0, devPtrX, Arow_d_size, &
     &           cudaMemcpyDeviceToDevice)
            ierr2 = ierr2 + ierr
        end if
      end do kspLoop
 9527 continue
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error synchronizing after iterations ", ierr2
          stop
      end if
      ! write(6,'(A)') ' Copy solution from GPU to CPU'
      if (.not. converged) then ! solution not found 
          KSPiter%niter=KSPiter%maxit
          ! return the best solution so far
          ierr = cudaMemcpy(xPtr,devPtrX0,Arow_d_size,cudaMemcpyDeviceToHost)
      else ! solution found
          KSPiter%niter=iter
          ! return the last solution
          ierr = cudaMemcpy(xPtr,devPtrX,Arow_d_size,cudaMemcpyDeviceToHost)
      end if
      if (ierr .ne. 0 ) then
          write(6, '(A, I2)') " cudaMemcpy back to host error: ", ierr
          stop
      end if
      ! \activiate lightsaber
      ! clear gpu mem
      ierr = cudaFree(devPtrArow)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrAcol)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrAval)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtriwus)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrLrow)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrLcol)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrLval)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrUrow)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrUcol)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrUval)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrX)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrX0)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrRHS)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrR)
      ierr2 = ierr2 + ierr 
      ierr = cudaFree(devPtrRT)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrP)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrPT)
      ierr2 = ierr2 + ierr 
      ierr = cudaFree(devPtrPH)
      ierr2 = ierr2 + ierr 
      ierr = cudaFree(devPtrS)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrST)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrSH)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtr32)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrT)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrV)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrAX)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(buffer)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(bufferL)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(bufferU)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error during cudafree: ',ierr2
          stop
      end if 
      ierr = cusparseSpSV_destroyDescr(LsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_destroyDescr(UsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroyDnVec(vecX)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroyDnVec(vecY)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroyDnVec(vecXs)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroyDnVec(vecYs)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroySpMat(matL)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroySpMat(matU)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroySpMat(matA)
      ierr2 = ierr2 + ierr
      ierr = cublasDestroy(cublasHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroy(cusparseHandle)
      ierr2 = ierr2 + ierr
      ierr = cudaEventDestroy(cuEvent1)
      ierr2 = ierr2 + ierr
      ierr = cudaEventDestroy(cuEvent2)
      ierr2 = ierr2 + ierr
      ierr = cudaStreamDestroy(cuStream)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error during cuda handle destruction: ',ierr2
          stop
      end if 
      return
end subroutine cuBiCGmix ! cuBiCGmix

#if defined(FG)
subroutine cuBiCGfg(b,x,KSPiter,comm_local,device_idx,adjt)
  ! fine-grained parallel version with GPU support 
  ! Stablized version of BiConjugate Gradient, set up for solving
  ! A x = b 
  ! this solves for the interior (edge) field
  ! essentially we divide the matrix A into several row blocks and 
  ! solve each block with one thread (GPU)
  !
  ! modified (no hell no) to call the CUDA lib to calculate with GPU 
  ! I kept the naming convention from my BiCGStab (see above CPU version)
  ! for instance X -> devPtrX, T -> devPtrT to manipulate the GPU memory
  ! also added the optional adjoint to solve adjoint system A^Tx = b
  ! 
  ! NOTE: BICG actually performs two sub (or half) line searches within a 
  !       iteration, but here we only store the relerr for the second sub
  !       just to be compatitive with QMR
  ! 
  ! interface...........
  ! redefining some of the interfaces for our convenience (locally)
  ! generic routines for vector operations for edge/ face nodes
  ! in a staggered grid
  !
  ! NOTE: this has not been extensively tested - I believe it feels a 
  ! little unstable (dispite the name)...
  ! if you have time reading this, test it!

      use modeloperator3d, only:  Alocal, Llocal, Ulocal, LHlocal, UHlocal,&
   &   vOmegaMuSigLoc
      implicit none
      !  b is right hand side
      complex (kind=prec),intent(in),target,dimension(:)    :: b
      !  solution vector is x ... on input is provided with the initial
      !  guess, on output is the iterate with smallest residual. 
      complex (kind=prec),intent(inout),target,dimension(:) :: x
      type (solverControl_t),intent(inout)                  :: KSPiter
      integer,intent(in)                                    :: comm_local
      integer,intent(in)                                    :: device_idx
      logical,intent(in),optional                           :: adjt
    
      ! local variables
      integer                                    :: ncol, nrow, nnz, iter
      integer                                    :: maxiter, i, j
      complex (kind=prec),pointer,dimension(:)   :: iwus

      real    (kind=prec), target                :: rnorm, bnorm, rnorm0, btol
      complex (kind=prec), target                :: bdot
      real    (kind=prec), target                :: xnorm
      complex (kind=prec)                        :: RHO1, ALPHA, BETA, OMEGA
      complex (kind=prec), target                :: RTV,TT,RHO
      real    (kind=SP), target                  :: ctime
      integer                                    :: TRANS
      integer                                    :: Lfillmode,Ldiagtype
      integer                                    :: Ufillmode,Udiagtype
      integer                                    :: ierr,ierr2, idx
      logical                                    :: converged, adjoint
      integer                                    :: interval
      logical                                    :: restart
      ! fine-grain parallel related    
      integer                                        :: rank_local, size_local
      integer,target                                 :: size_nccl, rank_nccl
      complex (kind=prec),pointer, dimension(:)      :: xbuff, rbuff, bbuff
      complex (kind=prec),target                     :: bdotLoc
      integer (c_size_t), allocatable, dimension(:)  :: isizes, displs
      integer (c_size_t)                             :: fsize, lsize

      ! size to determine corresponding buffer/memory
      integer*8                :: Arowp1_i_size,Arow_d_size
      integer*8                :: Acol_d_size,Annz_i_size,Annz_d_size
      integer*8                :: Mrow_d_size,Lnnz_i_size,Lnnz_d_size
      integer*8                :: Unnz_i_size,Unnz_d_size
      integer(c_int),target    :: zeroLoc
      ! integer(c_int)         :: mbsize
      integer(c_size_t)        :: bsize, lbsize, ubsize

      ! --------------------- pointers to *host* memory ------------------ !   
      type(c_ptr) :: cublasHandle 
      type(c_ptr) :: cusparseHandle
      type(c_ptr) :: cuStream
      type(c_ptr) :: cuEvent1
      type(c_ptr) :: cuEvent2
      type(c_ptr) :: matA
      type(c_ptr) :: matL
      type(c_ptr) :: matU
      type(c_ptr) :: vecX
      type(c_ptr) :: vecXloc
      type(c_ptr) :: vecY 
      type(c_ptr) :: LsolveHandle
      type(c_ptr) :: UsolveHandle
      type(c_ptr) :: ilu_info
      type(c_ptr) :: ArowPtr
      type(c_ptr) :: AcolPtr
      type(c_ptr) :: AvalPtr
      type(c_ptr) :: iwusPtr
      type(c_ptr) :: MrowPtr
      type(c_ptr) :: McolPtr
      type(c_ptr) :: MvalPtr
      type(c_ptr) :: xPtr  
      type(c_ptr) :: bPtr
      type(c_ptr) :: bdotLocPtr
      type(c_ptr) :: rnormPtr0
      type(c_ptr) :: rhoPtr
      type(c_ptr) :: rtvPtr
      type(c_ptr) :: ttPtr
      type(c_ptr) :: zeroLocPtr
      type(c_ptr) :: timePtr
      type(c_ptr) :: xBuffPtr
      type(c_ptr) :: bBuffPtr
      type(c_ptr) :: rBuffPtr
      ! type(ncclUniqueId)  :: uid          ! nccl id
      ! type(ncclComm)      :: comm_nccl    ! nccl communicator
      type(c_ptr)         :: rankPtr      
      type(c_ptr)         :: sizePtr      
      ! -------------------- pointers to *device* memory ----------------- ! 
      type(c_ptr) :: devPtrArow
      type(c_ptr) :: devPtrAcol
      type(c_ptr) :: devPtrAval
      type(c_ptr) :: devPtriwus
      type(c_ptr) :: devPtrLrow
      type(c_ptr) :: devPtrLcol
      type(c_ptr) :: devPtrLval
      type(c_ptr) :: devPtrUrow
      type(c_ptr) :: devPtrUcol
      type(c_ptr) :: devPtrUval
      type(c_ptr) :: devPtrX
      type(c_ptr) :: devPtrX0
      type(c_ptr) :: devPtrRHS
      type(c_ptr) :: devPtrR
      type(c_ptr) :: devPtrRT
      type(c_ptr) :: devPtrP
      type(c_ptr) :: devPtrPT
      type(c_ptr) :: devPtrPH
      type(c_ptr) :: devPtrS
      type(c_ptr) :: devPtrST
      type(c_ptr) :: devPtrSH
      type(c_ptr) :: devPtrV
      type(c_ptr) :: devPtrT
      type(c_ptr) :: devPtrAX
      type(c_ptr) :: buffer
      type(c_ptr) :: bufferL
      type(c_ptr) :: bufferU
      
      ! buffer for the full x, b and r
      type(c_ptr) :: devPtrXbuff ! full x
      ! note we don't check the input parameters compatibility 
      ierr2 = 0
      converged = .FALSE.
      zeroLoc = 0
      if (present(adjt)) then 
          ! write(6,'(A)') ' adjt = ', adjt
          adjoint = adjt
          if (adjt) then
              TRANS = CUSPARSE_OPERATION_TRANSPOSE
              Lfillmode = CUSPARSE_FILL_MODE_UPPER
              Ufillmode = CUSPARSE_FILL_MODE_LOWER
              Ldiagtype = CUSPARSE_DIAG_TYPE_NON_UNIT
              Udiagtype = CUSPARSE_DIAG_TYPE_UNIT
          else
              TRANS = CUSPARSE_OPERATION_NON_TRANSPOSE
              Lfillmode = CUSPARSE_FILL_MODE_LOWER
              Ufillmode = CUSPARSE_FILL_MODE_UPPER
              Ldiagtype = CUSPARSE_DIAG_TYPE_UNIT
              Udiagtype = CUSPARSE_DIAG_TYPE_NON_UNIT
          end if
      else
          adjoint = .FALSE.
          TRANS = CUSPARSE_OPERATION_NON_TRANSPOSE
          Lfillmode = CUSPARSE_FILL_MODE_LOWER
          Ufillmode = CUSPARSE_FILL_MODE_UPPER
          Ldiagtype = CUSPARSE_DIAG_TYPE_UNIT
          Udiagtype = CUSPARSE_DIAG_TYPE_NON_UNIT
      end if
      
      ! now see how many workers do we have
      call MPI_COMM_RANK(comm_local,rank_local,ierr)
      call MPI_COMM_SIZE(comm_local,size_local,ierr)
      ! firstly try to figure out how the data is distributed on
      ! each process
      allocate(isizes(size_local))
      allocate(displs(size_local))
      ! local array size
      lsize = size(x)
      ! try to get the local row sizes for each process, initialize by zero
      isizes = 0
      ! displacement
      displs = 0
      ! NOTE it can be used in-place, if the comm_local is a intra-communicator
      ! it should be - but I will refrain to do that as who knows our users
      ! will configure their MPI processes...
      call MPI_ALLGATHER(lsize,  2, MPI_INTEGER, isizes, 2, MPI_INTEGER &
 &            , comm_local, ierr)
      fsize = sum(isizes)
      ! and calculate the displacement
      do i=2,size_local
          displs(i) = sum(isizes(1:i-1))
      end do
      ! for debug
      ! write(6, *) 'local size =', lsize, ' displs = ', displs(rank_local+1), &
      !&       'full size =', fsize, 'rank =', rank_local
      ! note we still keep the Ae = CCGDe + iwuse 
      ! idea as in ModEM SP2, we will test if this will save some time 
      ! (probably not much)
      ! call CSR_R2Cdiag(AAii,VOmegaMuSig,Aii)
      ncol = fsize
      nrow = lsize
      nnz = size(Alocal%col)
      ArowPtr = c_loc(Alocal%row)
      AcolPtr = c_loc(Alocal%col)
      AvalPtr = c_loc(Alocal%val)
      ! other host pointers
      xPtr  = c_loc(x)  ! x = A \ b --> local x
      bPtr  = c_loc(b)  ! b = ones(n,1) --> local b
      ! buffer
      allocate(xbuff(fsize))
      allocate(bbuff(fsize))
      xBuffPtr  = c_loc(xbuff)  ! full x
      bBuffPtr  = c_loc(bbuff)  ! full b
      bdotLocPtr = c_loc(bdotLoc)
      allocate(iwus(nrow))
      ! only the local part of iwus
      iwus = VOmegaMuSigLoc*ISIGN*CMPLX(0.0,1.0,8)
      iwusPtr = c_loc(iwus) ! i V omega mu sigma (local)
      rnormPtr0 = c_loc(rnorm0)
      rhoPtr = c_loc(RHO)
      rtvPtr = c_loc(rtv)
      ttPtr = c_loc(tt)
      zeroLocPtr = c_loc(zeroLoc)
      ! pointer to store the time 
      ctime = 0.0 
      timePtr = c_loc(ctime)
      
      ! get the size of the vector and matrix (need to allocate on the device)
      Arowp1_i_size=sizeof(Alocal%row(1:nrow+1))
      Arow_d_size=sizeof(b(1:nrow))
      Acol_d_size=sizeof(x(1)) * ncol ! size of the (full) x vector
      Annz_i_size=sizeof(Alocal%col(1:nnz))
      Annz_d_size=sizeof(Alocal%val(1:nnz))

      ! select the current cuda device 
      ! note this only works for physical devices (not working for MIG devices)
      ! ierr = cudaSetDevice(device_idx);
      ! ierr2 = ierr2 + ierr
      ! firstly define the CUDA Stream and cuda handles
      ierr = cudaStreamCreateWithFlags(cuStream, cudaStreamNonBlocking) 
      ierr2 = ierr2 + ierr
      ! initialize the cusparse lib
      ierr = cusparseCreate(cusparseHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseSetStream(cusparseHandle,cuStream) 
      ierr2 = ierr2 + ierr
      ! now initialize the cublas lib
      ierr = cublasCreate(cublasHandle)
      ierr2 = ierr2 + ierr
      ierr = cublasSetStream(cublasHandle,cuStream)
      ierr2 = ierr2 + ierr
      ! and creates two cuda events to record the time
      ierr = cudaEventCreate(cuEvent1)
      ierr2 = ierr2 + ierr
      ierr = cudaEventCreate(cuEvent2)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I4)') 'Error during cuda initialize ',ierr2
          stop
      end if 
      ! if (ncclIsInit.eq.0) then
      !        ! and init the NCCL communicator
      !        if (rank_local .eq. 0) then
      !            ! leader generating the uniqueId
      !            ierr = ncclGetUniqueId(uid)
      !            ierr2 = ierr2 + ierr
      !            if (ierr2 .ne. 0) then
      !               write(0,*) 'error getting NCCL unique id on device:', ierr
      !               stop
      !            endif
      !        endif
      !        ! distribute the ID to all workers, using MPI
      !        call MPI_BCAST(uid%internal, NCCL_UNIQUE_ID_BYTES, MPI_CHAR, 0, &
      !      &         comm_local, ierr)
      !        ! for debug
      !        ! write(6,*) 'uid = ', uid%internal, ' @rank: ', rank_local
      sizePtr = c_loc(size_nccl)
      rankPtr = c_loc(rank_nccl)
      size_nccl = size_local
      rank_nccl = rank_local
      !        ! for debug
      !        write(0,*) 'initializing NCCL communicator... @ ', rank_nccl
      !        ierr = ncclCommInitRank(comm_nccl, size_nccl, uid, rank_nccl)
      !        ierr2 = ierr2 + ierr
      !        if (ierr2.ne.0) then
      !            write(0,'(A, I4)') 'Error initializing nccl ',ierr2
      !            stop
      !        end if 
      !        ncclIsInit = 1
      ! end if
      ierr = ncclCommUserRank(comm_nccl, rankPtr)
      ierr2 = ierr2 + ierr
      ierr = ncclCommCount(comm_nccl, sizePtr)
      ierr2 = ierr2 + ierr
      ! record the event before memory manipulation
      ! ierr = cudaEventRecord(cuEvent1, cuStream)
      ! ierr2 = ierr2 + ierr
      ! write(6,*) 'Allocating GPU memory'
      ierr = cudaMalloc(devPtrX,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrX0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrRHS,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrR,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrRT,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrP,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrPT,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrPH,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrS,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrST,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrSH,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrT,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrV,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrAX,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrAval,Annz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrAcol,Annz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrArow,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtriwus,Arow_d_size)
      ierr2 = ierr2 + ierr
      ! alloc the full vector buffer
      ierr = cudaMalloc(devPtrXbuff,Acol_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error during CUDA allocation: ',ierr2
          stop
      end if 
      ! write(6,*) 'reset GPU memory to all zeros'
      ierr = cudaMemset(devPtrX,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrX0,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrRHS,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrR,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrRT,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrP,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrPT,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrPH,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrS,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrST,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrSH,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrT,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrV,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrAX,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrAval,0,Annz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrAcol,0,Annz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrArow,0,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtriwus,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ! set the full buffer to zero
      ierr = cudaMemset(devPtrXbuff,0,Acol_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I3)') 'Error during device memory reseting : ',ierr2
          stop
      end if 
      ! transfer memory over to GPU
      ! write(6,*) 'Transferring (local) memory to GPU'
      ! initialize (local) A
      ierr = cudaMemcpyAsync(devPtrArow,ArowPtr,Arowp1_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrAcol,AcolPtr,Annz_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrAval,AvalPtr,Annz_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtriwus,iwusPtr,Arow_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ! finally need to create the SpMatDescr 
      ierr = cusparseCreateCsr(matA, nrow, ncol, nnz, devPtrArow, devPtrAcol, &
     &       devPtrAval, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, &
     &       CUSPARSE_INDEX_BASE_ONE, CUDA_R_64F) 
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, *) " local matrix size = ", nrow, ' * ', ncol
          write(6, '(A,I2)') " error assembling system matrix ", ierr2
          stop
      end if
      ! initialize rhs and x
      ierr = cudaMemcpyAsync(devPtrX,xPtr,Arow_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrRHS,bPtr,Arow_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error during cuda memcpy ", ierr2
          stop
      end if
      ! now deallocate temp array
      if ( c_associated(iwusPtr) ) then
          ! this is a little tricky as iwus is associated with
          ! iwusPtr in C
          nullify(iwus)
          call cf_free(iwusPtr)
      end if
      ! for debug
      ! write(6,*) 'nccl gather from all GPUs #', rank_nccl
      ! note that nccl doesn't really support the complex communication
      ierr =  ncclAllGatherV(devPtrX, lsize*2, ncclFloat64, devPtrXbuff, &
     &     isizes*2, displs*2, 0, size_nccl, comm_nccl, cuStream)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0) then 
          write(6,*) ' error gather data from devices: ', ierr
          stop
      endif
      ! firstly we need to create three dense vectors 
      ! as the user-hostile developers in Nvidia think of a new idea 
      ! to mess up the interfaces
      ! full x
      ierr = cusparseCreateDnVec(vecX, ncol, devPtrXbuff, CUDA_C_64F)
      ierr2 = ierr2 + ierr
      ! local x
      ierr = cusparseCreateDnVec(vecXloc, nrow, devPtrX, CUDA_C_64F)
      ierr2 = ierr2 + ierr
      ! local b/r 
      ierr = cusparseCreateDnVec(vecY, nrow, devPtrR, CUDA_C_64F)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, *) " error creating the dense Vecs "
          stop
      end if
      ierr = cusparseSpMV_bufferSize_cmplx(cusparseHandle,              &
     &        TRANS, C_ONE, matA, vecX, C_ONE, vecY,                    &
     &        CUDA_C_64F, CUSPARSE_SPMV_CSR_ALG2, bsize)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error estimating buffer for spMV ",  ierr2
          stop
      end if
      ! write(6,*) 'spmv buffersize = ', mbsize
      ! and finally (re)allocate the buffer
      ierr = cudaMalloc(buffer, bsize)
      ierr2 = ierr2 + ierr
      ! now let's deal with L and U
      ! write(6,'(A)') ' Setup L and U preconditioners on GPU'
      ! L first
      if (adjoint) then
          nnz = size(LHlocal%col)
          ! sizes (row size is the same as A
          Lnnz_i_size=sizeof(LHlocal%col(1:nnz))
          Lnnz_d_size=sizeof(LHlocal%val(1:nnz))
          MrowPtr = c_loc(LHlocal%row)
          McolPtr = c_loc(LHlocal%col)
          MvalPtr = c_loc(LHlocal%val)
      else
          nnz = size(Llocal%col)
          ! sizes (row size is the same as A
          Lnnz_i_size=sizeof(Llocal%col(1:nnz))
          Lnnz_d_size=sizeof(Llocal%val(1:nnz))
          MrowPtr = c_loc(Llocal%row)
          McolPtr = c_loc(Llocal%col)
          MvalPtr = c_loc(Llocal%val)
      end if
      ! initialize L in GPU memory
      ierr = cudaMalloc(devPtrLval,Lnnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrLcol,Lnnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrLrow,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! set everything to zeros
      ierr = cudaMemset(devPtrLval,0,Lnnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrLcol,0,Lnnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrLrow,0,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! now copy the values to device
      ierr = cudaMemcpyAsync(devPtrLval,MvalPtr,Lnnz_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrLcol,McolPtr,Lnnz_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrLrow,MrowPtr,Arowp1_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ! finish synchronizing
      ierr = cudaDeviceSynchronize()
      ierr2 = ierr2 + ierr
      ! finally need to create the SpMatDescr 
      ierr = cusparseCreateCsr(matL,nrow,nrow, nnz, devPtrLrow, devPtrLcol, &
     &       devPtrLval, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, &
     &       CUSPARSE_INDEX_BASE_ONE, CUDA_C_64F) 
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error assembling L matrix ", ierr2
          stop
      end if
      ! now start messing up with SpSV
      ! now estimate the buffersize needed by SpSV (Lsolve)
      ! solves y in L*y = a*x (if a=1)
      ierr = cusparseSpMatSetAttribute(matL, CUSPARSE_SPMAT_FILL_MODE, &
    &        Lfillmode, 4)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpMatSetAttribute(matL, CUSPARSE_SPMAT_DIAG_TYPE, &
    &        Ldiagtype, 4)
      ierr2 = ierr2 + ierr
      ! still need to establish a context handler for Lsolve
      ierr = cusparseSpSV_createDescr(LsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_bufferSize_dcmplx(cusparseHandle,             &
     &        CUSPARSE_OPERATION_NON_TRANSPOSE, C_ONE, matL, vecXloc, vecY,&
     &        CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle, lbsize)
      ierr2 = ierr2 + ierr

      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error estimating buffer for Lsolve ",  ierr2
          stop
      end if
      ierr = cudaMalloc(bufferL, lbsize)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_analysis_dcmplx(cusparseHandle,         &
     &        CUSPARSE_OPERATION_NON_TRANSPOSE, C_ONE, matL, vecXloc, vecY, &
     &        CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle, bufferL)
      ierr2 = ierr2 + ierr
      ! write(6,*) 'spsv Lsolve buffersize = ', lbsize
      ! then U 
      if (adjoint) then
          nnz = size(UHlocal%col)
          ! sizes (row size is the same as A
          Unnz_i_size=sizeof(UHlocal%col(1:nnz))
          Unnz_d_size=sizeof(UHlocal%val(1:nnz))
          MrowPtr = c_loc(UHlocal%row)
          McolPtr = c_loc(UHlocal%col)
          MvalPtr = c_loc(UHlocal%val)
      else
          nnz = size(Ulocal%col)
          ! sizes (row size is the same as A
          Unnz_i_size=sizeof(Ulocal%col(1:nnz))
          Unnz_d_size=sizeof(Ulocal%val(1:nnz))
          MrowPtr = c_loc(Ulocal%row)
          McolPtr = c_loc(Ulocal%col)
          MvalPtr = c_loc(Ulocal%val)
      end if
      ! initialize U in GPU memory
      ierr = cudaMalloc(devPtrUval,Unnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrUcol,Unnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMalloc(devPtrUrow,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! set everything to zeros
      ierr = cudaMemset(devPtrUval,0,Unnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrUcol,0,Unnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = cudaMemset(devPtrUrow,0,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! now copy the values to device
      ierr = cudaMemcpyAsync(devPtrUval,MvalPtr,Unnz_d_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrUcol,McolPtr,Unnz_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = cudaMemcpyAsync(devPtrUrow,MrowPtr,Arowp1_i_size, &
     &        cudaMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ! finally need to create the SpMatDescr 
      ierr = cusparseCreateCsr(matU, nrow, nrow, nnz, devPtrUrow, devPtrUcol, &
     &       devPtrUval, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, &
     &       CUSPARSE_INDEX_BASE_ONE, CUDA_C_64F) 
      ierr2 = ierr2 + ierr
      ! now estimate the buffersize needed by SpSV (Usolve)
      ! solves y in U*y = a*x (if a=1)
      ierr = cusparseSpMatSetAttribute(matU, CUSPARSE_SPMAT_FILL_MODE, &
    &        Ufillmode, 4)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpMatSetAttribute(matU, CUSPARSE_SPMAT_DIAG_TYPE, &
    &        Udiagtype, 4)
      ierr2 = ierr2 + ierr
      ! still need to establish a context handler for Usolve
      ierr = cusparseSpSV_createDescr(UsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_bufferSize_dcmplx(cusparseHandle,             &
     &        CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matU, vecXloc, vecY,&
     &        CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, UsolveHandle, ubsize)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error estimating buffer for Usolve ",  ierr2
          stop
      end if
      ! write(6,*) 'spsv Usolve buffersize = ', ubsize
      ierr = cudaMalloc(bufferU, ubsize)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_analysis_dcmplx(cusparseHandle,         &
     &        CUSPARSE_OPERATION_NON_TRANSPOSE, C_ONE, matU, vecXloc, vecY, &
     &        CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, UsolveHandle, bufferU)
      ierr2 = ierr2 + ierr
      !write(6, '(A, I8)') " allocated SpMV/SpSV buffer on GPU (MB) ",&
     !&            bsize/1024/1024
      ! see how long it takes for the intialization part
      ! ierr = cudaEventRecord(cuEvent2, cuStream)
      ! ierr2 = ierr2 + ierr
      ! ierr = cudaDeviceSynchronize()
      ! ierr = cudaEventElapsedTime(timePtr, cuEvent1, cuEvent2)
      ! ierr2 = ierr2 + ierr
      ! write(6,*) 'initial GPU memory cpy time is ', ctime, 'ms'
      ! if (ierr2 .ne. 0 ) then
      !     write(6, *) " error recording the time ", ierr2
      !     stop
      ! end if
      ! firstly pin the host memory using cudaHostRegister
      ! ierr = cudaHostRegister(bdotLocPtr, INT8(16), 1)
      ! ierr2 = ierr2 + ierr
      ! if (ierr2 .ne. 0 ) then
      !     write(0, *) " error pinning the host memory ", ierr
      !     stop
      ! end if
      ! Norm of rhs
      ! the idea is to calculate the dot product of blocal for all processes and
      ! sum the result
      ! bnorm = nrm2(b)
      ierr = cublasZdot(cublasHandle,nrow,devPtrRHS,1,devPtrRHS,1,bdotLocPtr)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(0, *) " error with Zdot operation ", ierr
          stop
      end if
      ! note that ALLREDUCE is equivalent to a REDUCE and a BCAST
      call MPI_ALLREDUCE(bdotLoc, bdot, 1, MPI_DOUBLE_COMPLEX, MPI_SUM,    &
 &        comm_local, ierr)
      ierr2 = ierr2 + ierr
      ! here every process stores a copy of bnorm
      bnorm = sqrt(abs(bdot))
      if (isnan(bnorm)) then
      ! this usually means an inadequate model, in which case Maxwell's failed
          write(0,*) 'Error: b in BICG contains NaNs; exiting...'
          stop
      else if ( bnorm .eq. 0.0) then ! zero rhs -> zero solution
          write(0,*) 'Warning: b in BICG has all zeros, returning zero &
     &        solution'
          x = b 
          KSPiter%niter=1
          KSPiter%failed=.false.
          KSPiter%rerr(1)=0.0
          converged = .true.
          goto 9527 
      endif
      ! now Compute the initial residual
      ! write(6,*) 'Computing initial residual'
      ! R = Diag(iwus)*X <-- diagonal
      call kernelc_hadac(devPtrX, devPtriwus, devPtrR, nrow, cuStream)
      ! setup the vecX and vecY
      ierr = cusparseDnVecSetValues(vecX, devPtrXbuff)
      ierr2 = ierr2 + ierr
      ierr = cusparseDnVecSetValues(vecY, devPtrR)
      ierr2 = ierr2 + ierr
      ! now calculate the y = Ax
      ! R = CC*X + Diag(iwus)*X
      ! thank God (or whatever deserves it) that SpMV does not 
      ! need to be analyzed
      ierr = cusparseSpMV_cmplx(cusparseHandle,TRANS,                  &
     &        C_ONE, matA, vecX, C_ONE, vecY, CUDA_C_64F,              &
     &        CUSPARSE_SPMV_CSR_ALG2, buffer)
      ierr2 = ierr2 + ierr
      ! and get the values back 
      ierr = cusparseDnVecGetValues(vecY, devPtrR)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, *) " error with SpMV operation ", ierr2
          stop
      end if
      ! now calculate the residual
      ! R = -Ax 
      ! ierr = cublasZscal(cublasHandle,nrow,C_MINUSONE,devPtrR,1)
      ! ierr2 = ierr2 + ierr
      ! R = b - Ax
      ! ierr = cublasZaxpy(cublasHandle,nrow,C_ONE,devPtrRHS,1,devPtrR,1)
      ! ierr2 = ierr2 + ierr
      ! TEST: R = b - Ax with an all-in-one kernel of xpby
      call kernelc_xpbyc(devPtrRHS, C_MINUSONE, devPtrR, nrow, cuStream)
      ! rnorm = nrm2(b - Ax)
      ierr = cublasZdot(cublasHandle,nrow,devPtrR,1,devPtrR,1,bdotLocPtr)
      ierr2 = ierr2 + ierr
      ! note that ALLREDUCE is equivalent to a REDUCE and a BCAST
      call MPI_ALLREDUCE(bdotLoc, bdot, 1, MPI_DOUBLE_COMPLEX, MPI_SUM,  &
 &        comm_local, ierr)
      ierr2 = ierr2 + ierr
      ! here every process stores a copy of rnorm
      rnorm = sqrt(abs(bdot))
      rnorm0 = rnorm
      if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error during residual estimation", ierr2
          stop
      end if
      ! write(6,'(A, ES12.6)') ' initial relative residual = ', rnorm/bnorm
      !================= Now start configuring the iteration ================!
      btol = KSPiter%tol*bnorm
      maxIter = KSPiter%maxIt
      KSPiter%rerr(1) = real(rnorm/bnorm)
      if ( rnorm .le. btol ) then ! the first guess is already good enough
         ! returning
          write(6, *) " The first guess is good enough, exiting..."
          iter=1
          KSPiter%failed=.false.
          converged = .true.
          goto 9527 
      end if
      ! intial the parameters for restarting
      restart = .FALSE.
      ! hard coded here, need a more elegant way to deal with it
      interval = 120
      ! write(6,'(A,I4)') ' maxiter = ', maxiter
      ! RHO = C_ONE
      kspLoop: do iter = 1,maxIter
        ! write(6,'(A, I4)') ' KSP iteration #', iter
        if (mod(iter, interval) .eq. 1) then ! hard coded here
            restart = .TRUE.
        end if
        if (restart) then 
            ! restart the iteration (to steepest decend) every interval times
            ! current and previous RHO, RHO0 should be one
            RHO1 = C_ONE
            OMEGA = C_ONE
            ! RT = R (local)
            ierr = cublasZcopy(cublasHandle,nrow,devPtrR,1,devPtrRT,1)
            ierr2 = ierr2 + ierr
            ! RHO = dot(RT,R)
            ierr = cublasZdot(cublasHandle,nrow,devPtrRT,1,devPtrR,1,bdotLocPtr)
            ierr2 = ierr2 + ierr
            call MPI_ALLREDUCE(bdotLoc, bdot, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, &
    &           comm_local, ierr)
            ierr2 = ierr2 + ierr
            RHO = bdot
            if (RHO .eq. 0.0) then !bad beta
                KSPiter%failed = .true.
                exit
            endif
            ! P = R (local)
            ierr = cublasZcopy(cublasHandle,nrow,devPtrR,1,devPtrP,1)
            ierr2 = ierr2 + ierr
            restart = .FALSE.
            if (ierr2 .ne. 0 ) then
              write(6, '(A, I2)') " Error steering search direction ", ierr2
              stop
            end if
        else
            ! save the previous RHO
            RHO1 = RHO
            ! RHO = dot(RT,R)
            ierr = cublasZdot(cublasHandle,nrow,devPtrRT,1,devPtrR,1,bdotLocPtr)
            ierr2 = ierr2 + ierr
            call MPI_ALLREDUCE(bdotLoc, bdot, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, &
    &           comm_local, ierr)
            ierr2 = ierr2 + ierr
            RHO = bdot
            BETA=(RHO/RHO1)*(ALPHA/OMEGA)
            if (BETA .eq. 0.0) then !bad beta
                converged=.FALSE.
                KSPiter%failed = .true.
                exit kspLoop
            end if
            ! P = R + BETA * (P - OMEGA * V)
            ! P = P - OMEGA * V
            ! ierr = cublasZaxpy(cublasHandle,nrow,-(OMEGA),devPtrV,1,devPtrP,1)
            ! ierr2 = ierr2 + ierr
            ! P = BETA * P
            ! ierr = cublasZscal(cublasHandle,nrow,BETA,devPtrP,1)
            ! ierr2 = ierr2 + ierr
            ! P = R + P
            ! ierr = cublasZaxpy(cublasHandle,nrow,C_ONE,devPtrR,1,devPtrP,1)
            ! ierr2 = ierr2 + ierr
            ! TEST: P = R + BETA * P  with an all-in-one kernel 
            ! call kernelc_xpbyc(devPtrR, BETA, devPtrP, n)
            ! TEST: update (local) P with an all-in-one kernel 
            call kernelc_update_pc(devPtrR, devPtrV,BETA,OMEGA, devPtrP,nrow,&
      &         cuStream)
            if (ierr2 .ne. 0 ) then
              write(6, '(A, I2)') " Error steering search direction ", ierr2
              stop
            end if
        end if
        ! record - start of two SPSVs
        ! ierr = cudaEventRecord(cuEvent1, cuStream)
        ! ierr2 = ierr2 + ierr
        ! ============== first half of the conjugate iteration ============= !
        ! L solve --> L*PT = P
        ! write(6,'(A)') ' Lsolve '
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecXloc, devPtrP)
        ierr2 = ierr2 + ierr
        ierr = cusparseDnVecSetValues(vecY, devPtrPT)
        ierr2 = ierr2 + ierr
        ierr = cusparseSpSV_solve_dcmplx(cusparseHandle,            &
     &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matL, vecXloc, vecY, &
     &         CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecY, devPtrPT)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I2)') " Error during Lsolve ", ierr2
          stop
        end if
        ! U solve --> U*PH = PT
        ! write(6,'(A)') ' Usolve '
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecXloc, devPtrPT)
        ierr2 = ierr2 + ierr
        ierr = cusparseDnVecSetValues(vecY, devPtrPH)
        ierr2 = ierr2 + ierr
        ierr = cusparseSpSV_solve_dcmplx(cusparseHandle,            &
     &         CUSPARSE_OPERATION_NON_TRANSPOSE, C_ONE, matU, vecXloc, vecY, &
     &         CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, UsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecY, devPtrPH)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I2)') " Error during Usolve ", ierr2
          stop
        end if
        ! record - end of two SPSVs
        ! ierr = cudaEventRecord(cuEvent2, cuStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = cudaDeviceSynchronize()
        ! ierr = cudaEventElapsedTime(timePtr, cuEvent1, cuEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '2 SPSV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !    write(6,'(A,I4)') " error recording the time ", ierr2
        !     stop
        ! end if
        ! record - start of one SPMV
        ! ierr = cudaEventRecord(cuEvent1, cuStream)
        ! ierr2 = ierr2 + ierr
        ! V = A*PH
        ! write(6,'(A)') ' Axpy '
        ! for debug
        ! write(6,*) 'nccl gather from all GPUs #', rank_nccl
        ! need to gather all local PH results to xbuff
        ierr =  ncclAllGatherV(devPtrPH, lsize*2, ncclFloat64, devPtrXbuff, &
     &     isizes*2, displs*2, 0, size_nccl, comm_nccl, cuStream)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
            write(6,'(A,I4)') " Error gathering PH ", ierr2
            stop
        end if
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecX, devPtrXbuff)
        ierr2 = ierr2 + ierr
        ! now calculate V = A * PH (local)
        ! V = Diag(iwus)*PH (local)
        call kernelc_hadac(devPtrPH, devPtriwus, devPtrV, nrow, cuStream)
        ierr = cusparseDnVecSetValues(vecY, devPtrV)
        ierr2 = ierr2 + ierr
        ! V = CC*PH + Diag(iwus)*PH (local)
        ierr = cusparseSpMV_cmplx(cusparseHandle,TRANS, &
     &         C_ONE, matA, vecX, C_ONE, vecY, CUDA_C_64F,  &
     &         CUSPARSE_SPMV_CSR_ALG2, buffer)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecY, devPtrV)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
            write(6,'(A,I4)') " Error calculating A*PH", ierr2
            stop
        end if
        ! record - end of one SPMV
        ! ierr = cudaEventRecord(cuEvent2, cuStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = cudaDeviceSynchronize()
        ! ierr = cudaEventElapsedTime(timePtr, cuEvent1,cuEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '1 SPMV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !    write(6, *) " error recording the time ", ierr2
        !    stop
        ! end if
        ! now check the residual of first half of iteration
        ! rtv = dot(RT,V)
        ierr = cublasZdot(cublasHandle,nrow,devPtrRT,1,devPtrV,1,bdotLocPtr)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
            write(6,'(A,I4)') " Error calculating dot(RT,V)", ierr2
            stop
        end if
        ! note that ALLREDUCE is equivalent to a REDUCE and a BCAST
        call MPI_ALLREDUCE(bdotLoc, bdot, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, &
 &          comm_local, ierr)
        ierr2 = ierr2 + ierr
        RTV = bdot
        ! wait for the others
        call MPI_BARRIER(comm_local, ierr)
        ierr2 = ierr2 + ierr
        ALPHA = RHO/RTV
        ! write(6,*) 'ALPHA = ', ALPHA
        if (ALPHA .eq. 0.0) then !bad alpha
            KSPiter%failed = .true.
            converged=.FALSE.
            exit kspLoop
        end if
        ! x = x + ALPHA * PH (local)
        ierr = cublasZaxpy(cublasHandle,nrow,ALPHA,devPtrPH,1,devPtrX,1)
        ierr2 = ierr2 + ierr
        ! S = R (local)
        ierr = cudaMemcpyAsync(devPtrS, devPtrR, Arow_d_size, &
       &       cudaMemcpyDeviceToDevice)
        ierr2 = ierr2 + ierr
        ! S = R - ALPHA * V (local)
        ierr = cublasZaxpy(cublasHandle,nrow,-(ALPHA),devPtrV,1,devPtrS,1)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
            write(6, '(A, I4)') " Error estimating the first Conj res. ", ierr2
            stop
        end if
        ! ============== second half of the conjugate iteration ============= !
        ! record - start of two SPSV
        ! ierr = cudaEventRecord(cuEvent1, cuStream)
        ! ierr2 = ierr2 + ierr
        ! L solve --> L*ST = S
        ! write(6,'(A)') ' Lsolve '
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecXloc, devPtrS)
        ierr2 = ierr2 + ierr
        ierr = cusparseDnVecSetValues(vecY, devPtrST)
        ierr2 = ierr2 + ierr
    !    ierr = cusparseSpSV_analysis_dcmplx(cusparseHandle,        &
    ! &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matL, vecXloc, vecY, &
    ! &         CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle, bufferL)
    !    ierr2 = ierr2 + ierr
        ierr = cusparseSpSV_solve_dcmplx(cusparseHandle,           &
     &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matL, vecXloc, vecY, &
     &         CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, LsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecY, devPtrST)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
            write(6, '(A, I2)') " Error during Lsolve ", ierr2
            stop
        end if
        ! U solve --> U*SH = ST
        ! write(6,'(A)') ' Usolve '
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecXloc, devPtrST)
        ierr2 = ierr2 + ierr
        ierr = cusparseDnVecSetValues(vecY, devPtrSH)
        ierr2 = ierr2 + ierr
    !    ierr = cusparseSpSV_analysis_dcmplx(cusparseHandle,         &
    ! &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matU, vecXloc, vecY, &
    ! &         CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, UsolveHandle, bufferU)
    !    ierr2 = ierr2 + ierr
        ierr = cusparseSpSV_solve_dcmplx(cusparseHandle,            &
     &         CUSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matU, vecXloc, vecY, &
     &         CUDA_C_64F,CUSPARSE_SPSV_ALG_DEFAULT, UsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = cusparseDnVecGetValues(vecY, devPtrSH)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
            write(6, '(A, I2)') " Error during Usolve ", ierr2
            stop
        end if

        ! record - end of 2 SPSVs
        ! ierr = cudaEventRecord(cuEvent2, cuStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = cudaDeviceSynchronize()
        ! ierr = cudaEventElapsedTime(timePtr, cuEvent1,cuEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '2 SPSV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !     write(6, *) " error recording the time ", ierr2
        !     stop
        ! end if
        ! record - start of 1 SPMV
        ! ierr = cudaEventRecord(cuEvent1, cuStream)
        ! ierr2 = ierr2 + ierr
        ! T = A*SH 
        ! write(6,'(A)') ' Axpy '
        ! for debug
        ! write(6,*) 'nccl gather from all GPUs #', rank_nccl
        ! need to gather all local PH results to xbuff
        ierr =  ncclAllGatherV(devPtrSH, lsize*2, ncclFloat64, devPtrXbuff, &
     &     isizes*2, displs*2, 0, size_nccl, comm_nccl, cuStream)
        ierr2 = ierr2 + ierr
        ! still need to use the vecX/Y types
        ierr = cusparseDnVecSetValues(vecX, devPtrXbuff)
        ierr2 = ierr2 + ierr
        ! T = Diag(iwus)*SH (local)
        call kernelc_hadac(devPtrSH, devPtriwus, devPtrT, nrow, cuStream)
        ierr = cusparseDnVecSetValues(vecY, devPtrT)
        ierr2 = ierr2 + ierr
        ! T = A*SH + Diag(iwus)*SH (local)
        ierr = cusparseSpMV_cmplx(cusparseHandle,TRANS, &
       &        C_ONE, matA, vecX, C_ONE, vecY, CUDA_C_64F,                 & 
       &        CUSPARSE_SPMV_CSR_ALG2, buffer)
        ierr2 = ierr2 + ierr
        ! and get the (local) values back 
        ierr = cusparseDnVecGetValues(vecY, devPtrT)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
            write(6,'(A,I4)') " Error calculating A*SH", ierr2
            stop
        end if

        ! record - end of 1 SPMV 
        ! ierr = cudaEventRecord(cuEvent2, cuStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = cudaDeviceSynchronize()
        ! ierr = cudaEventElapsedTime(timePtr, cuEvent1,cuEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '2 SPMV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !     write(6, *) " error recording the time ", ierr2
        !     stop
        ! end if

        ! now check the second half of iteration
        ! calculate the residual norm for the second half of iteration 
        ! tt = dot(T,T)
        ierr = cublasZdot(cublasHandle,nrow,devPtrT,1,devPtrT,1,bdotLocPtr)
        ierr2 = ierr2 + ierr
        ! note that ALLREDUCE is equivalent to a REDUCE and a BCAST
        call MPI_ALLREDUCE(bdotLoc, bdot, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, &
 &          comm_local, ierr)
        ierr2 = ierr2 + ierr
        TT = bdot
        ! tts = dot(T,S)
        ierr = cublasZdot(cublasHandle,nrow,devPtrT,1,devPtrS,1,bdotLocPtr)
        ierr2 = ierr2 + ierr
        ! note that ALLREDUCE is equivalent to a REDUCE and a BCAST
        call MPI_ALLREDUCE(bdotLoc, bdot, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, &
 &          comm_local, ierr)
        ierr2 = ierr2 + ierr
        OMEGA = bdot
        OMEGA = OMEGA/TT
        ! wait for the others
        call MPI_BARRIER(comm_local, ierr)
        ierr2 = ierr2 + ierr
        ! write(6,*) 'OMEGA = ', OMEGA
        if (OMEGA .eq. 0.0) then !bad omega
            KSPiter%failed = .true.
            converged=.false.
            exit kspLoop
        end if
        ! x = x + OMEGA * SH (local)
        ierr = cublasZaxpy(cublasHandle,nrow,OMEGA,devPtrSH,1,devPtrX,1)
        ierr2 = ierr2 + ierr
        ! R = S
        ierr = cudaMemcpyAsync(devPtrR, devPtrS, Arow_d_size, &
       &       cudaMemcpyDeviceToDevice)
        ierr2 = ierr2 + ierr
        ! R = S - OMEGA * T (local)
        ierr = cublasZaxpy(cublasHandle,nrow,-(OMEGA),devPtrT,1,devPtrR,1)
        ierr2 = ierr2 + ierr
        ! early second half convergence check (norm of R)
        ! rnorm = norm(R)
        ierr = cublasZdot(cublasHandle,nrow,devPtrR,1,devPtrR,1,bdotLocPtr)
        ierr2 = ierr2 + ierr
        call MPI_ALLREDUCE(bdotLoc, bdot, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, &
    &          comm_local, ierr)
        rnorm = sqrt(bdot)
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error estimating the second Conj res. ", ierr2
          stop
        end if
        KSPiter%rerr(iter) = rnorm/bnorm
        ! for debug
        ! write(6,'(A12,I4,A10,ES12.6)') 'iteration #', iter, ' relres= ', &
        !&           KSPiter%rerr(iter)
        ! now check the second half of iteration
        ! NOTE WE TEST AN IDEA THAT OMIT THE X, WHICH SAVES ANOTHER 
        ! SPMV FOR US
        if (rnorm.lt.btol) then
            KSPiter%rerr(iter) = real(rnorm/bnorm)
            converged=.TRUE.
            KSPiter%failed = .false.
            KSPiter%niter = iter
            exit kspLoop
        end if
        if (rnorm .lt. rnorm0) then
            rnorm0 = rnorm
            ierr = cudaMemCpyAsync(devPtrX0, devPtrX, Arow_d_size, &
                cudaMemcpyDeviceToDevice)
            ierr2 = ierr2 + ierr
        end if
        ierr = cudaDeviceSynchronize()
        ierr2 = ierr2 + ierr
      end do kspLoop
 9527 continue
      if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error synchronizing after iterations ", ierr2
          stop
      end if
      ! for debug
      ! write(6,'(A)') ' Copy solution from GPU to CPU'
      if (.not. converged) then ! solution not found 
          KSPiter%niter=KSPiter%maxit
          ! return the best solution so far
          ierr = cudaMemcpy(xPtr,devPtrX0,Arow_d_size,cudaMemcpyDeviceToHost)
      else ! solution found
          KSPiter%niter=iter
          ! return the last solution
          ierr = cudaMemcpy(xPtr,devPtrX,Arow_d_size,cudaMemcpyDeviceToHost)
      end if
      if (ierr .ne. 0 ) then
          write(6, '(A, I2)') " cudaMemcpy back to host error: ", ierr
          stop
      end if
      ! ierr = cudaHostUnregister(bdotLocPtr)
      ! ierr2 = ierr2 + ierr
      ! if (ierr2 .ne. 0 ) then
      !     write(0, *) " error unpinning the host memory ", ierr2
      !     stop
      ! end if
      ! \activiate lightsaber
      deallocate(xbuff)
      deallocate(bbuff)
      ! clear gpu mem
      ierr = cudaFree(devPtrArow)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrAcol)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrAval)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtriwus)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrLrow)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrLcol)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrLval)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrUrow)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrUcol)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrUval)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrX)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrX0)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrRHS)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrR)
      ierr2 = ierr2 + ierr 
      ierr = cudaFree(devPtrRT)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrP)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrPT)
      ierr2 = ierr2 + ierr 
      ierr = cudaFree(devPtrPH)
      ierr2 = ierr2 + ierr 
      ierr = cudaFree(devPtrS)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrST)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrSH)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrT)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrV)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrAX)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(buffer)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(bufferL)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(bufferU)
      ierr2 = ierr2 + ierr
      ierr = cudaFree(devPtrXbuff)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error during cudafree: ',ierr2
          stop
      end if 
      ierr = cusparseSpSV_destroyDescr(LsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseSpSV_destroyDescr(UsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroyDnVec(vecX)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroyDnVec(vecXloc)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroyDnVec(vecY)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroySpMat(matL)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroySpMat(matU)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroySpMat(matA)
      ierr2 = ierr2 + ierr
      ierr = cublasDestroy(cublasHandle)
      ierr2 = ierr2 + ierr
      ierr = cusparseDestroy(cusparseHandle)
      ierr2 = ierr2 + ierr
      ierr = cudaStreamDestroy(cuStream)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error during cuda handle destruction: ',ierr2
          stop
      end if 
      ierr = cf_resetFlag(device_idx)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error setting device flags: ',ierr2
          stop
      end if 
      return
end subroutine cuBiCGfg ! cuBiCGfg
#endif

#elif defined(HIP)
subroutine hipBiCG(b,x,KSPiter,device_idx,adjt)
  ! Stablized version of BiConjugate Gradient, set up for solving
  ! A x = b 
  ! solves for the interior (edge) field
  !
  ! modified (no hell no) to call the HIP lib to calculate with GPU 
  ! I kept the naming convention from my BiCGStab (see above CPU version)
  ! for instance X -> devPtrX, T -> devPtrT to manipulate the GPU memory
  ! also added the optional adjoint to solve adjoint system A^Tx = b
  ! 
  ! NOTE: BICG actually performs two sub (or half) line searches within a 
  !       iteration, but here we only store the relerr for the second sub
  !       just to be compatitive with QMR
  ! 
  ! interface...........
  ! redefining some of the interfaces for our convenience (locally)
  ! generic routines for vector operations for edge/ face nodes
  ! in a staggered grid
  !
  ! NOTE: this has not been extensively tested - I believe it feels a 
  ! little unstable (dispite the name)...
  ! if you have time reading this, test it!

      use modeloperator3d, only:  AAii, L, U, LH, UH,vOmegaMuSig
      use spoptools
      implicit none
      !  b is right hand side
      complex (kind=prec),intent(in),target,dimension(:)    :: b
      !  solution vector is x ... on input is provided with the initial
      !  guess, on output is the iterate with smallest residual. 
      complex (kind=prec),intent(inout),target,dimension(:) :: x
      type (solverControl_t),intent(inout)                  :: KSPiter
      integer,intent(in)                                    :: device_idx
      logical,intent(in),optional                           :: adjt
    
      ! local variables
      integer                                    :: n, nnz, iter, maxIter
      complex (kind=prec),pointer,dimension(:)   :: iwus

      real    (kind=prec), target                :: rnorm, bnorm, rnorm0, btol
      real    (kind=prec), target                :: xnorm
      complex (kind=prec)                        :: RHO1, ALPHA, BETA, OMEGA
      complex (kind=prec), target                :: RTV,TT,RHO
      real    (kind=SP), target                  :: ctime
      integer                                    :: TRANS
      integer                                    :: Lfillmode,Ldiagtype
      integer                                    :: Ufillmode,Udiagtype
      integer                                    :: ierr,ierr2, idx
      logical                                    :: converged, adjoint
      integer                                    :: interval
      logical                                    :: restart

      ! size to determine corresponding buffer/memory
      integer*8                :: Arowp1_i_size,Arow_d_size
      integer*8                :: Acol_d_size,Annz_i_size,Annz_d_size
      integer*8                :: Mrow_d_size,Lnnz_i_size,Lnnz_d_size
      integer*8                :: Unnz_i_size,Unnz_d_size
      integer(c_int),target    :: zeroLoc
      ! integer(c_int)           :: mbsize
      integer(c_size_t)        :: bsize, lbsize, ubsize

      ! --------------------- pointers to *host* memory ------------------ !   
      type(c_ptr) :: hipblasHandle 
      type(c_ptr) :: hipsparseHandle
      type(c_ptr) :: hipStream
      type(c_ptr) :: hipEvent1
      type(c_ptr) :: hipEvent2
      type(c_ptr) :: matA
      type(c_ptr) :: matL
      type(c_ptr) :: matU
      type(c_ptr) :: vecX
      type(c_ptr) :: vecY 
      type(c_ptr) :: LsolveHandle
      type(c_ptr) :: UsolveHandle
      type(c_ptr) :: ilu_info
      type(c_ptr) :: ArowPtr
      type(c_ptr) :: AcolPtr
      type(c_ptr) :: AvalPtr
      type(c_ptr) :: iwusPtr
      type(c_ptr) :: MrowPtr
      type(c_ptr) :: McolPtr
      type(c_ptr) :: MvalPtr
      type(c_ptr) :: xPtr  
      type(c_ptr) :: bPtr
      type(c_ptr) :: bnormPtr
      type(c_ptr) :: rnormPtr
      type(c_ptr) :: rnormPtr0
      type(c_ptr) :: rhoPtr
      type(c_ptr) :: rtvPtr
      type(c_ptr) :: ttPtr
      type(c_ptr) :: zeroLocPtr
      type(c_ptr) :: timePtr
      ! -------------------- pointers to *device* memory ----------------- ! 
      type(c_ptr) :: devPtrArow
      type(c_ptr) :: devPtrAcol
      type(c_ptr) :: devPtrAval
      type(c_ptr) :: devPtriwus
      type(c_ptr) :: devPtrLrow
      type(c_ptr) :: devPtrLcol
      type(c_ptr) :: devPtrLval
      type(c_ptr) :: devPtrUrow
      type(c_ptr) :: devPtrUcol
      type(c_ptr) :: devPtrUval
      type(c_ptr) :: devPtrX
      type(c_ptr) :: devPtrX0
      type(c_ptr) :: devPtrRHS
      type(c_ptr) :: devPtrR
      type(c_ptr) :: devPtrRT
      type(c_ptr) :: devPtrP
      type(c_ptr) :: devPtrPT
      type(c_ptr) :: devPtrPH
      type(c_ptr) :: devPtrS
      type(c_ptr) :: devPtrST
      type(c_ptr) :: devPtrSH
      type(c_ptr) :: devPtrV
      type(c_ptr) :: devPtrT
      type(c_ptr) :: devPtrAX
      type(c_ptr) :: buffer
      type(c_ptr) :: bufferL
      type(c_ptr) :: bufferU
      
      ! note we don't check the input parameters compatibility 
      ierr2 = 0
      converged = .FALSE.
      zeroLoc = 0
      if (present(adjt)) then 
          ! write(6,'(A)') ' adjt = ', adjt
          adjoint = adjt
          if (adjt) then
              TRANS = HIPSPARSE_OPERATION_TRANSPOSE
              Lfillmode = HIPSPARSE_FILL_MODE_UPPER
              Ufillmode = HIPSPARSE_FILL_MODE_LOWER
              Ldiagtype = HIPSPARSE_DIAG_TYPE_NON_UNIT
              Udiagtype = HIPSPARSE_DIAG_TYPE_UNIT
          else
              TRANS = HIPSPARSE_OPERATION_NON_TRANSPOSE
              Lfillmode = HIPSPARSE_FILL_MODE_LOWER
              Ufillmode = HIPSPARSE_FILL_MODE_UPPER
              Ldiagtype = HIPSPARSE_DIAG_TYPE_UNIT
              Udiagtype = HIPSPARSE_DIAG_TYPE_NON_UNIT
          end if
      else
          adjoint = .FALSE.
          TRANS = HIPSPARSE_OPERATION_NON_TRANSPOSE
          Lfillmode = HIPSPARSE_FILL_MODE_LOWER
          Ufillmode = HIPSPARSE_FILL_MODE_UPPER
          Ldiagtype = HIPSPARSE_DIAG_TYPE_UNIT
          Udiagtype = HIPSPARSE_DIAG_TYPE_NON_UNIT
      end if
      
      ! firstly need to translate the ModEM SP datatypes to simple
      ! CSR plain vectors
      ! note we still keep the Ae = CCGDe + iwuse 
      ! idea as in ModEM SP2, we will test if this will save some time 
      ! (probably not much)
      ! call CSR_R2Cdiag(AAii,VOmegaMuSig,Aii)
      n = AAii%nrow
      nnz = size(AAii%col)
      ArowPtr = c_loc(AAii%row)
      AcolPtr = c_loc(AAii%col)
      AvalPtr = c_loc(AAii%val)
      ! other host pointers
      xPtr  = c_loc(x)  ! x = A \ b
      bPtr  = c_loc(b)  ! b = ones(n,1)
      allocate(iwus(n))
      iwus = VOmegaMuSig*ISIGN*CMPLX(0.0,1.0,8)
      iwusPtr = c_loc(iwus) ! i V omega mu sigma
      bnormPtr = c_loc(bnorm)
      rnormPtr = c_loc(rnorm)
      rnormPtr0 = c_loc(rnorm0)
      rhoPtr = c_loc(RHO)
      rtvPtr = c_loc(rtv)
      ttPtr = c_loc(tt)
      zeroLocPtr = c_loc(zeroLoc)
      ! pointer to store the time 
      ctime = 0.0 
      timePtr = c_loc(ctime)
      ! now remove the Aii matrix structure
      ! call deall_spMatCSR(Aii)
      
      ! get the size of the vector and matrix (need to allocate on the device)
      Arowp1_i_size=sizeof(AAii%row(1:n+1))
      Arow_d_size=sizeof(b(1:n))
      Acol_d_size=sizeof(x(1:n)) ! this is useful if A is not square
      Annz_i_size=sizeof(AAii%col(1:nnz))
      Annz_d_size=sizeof(AAii%val(1:nnz))

      ! select the current gpu device 
      ! note this only works for physical devices (not working for MIG devices)
      ierr = hipSetDevice(device_idx);
      ierr2 = ierr2 + ierr
      ! firstly define the HIP Stream and handles
      ierr = hipStreamCreateWithFlags(hipStream, hipStreamNonBlocking) 
      ierr2 = ierr2 + ierr
      ! initialize the cusparse lib
      ierr = hipsparseCreate(hipsparseHandle)
      ierr2 = ierr2 + ierr
      ierr = hipsparseSetStream(hipsparseHandle,hipStream) 
      ierr2 = ierr2 + ierr
      ! now initialize the cublas lib
      ierr = hipblasCreate(hipblasHandle)
      ierr2 = ierr2 + ierr
      ierr = hipblasSetStream(hipblasHandle,hipStream)
      ierr2 = ierr2 + ierr
      ! and creates two hip events to record the time
      ierr = hipEventCreate(hipEvent1)
      ierr2 = ierr2 + ierr
      ierr = hipEventCreate(hipEvent2)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I4)') 'Error during HIP initialize ',ierr2
          stop
      end if 
      ! record the event before memory manipulation
      ! ierr = hipEventRecord(hipEvent1, hipStream)
      ! ierr2 = ierr2 + ierr
      ! write(6,*) 'Allocating GPU memory'
      ierr = hipMalloc(devPtrX,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrX0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrRHS,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrR,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrRT,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrP,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrPT,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrPH,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrS,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrST,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrSH,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrT,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrV,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrAX,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrAval,Annz_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrAcol,Annz_i_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrArow,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtriwus,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error during GPU MEM allocation: ',ierr2
          stop
      end if 
      ! write(6,*) 'reset GPU memory to all zeros'
      ierr = hipMemset(devPtrX,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrX0,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrRHS,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrR,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrRT,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrP,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrPT,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrPH,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrS,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrST,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrSH,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrT,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrV,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrAX,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrAval,0,Annz_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrAcol,0,Annz_i_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrArow,0,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtriwus,0,Arow_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I3)') 'Error during device memory reseting : ',ierr2
          stop
      end if 
      ! transfer memory over to GPU
      ! write(6,*) 'Transferring memory to GPU'
      ! initialize A
      ierr = hipMemcpyAsync(devPtrArow,ArowPtr,Arowp1_i_size, &
     &        hipMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = hipMemcpyAsync(devPtrAcol,AcolPtr,Annz_i_size, &
     &        hipMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = hipMemcpyAsync(devPtrAval,AvalPtr,Annz_d_size, &
     &        hipMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = hipMemcpyAsync(devPtriwus,iwusPtr,Arow_d_size, &
     &        hipMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ! finally need to create the SpMatDescr 
      ierr = hipsparseCreateCsr(matA, n, n, nnz, devPtrArow, devPtrAcol, &
     &       devPtrAval, HIPSPARSE_INDEX_32I, HIPSPARSE_INDEX_32I, &
     &       HIPSPARSE_INDEX_BASE_ONE, HIP_R_64F) 
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error assembling system matrix ", ierr2
          stop
      end if
      ! initialize rhs and x
      ierr = hipMemcpyAsync(devPtrX,xPtr,Arow_d_size, &
     &        hipMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = hipMemcpyAsync(devPtrRHS,bPtr,Arow_d_size, &
     &        hipMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = hipDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error during hip memcpy ", ierr2
          stop
      end if
      ! now deallocate temp array
      if ( c_associated(iwusPtr) ) then
          ! this is a little tricky as iwus is associated with
          ! iwusPtr in C
          nullify(iwus)
          call cf_free(iwusPtr)
      end if
      ! firstly we need to create two dense vectors 
      ! as the user-hostile developers in Nvidia think of a new idea 
      ! to mess up the interfaces
      ierr = hipsparseCreateDnVec(vecX, n, devPtrX, HIP_C_64F)
      ierr2 = ierr2 + ierr
      ierr = hipsparseCreateDnVec(vecY, n, devPtrR, HIP_C_64F)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, *) " error creating the dense Vecs "
          stop
      end if
      ierr = hipsparseSpMV_bufferSize_cmplx(hipsparseHandle,              &
     &        TRANS, C_ONE, matA, vecX, C_ONE, vecY,                   &
     &        HIP_C_64F, HIPSPARSE_SPMV_CSR_ALG1, bsize)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error estimating buffer for spMV ",  ierr2
          stop
      end if
      ! and finally (re)allocate the buffer
      ierr = hipMalloc(buffer, bsize)
      ierr2 = ierr2 + ierr
      ! now let's deal with L and U
      ! write(6,'(A)') ' Setup L and U preconditioners on GPU'
      ! L first
      if (adjoint) then
          nnz = size(LH%col)
          ! sizes (row size is the same as A
          Lnnz_i_size=sizeof(LH%col(1:nnz))
          Lnnz_d_size=sizeof(LH%val(1:nnz))
          MrowPtr = c_loc(LH%row)
          McolPtr = c_loc(LH%col)
          MvalPtr = c_loc(LH%val)
      else
          nnz = size(L%col)
          ! sizes (row size is the same as A
          Lnnz_i_size=sizeof(L%col(1:nnz))
          Lnnz_d_size=sizeof(L%val(1:nnz))
          MrowPtr = c_loc(L%row)
          McolPtr = c_loc(L%col)
          MvalPtr = c_loc(L%val)
      end if
      ! initialize L in GPU memory
      ierr = hipMalloc(devPtrLval,Lnnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrLcol,Lnnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrLrow,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! set everything to zeros
      ierr = hipMemset(devPtrLval,0,Lnnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrLcol,0,Lnnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrLrow,0,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! now copy the values to device
      ierr = hipMemcpyAsync(devPtrLval,MvalPtr,Lnnz_d_size, &
     &        hipMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = hipMemcpyAsync(devPtrLcol,McolPtr,Lnnz_i_size, &
     &        hipMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = hipMemcpyAsync(devPtrLrow,MrowPtr,Arowp1_i_size, &
     &        hipMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ! finish synchronizing
      ierr = hipDeviceSynchronize()
      ierr2 = ierr2 + ierr
      ! finally need to create the SpMatDescr 
      ierr = hipsparseCreateCsr(matL, n, n, nnz, devPtrLrow, devPtrLcol, &
     &       devPtrLval, HIPSPARSE_INDEX_32I, HIPSPARSE_INDEX_32I, &
     &       HIPSPARSE_INDEX_BASE_ONE, HIP_C_64F) 
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error assembling L matrix ", ierr2
          stop
      end if
      ! now start messing up with SpSV
      ! now estimate the buffersize needed by SpSV (Lsolve)
      ! solves y in L*y = a*x (if a=1)
      ierr = hipsparseSpMatSetAttribute(matL, HIPSPARSE_SPMAT_FILL_MODE, &
    &        Lfillmode, 4)
      ierr2 = ierr2 + ierr
      ierr = hipsparseSpMatSetAttribute(matL, HIPSPARSE_SPMAT_DIAG_TYPE, &
    &        Ldiagtype, 4)
      ierr2 = ierr2 + ierr
      ! still need to establish a context handler for Lsolve
      ierr = hipsparseSpSV_createDescr(LsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = hipsparseSpSV_bufferSize_dcmplx(hipsparseHandle,             &
     &        HIPSPARSE_OPERATION_NON_TRANSPOSE, C_ONE, matL, vecX, vecY,&
     &        HIP_C_64F,HIPSPARSE_SPSV_ALG_DEFAULT, LsolveHandle, lbsize)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error estimating buffer for Lsolve ",  ierr2
          stop
      end if
      ierr = hipMalloc(bufferL, lbsize)
      ierr2 = ierr2 + ierr
      ierr = hipsparseSpSV_analysis_dcmplx(hipsparseHandle,         &
     &        HIPSPARSE_OPERATION_NON_TRANSPOSE, C_ONE, matL, vecX, vecY, &
     &        HIP_C_64F,HIPSPARSE_SPSV_ALG_DEFAULT, LsolveHandle, bufferL)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error analysing L matrix",  ierr2
          stop
      end if
      ! write(6,*) 'spsv Lsolve buffersize = ', lbsize
      ! then U 
      if (adjoint) then
          nnz = size(UH%col)
          ! sizes (row size is the same as A
          Unnz_i_size=sizeof(UH%col(1:nnz))
          Unnz_d_size=sizeof(UH%val(1:nnz))
          MrowPtr = c_loc(UH%row)
          McolPtr = c_loc(UH%col)
          MvalPtr = c_loc(UH%val)
      else
          nnz = size(U%col)
          ! sizes (row size is the same as A
          Unnz_i_size=sizeof(U%col(1:nnz))
          Unnz_d_size=sizeof(U%val(1:nnz))
          MrowPtr = c_loc(U%row)
          McolPtr = c_loc(U%col)
          MvalPtr = c_loc(U%val)
      end if
      ! initialize U in GPU memory
      ierr = hipMalloc(devPtrUval,Unnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrUcol,Unnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = hipMalloc(devPtrUrow,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! set everything to zeros
      ierr = hipMemset(devPtrUval,0,Unnz_d_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrUcol,0,Unnz_i_size)
      ierr2 = ierr2 + ierr
      ierr = hipMemset(devPtrUrow,0,Arowp1_i_size)
      ierr2 = ierr2 + ierr
      ! now copy the values to device
      ierr = hipMemcpyAsync(devPtrUval,MvalPtr,Unnz_d_size, &
     &        hipMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = hipMemcpyAsync(devPtrUcol,McolPtr,Unnz_i_size, &
     &        hipMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ierr = hipMemcpyAsync(devPtrUrow,MrowPtr,Arowp1_i_size, &
     &        hipMemcpyHostToDevice)
      ierr2 = ierr2 + ierr
      ! finally need to create the SpMatDescr 
      ierr = hipsparseCreateCsr(matU, n, n, nnz, devPtrUrow, devPtrUcol, &
     &       devPtrUval, HIPSPARSE_INDEX_32I, HIPSPARSE_INDEX_32I, &
     &       HIPSPARSE_INDEX_BASE_ONE, HIP_C_64F) 
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error creating U matrix",  ierr2
          stop
      end if
      ! now estimate the buffersize needed by SpSV (Usolve)
      ! solves y in U*y = a*x (if a=1)
      ierr = hipsparseSpMatSetAttribute(matU, HIPSPARSE_SPMAT_FILL_MODE, &
    &        Ufillmode, 4)
      ierr2 = ierr2 + ierr
      ierr = hipsparseSpMatSetAttribute(matU, HIPSPARSE_SPMAT_DIAG_TYPE, &
    &        Udiagtype, 4)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error setting attribute for U ",  ierr2
          stop
      end if
      ! still need to establish a context handler for Usolve
      ierr = hipsparseSpSV_createDescr(UsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = hipsparseSpSV_bufferSize_dcmplx(hipsparseHandle,             &
     &        HIPSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matU, vecX, vecY,&
     &        HIP_C_64F,HIPSPARSE_SPSV_ALG_DEFAULT, UsolveHandle, ubsize)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A,I2)') " error estimating buffer for Usolve ",  ierr2
          stop
      end if
      ! write(6,*) 'spsv Usolve buffersize = ', ubsize
      ierr = hipMalloc(bufferU, ubsize)
      ierr2 = ierr2 + ierr
      ierr = hipsparseSpSV_analysis_dcmplx(hipsparseHandle,         &
     &        HIPSPARSE_OPERATION_NON_TRANSPOSE, C_ONE, matU, vecX, vecY, &
     &        HIP_C_64F,HIPSPARSE_SPSV_ALG_DEFAULT, UsolveHandle, bufferU)
      ierr2 = ierr2 + ierr
      !write(6, '(A, I8)') " allocated SpMV/SpSV buffer on GPU (MB) ",&
     !&            bsize/1024/1024
      ! see how long it takes for the intialization part
      ! ierr = hipEventRecord(hipEvent2, hipStream)
      ! ierr2 = ierr2 + ierr
      ! ierr = hipDeviceSynchronize()
      ! ierr = hipEventElapsedTime(timePtr, hipEvent1, hipEvent2)
      ! ierr2 = ierr2 + ierr
      ! write(6,*) 'initial GPU memory cpy time is ', ctime, 'ms'
      ! if (ierr2 .ne. 0 ) then
      !     write(6, *) " error recording the time ", ierr2
      !     stop
      ! end if
      ! now Compute the initial residual
      ! write(6,*) 'Computing initial residual'
      ! setup the vecX and vecY
      ierr = hipsparseDnVecSetValues(vecX, devPtrX)
      ierr2 = ierr2 + ierr
      ! R = Diag(iwus)*X <-- diagonal
      call kernelc_hadac(devPtrX, devPtriwus, devPtrR, n, cuStream)
      ierr = hipsparseDnVecSetValues(vecY, devPtrR)
      ierr2 = ierr2 + ierr
      ! now calculate the y = Ax
      ! R = CC*X + Diag(iwus)*X
      ! thank God (or whatever deserves it) that SpMV does not 
      ! need to be analyzed
      ierr = hipsparseSpMV_cmplx(hipsparseHandle,TRANS,                  &
     &        C_ONE, matA, vecX, C_ONE, vecY, HIP_C_64F,             &
     &        HIPSPARSE_SPMV_CSR_ALG1, buffer)
      ierr2 = ierr2 + ierr
      ! and get the values back 
      ierr = hipsparseDnVecGetValues(vecY, devPtrR)
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, *) " error with SpMV operation ", ierr2
          stop
      end if
      ! bnorm = nrm2(b)
      ierr = hipblasZnrm2(hipblasHandle,n,devPtrRHS,1,bnormPtr)
      ierr2 = ierr2 + ierr
      if (isnan(bnorm)) then
      ! this usually means an inadequate model, in which case Maxwell's fails
          write(6,*) 'Error: b in BICG contains NaNs; exiting...'
          stop
      else if ( bnorm .eq. 0.0) then ! zero rhs -> zero solution
          write(6,*) 'Warning: b in BICG has all zeros, returning zero &
     &        solution'
          x = b 
          KSPiter%niter=1
          KSPiter%failed=.false.
          KSPiter%rerr(1)=0.0
          converged = .true.
          goto 9527 
      endif
      ! R = -Ax 
      ! ierr = hipblasZscal(hipblasHandle,n,C_MINUSONE,devPtrR,1)
      ! ierr2 = ierr2 + ierr
      ! R = b - Ax
      ! ierr = hipblasZaxpy(hipblasHandle,n,C_ONE,devPtrRHS,1,devPtrR,1)
      ! ierr2 = ierr2 + ierr
      ! TEST: R = b - Ax with an all-in-one kernel of xpby
      call kernelc_xpbyc(devPtrRHS, C_MINUSONE, devPtrR, n, cuStream)
      ! rnorm = nrm2(b - Ax)
      ierr = hipblasZnrm2(hipblasHandle,n,devPtrR,1,rnormPtr)
      ierr2 = ierr2 + ierr
      rnorm0 = rnorm
      if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error during residual estimation", ierr2
          stop
      end if
      ! write(6,'(A, ES12.6)') ' initial relative residual = ', rnorm/bnorm
      !================= Now start configuring the iteration ================!
      btol = KSPiter%tol*bnorm
      maxIter = KSPiter%maxIt
      KSPiter%rerr(1) = real(rnorm/bnorm)
      if ( rnorm .le. btol ) then ! the first guess is already good enough
         ! returning
          write(6, *) " The first guess is good enough, exiting..."
          iter=1
          KSPiter%failed=.false.
          converged = .true.
          goto 9527 
      end if
      ! intial the parameters for restarting
      restart = .FALSE.
      interval = 120
      ! write(6,'(A,I4)') ' maxiter = ', maxiter
      ! RHO = C_ONE
      kspLoop: do iter = 1,maxIter
        ! write(6,'(A, I4)') ' KSP iteration #', iter
        if (mod(iter, interval) .eq. 1) then ! hard coded here
            restart = .TRUE.
        end if
        if (restart) then 
            ! restart the iteration (to steepest decend) every interval times
            ! RT = R
            ierr = hipblasZcopy(hipblasHandle,n,devPtrR,1,devPtrRT,1)
            ierr2 = ierr2 + ierr
            ! current and previous RHO, RHO0 should be one
            RHO1 = C_ONE
            ! RHO = dot(RT,R)
            ierr = hipblasZdot(hipblasHandle,n,devPtrRT,1,devPtrR,1,rhoPtr)
            ierr2 = ierr2 + ierr
            ! P = R
            ierr = hipblasZcopy(hipblasHandle,n,devPtrR,1,devPtrP,1)
            ierr2 = ierr2 + ierr
            OMEGA = C_ONE
            restart = .FALSE.
            if (ierr2 .ne. 0 ) then
              write(6, '(A, I2)') " Error steering search direction ", ierr2
              stop
            end if
        else
            ! current and previous RHO, RHO0 should be one
            RHO1 = RHO
            ! RHO = dot(RT,R)
            ierr = hipblasZdot(hipblasHandle,n,devPtrRT,1,devPtrR,1,rhoPtr)
            ierr2 = ierr2 + ierr
            BETA=(RHO/RHO1)*(ALPHA/OMEGA)
            if (BETA .eq. 0.0) then !bad beta
                converged=.FALSE.
                KSPiter%failed = .true.
                exit kspLoop
            end if
            ! P = R + BETA * (P - OMEGA * V)
            ! P = P - OMEGA * V
            ! ierr = hipblasZaxpy(hipblasHandle,n,-(OMEGA),devPtrV,1,devPtrP,1)
            ! ierr2 = ierr2 + ierr
            ! P = BETA * P
            ! ierr = hipblasZscal(hipblasHandle,n,BETA,devPtrP,1)
            ! ierr2 = ierr2 + ierr
            ! P = R + P
            ! ierr = hipblasZaxpy(hipblasHandle,n,C_ONE,devPtrR,1,devPtrP,1)
            ! ierr2 = ierr2 + ierr
            ! TEST: P = R + BETA * P  with an all-in-one kernel 
            ! call kernelc_xpbyc(devPtrR, BETA, devPtrP, n, cuStream)
            ! TEST: update P with an all-in-one kernel 
            call kernelc_update_pc(devPtrR, devPtrV, BETA, OMEGA, devPtrP, n,&
       &            cuStream)
            if (ierr2 .ne. 0 ) then
              write(6, '(A, I2)') " Error steering search direction ", ierr2
              stop
            end if
        end if
        ! record - start of two SPSVs
        ! ierr = hipEventRecord(hipEvent1, hipStream)
        ! ierr2 = ierr2 + ierr
        ! ============== first half of the conjugate iteration ============= !
        ! L solve --> L*PT = P
        ! write(6,'(A)') ' Lsolve '
        ! still need to use the vecX/Y types
        ierr = hipsparseDnVecSetValues(vecX, devPtrP)
        ierr2 = ierr2 + ierr
        ierr = hipsparseDnVecSetValues(vecY, devPtrPT)
        ierr2 = ierr2 + ierr
        ierr = hipsparseSpSV_solve_dcmplx(hipsparseHandle,            &
     &         HIPSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matL, vecX, vecY, &
     &         HIP_C_64F,HIPSPARSE_SPSV_ALG_DEFAULT, LsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = hipsparseDnVecGetValues(vecY, devPtrPT)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I2)') " Error during Lsolve ", ierr2
          stop
        end if
        ! U solve --> U*PH = PT
        ! write(6,'(A)') ' Usolve '
        ! still need to use the vecX/Y types
        ierr = hipsparseDnVecSetValues(vecX, devPtrPT)
        ierr2 = ierr2 + ierr
        ierr = hipsparseDnVecSetValues(vecY, devPtrPH)
        ierr2 = ierr2 + ierr
        ierr = hipsparseSpSV_solve_dcmplx(hipsparseHandle,            &
     &         HIPSPARSE_OPERATION_NON_TRANSPOSE, C_ONE, matU, vecX, vecY, &
     &         HIP_C_64F,HIPSPARSE_SPSV_ALG_DEFAULT, UsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = hipsparseDnVecGetValues(vecY, devPtrPH)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I2)') " Error during Usolve ", ierr2
          stop
        end if
        ! record - end of two SPSVs
        ! ierr = hipEventRecord(hipEvent2, hipStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = hipDeviceSynchronize()
        ! ierr = hipEventElapsedTime(timePtr, hipEvent1, hipEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '2 SPSV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !    write(6,'(A,I4)') " error recording the time ", ierr2
        !     stop
        ! end if
        ! record - start of one SPMV
        ! ierr = hipEventRecord(hipEvent1, hipStream)
        ! ierr2 = ierr2 + ierr
        ! V = A*PH
        ! write(6,'(A)') ' Axpy '
        ! still need to use the vecX/Y types
        ierr = hipsparseDnVecSetValues(vecX, devPtrPH)
        ierr2 = ierr2 + ierr
        ! V = Diag(iwus)*PH
        call kernelc_hadac(devPtrPH, devPtriwus, devPtrV, n, cuStream)
        ierr = hipsparseDnVecSetValues(vecY, devPtrV)
        ierr2 = ierr2 + ierr
        ! V = CC*PH + Diag(iwus)*PH
        ierr = hipsparseSpMV_cmplx(hipsparseHandle,TRANS, &
     &         C_ONE, matA, vecX, C_ONE, vecY, HIP_C_64F,  &
     &         HIPSPARSE_SPMV_CSR_ALG1, buffer)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = hipsparseDnVecGetValues(vecY, devPtrV)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
            write(6,'(A,I4)') " Error calculating A*PH", ierr2
            stop
        end if
        ! record - end of one SPMV
        ! ierr = hipEventRecord(hipEvent2, hipStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = hipDeviceSynchronize()
        ! ierr = hipEventElapsedTime(timePtr, hipEvent1,hipEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '1 SPMV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !    write(6, *) " error recording the time ", ierr2
        !    stop
        ! end if
        ! now check the residual of first half of iteration
        ! rtv = dot(RT,V)
        ierr = hipblasZdot(hipblasHandle,n,devPtrRT,1,devPtrV,1,rtvPtr)
        ierr2 = ierr2 + ierr
        ALPHA = RHO/rtv
        ! write(6,*) 'alpha = ', ALPHA
        if (ALPHA .eq. 0.0) then !bad alpha
            KSPiter%failed = .true.
            converged=.FALSE.
            exit kspLoop
        end if
        ! x = x + ALPHA * PH
        ierr = hipblasZaxpy(hipblasHandle,n,ALPHA,devPtrPH,1,devPtrX,1)
        ierr2 = ierr2 + ierr
        ! S = R
        ierr = hipMemcpyAsync(devPtrS, devPtrR, Arow_d_size, &
       &       hipMemcpyDeviceToDevice)
        ierr2 = ierr2 + ierr
        ! S = R - ALPHA * V
        ierr = hipblasZaxpy(hipblasHandle,n,-(ALPHA),devPtrV,1,devPtrS,1)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error estimating the first Conj res. ", ierr2
          stop
        end if
        ! ============== second half of the conjugate iteration ============= !
        ! record - start of two SPSV
        ! ierr = hipEventRecord(hipEvent1, hipStream)
        ! ierr2 = ierr2 + ierr
        ! L solve --> L*ST = S
        ! write(6,'(A)') ' Lsolve '
        ! still need to use the vecX/Y types
        ierr = hipsparseDnVecSetValues(vecX, devPtrS)
        ierr2 = ierr2 + ierr
        ierr = hipsparseDnVecSetValues(vecY, devPtrST)
        ierr2 = ierr2 + ierr
    !    ierr = hipsparseSpSV_analysis_dcmplx(hipsparseHandle,        &
    ! &         HIPSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matL, vecX, vecY, &
    ! &         HIP_C_64F,HIPSPARSE_SPSV_ALG_DEFAULT, LsolveHandle, bufferL)
    !    ierr2 = ierr2 + ierr
        ierr = hipsparseSpSV_solve_dcmplx(hipsparseHandle,           &
     &         HIPSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matL, vecX, vecY, &
     &         HIP_C_64F,HIPSPARSE_SPSV_ALG_DEFAULT, LsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = hipsparseDnVecGetValues(vecY, devPtrST)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I2)') " Error during Lsolve ", ierr2
          stop
        end if
        
        ! U solve --> U*SH = ST
        ! write(6,'(A)') ' Usolve '
        ! still need to use the vecX/Y types
        ierr = hipsparseDnVecSetValues(vecX, devPtrST)
        ierr2 = ierr2 + ierr
        ierr = hipsparseDnVecSetValues(vecY, devPtrSH)
        ierr2 = ierr2 + ierr
    !    ierr = hipsparseSpSV_analysis_dcmplx(hipsparseHandle,         &
    ! &         HIPSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matU, vecX, vecY, &
    ! &         HIP_C_64F,HIPSPARSE_SPSV_ALG_DEFAULT, UsolveHandle, bufferU)
    !    ierr2 = ierr2 + ierr
        ierr = hipsparseSpSV_solve_dcmplx(hipsparseHandle,            &
     &         HIPSPARSE_OPERATION_NON_TRANSPOSE,C_ONE, matU, vecX, vecY, &
     &         HIP_C_64F,HIPSPARSE_SPSV_ALG_DEFAULT, UsolveHandle)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = hipsparseDnVecGetValues(vecY, devPtrSH)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I2)') " Error during Usolve ", ierr2
          stop
        end if

        ! record - end of 2 SPSVs
        ! ierr = hipEventRecord(hipEvent2, hipStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = hipDeviceSynchronize()
        ! ierr = hipEventElapsedTime(timePtr, hipEvent1,hipEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '2 SPSV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !     write(6, *) " error recording the time ", ierr2
        !     stop
        ! end if
        ! record - start of 1 SPMV
        ! ierr = hipEventRecord(hipEvent1, hipStream)
        ! ierr2 = ierr2 + ierr
        ! T = A*SH 
        ! write(6,'(A)') ' Axpy '
        ! still need to use the vecX/Y types
        ierr = hipsparseDnVecSetValues(vecX, devPtrSH)
        ierr2 = ierr2 + ierr
        ! T = Diag(iwus)*SH
        call kernelc_hadac(devPtrSH, devPtriwus, devPtrT, n, cuStream)
        ierr = hipsparseDnVecSetValues(vecY, devPtrT)
        ierr2 = ierr2 + ierr
        ! T = A*SH + Diag(iwus)*SH
        ierr = hipsparseSpMV_cmplx(hipsparseHandle,TRANS, &
       &        C_ONE, matA, vecX, C_ONE, vecY, HIP_C_64F,                 & 
       &        HIPSPARSE_SPMV_CSR_ALG1, buffer)
        ierr2 = ierr2 + ierr
        ! and get the values back 
        ierr = hipsparseDnVecGetValues(vecY, devPtrT)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
            write(6,'(A,I4)') " Error calculating A*SH", ierr2
            stop
        end if

        ! record - end of 1 SPMV 
        ! ierr = hipEventRecord(hipEvent2, hipStream)
        ! ierr2 = ierr2 + ierr
        ! ierr = hipDeviceSynchronize()
        ! ierr = hipEventElapsedTime(timePtr, hipEvent1,hipEvent2)
        ! ierr2 = ierr2 + ierr
        ! write(6,*) '2 SPMV time = ', ctime, 'ms'
        ! if (ierr2 .ne. 0 ) then
        !     write(6, *) " error recording the time ", ierr2
        !     stop
        ! end if

        ! now check the second half of iteration
        ! calculate the residual norm for the second half of iteration 
        ! tt = dot(T,T)
        ierr = hipblasZdot(hipblasHandle,n,devPtrT,1,devPtrT,1,ttPtr)
        ierr2 = ierr2 + ierr
        ! rtv = dot(T,S)
        ierr = hipblasZdot(hipblasHandle,n,devPtrT,1,devPtrS,1,rtvPtr)
        ierr2 = ierr2 + ierr
        OMEGA = rtv/tt
        ! write(6,*) 'omega = ', OMEGA
        if (OMEGA .eq. 0.0) then !bad omega
            KSPiter%failed = .true.
            converged=.false.
            exit kspLoop
        end if
        ! x = x + OMEGA * SH
        ierr = hipblasZaxpy(hipblasHandle,n,OMEGA,devPtrSH,1,devPtrX,1)
        ierr2 = ierr2 + ierr
        ! R = S
        ierr = hipMemcpyAsync(devPtrR, devPtrS, Arow_d_size, &
       &       hipMemcpyDeviceToDevice)
        ierr2 = ierr2 + ierr
        ! R = S - OMEGA * T
        ierr = hipblasZaxpy(hipblasHandle,n,-(OMEGA),devPtrT,1,devPtrR,1)
        ierr2 = ierr2 + ierr
        ierr = hipDeviceSynchronize()
        ierr2 = ierr2 + ierr
        ! early second half convergence check (norm of R)
        ! rnorm = norm(R)
        ierr = hipblasZnrm2(hipblasHandle,n,devPtrR,1,rnormPtr)
        ierr2 = ierr2 + ierr
        if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error estimating the second Conj res. ", ierr2
          stop
        end if
        KSPiter%rerr(iter) = rnorm/bnorm
        ! write(6,'(A12,I4,A10,ES12.6)') 'iteration #', iter, ' relres= ', &
        !&           KSPiter%rerr(iter)
        ! now check the second half of iteration
        ! NOTE WE TEST AN IDEA THAT OMIT THE X, WHICH SAVES ANOTHER 
        ! SPMV FOR US
        if (rnorm.lt.btol) then
            ! double check if the residual is really less than tol
            ! R = A*X
            ! still need to use the vecX/Y types
            ierr = hipsparseDnVecSetValues(vecX, devPtrX)
            ierr2 = ierr2 + ierr
            ! R = Diag(iwus)*X
            call kernelc_hadac(devPtrX, devPtriwus, devPtrR, n, cuStream)
            ierr = hipsparseDnVecSetValues(vecY, devPtrR)
            ierr2 = ierr2 + ierr
            ! R = CC*X + Diag(iwus)*X
            ierr = hipsparseSpMV_cmplx(hipsparseHandle,TRANS, &
     &          C_ONE, matA, vecX, C_ONE, vecY, HIP_C_64F,  &
     &          HIPSPARSE_SPMV_CSR_ALG1, buffer)
            ierr2 = ierr2 + ierr
            ! and get the values back
            ierr = hipsparseDnVecGetValues(vecY, devPtrR)
            ierr2 = ierr2 + ierr
            if (ierr2 .ne. 0 ) then
                write(6,'(A,I4)') " Error calculating A*X", ierr2
                stop
            end if
            ! R = -Ax
            ! ierr = hipblasZscal(hipblasHandle,n,C_MINUSONE,devPtrR,1)
            ! ierr2 = ierr2 + ierr
            ! R = b - Ax
            ! ierr = hipblasZaxpy(hipblasHandle,n,C_ONE,devPtrRHS,1,devPtrR,1)
            ! ierr2 = ierr2 + ierr
            ! TEST: R = b - Ax with an all-in-one kernel of xpby
            call kernelc_xpbyc(devPtrRHS, C_MINUSONE, devPtrR, n, cuStream)
            ! rnorm = norm(R)
            ierr = hipblasZnrm2(hipblasHandle,n,devPtrR,1,rnormPtr)
            ierr2 = ierr2 + ierr
            if (rnorm.lt.btol) then
                KSPiter%rerr(iter) = real(rnorm/bnorm)
                converged=.TRUE.
                KSPiter%failed = .false.
                KSPiter%niter = iter
                exit kspLoop
            end if
        end if
        if (rnorm .lt. rnorm0) then
            rnorm0 = rnorm
            ierr = hipMemcpyAsync(devPtrX0, devPtrX, Arow_d_size, &
                hipMemcpyDeviceToDevice)
            ierr2 = ierr2 + ierr
        end if
      end do kspLoop
 9527 continue
      ierr = hipDeviceSynchronize()
      ierr2 = ierr2 + ierr
      if (ierr2 .ne. 0 ) then
          write(6, '(A, I4)') " Error synchronizing after iterations ", ierr2
          stop
      end if
      ! write(6,'(A)') ' Copy solution from GPU to CPU'
      if (.not. converged) then ! solution not found 
          KSPiter%niter=KSPiter%maxit
          ! return the best solution so far
          ierr = hipMemcpy(xPtr,devPtrX0,Arow_d_size,hipMemcpyDeviceToHost)
      else ! solution found
          KSPiter%niter=iter
          ! return the last solution
          ierr = hipMemcpy(xPtr,devPtrX,Arow_d_size,hipMemcpyDeviceToHost)
      end if
      if (ierr .ne. 0 ) then
          write(6, '(A, I2)') " hipMemcpy back to host error: ", ierr
          stop
      end if
      ! \activiate lightsaber
      ! clear gpu mem
      ierr = hipFree(devPtrArow)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrAcol)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrAval)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtriwus)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrLrow)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrLcol)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrLval)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrUrow)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrUcol)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrUval)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrX)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrX0)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrRHS)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrR)
      ierr2 = ierr2 + ierr 
      ierr = hipFree(devPtrRT)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrP)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrPT)
      ierr2 = ierr2 + ierr 
      ierr = hipFree(devPtrPH)
      ierr2 = ierr2 + ierr 
      ierr = hipFree(devPtrS)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrST)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrSH)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrT)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrV)
      ierr2 = ierr2 + ierr
      ierr = hipFree(devPtrAX)
      ierr2 = ierr2 + ierr
      ierr = hipFree(buffer)
      ierr2 = ierr2 + ierr
      ierr = hipFree(bufferL)
      ierr2 = ierr2 + ierr
      ierr = hipFree(bufferU)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error during hip free: ',ierr2
          stop
      end if 
      ierr = hipsparseSpSV_destroyDescr(LsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = hipsparseSpSV_destroyDescr(UsolveHandle)
      ierr2 = ierr2 + ierr
      ierr = hipsparseDestroyDnVec(vecX)
      ierr2 = ierr2 + ierr
      ierr = hipsparseDestroyDnVec(vecY)
      ierr2 = ierr2 + ierr
      ierr = hipsparseDestroySpMat(matL)
      ierr2 = ierr2 + ierr
      ierr = hipsparseDestroySpMat(matU)
      ierr2 = ierr2 + ierr
      ierr = hipsparseDestroySpMat(matA)
      ierr2 = ierr2 + ierr
      ierr = hipblasDestroy(hipblasHandle)
      ierr2 = ierr2 + ierr
      ierr = hipsparseDestroy(hipsparseHandle)
      ierr2 = ierr2 + ierr
      ierr = hipStreamDestroy(hipStream)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error during hip handle destruction: ',ierr2
          stop
      end if 
      ierr = cf_resetFlag(device_idx)
      ierr2 = ierr2 + ierr
      if (ierr2.ne.0) then
          write(6,'(A, I2)') 'Error setting device flags: ',ierr2
          stop
      end if 
      return
end subroutine hipBiCG ! hipBiCG

#endif

end module solver ! spsolver
