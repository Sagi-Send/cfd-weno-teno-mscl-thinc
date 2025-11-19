MODULE types_vars
    
    ! This module defines the KIND types of all the variables used in the code:
    ! I4B, I2B and I1B for integer variables, SP and DP for real variables (and
    ! SPC and DPC for corresponding complex cases), and LGT for the default
    ! logical type. This follows the convention used the Numerical Recipes for
    ! Fortran 90 types module 'nrtype', pp. 1361
    !
    ! Symbolic names for kind types of 4-, 2- and 1-byte integers:
    
    INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
    INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
    INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
    
    ! Symbolic names for kind types of single- and double-precison reals
    
    INTEGER, PARAMETER :: SP = KIND(1.0)
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
    
    ! Symbolic names for kind types of single- and double-precison complex
    
    INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
    INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
    
    ! Symbolic name for kind type of default logical
    
    INTEGER, PARAMETER :: LOGIC = KIND(.true.)
    
    ! Frequently used mathematical constants (with precision to spare)
    
    REAL(SP), PARAMETER :: one   = 1.0_sp
    REAL(SP), PARAMETER :: zero  = 0.0_sp
    REAL(SP), PARAMETER :: pi    = 3.141592653589793238462643383279502884197_sp
    REAL(SP), PARAMETER :: pio2  = 1.57079632679489661923132169163975144209858_sp
    REAL(SP), PARAMETER :: twopi = 6.283185307179586476925286766559005768394_sp
    
    REAL(DP), PARAMETER :: one_d   = 1.0_dp
    REAL(DP), PARAMETER :: two_d   = 2.0_dp
    REAL(DP), PARAMETER :: zero_d  = 0.0_dp
    REAL(DP), PARAMETER :: pi_d    = 3.141592653589793238462643383279502884197_dp
    REAL(DP), PARAMETER :: pio2_d  = 1.57079632679489661923132169163975144209858_dp
    REAL(DP), PARAMETER :: twopi_d = 6.283185307179586476925286766559005768394_dp
    REAL(DP), PARAMETER :: onethird_d = 1.0_dp/3.0_dp
    REAL(DP), PARAMETER :: twothird_d = 2.0_dp/3.0_dp
    REAL(DP), PARAMETER :: onefourth_d = 1.0_dp/4.0_dp
    REAL(DP), PARAMETER :: threefourth_d = 3.0_dp/4.0_dp
    REAL(DP), PARAMETER :: onefifth_d= 1.0_dp/5.0_dp
    REAL(DP), PARAMETER :: onesixth_d= 1.0_dp/6.0_dp
    REAL(DP), PARAMETER :: half_d = 0.5_dp
    
END MODULE types_vars

!**************************************************************

MODULE subroutines
    USE types_vars
    IMPLICIT NONE
    
    INTEGER     :: nx, ny, nc_x, nc_y, ng, nv, nt, ntsteps, nf, kf
    REAL(DP)    :: xl, xr, yl, yr, t, tfinal, dt, dx, dy, dtodx, dxodt, dtody, dyodt, gam=1.4d0, gm1, gm1i, wavespeed
    REAL(DP)    :: Mach, Pr, cp, ka, Re, mu, R
    
    CONTAINS
    
    !**************************************************************
    
    SUBROUTINE grid(x, y)
        IMPLICIT NONE
        INTEGER :: i
        REAL(DP),INTENT(OUT) :: x(:), y(:)
        
        ! This subroutine generates a finite-volume grid in 1D
        ! Specifically, x(i) contains the coordinates of the 
        ! cell centers starting at h/2
        ! |---x---|---x---|.....|---x---|---x---|
        
        do i = ng + 1, nx - ng
            x(i) = xl + float(i-1)*dx + 0.5d0*dx
        end do
        do i = ng + 1, ny - ng
            y(i) = yl + float(i-1)*dy + 0.5d0*dy
        end do
        
    END SUBROUTINE grid
    
    !**************************************************************
    
    SUBROUTINE init_shock_tube(x, rho, u_x, u_y, p, T)
        IMPLICIT NONE
        INTEGER :: i, j
        REAL(DP), INTENT(IN) :: x(:)
        REAL(DP), INTENT(OUT) :: rho(:,:), u_x(:,:), u_y(:,:), p(:,:),T(:,:)
    
        DO j = 1, ny
          DO i = 1, nx              
              IF (x(i) <= 0.5_dp) THEN
                  rho(i,j) = 120.0_dp
                  u_x(i,j) = 0.0_dp
                  u_y(i,j) = 0.0_dp
                  p(i,j) = 120.0_dp / gam
              ELSE
                  rho(i,j) = 1.2_dp
                  u_x(i,j) = 0.0_dp
                  u_y(i,j) = 0.0_dp
                  p(i,j) = 1.2_dp / gam
              END IF
              T(i,j) = P(i,j) / (rho(i,j) * R)
          END DO
      END DO
  
    END SUBROUTINE init_shock_tube
    
    !**************************************************************
    
    SUBROUTINE solvec(rho, u_x, u_y, p, eu)
        IMPLICIT NONE
        REAL(DP), INTENT(IN)  :: rho(:,:), u_x(:,:), u_y(:,:), p(:,:)
        REAL(DP), INTENT(OUT) :: eu(:,:,:)
        INTEGER               :: i, j
        
        ! This subroutine constructs the conserved variable
        ! or solution vector from the primitive variables
        !
        ! U = [rho, rho*u, E]^T
        ! p = (gam-1)*(E-0.5*rho*u^2)
        ! E = p/(gam-1) + 0.5*rho*u^2
        
        DO i = 1, nx
            DO j = 1, ny
                eu(i, j, 1) = rho(i, j)
                eu(i, j, 2) = rho(i, j) * u_x(i, j)
                eu(i, j, 3) = rho(i, j) * u_y(i, j)
                eu(i, j, 4) = p(i, j) * gm1i + 0.5d0 * rho(i, j) * (u_x(i, j)**2 + u_y(i, j)**2)
            END DO
        END DO
        
    END SUBROUTINE solvec
    
    !**************************************************************
    
    SUBROUTINE fluxvec(u,f)
        IMPLICIT NONE
        REAL(DP), INTENT(in), DIMENSION(:,:) :: u
        REAL(DP), INTENT(out), DIMENSION(:,:) :: f
        INTEGER :: i
        
        ! This subroutine constructs the flux vector F from
        ! the conserved variable vector U based on the 1D
        ! Euler equations
        !
        ! u(1) = rho, u(2) = rho*u_x, u(3) = E
        ! f(1) = rho*u_x, f(2) = rho*u_x^2+p, f(3) = u_x*(E+p)
        ! p = (gam-1)*(E-0.5*rho*u^2)
        ! E = p/(gam-1) + 0.5*rho*u^2
        
        do i = 1, nx
            f(i,1) = u(i,2)
            f(i,2) = u(i,2)**2/u(i,1) + gm1*(u(i,3)-0.5d0*u(i,2)**2/u(i,1))
            f(i,3) = u(i,2)/u(i,1)*(u(i,3) + gm1*(u(i,3)-0.5d0*u(i,2)**2/u(i,1)))
        end do
        
    END SUBROUTINE fluxvec
    
    !**************************************************************
    
    SUBROUTINE decomp(rho, u_x, u_y, p, eu, c, Temp)
        IMPLICIT NONE
        REAL(DP), INTENT(IN), DIMENSION(:,:,:) :: eu
        REAL(DP), INTENT(OUT), DIMENSION(:,:) :: rho, u_x, u_y, p, c, Temp
        INTEGER :: i, j
        
        ! This suborutine decomposes the solution vector U back
        ! to the primitive variables according to the defintion
        ! for 1D Euler equations
        !
        ! rho, rho*u, E
        ! p = (gam-1)*(E-0.5*rho*u^2)
        ! E = p/(gam-1) + 0.5*rho*u^2
        
        DO i = 1, nx
            DO j = 1, ny
                rho(i, j) = MAX(eu(i,j,1), 1.0e-6_dp)
                u_x(i, j) = eu(i,j,2) / rho(i, j)
                u_y(i, j) = eu(i,j,3) / rho(i, j)
                p(i, j)   = MAX(gm1*(eu(i, j, 4) - 0.5d0 * rho(i, j) * (u_x(i, j)**2 + u_y(i,j)**2)), 1.0e-6_dp)
                c(i, j)   = dsqrt(gam * p(i, j) / rho(i, j))
                Temp(i,j) = p(i,j) / (rho(i,j) * R)
            END DO
        END DO
        
    END SUBROUTINE decomp
    
    !**************************************************************
    SUBROUTINE TENOFlux(fh_x, fh_y, rho, u_x, u_y, p, ct, solver)
        IMPLICIT NONE
        INTEGER                      :: k, i, j
        REAL(DP), INTENT(IN)         :: rho(:,:), u_x(:,:), u_y(:,:), p(:,:), ct
        CHARACTER(LEN=4), INTENT(IN) :: solver
        REAL(DP), INTENT(OUT)        :: fh_x(:,:,:), fh_y(:,:,:)
        REAL(DP), ALLOCATABLE        :: qL_x(:,:,:), qR_x(:,:,:), qL_y(:,:,:),&
         qR_y(:,:,:), q(:,:,:),dqL(:,:,:), dqR(:,:,:),dq(:,:,:)
        REAL(DP)                     :: beta(3), omega(3), cell_values(5), p_stencils(5)
        ALLOCATE(q(nx,ny,nv),qL_x(nx,ny,nv),qR_x(nx,ny,nv), qL_y(nx,ny,nv),qR_y(nx,ny,nv), dqL(nx,ny,nv),dqR(nx,ny,nv),dq(nx,ny,nv))
        q(:,:,:) = 0.0_dp
        qL_x(:,:,:) = 0.0_dp
        qR_x(:,:,:) = 0.0_dp
        qL_y(:,:,:) = 0.0_dp
        qR_y(:,:,:) = 0.0_dp
        dqL(:,:,:) = 0.0_dp
        dqR(:,:,:) = 0.0_dp

        !Build primitive variables array
        DO i = 1, nx
            DO j = 1, ny
                q(i, j, 1) = rho(i, j)
                q(i, j, 2) = u_x(i, j)
                q(i, j, 3) = u_y(i, j)
                q(i, j, 4) = p(i, j)
            END DO
        END DO
        
        DO k = 1, nv                ! Loop over rho, u, v and p components
            DO j = 1, ny
                DO i = 3, nx - 3
                    ! Assemble left-biased stencil cell primitive values
                    cell_values = (/q(i-2, j, k), q(i-1, j, k), q(i, j, k), q(i+1, j, k), q(i+2, j, k)/)
                    
                    ! Left-biased stencil primitive values before weights
                    CALL stencilsValues(cell_values, p_stencils)
                    
                    ! Smoothness indicators for left-biased stencil
                    CALL smoothnessIndicators(cell_values, beta)
                    
                    ! Compute non-linear weights
                    CALL TENOWeights(beta, ct, omega)
                    
                    ! Compute the left-biased primitive values at the interface
                    qL_x(i, j, k) = omega(1) * p_stencils(1) + omega(2) * p_stencils(2) + omega(3) * p_stencils(3)
                    
                    ! Assemble left-biased stencil cell primitive values
                    cell_values = (/q(i+3, j, k), q(i+2, j, k), q(i+1, j, k), q(i, j, k), q(i-1, j, k)/)
                    
                    ! Right-biased stencil primitive values before weights
                    CALL stencilsValues(cell_values, p_stencils)
                    
                    ! Smoothness indicators for right-biased stenci
                    CALL smoothnessIndicators(cell_values, beta)
                    
                    ! Compute non-linear weights
                    CALL TENOWeights(beta, ct, omega)
                    
                    ! Compute the right-biased primitive values at the interface
                    qR_x(i, j, k) = omega(1) * p_stencils(1) + omega(2) * p_stencils(2) + omega(3) * p_stencils(3)  
                END DO    
            END DO
        END DO

        DO k = 1, nv                ! Loop over rho, u, and p components
            DO i = 1, nx
                DO j = 3, ny - 3
                    ! Assemble left-biased stencil cell primitive values
                    cell_values = (/q(i, j-2, k), q(i, j-1, k), q(i, j, k), q(i, j+1, k), q(i, j+2, k)/)
                    
                    ! Left-biased stencil primitive values before weights
                    CALL stencilsValues(cell_values, p_stencils)
                    
                    ! Smoothness indicators for left-biased stencil
                    CALL smoothnessIndicators(cell_values, beta)
                    
                    ! Compute non-linear weights
                    CALL TENOWeights(beta, ct, omega)
                    
                    ! Compute the left-biased primitive values at the interface
                    qL_y(i,j,k) = omega(1) * p_stencils(1) + omega(2) * p_stencils(2) + omega(3) * p_stencils(3)
                    
                    ! Assemble right-biased stencil cell primitive values
                    cell_values = (/q(i, j+3, k), q(i, j+2, k), q(i, j+1, k), q(i, j, k), q(i, j-1, k)/)
                    
                    ! Right-biased stencil primitive values before weights
                    CALL stencilsValues(cell_values, p_stencils)
                    
                    ! Smoothness indicators for right-biased stenci
                    CALL smoothnessIndicators(cell_values, beta)
                    
                    ! Compute non-linear weights
                    CALL TENOWeights(beta, ct, omega)
                    
                    ! Compute the right-biased primitive values at the interface
                    qR_y(i, j, k) = omega(1) * p_stencils(1) + omega(2) * p_stencils(2) + omega(3) * p_stencils(3)  
                END DO    
            END DO
        END DO

        ! Perform interpolation over ghost cells
        CALL WTENOghostCellsInterpolation(qL_x, qL_y, qR_x, qR_y)

        ! Call approximate Riemann solver
        CALL riemann_solver(qL_x, qR_x, fh_x, qL_y, qR_y, fh_y, solver)
    END SUBROUTINE TENOFlux
    
    !**************************************************************
    SUBROUTINE TENOWeights(beta, ct, omega)
        
        IMPLICIT NONE
        INTEGER                             :: stencil_i
        REAL(DP)                            :: eps
        REAL(DP), INTENT(IN)                :: beta(3), ct
        REAL(DP), DIMENSION(3), INTENT(OUT) :: omega
        REAL(DP), DIMENSION(3)              :: gamma, chi, delta, d, ddelta
        
        eps = 1.0e-6_dp
        d(1) = 1.d0 / 10.d0
        d(2) = 6.d0 / 10.d0
        d(3) = 3.d0 / 10.d0
        
        ! Initial non-linear weights
        DO stencil_i = 1,3
            gamma(stencil_i) = (1.d0 + ABS(beta(3) - beta(1)) / (eps + beta(stencil_i)))**6
        END DO
        
        ! Normalize non-linear weights and determine ENO selection
        DO stencil_i = 1,3
            chi(stencil_i) = gamma(stencil_i) / SUM(gamma)
            
            IF (chi(stencil_i) < ct) THEN
                delta(stencil_i) = 0
            ELSE
                delta(stencil_i) = 1
            END IF
        END DO
        
        ! New non-linear weights
        ddelta = d * delta
        DO stencil_i = 1,3
            omega(stencil_i) = (d(stencil_i) * delta(stencil_i)) / SUM(ddelta)
        END DO
        
    END SUBROUTINE TENOWeights
    
    !**************************************************************
    SUBROUTINE WENOFlux(fh_x, fh_y, rho, u_x, u_y, p, solver)
        IMPLICIT NONE
        INTEGER                      :: k, i, j
        REAL(DP), INTENT(IN)         :: rho(:,:), u_x(:,:), u_y(:,:), p(:,:)
        CHARACTER(LEN=4), INTENT(IN) :: solver
        REAL(DP), INTENT(OUT)        :: fh_x(:,:,:), fh_y(:,:,:)
        REAL(DP), ALLOCATABLE        :: qL_x(:,:,:), qR_x(:,:,:), qL_y(:,:,:),&
         qR_y(:,:,:), q(:,:,:),dqL(:,:,:), dqR(:,:,:),dq(:,:,:)
        REAL(DP)                     :: beta(3), omega(3), cell_values(5), p_stencils(5)
        ALLOCATE(q(nx,ny,nv),qL_x(nx,ny,nv),qR_x(nx,ny,nv), qL_y(nx,ny,nv),qR_y(nx,ny,nv), dqL(nx,ny,nv),dqR(nx,ny,nv),dq(nx,ny,nv))
        q(:,:,:) = 0.0_dp
        qL_x(:,:,:) = 0.0_dp
        qR_x(:,:,:) = 0.0_dp
        qL_y(:,:,:) = 0.0_dp
        qR_y(:,:,:) = 0.0_dp
        dqL(:,:,:) = 0.0_dp
        dqR(:,:,:) = 0.0_dp

        !Build primitive variables array
        DO i = 1, nx
            DO j = 1, ny
                q(i, j, 1) = rho(i, j)
                q(i, j, 2) = u_x(i, j)
                q(i, j, 3) = u_y(i, j)
                q(i, j, 4) = p(i, j)
            END DO
        END DO
        
        DO k = 1, nv                ! Loop over rho, u, v and p components
            DO j = 1, ny
                DO i = 3, nx - 3
                    ! Assemble left-biased stencil cell primitive values
                    cell_values = (/q(i-2, j, k), q(i-1, j, k), q(i, j, k), q(i+1, j, k), q(i+2, j, k)/)
                    
                    ! Left-biased stencil primitive values before weights
                    CALL stencilsValues(cell_values, p_stencils)
                    
                    ! Smoothness indicators for left-biased stencil
                    CALL smoothnessIndicators(cell_values, beta)
                    
                    ! Compute non-linear weights
                    CALL WENOWeights(beta, omega)
                    
                    ! Compute the left-biased primitive values at the interface
                    qL_x(i, j, k) = omega(1) * p_stencils(1) + omega(2) * p_stencils(2) + omega(3) * p_stencils(3)
                    
                    ! Assemble left-biased stencil cell primitive values
                    cell_values = (/q(i+3, j, k), q(i+2, j, k), q(i+1, j, k), q(i, j, k), q(i-1, j, k)/)
                    
                    ! Right-biased stencil primitive values before weights
                    CALL stencilsValues(cell_values, p_stencils)
                    
                    ! Smoothness indicators for right-biased stenci
                    CALL smoothnessIndicators(cell_values, beta)
                    
                    ! Compute non-linear weights
                    CALL WENOWeights(beta, omega)
                    
                    ! Compute the right-biased primitive values at the interface
                    qR_x(i, j, k) = omega(1) * p_stencils(1) + omega(2) * p_stencils(2) + omega(3) * p_stencils(3)  
                END DO    
            END DO
        END DO

        DO k = 1, nv                ! Loop over rho, u, and p components
            DO i = 1, nx
                DO j = 3, ny - 3
                    ! Assemble left-biased stencil cell primitive values
                    cell_values = (/q(i, j-2, k), q(i, j-1, k), q(i, j, k), q(i, j+1, k), q(i, j+2, k)/)
                    
                    ! Left-biased stencil primitive values before weights
                    CALL stencilsValues(cell_values, p_stencils)
                    
                    ! Smoothness indicators for left-biased stencil
                    CALL smoothnessIndicators(cell_values, beta)
                    
                    ! Compute non-linear weights
                    CALL WENOWeights(beta, omega)
                    
                    ! Compute the left-biased primitive values at the interface
                    qL_y(i,j,k) = omega(1) * p_stencils(1) + omega(2) * p_stencils(2) + omega(3) * p_stencils(3)
                    
                    ! Assemble right-biased stencil cell primitive values
                    cell_values = (/q(i, j+3, k), q(i, j+2, k), q(i, j+1, k), q(i, j, k), q(i, j-1, k)/)
                    
                    ! Right-biased stencil primitive values before weights
                    CALL stencilsValues(cell_values, p_stencils)
                    
                    ! Smoothness indicators for right-biased stenci
                    CALL smoothnessIndicators(cell_values, beta)
                    
                    ! Compute non-linear weights
                    CALL WENOWeights(beta, omega)
                    
                    ! Compute the right-biased primitive values at the interface
                    qR_y(i, j, k) = omega(1) * p_stencils(1) + omega(2) * p_stencils(2) + omega(3) * p_stencils(3)  
                END DO    
            END DO
        END DO

        ! Perform interpolation over ghost cells
        CALL WTENOghostCellsInterpolation(qL_x, qL_y, qR_x, qR_y)

        ! Call approximate Riemann solver
        CALL riemann_solver(qL_x,qR_x,fh_x,qL_y,qR_y,fh_y, solver)
    END SUBROUTINE WENOFlux
    
    !**************************************************************
    SUBROUTINE WENOWeights(beta, omega)
        
        IMPLICIT NONE
        INTEGER                             :: stencil_i
        REAL(DP)                            :: alpha(3), d(3), eps
        REAL(DP), DIMENSION(3), INTENT(IN)  :: beta
        REAL(DP), DIMENSION(3), INTENT(OUT) :: omega
        
        eps = 1.0e-6_dp
        d(1) = 1.d0 / 10.d0
        d(2) = 6.d0/ 10.d0
        d(3) = 3.d0/ 10.d0
        ! Non-linear weights
        DO stencil_i = 1,3
            alpha(stencil_i) = d(stencil_i) / (eps + beta(stencil_i))**2
        END DO
        
        ! Normalize non-linear weights
        DO stencil_i = 1,3
            omega(stencil_i) = alpha(stencil_i) / SUM(alpha)
        END DO
    END SUBROUTINE WENOWeights
    
    !**************************************************************
    SUBROUTINE smoothnessIndicators(q, beta)
        IMPLICIT NONE
        REAL(DP), INTENT(IN)   :: q(:)
        REAL(DP), INTENT(OUT)  :: beta(:)
        
        ! Smoothness indicators
        beta(1) = (13.0_dp / 12.0_dp) * (q(1) - 2.0_dp * q(2) +&
        q(3))**2.0_dp + (1.0_dp / 4.0_dp) * (q(1) - 4.0_dp * q(2) + 3.0_dp * q(3))**2.0_dp
        beta(2) = (13.0_dp / 12.0_dp) * (q(2) - 2.0_dp * q(3) +&
        q(4))**2.0_dp + (1.0_dp / 4.0_dp) * (q(2) - q(4))**2.0_dp
        beta(3) = (13.0_dp / 12.0_dp) * (q(3) - 2.0_dp * q(4) +&
        q(5))**2.0_dp + (1.0 / 4.0) * (3.0_dp * q(3) - 4.0_dp * q(4) + q(5))**2.0_dp
        
    END SUBROUTINE smoothnessIndicators
    
    !**************************************************************
    SUBROUTINE stencilsValues(q, p)
        IMPLICIT NONE
        REAL(DP), INTENT(IN)   :: q(:)
        REAL(DP), INTENT(OUT)  :: p(:)
        
        ! Stencils values
        p(1) = (2.0_dp * q(1) - 7.0_dp * q(2) + 11.0_dp * q(3)) / 6.0_dp
        p(2) = (-q(2) + 5.0_dp * q(3) + 2.0_dp * q(4)) / 6.0_dp
        p(3) = (2.0_dp * q(3) + 5.0_dp * q(4) - q(5)) / 6.0_dp
        
    END SUBROUTINE stencilsValues

    !**************************************************************
    SUBROUTINE WTENOghostCellsInterpolation(qL_x, qL_y, qR_x, qR_y)
        IMPLICIT NONE
        REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)  :: qL_x, qL_y, qR_x, qR_y
        
        qL_x(2,:,:)     = qL_x(3,:,:)
        qL_x(1,:,:)     = qL_x(2,:,:)
        qR_x(2,:,:)     = qR_x(3,:,:)
        qR_x(1,:,:)     = qR_x(2,:,:)
        qL_x(nx-2,:,:)  = qL_x(nx-3,:,:)
        qL_x(nx-1,:,:)  = qL_x(nx-2,:,:)
        qL_x(nx,:,:)    = qL_x(nx-1,:,:)
        qR_x(nx-2,:,:)  = qR_x(nx-3,:,:)
        qR_x(nx-1,:,:)  = qR_x(nx-2,:,:)
        qR_x(nx,:,:)    = qR_x(nx-1,:,:)

        qL_y(:,2,:)     = qL_y(:,3,:)
        qL_y(:,1,:)     = qL_y(:,2,:)
        qR_y(:,2,:)     = qR_y(:,3,:)
        qR_y(:,1,:)     = qR_y(:,2,:)
        qL_y(:,ny-2,:)  = qL_y(:,ny-3,:)
        qL_y(:,ny-1,:)  = qL_y(:,ny-2,:)
        qL_y(:,ny,:)    = qL_y(:,ny-1,:)
        qR_y(:,ny-2,:)  = qR_y(:,ny-3,:)
        qR_y(:,ny-1,:)  = qR_y(:,ny-2,:)
        qR_y(:,ny,:)    = qR_y(:,ny-1,:)

    END SUBROUTINE WTENOghostCellsInterpolation
    
    !**************************************************************
    SUBROUTINE MTghostCellsInterpolation(qL_x, qL_y, qR_x, qR_y)
        IMPLICIT NONE
        REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)  :: qL_x, qL_y, qR_x, qR_y
        
        qL_x(1,:,:) = qL_x(2,:,:)
        qL_x(nx,:,:) = qL_x(nx-1,:,:)

        qR_x(1,:,:) = qR_x(2,:,:)
        qR_x(nx,:,:) = qR_x(nx-1,:,:)
      
        qL_y(:,ny-1,:) = qL_y(:,ny-2,:)
        qL_y(:,ny,:) = qL_y(:,ny-1,:)
        
        qR_y(:,ny-1,:) = qR_y(:,ny-2,:)
        qR_y(:,ny,:) = qR_y(:,ny-1,:)

    END SUBROUTINE MTghostCellsInterpolation
    
    !**************************************************************
    SUBROUTINE MUSCL_THINC(fh_x, fh_y, rho, u_x, u_y, p, beta, solver)
      
        IMPLICIT NONE
        INTEGER :: k, i, j
        REAL(DP), INTENT(IN) :: rho(:,:), u_x(:,:), u_y(:,:), p(:,:), beta
        CHARACTER(4), INTENT(IN):: solver
        REAL(DP), INTENT(INOUT) :: fh_x(:,:,:), fh_y(:,:,:)
        REAL(DP) :: qi,qip1,qim1
        REAL(DP) :: num,den,r
        REAL(DP) :: qLiphM,qRimhM
        REAL(DP) :: qLiphT,qRimhT
        REAL(DP) :: qLiphMT,qRimhMT
        REAL(DP) :: qmin,delq,theta,A,B,term
        REAL(DP) :: dpdx,drhodx,phip,phirho,phiprho,xsiph
        REAL(DP) :: eps, eta, arg1, arg2, xi
        REAL(DP), ALLOCATABLE :: q(:,:,:), qL_x(:,:,:), qR_x(:,:,:), qL_y(:,:,:), qR_y(:,:,:)
      
        ALLOCATE(q(nx,ny,nv),qL_x(nx,ny,nv),qR_x(nx,ny,nv), qL_y(nx,ny,nv),qR_y(nx,ny,nv))
      
        eta = 1.d0/3.d0
        eps = 1e-20
      
        !Build primitive variables array
        DO i = 1, nx
            DO j = 1, ny
                q(i, j, 1) = rho(i, j)
                q(i, j, 2) = u_x(i, j)
                q(i, j, 3) = u_y(i, j)
                q(i, j, 4) = p(i, j)
            END DO
        END DO
      
        ! X reconstruction
        DO k = 1, nv
            DO j = 1, ny
                ! Loop over each control volume
                DO i = 2, nx - 1
            
                    qi   = q(i,j,k)
                    qip1 = q(i+1,j,k)
                    qim1 = q(i-1,j,k)
            
                    num = qip1 - qi
                    den = qi - qim1
                    r = num/(den + 1e-08)
            
                    ! MUSCL cell edge predictors
                    qLiphM = qi + 0.25d0*((1.d0-eta)*psi(r)*(qi-qim1) &
                                    +         (1.d0+eta)*psi(1.d0/(r+1e-08))*(qip1-qi))
                    qRimhM = qi - 0.25d0*((1.d0+eta)*psi(r)*(qi-qim1) &
                                    -         (1.d0-eta)*psi(1.d0/(r+1e-08))*(qip1-qi))
                    qL_x(i,j,k)   = qLiphM
                    qR_x(i-1,j,k) = qRimhM
            
                    if ((qip1 - qi)*(qi - qim1) > 1e-30) then
            
                    ! THINC cell edge predictors
                    qmin = min(qim1,qip1)
                    delq = abs(qip1-qim1)
                    theta = sign(1.d0,qip1-qim1)
                    B = exp(theta*beta*((2.d0*(qi-qmin + eps))/(delq+eps) - 1.d0))
                    A = (B/cosh(beta)-1.d0)/tanh(beta)
                    term = (tanh(beta) + A)/(1.d0+A*tanh(beta))
                    qLiphT = qmin + delq*0.5d0*(1.d0+theta*term)
                    qRimhT = qmin + delq*0.5d0*(1.d0+theta*A)
            
                    arg1 = (qLiphM - qRimhM)/(qip1 - qi + eps)
                    arg2 = (qLiphM - qRimhM)/(qi - qim1 + eps)
                    xi = 1.d0 - min(arg1,arg2)
            
                    ! New hybrid scheme
                    dpdx   = (p(i+1,j)-p(i,j))/dx
                    drhodx = (rho(i+1,j)-rho(i,j))/dx
                    phip     = max(p(i,j),p(i+1,j))/min(p(i,j),p(i+1,j))
                    phirho   = max(rho(i,j),rho(i+1,j))/min(rho(i,j),rho(i+1,j))
                    if (dpdx*drhodx >= 0) then
                        phiprho = phip/(phirho+1e-08)
                        xsiph = exp(-25.d0*(max(1.d0,phiprho)-1.d0))
                    else
                        xsiph = 1.d0
                    end if
            
                ! MUSCL-THINC hybrid cell edge predictor
                    qLiphMT = (1.d0-xi*xsiph)*qLiphM + xi*xsiph*qLiphT
                    qRimhMT = (1.d0-xi*xsiph)*qRimhM + xi*xsiph*qRimhT
            
                    qL_x(i,j,k)   = qLiphMT
                    qR_x(i-1,j,k) = qRimhMT
                    end if
                END DO
            END DO
        END DO

        ! Y reconstruction
        DO k = 1, nv
            DO i = 1, nx
                ! Loop over each control volume
                DO j = 2, ny-  1
            
                    qi   = q(i,j,k)
                    qip1 = q(i,j+1,k)
                    qim1 = q(i,j-1,k)
            
                    num = qip1 - qi
                    den = qi - qim1
                    r = num/(den + 1e-08)
            
                    ! MUSCL cell edge predictors
                    qLiphM = qi + 0.25d0*((1.d0-eta)*psi(r)*(qi-qim1) &
                                    +         (1.d0+eta)*psi(1.d0/(r+1e-08))*(qip1-qi))
                    qRimhM = qi - 0.25d0*((1.d0+eta)*psi(r)*(qi-qim1) &
                                    -         (1.d0-eta)*psi(1.d0/(r+1e-08))*(qip1-qi))
                    qL_y(i,j,k)   = qLiphM
                    qR_y(i,j-1,k) = qRimhM
            
                    if ((qip1 - qi)*(qi - qim1) > 1e-30) then
            
                    ! THINC cell edge predictors
                    qmin = min(qim1,qip1)
                    delq = abs(qip1-qim1)
                    theta = sign(1.d0,qip1-qim1)
                    B = exp(theta*beta*((2.d0*(qi-qmin + eps))/(delq+eps) - 1.d0))
                    A = (B/cosh(beta)-1.d0)/tanh(beta)
                    term = (tanh(beta) + A)/(1.d0+A*tanh(beta))
                    qLiphT = qmin + delq*0.5d0*(1.d0+theta*term)
                    qRimhT = qmin + delq*0.5d0*(1.d0+theta*A)
            
                    arg1 = (qLiphM - qRimhM)/(qip1 - qi + eps)
                    arg2 = (qLiphM - qRimhM)/(qi - qim1 + eps)
                    xi = 1.d0 - min(arg1,arg2)
            
                    ! New hybrid scheme
                    dpdx   = (p(i,j+1)-p(i,j)) / dy
                    drhodx = (rho(i,j+1)-rho(i,j)) / dy
                    phip     = max(p(i,j),p(i,j+1))/min(p(i,j),p(i,j+1))
                    phirho   = max(rho(i,j),rho(i,j+1))/min(rho(i,j),rho(i,j+1))
                    if (dpdx*drhodx >= 0) then
                        phiprho = phip/(phirho+1e-08)
                        xsiph = exp(-25.d0*(max(1.d0,phiprho)-1.d0))
                    else
                        xsiph = 1.d0
                    end if
            
                ! MUSCL-THINC hybrid cell edge predictor
                    qLiphMT = (1.d0-xi*xsiph)*qLiphM + xi*xsiph*qLiphT
                    qRimhMT = (1.d0-xi*xsiph)*qRimhM + xi*xsiph*qRimhT
            
                    qL_y(i,j,k)   = qLiphMT
                    qR_y(i,j-1,k) = qRimhMT
                    end if
                END DO
            END DO
        END DO

        CALL MTghostCellsInterpolation(qL_x, qL_y, qR_x, qR_y)
      
        ! Call approximate Riemann solver
        CALL riemann_solver(qL_x, qR_x, fh_x, qL_y, qR_y, fh_y, solver)
      
    END SUBROUTINE MUSCL_THINC

    REAL FUNCTION psi(arg)
        IMPLICIT NONE
        REAL(DP) :: arg
    
        ! This is needed for the MUSCL-THINC scheme
        psi = max(0.d0,min(1.d0,arg))
    
        RETURN
    END FUNCTION psi
    !**************************************************************
    SUBROUTINE riemann_HLLC(qL_x, qR_x, fh_x, qL_y, qR_y, fh_y)
        IMPLICIT NONE
        INTEGER :: i,j
        REAL(DP), INTENT(IN), DIMENSION(:,:,:)  :: qL_x,qR_x,qL_y,qR_y
        REAL(DP), INTENT(INOUT), DIMENSION(:,:,:)  :: fh_x,fh_y
        REAL(DP) :: pR, pL, cL, cR, eps, hL, hR
        REAL(DP) :: uL1,uL2,uL3,uL4
        REAL(DP) :: uR1,uR2,uR3,uR4
        REAL(DP) :: fL1,fL2,fL3,fL4
        REAL(DP) :: fR1,fR2,fR3,fR4
        REAL(DP) :: SL, SR, sStar, delta_p,pStar
        REAL(DP) :: uL, uR, vL, vR
        REAL(DP) :: uL_cof, uR_cof
        REAL(DP) :: uL_s1, uL_s2, uL_s3, uL_s4, fL_s1, fL_s2, fL_s3, fL_s4
        REAL(DP) :: uR_s1, uR_s2, uR_s3, uR_s4, fR_s1, fR_s2, fR_s3, fR_s4
        
        eps = 1.0e-6_dp
        
        !========
        !X Fluxes
        !========
        DO i = 1, nx
            DO j=1, ny
                uL1 = max(qL_x(i,j,1), eps)
                uL2 = uL1 * qL_x(i,j,2) 
                uL3 = uL1 * qL_x(i,j,3) 
                uL4 = qL_x(i,j,4) * gm1i + 0.5d0 * uL1 * (qL_x(i,j,2)**2 + qL_x(i,j,3)**2)
                
                uR1 = max(qR_x(i,j,1), eps)
                uR2 = uR1 * qR_x(i,j,2)
                uR3 = uR1 * qR_x(i,j,3)
                uR4 = qR_x(i,j,4) * gm1i + 0.5d0 * uR1 * (qR_x(i,j,2)**2 + qR_x(i,j,3)**2) 
                
                uL = qL_x(i,j,2)
                uR = qR_x(i,j,2)
                vL = qL_x(i,j,3)
                vR = qR_x(i,j,3)
                
                pL = MAX(gm1 * (uL4 - 0.5d0 * (uL2**2 + uL3**2) / (uL1 + eps)), eps)
                pR = MAX(gm1 * (uR4 - 0.5d0 * (uR2**2 + uR3**2) / (uR1 + eps)), eps)
                
                ! Compute left and right enthalpy
                hL = (uL4 + pL) / uL1
                hR = (uR4 + pR) / uR1
                
                cL = SQRT(gam * pL / uL1)
                cR = SQRT(gam * pR / uR1)
                
                ! Pressure-based entropy fix
                delta_p = 0.1_dp * ABS(pR - pL) / MAX(pL, pR, eps)
                SL = MIN(uL - cL, uR - cR) - delta_p
                SR = MAX(uL + cL, uR + cR) + delta_p
                
                fL1 = uL2                               ! rho_L * uL
                fL2 = uL2**2 / uL1 + pL                 ! rho_L * uL^2 + p_L
                fL3 = uL2*uL3/uL1                       ! rho_L * uL * vL
                fL4 = uL2 / uL1 * (uL4 + pL)            ! uL * (E_L + p_L)
                
                fR1 = uR2                               ! rho_R * uR
                fR2 = uR2**2 / uR1 + pR                 ! rho_R * uR^2 + p_R
                fR3 = uR2*uR3/uR1                       ! rho_R * uR * vR
                fR4 = uR2 / uR1 * (uR4 + pR)            ! uR * (E_R + p_R)
                
                IF (SL > 0.0_dp) THEN
                    fh_x(i,j,1) = fL1
                    fh_x(i,j,2) = fL2
                    fh_x(i,j,3) = fL3
                    fh_x(i,j,4) = fL4
                ELSE IF (SR < 0.0_dp) THEN
                    fh_x(i,j,1) = fR1
                    fh_x(i,j,2) = fR2
                    fh_x(i,j,3) = fR3
                    fh_x(i,j,4) = fR4
                ELSE
                    ! Contact wave velocity
                    sStar = (pR - pL + uL2 * (SL - uL) - uR2 * (SR - uR)) / &
                    (uL1 * (SL - uL) - uR1 * (SR - uR) + eps)
                    
                    ! Star pressure stabilization
                    pStar = 0.5d0 * (pL + pR) - 0.5d0 * (SR - SL) * (uR1 - uL1)
                    
                    ! Left star region
                    uL_cof = uL1 * (SL - uL) / (SL - sStar + eps)
                    uL_s1 = uL_cof
                    uL_s2 = uL_cof * sStar
                    uL_s3 = uL_cof * vL
                    uL_s4 = uL_cof * (uL4 / uL1 + (sStar - uL) * (sStar + pStar / (uL1 * (SL - uL) + eps)))
                    
                    fL_s1 = fL1 + SL * (uL_s1 - uL1)
                    fL_s2 = fL2 + SL * (uL_s2 - uL2)
                    fL_s3 = fL3 + SL * (uL_s3 - uL3)
                    fL_s4 = fL4 + SL * (uL_s4 - uL4)
                    
                    ! Right star region
                    uR_cof = uR1 * (SR - uR) / (SR - sStar + eps)
                    uR_s1 = uR_cof
                    uR_s2 = uR_cof * sStar
                    uR_s3 = uR_cof * vR
                    uR_s4 = uR_cof * (uR4 / uR1 + (sStar - uR) * (sStar + pStar / (uR1 * (SR - uR) + eps)))
                    
                    fR_s1 = fR1 + SR * (uR_s1 - uR1)
                    fR_s2 = fR2 + SR * (uR_s2 - uR2)
                    fR_s3 = fR3 + SR * (uR_s3 - uR3)
                    fR_s4 = fR4 + SR * (uR_s4 - uR4)
                    
                    ! If 0 <= sStar => use left star
                    IF (sStar >= 0.0_dp) THEN
                        fh_x(i,j,1) = fL_s1
                        fh_x(i,j,2) = fL_s2
                        fh_x(i,j,3) = fL_s3
                        fh_x(i,j,4) = fL_s4
                    ELSE
                        fh_x(i,j,1) = fR_s1
                        fh_x(i,j,2) = fR_s2
                        fh_x(i,j,3) = fR_s3
                        fh_x(i,j,4) = fR_s4
                    END IF
                END IF
                
            END DO
        END DO
        
        !========
        !Y Fluxes
        !========
        DO i = 1, nx
            DO j=1, ny
                uL1 = max(qL_y(i,j,1), eps)
                uL2 = uL1 * qL_y(i,j,2) 
                uL3 = uL1 * qL_y(i,j,3) 
                uL4 = qL_y(i,j,4) * gm1i + 0.5d0 * uL1 * (qL_y(i,j,2)**2 + qL_y(i,j,3)**2)
                
                uR1 = max(qR_y(i,j,1), eps)
                uR2 = uR1 * qR_y(i,j,2)
                uR3 = uR1 * qR_y(i,j,3)
                uR4 = qR_y(i,j,4) * gm1i + 0.5d0 * uR1 * (qR_y(i,j,2)**2 + qR_y(i,j,3)**2) 
                
                uL = qL_y(i,j,2)
                uR = qR_y(i,j,2)
                vL = qL_y(i,j,3)
                vR = qR_y(i,j,3)
                
                pL = MAX(gm1 * (uL4 - 0.5d0 * (uL2**2 + uL3**2) / (uL1 + eps)), eps)
                pR = MAX(gm1 * (uR4 - 0.5d0 * (uR2**2 + uR3**2) / (uR1 + eps)), eps)
                
                ! Compute left and right enthalpy
                hL = (uL4 + pL) / uL1
                hR = (uR4 + pR) / uR1
                
                cL = SQRT(gam * pL / uL1)
                cR = SQRT(gam * pR / uR1)
                
                ! Pressure-based entropy fix
                delta_p = 0.1_dp * ABS(pR - pL) / MAX(pL, pR, eps)
                SL = MIN(vL - cL, vR - cR) - delta_p
                SR = MAX(vL + cL, vR + cR) + delta_p
                
                fL1 = uL3                           ! rho_L * vL
                fL2 = uL3 * uL2 / uL1               ! rho_L * uL * vL
                fL3 = uL3**2 / uL1 + pL             ! rho_L * vL^2 + p_L
                fL4 = uL3 / uL1 * (uL4 + pL)        ! vL * (E_L + p_L)
                
                fR1 = uR3                              ! rho_R * uR
                fR2 =  uR2*uR3/uR1                ! rho_R * uR * vR
                fR3 =  uR3**2 / uR1 + pR                      ! rho_R * uR^2 + p_R
                fR4 = uR3 / uR1 * (uR4 + pR)            ! uR * (E_R + p_R)
                
                IF (SL > 0.0_dp) THEN
                    fh_y(i,j,1) = fL1
                    fh_y(i,j,2) = fL2
                    fh_y(i,j,3) = fL3
                    fh_y(i,j,4) = fL4
                ELSE IF (SR < 0.0_dp) THEN
                    fh_y(i,j,1) = fR1
                    fh_y(i,j,2) = fR2
                    fh_y(i,j,3) = fR3
                    fh_y(i,j,4) = fR4
                ELSE
                    ! Contact wave velocity
                    sStar = (pR - pL + uL3 * (SL - vL) - uR3 * (SR - vR)) / &
                    (uL1 * (SL - vL) - uR1 * (SR - vR) + eps)
                    
                    ! Star pressure stabilization
                    pStar = 0.5d0 * (pL + pR) - 0.5d0 * (SR - SL) * (uR1 - uL1)
                    
                    ! Left star region
                    uL_cof = uL1 * (SL - vL) / (SL - sStar + eps)
                    uL_s1 = uL_cof
                    uL_s2 = uL_cof * uL
                    uL_s3 = uL_cof * sStar
                    uL_s4 = uL_cof * (uL4 / uL1 + (sStar - vL) * (sStar + pStar / (uL1 * (SL - vL) + eps)))
                    
                    fL_s1 = fL1 + SL * (uL_s1 - uL1)
                    fL_s2 = fL2 + SL * (uL_s2 - uL2)
                    fL_s3 = fL3 + SL * (uL_s3 - uL3)
                    fL_s4 = fL4 + SL * (uL_s4 - uL4)
                    
                    ! Right star region
                    uR_cof = uR1 * (SR - vR) / (SR - sStar + eps)
                    uR_s1 = uR_cof
                    uR_s2 = uR_cof * uR
                    uR_s3 = uR_cof * sStar
                    uR_s4 = uR_cof * (uR4 / uR1 + (sStar - vR) * (sStar + pStar / (uR1 * (SR - vR) + eps)))
                    
                    fR_s1 = fR1 + SR * (uR_s1 - uR1)
                    fR_s2 = fR2 + SR * (uR_s2 - uR2)
                    fR_s3 = fR3 + SR * (uR_s3 - uR3)
                    fR_s4 = fR4 + SR * (uR_s4 - uR4)
                    
                    ! If 0 <= sStar => use left star
                    IF (sStar >= 0.0_dp) THEN
                        fh_y(i,j,1) = fL_s1
                        fh_y(i,j,2) = fL_s2
                        fh_y(i,j,3) = fL_s3
                        fh_y(i,j,4) = fL_s4
                    ELSE
                        fh_y(i,j,1) = fR_s1
                        fh_y(i,j,2) = fR_s2
                        fh_y(i,j,3) = fR_s3
                        fh_y(i,j,4) = fR_s4
                    END IF
                END IF
                
            END DO
        END DO
        
        
    END SUBROUTINE riemann_HLLC
        
    !**************************************************************
    SUBROUTINE ApplyBoundaryConditions(rho, ux, uy, p)
        IMPLICIT NONE
        REAL(DP), INTENT(INOUT) :: rho(:,:), ux(:,:), uy(:,:), p(:,:)
        INTEGER :: i, j, n, m
        REAL(DP), PARAMETER :: EPSILON = 1.0E-12_DP
      
        ! Left boundary (x = 1)

        DO j = 1, ny
          DO n = 1, ng
            rho(ng+1-n, j) = rho(ng+n+1, j)
            ux(ng+1-n, j) =  -ux(ng+n+1, j)
            uy(ng+1-n, j) =  -uy(ng+n+1, j)
            p(ng+1-n, j) =     p(ng+n+1, j)
          END DO
        END DO
      
        ! Right boundary (x = nx)
        DO j = 1, ny
          DO n = 1, ng
            rho(nx-ng+n, j) = rho(nx-ng-n, j)
            ux(nx-ng+n, j) =  -ux(nx-ng-n, j)
            uy(nx-ng+n, j) =  -uy(nx-ng-n, j)
            p(nx-ng+n, j) =     p(nx-ng-n, j)
          END DO
        END DO
      
        
      
        DO i = 1, nx
          DO n = 1, ng
            ! bottom boundary
            rho(i, ng+1-n) = rho(i, ng+n+1)
            ux(i, ng+1-n) =  -ux(i, ng+n+1)
            uy(i, ng+1-n) =  -uy(i, ng+n+1)
            p(i, ng+1-n) =     p(i, ng+n+1)

            ! top boundary - symmetry
            rho(i, ny+n-ng) =  rho(i, ny-ng-n)
            ux(i, ny+n-ng) =    ux(i, ny-ng-n)
            uy(i, ny+n-ng) =    uy(i, ny-ng-n)
            p(i, ny+n-ng) =      p(i, ny-ng-n)
          END DO
        END DO


        ! hemanth: fill the corner points
      ! Corner Treatment
          DO n = 1, ng
            DO m = 1, ng
              ! Bottom-left corner (x = 1, y = 1)
              rho(ng+1-n, ng+1-m) = rho(ng+1, ng+1)
              ux(ng+1-n, ng+1-m) = -ux(ng+1, ng+1)
              uy(ng+1-n, ng+1-m) = -uy(ng+1, ng+1)
              p(ng+1-n, ng+1-m) = p(ng+1, ng+1)

              ! Bottom-right corner (x = nx, y = 1)
              rho(nx-ng+n, ng+1-m) = rho(nx-ng, ng+1)
              ux(nx-ng+n, ng+1-m) = -ux(nx-ng, ng+1)
              uy(nx-ng+n, ng+1-m) = -uy(nx-ng, ng+1)
              p(nx-ng+n, ng+1-m) = p(nx-ng, ng+1)

              ! Top-left corner (x = 1, y = ny)
              rho(ng+1-n, ny+m-ng) = rho(ng+1, ny-ng)
              ux(ng+1-n, ny+m-ng) = -ux(ng+1, ny-ng)
              uy(ng+1-n, ny+m-ng) = -uy(ng+1, ny-ng)
              p(ng+1-n, ny+m-ng) = p(ng+1, ny-ng)

              ! Top-right corner (x = nx, y = ny)
              rho(nx+n-ng, ny+m-ng) = rho(nx-ng, ny-ng)
              ux(nx+n-ng, ny+m-ng) = -ux(nx-ng, ny-ng)
              uy(nx+n-ng, ny+m-ng) = -uy(nx-ng, ny-ng)
              p(nx+n-ng, ny+m-ng) = p(nx-ng, ny-ng)
            END DO
          END DO
    END SUBROUTINE ApplyBoundaryConditions
    
    !**************************************************************
    SUBROUTINE viscousFluxes(u, v, dFvx_dx, dFvy_dy, T)
        IMPLICIT NONE
        INTEGER :: i, j, k
        REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)  :: dFvx_dx, dFvy_dy
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: u, v, T
        REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Fvx_L, Fvx_R, Fvy_L, Fvy_R
        REAL(DP), ALLOCATABLE :: Tau_xx_L(:,:), Tau_yy_L(:,:), Tau_xy_L(:,:), Tau_yx_L(:,:)
        REAL(DP), ALLOCATABLE :: Tau_xx_R(:,:), Tau_yy_R(:,:), Tau_xy_R(:,:), Tau_yx_R(:,:)
        REAL(DP), ALLOCATABLE :: qx_L(:,:), qy_L(:,:), qx_R(:,:), qy_R(:,:)
        REAL(DP) :: u_L, u_R, v_L, v_R
    
        ALLOCATE(Fvx_L(nx,ny,nv), Fvx_R(nx,ny,nv), Fvy_L(nx,ny,nv), Fvy_R(nx,ny,nv))
        ALLOCATE(Tau_xx_L(nx,ny), Tau_yy_L(nx,ny), Tau_xy_L(nx,ny),Tau_yx_L(nx,ny))
        ALLOCATE(Tau_xx_R(nx,ny), Tau_yy_R(nx,ny), Tau_xy_R(nx,ny),Tau_yx_R(nx,ny))
        ALLOCATE(qx_L(nx,ny), qy_L(nx,ny), qx_R(nx,ny), qy_R(nx,ny))
    
        ! Initialize fluxes
        Fvx_L(:,:,:) = 0.0_dp
        Fvx_R(:,:,:) = 0.0_dp
        Fvy_L(:,:,:) = 0.0_dp
        Fvy_R(:,:,:) = 0.0_dp
    
        ! Viscous stresses
        CALL viscousStress(Tau_xx_L, Tau_yy_L, Tau_xy_L, Tau_xx_R, Tau_yy_R, Tau_xy_R&
        , Tau_yx_L, Tau_yx_R, u, v)
    
        ! Heat fluxes
        CALL HeatFluxes(qx_L, qx_R, qy_L, qy_R, T)
    
        ! X viscous fluxes
        DO j = 3, ny-2
            DO i = 3, nx-2
                ! Left and right velocities at cell interface
                u_L = (-u(i-2,j) + 9.0_dp * u(i-1,j) + 9.0_dp * u(i,j) - u(i+1,j)) / 16.0_dp
                u_R = (-u(i-1,j) + 9.0_dp * u(i,j) + 9.0_dp * u(i+1,j) - u(i+2,j)) / 16.0_dp
    
                v_L = (-v(i-2,j) + 9.0_dp * v(i-1,j) + 9.0_dp * v(i+0,j) - v(i+1,j)) / 16.0_dp
                v_R = (-v(i-1,j) + 9.0_dp * v(i+0,j) + 9.0_dp * v(i+1,j) - v(i+2,j)) / 16.0_dp
    
                ! Left and right viscous fluxes at (i+1/2, j)
                Fvx_L(i,j,1) = 0.0_dp
                Fvx_L(i,j,2) = Tau_xx_L(i,j)
                Fvx_L(i,j,3) = Tau_xy_L(i,j)
                Fvx_L(i,j,4) = u_L * Tau_xx_L(i,j) + v_L * Tau_xy_L(i,j) - qx_L(i,j)
    
                Fvx_R(i,j,1) = 0.0_dp
                Fvx_R(i,j,2) = Tau_xx_R(i,j)
                Fvx_R(i,j,3) = Tau_xy_R(i,j)
                Fvx_R(i,j,4) = u_R * Tau_xx_R(i,j) + v_R * Tau_xy_R(i,j) - qx_R(i,j)
            END DO
        END DO
    
        ! Y viscous fluxes
        DO j = 3, ny-2
            DO i = 3, nx-2
    
                u_L = (-u(i,j-2) + 9.0_dp * u(i,j-1) + 9.0_dp * u(i,j+0) - u(i,j+1)) / 16.0_dp
                u_R = (-u(i,j-1) + 9.0_dp * u(i,j+0) + 9.0_dp * u(i,j+1) - u(i,j+2)) / 16.0_dp
    
                v_L = (-v(i,j-2) + 9.0_dp * v(i,j-1) + 9.0_dp * v(i,j) - v(i,j+1)) / 16.0_dp
                v_R = (-v(i,j-1) + 9.0_dp * v(i,j) + 9.0_dp * v(i,j+1) - v(i,j+2)) / 16.0_dp
                
                ! Left and right viscous fluxes at (i, j+1/2)
                Fvy_L(i,j,1) = 0.0_dp
                Fvy_L(i,j,2) = Tau_yx_L(i,j)
                Fvy_L(i,j,3) = Tau_yy_L(i,j)
                Fvy_L(i,j,4) = u_L * Tau_yx_L(i,j) + v_L * Tau_yy_L(i,j) - qy_L(i,j)
    
                Fvy_R(i,j,1) = 0.0_dp!
                Fvy_R(i,j,2) = Tau_yx_R(i,j)
                Fvy_R(i,j,3) = Tau_yy_R(i,j)
                Fvy_R(i,j,4) = u_R * Tau_yx_R(i,j) + v_R * Tau_yy_R(i,j) - qy_R(i,j)
            END DO
        END DO
    
        ! dFvx/dx using finite volume method
        DO j = 1, ny
            DO i = 1, nx
                DO k = 1, 4
                    dFvx_dx(i,j,k) = (Fvx_R(i,j,k) - Fvx_L(i,j,k)) / dx
                END DO
            END DO
        END DO
    
        ! dFvy/dy using finite volume method
        DO j = 1, ny
            DO i = 1, nx
                DO k = 1, 4
                    dFvy_dy(i,j,k) = (Fvy_R(i,j,k) - Fvy_L(i,j,k)) / dy
                END DO
            END DO
        END DO
        
        ! Deallocate
        DEALLOCATE(Fvx_L, Fvx_R, Fvy_L, Fvy_R)
        DEALLOCATE(Tau_xx_L, Tau_yy_L, Tau_xy_L, Tau_yx_L)
        DEALLOCATE(Tau_xx_R, Tau_yy_R, Tau_xy_R, Tau_yx_R)
        DEALLOCATE(qx_L, qy_L, qx_R, qy_R)
    
    END SUBROUTINE viscousFluxes
    
    !**************************************************************
    SUBROUTINE viscousStress(Tau_xx_L, Tau_yy_L, Tau_xy_L,Tau_xx_R, Tau_yy_R, Tau_xy_R, Tau_yx_L, Tau_yx_R, u, v)
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: u(nx,ny), v(nx,ny)
        REAL(DP), INTENT(INOUT) :: Tau_xx_L(nx,ny), Tau_xx_R(nx,ny)
        REAL(DP), INTENT(INOUT) :: Tau_yy_L(nx,ny), Tau_yy_R(nx,ny)
        REAL(DP), INTENT(INOUT) :: Tau_xy_L(nx,ny), Tau_xy_R(nx,ny)
        REAL(DP), INTENT(INOUT) :: Tau_yx_L(nx,ny), Tau_yx_R(nx,ny)
        REAL(DP) :: dudx_L(nx,ny), dudx_R(nx,ny), dvdy_L(nx,ny), dvdy_R(nx,ny)
        REAL(DP) :: dudy_L(nx,ny), dudy_R(nx,ny), dvdx_L(nx,ny), dvdx_R(nx,ny)
        REAL(DP) :: dudx(nx,ny), dvdy(nx,ny), dudy(nx,ny), dvdx(nx,ny)
        INTEGER :: i, j

        Tau_xx_L(:,:) = 0.0_dp
        Tau_xx_R(:,:) = 0.0_dp
        Tau_yy_L(:,:) = 0.0_dp
        Tau_yy_R(:,:) = 0.0_dp
        Tau_xy_L(:,:) = 0.0_dp
        Tau_xy_R(:,:) = 0.0_dp
        dudx(:,:) = 0.0_dp
        dvdy(:,:) = 0.0_dp
        dudx_L(:,:) = 0.0_dp
        dvdy_L(:,:) = 0.0_dp
        dudx_R(:,:) = 0.0_dp
        dvdy_R(:,:) = 0.0_dp

        ! Central differences at cell centers
        DO j = 3, ny-2
            DO i = 3, nx-2
                dudx(i,j) = (2.0_dp / 3.0_dp) * ((u(i+1,j) - u(i-1,j)) / dx) - &
                            (1.0_dp / 12.0_dp) * ((u(i+2,j) - u(i-2,j)) / dx)

                dvdy(i,j) = (2.0_dp / 3.0_dp) * ((v(i,j+1) - v(i,j-1)) / dy) - &
                            (1.0_dp / 12.0_dp) * ((v(i,j+2) - v(i,j-2)) / dy)

                dudy(i,j) = (2.0_dp / 3.0_dp) * ((u(i,j+1) - u(i,j-1)) / dy) - &
                            (1.0_dp / 12.0_dp) * ((u(i,j+2) - u(i,j-2)) / dy)

                dvdx(i,j) = (2.0_dp / 3.0_dp) * ((v(i+1,j) - v(i-1,j)) / dx) - &
                            (1.0_dp / 12.0_dp) * ((v(i+2,j) - v(i-2,j)) / dx)

            END DO
        END DO

        ! Upwind scheme near boundaries
        DO j = 1, ny
            ! Left boundary
            dudx(1,j) = (-25.0_dp * u(1,j) + 48.0_dp * u(2,j) - 36.0_dp * u(3,j) + &
                        16.0_dp * u(4,j) - 3.0_dp * u(5,j)) / (12.0_dp * dx)
            dudx(2,j) = (-3.0_dp * u(1,j) - 10.0_dp * u(2,j) +18.0_dp * u(3,j) - &
                        6.0_dp * u(4,j) +  u(5,j)) / (12.0_dp * dx)
        
            dvdx(1,j) = (-25.0_dp * v(1,j) + 48.0_dp * v(2,j) - 36.0_dp * v(3,j) + &
                        16.0_dp * v(4,j) - 3.0_dp * v(5,j)) / (12.0_dp * dx)
            dvdx(2,j) = (-3.0_dp * v(1,j) - 10.0_dp * v(2,j) +18.0_dp * v(3,j) - &
                        6.0_dp * v(4,j) +  v(5,j)) / (12.0_dp * dx)
        
            ! Right boundary
            dudx(nx-1,j) = (-3.0_dp * u(nx,j) - 10.0_dp * u(nx-1,j) +18.0_dp * u(nx-2,j) - &
                            6.0_dp * u(nx-3,j) +  u(nx-4,j)) / (12.0_dp * dx)
            dudx(nx,j) = (25.0_dp * u(nx,j) - 48.0_dp * u(nx-1,j) + 36.0_dp * u(nx-2,j) - &
                        16.0_dp * u(nx-3,j) + 3.0_dp * u(nx-4,j)) / (12.0_dp * dx)
        
            dvdx(nx-1,j) = (-3.0_dp * v(nx,j) - 10.0_dp * v(nx-1,j) +18.0_dp * v(nx-2,j) - &
                        6.0_dp * v(nx-3,j) +  v(nx-4,j)) / (12.0_dp * dx)
            dvdx(nx,j) = (25.0_dp * v(nx,j) - 48.0_dp * v(nx-1,j) + 36.0_dp * v(nx-2,j) - &
                        16.0_dp * v(nx-3,j) + 3.0_dp * v(nx-4,j)) / (12.0_dp * dx)
        END DO

        DO i = 1, nx
            ! Bottom boundary
            dvdy(i,1) = (-25.0_dp * v(i,1) + 48.0_dp * v(i,2) - 36.0_dp * v(i,3) + &
                        16.0_dp * v(i,4) - 3.0_dp * v(i,5)) / (12.0_dp * dy)
            dvdy(i,2) = (-3.0_dp * v(i,1) - 10.0_dp * v(i,2) +18.0_dp * v(i,3) - &
            6.0_dp * v(i,4) +  v(i,5)) / (12.0_dp * dy)

            dudy(i,1) = (-25.0_dp * u(i,1) + 48.0_dp * u(i,2) - 36.0_dp * u(i,3) + &
                        16.0_dp * u(i,4) - 3.0_dp * u(i,5)) / (12.0_dp * dy)
            dudy(i,2) = (-3.0_dp * u(i,1) - 10.0_dp * u(i,2) +18.0_dp * u(i,3) - &
            6.0_dp * u(i,4) +  u(i,5)) / (12.0_dp * dy)

            ! Top boundary
            dvdy(i,ny-1) = (-3.0_dp * v(i,ny) - 10.0_dp * v(i,ny-1) +18.0_dp * v(i,ny-2) - &
                            6.0_dp * v(i,ny-3) +  v(i,ny-4)) / (12.0_dp * dy) 
            dvdy(i,ny) = (25.0_dp * v(i,ny) - 48.0_dp * v(i,ny-1) + 36.0_dp * v(i,ny-2) - &
                        16.0_dp * v(i,ny-3) + 3.0_dp * v(i,ny-4)) / (12.0_dp * dy)

            dudy(i,ny-1) = (-3.0_dp * u(i,ny) - 10.0_dp * u(i,ny-1) +18.0_dp * u(i,ny-2) - &
                        6.0_dp * u(i,ny-3) +  u(i,ny-4)) / (12.0_dp * dy) 
            dudy(i,ny) = (25.0_dp * u(i,ny) - 48.0_dp * u(i,ny-1) + 36.0_dp * u(i,ny-2) - &
                        16.0_dp * u(i,ny-3) + 3.0_dp * u(i,ny-4)) / (12.0_dp * dy)                  
        END DO

        ! X staggered differences
        DO j = 3, ny-2
            DO i = 3, nx-2
                dudx_R(i,j) = (9.0_dp / 8.0_dp) * ((u(i+1,j) - u(i,j)) / dx) - &
                            (1.0_dp / 8.0_dp) * ((u(i+2,j) - u(i-1,j)) / (3.0_dp * dx))

                dudx_L(i,j) = (9.0_dp / 8.0_dp) * ((u(i,j) - u(i-1,j)) / dx) - &
                            (1.0_dp / 8.0_dp) * ((u(i+1,j) - u(i-2,j)) / (3.0_dp * dx))

                dvdx_R(i,j) = (9.0_dp / 8.0_dp) * ((v(i+1,j) - v(i,j)) / dx) - &
                            (1.0_dp / 8.0_dp) * ((v(i+2,j) - v(i-1,j)) / (3.0_dp * dx))

                dvdx_L(i,j) = (9.0_dp / 8.0_dp) * ((v(i,j) - v(i-1,j)) / dx) - &
                            (1.0_dp / 8.0_dp) * ((v(i+1,j) - v(i-2,j)) / (3.0_dp * dx))                          

                dvdy_R(i,j) = (-dvdy(i-1,j) + 7.0_dp * dvdy(i,j) + 7.0_dp * dvdy(i+1,j) - dvdy(i+2,j)) / 16.0_dp
                dvdy_L(i,j) = (-dvdy(i-2,j) + 7.0_dp * dvdy(i-1,j) + 7.0_dp * dvdy(i,j) - dvdy(i+1,j)) / 16.0_dp

                dudy_R(i,j) = (-dudy(i-1,j) + 7.0_dp * dudy(i,j) + 7.0_dp * dudy(i+1,j) - dudy(i+2,j)) / 16.0_dp
                dudy_L(i,j) = (-dudy(i-2,j) + 7.0_dp * dudy(i-1,j) + 7.0_dp * dudy(i,j) - dudy(i+1,j)) / 16.0_dp
            END DO
        END DO

        ! X viscous stresses
        DO j = 3, ny-2
            DO i = 3, nx-2
                ! Tau_xx at (i+1/2, j) and (i-1/2, j)
                Tau_xx_R(i,j) = (4.0_dp / 3.0_dp) * mu * dudx_R(i,j) - (2.0_dp / 3.0_dp) * mu * dvdy_R(i,j)
                Tau_xx_L(i,j) = (4.0_dp / 3.0_dp) * mu * dudx_L(i,j) - (2.0_dp / 3.0_dp) * mu * dvdy_L(i,j)

                Tau_xy_R(i,j) = mu * (dudy_R(i,j)+dvdx_R(i,j))
                Tau_xy_L(i,j) = mu * (dudy_L(i,j)+dvdx_L(i,j))
            END DO
        END DO

        ! Y staggered differences
        DO j = 3, ny-2
            DO i = 3, nx-2
                dvdy_R(i,j) = (9.0_dp / 8.0_dp) * ((v(i,j+1) - v(i,j)) / dy) - &
                            (1.0_dp / 8.0_dp) * ((v(i,j+2) - v(i,j-1)) / (3.0_dp * dy))

                dvdy_L(i,j) = (9.0_dp / 8.0_dp) * ((v(i,j) - v(i,j-1)) / dy) - &
                            (1.0_dp / 8.0_dp) * ((v(i,j+1) - v(i,j-2)) / (3.0_dp * dy))

                dudy_R(i,j) = (9.0_dp / 8.0_dp) * ((u(i,j+1) - u(i,j)) / dy) - &
                            (1.0_dp / 8.0_dp) * ((u(i,j+2) - u(i,j-1)) / (3.0_dp * dy))

                dudy_L(i,j) = (9.0_dp / 8.0_dp) * ((u(i,j) - u(i,j-1)) / dy) - &
                            (1.0_dp / 8.0_dp) * ((u(i,j+1) - u(i,j-2)) / (3.0_dp * dy))

                dudx_R(i,j) = (-dudx(i,j-1) + 7.0_dp * dudx(i,j) + 7.0_dp * dudx(i,j+1) - dudx(i,j+2)) / 16.0_dp
                dudx_L(i,j) = (-dudx(i,j-2) + 7.0_dp * dudx(i,j-1) + 7.0_dp * dudx(i,j) - dudx(i,j+1)) / 16.0_dp

                dvdx_R(i,j) = (-dvdx(i,j-1) + 7.0_dp * dvdx(i,j) + 7.0_dp * dvdx(i,j+1) - dvdx(i,j+2)) / 16.0_dp
                dvdx_L(i,j) = (-dvdx(i,j-2) + 7.0_dp * dvdx(i,j-1) + 7.0_dp * dvdx(i,j) - dvdx(i,j+1)) / 16.0_dp
            END DO
        END DO

        ! Y viscous stresses
        DO j = 3, ny-2
            DO i = 3, nx-2
                ! Compute Tau_yy at (i, j+1/2) and (i, j-1/2)
                Tau_yy_R(i,j) = (4.0_dp / 3.0_dp) * mu * dvdy_R(i,j) - (2.0_dp / 3.0_dp) * mu * dudx_R(i,j)
                Tau_yy_L(i,j) = (4.0_dp / 3.0_dp) * mu * dvdy_L(i,j) - (2.0_dp / 3.0_dp) * mu * dudx_L(i,j)

                ! Compute Tau_xy at (i+1/2, j) and (i-1/2, j)
                Tau_yx_R(i,j) = mu * (dudy_R(i,j)+dvdx_R(i,j))
                Tau_yx_L(i,j) = mu * (dudy_L(i,j)+dvdx_L(i,j))
            END DO
        END DO
            
    END SUBROUTINE viscousStress
    
    !**************************************************************
    SUBROUTINE heatFluxes(qx_L, qx_R, qy_L, qy_R, T)
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: T(:,:)
        REAL(DP), INTENT(OUT) :: qx_L(:,:), qx_R(:,:), qy_L(:,:), qy_R(:,:)
        REAL(DP) :: dTdx_L(nx,ny), dTdx_R(nx,ny), dTdy_L(nx,ny), dTdy_R(nx,ny)
        INTEGER :: i, j
    
        ! Initialize arrays
        qx_L(:,:) = 0.0_dp
        qx_R(:,:) = 0.0_dp
        qy_L(:,:) = 0.0_dp
        qy_R(:,:) = 0.0_dp
    
        ! Compute dT/dx at left and right staggered interfaces (i-1/2, j) and (i+1/2, j)
        DO j = 3, ny-2
            DO i = 3, nx-2
                dTdx_L(i,j) = (T(i-2,j) - 9.0_dp * T(i-1,j) + 9.0_dp * T(i,j) - T(i+1,j)) / (8.0_dp * dx)
                dTdx_R(i,j) = (T(i-1,j) - 9.0_dp * T(i,j) + 9.0_dp * T(i+1,j) - T(i+2,j)) / (8.0_dp * dx)
            END DO
        END DO
    
        ! Compute heat flux qx at left and right staggered interfaces
        DO j = 3, ny-2
            DO i = 3, nx-2
                qx_L(i,j) = -ka * dTdx_L(i,j)
                qx_R(i,j) = -ka * dTdx_R(i,j)
            END DO
        END DO
    
        ! Compute dT/dy at bottom and top staggered interfaces (i, j-1/2) and (i, j+1/2)
        DO j = 3, ny-2
            DO i = 3, nx-2
                dTdy_L(i,j) = (T(i,j-2) - 9.0_dp * T(i,j-1) + 9.0_dp * T(i,j)   - T(i,j+1)) / (8.0_dp * dy)
                dTdy_R(i,j) = (T(i,j-1) - 9.0_dp * T(i,j)   + 9.0_dp * T(i,j+1) - T(i,j+2)) / (8.0_dp * dy)
            END DO
        END DO
    
        ! Compute heat flux qy at bottom and top staggered interfaces
        DO j = 3, ny-2
            DO i = 3, nx-2
                qy_L(i,j) = -ka * dTdy_L(i,j)
                qy_R(i,j) = -ka * dTdy_R(i,j)
            END DO
        END DO
    
        ! Apply boundary conditions
        DO j = 1, ny
            qx_L(1,j) = qx_L(2,j)
            qx_L(nx,j) = qx_L(nx-1,j)
            qx_R(1,j) = qx_R(2,j)
            qx_R(nx,j) = qx_R(nx-1,j)
    
            qy_L(1,j) = qy_L(2,j)
            qy_L(nx,j) = qy_L(nx-1,j)
            qy_R(1,j) = qy_R(2,j)
            qy_R(nx,j) = qy_R(nx-1,j)
        END DO
    
        DO i = 1, nx
            qx_L(i,1) = qx_L(i,2)
            qx_L(i,ny) = qx_L(i,ny-1)
            qx_R(i,1) = qx_R(i,2)
            qx_R(i,ny) = qx_R(i,ny-1)
    
            qy_L(i,1) = qy_L(i,2)
            qy_L(i,ny) = qy_L(i,ny-1)
            qy_R(i,1) = qy_R(i,2)
            qy_R(i,ny) = qy_R(i,ny-1)
        END DO
    
    END SUBROUTINE heatFluxes

    !**************************************************************
    SUBROUTINE RK3TVD(eu, cfl, ct, sharpness, rho, u_x, u_y, p, c, fh_x, fh_y, Temp, scheme, solver)
        IMPLICIT NONE
        INTEGER i, j, state_count
        REAL(DP)                  , INTENT(IN)    :: cfl, ct, sharpness
        CHARACTER(LEN=4)          , INTENT(IN)    :: scheme, solver
        REAL(DP), DIMENSION(:,:)  , INTENT(INOUT) :: rho, u_x, u_y, p, c, Temp
        REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: fh_x, fh_y, eu

        REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: dFvx_dx, dFvy_dy
        REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: eu0, eu1, eu2
        
        ALLOCATE(dFvx_dx(nx, ny, nv), dFvy_dy(nx, ny, nv))
        ALLOCATE(eu0(nx,ny,nv), eu1(nx,ny,nv), eu2(nx,ny,nv))

        t = 0.d0
        ! Count time steps for printing
        state_count = 0
        DO while (t<tfinal)
            state_count = state_count + 1
            ! Compute interative wave speed
            wavespeed = maxval(SQRT(u_x**2 + u_y**2) + sqrt(gam * p / rho))
            ! Print state of analysis
            CALL printState(scheme, state_count)
            
            dt = cfl * MIN(dx, dy) / abs(wavespeed)
            
            if(t+dt > tfinal) then
                dt = tfinal-t
            end if
            dtodx = dt/dx
            dxodt = dx/dt
            dtody = dt/dy
            dyodt = dy/dt
            
            ! Stage 1
            eu0 = eu               ! Save old solution
            CALL activateScheme(scheme, rho, u_x, u_y, p, fh_x, fh_y, ct, sharpness, solver)
            CALL viscousFluxes(u_x, u_y, dFvx_dx, dFvy_dy, Temp)
            DO i = 2, nx - 1
                DO j = 2, ny - 1
                    eu1(i,j,:) = eu0(i,j,:) - dtodx*(fh_x(i,j,:) - fh_x(i-1,j,:)) - dtody*(fh_y(i,j,:) - fh_y(i,j-1,:))&
                    + dt * (1.0_dp / Re) * (dFvx_dx(i,j,:) + dFvy_dy(i,j,:))
                END DO
            END DO

            
            ! Stage 2
            CALL decomp(rho, u_x, u_y, p, eu1, c, Temp)
            CALL ApplyBoundaryConditions(rho, u_x, u_y, p)
            CALL solvec(rho, u_x, u_y, p, eu1)

            CALL activateScheme(scheme, rho, u_x, u_y, p, fh_x, fh_y, ct, sharpness, solver)
            CALL viscousFluxes(u_x, u_y, dFvx_dx, dFvy_dy, Temp)
            DO i = 2, nx - 1
                DO j = 2, ny - 1
                    eu2(i,j,:) = (3.0_dp/4.0_dp) * eu0(i,j,:) + (1.0_dp/4.0_dp)*&
                    (eu1(i,j,:) - dtodx*(fh_x(i,j,:) - fh_x(i-1,j,:)) - dtody*(fh_y(i,j,:) - fh_y(i,j-1,:))&
                    + dt * (1.0_dp / Re) * (dFvx_dx(i,j,:) + dFvy_dy(i,j,:)))
                END DO
            END DO
            
            ! Stage 3
            CALL decomp(rho, u_x, u_y, p, eu2, c, Temp)
            CALL ApplyBoundaryConditions(rho, u_x, u_y, p)
            CALL solvec(rho, u_x, u_y, p, eu2)

            CALL activateScheme(scheme, rho, u_x, u_y, p, fh_x, fh_y, ct, sharpness, solver)
            CALL viscousFluxes(u_x, u_y, dFvx_dx, dFvy_dy, Temp)
            DO i = 2, nx - 1
                DO j = 2, ny - 1
                    eu(i,j,:) = (1.0_dp / 3.0_dp)*eu0(i,j,:) + &
                    (2.0_dp/3.0_dp)*(eu2(i,j,:) - dtodx*(fh_x(i,j,:) - fh_x(i-1,j,:)) - dtody*(fh_y(i,j,:) - fh_y(i,j-1,:))&
                    + dt * (1.0_dp / Re) * (dFvx_dx(i,j,:) + dFvy_dy(i,j,:)))
                END DO
            END DO
            
            CALL decomp(rho, u_x, u_y, p, eu, c, Temp)
            CALL ApplyBoundaryConditions(rho, u_x, u_y, p)
            CALL solvec(rho, u_x, u_y, p, eu)

            t = t + dt
            if (t>=tfinal) then
                CALL ApplyBoundaryConditions(rho, u_x, u_y, p)
                EXIT
            end if
        END DO
        DEALLOCATE(eu0, eu1, eu2)
    END SUBROUTINE RK3TVD
        
    !**************************************************************
    SUBROUTINE activateScheme(scheme, rho, u_x, u_y, p, fh_x, fh_y, ct, sharpness, solver)
            IMPLICIT NONE
            REAL(DP)        , INTENT(IN) :: ct, sharpness, rho(:,:), u_x(:,:), u_y(:,:), p(:,:)
            CHARACTER(LEN=4), INTENT(IN) :: scheme, solver
            REAL(DP)        , INTENT(OUT):: fh_x(:,:,:), fh_y(:,:,:)
            
            ! Determine which scheme to activate.
            IF (scheme == "TENO") THEN
                CALL TENOFlux(fh_x, fh_y, rho, u_x, u_y, p, ct, solver)

             ELSEIF (scheme == "WENO") THEN
                CALL WENOFlux(fh_x, fh_y, rho, u_x, u_y, p, solver)

             ELSEIF (scheme == "MUSC") THEN
                 CALL MUSCL_THINC(fh_x, fh_y, rho, u_x, u_y, p, sharpness, solver)
            END IF    
    END SUBROUTINE activateScheme
                
    !**************************************************************
    SUBROUTINE riemann_solver(qL_x, qR_x, fh_x, qL_y, qR_y, fh_y, solver)
        IMPLICIT NONE
        REAL(DP)        , INTENT(IN) :: qL_x(:,:,:), qR_x(:,:,:), qL_y(:,:,:), qR_y(:,:,:)
        CHARACTER(LEN=4), INTENT(IN) :: solver
        REAL(DP)        , INTENT(OUT):: fh_x(:,:,:), fh_y(:,:,:)
        
        ! Determine which solver to activate.
        IF (solver == "HLLC") THEN
            CALL riemann_HLLC(qL_x, qR_x, fh_x, qL_y, qR_y, fh_y)
        END IF  
    END SUBROUTINE riemann_solver
            
    !**************************************************************
    SUBROUTINE output(x, y, rho, u_x, u_y, p, scheme, solver, Re, Temp)
        IMPLICIT NONE
        INTEGER i, j
        CHARACTER(LEN=50)            :: filename
        CHARACTER(LEN=4), INTENT(IN) :: scheme, solver
        REAL(DP)        , INTENT(IN) :: x(:), y(:), rho(:,:), u_x(:,:), u_y(:,:), p(:,:), Temp(:,:), Re
        
        WRITE(filename, '(A,A,A,A,I0,A,F7.2,A,A)') scheme,"_",solver,"_grid_", nx, "_Re_", Re,"_0.6_",".dat"
        OPEN(UNIT=18, FILE=TRIM(filename), STATUS="REPLACE", ACTION="WRITE", FORM="FORMATTED")
        DO i = 1, nx
            DO j = 1, ny
                WRITE(18,*) x(i), y(j), rho(i,j), Temp(i,j), u_x(i,j), u_y(i,j), p(i,j)
            END DO
        END DO
        CLOSE(18)
    END SUBROUTINE output

    !**************************************************************
    SUBROUTINE printState(scheme, state_count)
        IMPLICIT NONE
        CHARACTER(LEN=4)          , INTENT(IN) :: scheme
        INTEGER, INTENT(IN)                    :: state_count
        IF (MOD(state_count, 100) == 0) THEN
            print *, "====================="
            print *, "Reynolds Number: ", Re," Scheme: ", scheme, " Grid: ", nx 
            print *, "Time: ", t
            ! Output for debugging
            ! print *, "wavespeed: ", wavespeed, " u_x: ", maxval(u_x), " u_y: ", maxval(u_y),&
            ! " p: ", maxval(p), " rho: ", maxval(rho)
        END IF
    END SUBROUTINE printState
    END MODULE subroutines
            
    !**************************************************************
    PROGRAM myeuler1d 
                USE types_vars
                USE subroutines
                
                IMPLICIT NONE
                INTEGER                                 :: scheme_index, solver_index, grid_index, Re_index
                REAL(DP)                                :: cfl, ct, sharpness
                REAL(DP), DIMENSION(:), ALLOCATABLE     :: x, y
                REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: rho, u_x, u_y, p, c, Temp
                REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: eu, fh_x, fh_y
                INTEGER                                 :: grid_size(1) = (/600/)
                REAL(DP)                                :: Res(2) = (/1000, 200/)
                CHARACTER(LEN=4)                        :: scheme(3) = (/"WENO", "MUSC", "TENO"/), solver(1) = (/"HLLC"/)

                
                ! Run TENO, WENO and MUSCL schemes.
                DO Re_index = 1, SIZE(Res)
                    Re = Res(Re_index)
                    DO scheme_index = 1, SIZE(scheme)
                        ! Number of ghost cells
                        IF (scheme(scheme_index) == "WENO" .OR. scheme(scheme_index) == "TENO") THEN
                            ng = 5
                        ELSE IF (scheme(scheme_index) == "MUSC") THEN
                            ng = 2
                        END IF
        
                        DO grid_index = 1, SIZE(grid_size)
                            nx = grid_size(grid_index) + 2 * ng
                            ny = grid_size(grid_index) / 2 + 2 * ng
                            
                            ! Adiabatic index expressions
                            gm1 = gam-1.d0
                            gm1i = 1.d0/gm1
                            ! Viscosity parameters
                            Mach=1.0d0
                            Pr=0.72d0
                            cp=1005.0d0
                            mu=1.0d0
                            ka=cp * mu / Pr
                            R=287.05d0
                            ! Spatial and temporal parameters
                            nc_x = nx + 1
                            nc_y = ny + 1
                            nv = 4
                            xl = -0.0d0
                            xr = 1.0d0
                            yl = -0.0d0
                            yr = 0.5d0
                            dx = (xr - xl) / (nx - 1)
                            dy = (yr - yl) / (ny - 1)
                            cfl = 0.1
                            tfinal = 0.6d0
                            ct = 10.d0**(-7.d0)
                            sharpness = 2.d0
        
        
                                DO solver_index = 1, SIZE(solver)
                                    ALLOCATE(x(nx), y(ny))
                                    ! Velocities array is defined as (x value | y values)
                                    ALLOCATE(rho(nx,ny),u_x(nx,ny), u_y(nx,ny), p(nx,ny), c(nx,ny), Temp(nx, ny))
                                    ! Converved values array is defined as (x co. | y co. | x values | y values)
                                    ALLOCATE(eu(nx,ny,nv), fh_x(nc_x,nc_y,nv), fh_y(nc_x,nc_y,nv))
                                    
                                    CALL grid(x, y)
                                    CALL init_shock_tube(x, rho, u_x, u_y, p, Temp)
                                    CALL solvec(rho, u_x, u_y, p, eu)
                                    
                                    
                                    CALL RK3TVD(eu, cfl, ct, sharpness, rho, u_x, u_y, p, c, fh_x, fh_y,&
                                    Temp, scheme(scheme_index), solver(solver_index))
                                    CALL output(x, y, rho, u_x, u_y, p, scheme(scheme_index), solver(solver_index), Re, Temp)
                                    
                                    DEALLOCATE(x, y)
                                    DEALLOCATE(rho, u_x, u_y, p, c, Temp)
                                    DEALLOCATE(eu, fh_x, fh_y)
                                END DO
                            END DO
                        END DO
                END DO
                    
    END PROGRAM myeuler1d