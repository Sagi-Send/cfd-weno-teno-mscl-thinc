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
    
    INTEGER :: nx, ny, nc_x, nc_y, ng, nv, nt, ntsteps, nf, kf
    REAL(DP) :: xl, xr, yl, yr, t, tfinal, dt, dx, dy, dtodx, dxodt, dtody, dyodt, gam=1.4d0, gm1, gm1i, wavespeed
    
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
        
        !$omp parallel do default(shared) private(i) schedule(static)
        do i = 1, nx
            x(i) = xl + REAL(i-ng-1, DP)*dx + 0.5d0*dx
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(i) schedule(static)
        do i = 1, ny
            y(i) = yl + REAL(i-ng-1, DP)*dy + 0.5d0*dy
        end do
        !$omp end parallel do
        
    END SUBROUTINE grid
    
    !**************************************************************
    
    SUBROUTINE init_riemann(x, y, rho, u_x, u_y, p)
        IMPLICIT NONE
        REAL(DP), INTENT(IN), DIMENSION(:) :: x, y
        REAL(DP), INTENT(OUT)              :: rho(:,:), u_x(:,:), u_y(:,:), p(:,:)
        INTEGER                            :: i, j
        ! 2D Riemann problem
        !$omp parallel do collapse(2) default(shared) private(i,j) schedule(static)
        DO i = 1, SIZE(x)
            DO j = 1, SIZE(y)
                IF (x(i) <= 0.5d0 .AND. y(j) >= 0.5d0) THEN
                    rho(i, j) = 1.d0
                    u_x(i, j) = 0.726d0
                    u_y(i, j) = 0.d0
                    p(i, j)   = 1.d0
                ELSE IF (x(i) <= 0.5d0 .AND. y(j) < 0.5d0) THEN
                    rho(i, j) = 0.8d0
                    u_x(i, j) = 0.d0
                    u_y(i, j) = 0.d0
                    p(i, j)   = 1.d0
                ELSE IF (x(i) > 0.5d0 .AND. y(j) >= 0.5d0) THEN
                    rho(i, j) = 0.5313d0
                    u_x(i, j) = 0.d0
                    u_y(i, j)   = 0.d0
                    p(i, j)  = 0.4d0
                ELSE IF (x(i) > 0.5d0 .AND. y(j) < 0.5d0) THEN
                    rho(i, j) = 1.d0
                    u_x(i, j) = 0.d0
                    u_y(i, j) = 0.7276d0
                    p(i, j)   = 1.d0
                END IF
            END DO
        END DO
        !$omp end parallel do
    END SUBROUTINE init_riemann
    
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
        
        !$omp parallel do collapse(2) default(shared) private(i,j) schedule(static)
        DO i = 1, nx
            DO j = 1, ny
                eu(i, j, 1) = rho(i, j)
                eu(i, j, 2) = rho(i, j) * u_x(i, j)
                eu(i, j, 3) = rho(i, j) * u_y(i, j)
                eu(i, j, 4) = p(i, j) * gm1i + 0.5d0 * rho(i, j) * (u_x(i, j)**2 + u_y(i, j)**2)
            END DO
        END DO
        !$omp end parallel do
        
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
        ! u(1) = rho, u(2) = rho*ux, u(3) = E
        ! f(1) = rho*ux, f(2) = rho*ux^2+p, f(3) = ux*(E+p)
        ! p = (gam-1)*(E-0.5*rho*u^2)
        ! E = p/(gam-1) + 0.5*rho*u^2
        
        do i = 1, nx
            f(i,1) = u(i,2)
            f(i,2) = u(i,2)**2/u(i,1) + gm1*(u(i,3)-0.5d0*u(i,2)**2/u(i,1))
            f(i,3) = u(i,2)/u(i,1)*(u(i,3) + gm1*(u(i,3)-0.5d0*u(i,2)**2/u(i,1)))
        end do
        
    END SUBROUTINE fluxvec
    
    !**************************************************************
    
    SUBROUTINE decomp(rho, u_x, u_y, p, eu, c)
        IMPLICIT NONE
        REAL(DP), INTENT(IN), DIMENSION(:,:,:) :: eu
        REAL(DP), INTENT(OUT), DIMENSION(:,:) :: rho, u_x, u_y, p, c
        INTEGER :: i, j
        
        ! This suborutine decomposes the solution vector U back
        ! to the primitive variables according to the defintion
        ! for 1D Euler equations
        !
        ! rho, rho*u, E
        ! p = (gam-1)*(E-0.5*rho*u^2)
        ! E = p/(gam-1) + 0.5*rho*u^2
        
        !$omp parallel do collapse(2) default(shared) private(i,j) schedule(static)
        DO i = 1, nx
            DO j = 1, ny
                rho(i, j) = eu(i,j,1)
                u_x(i, j) = eu(i,j,2) / rho(i, j)
                u_y(i, j) = eu(i,j,3) / rho(i, j)
                p(i, j)   = gm1*(eu(i, j, 4) - 0.5d0 * rho(i, j) * (u_x(i, j)**2 + u_y(i,j)**2))
                c(i, j)   = dsqrt(gam * p(i, j) / rho(i, j))
            END DO
        END DO
        !$omp end parallel do
        
    END SUBROUTINE decomp
    
    !**************************************************************
    SUBROUTINE TENOFlux(fh_x, fh_y, rho, u_x, u_y, p, ct, solver)
        IMPLICIT NONE
        INTEGER                      :: k, i, j
        REAL(DP), INTENT(IN)         :: rho(:,:), u_x(:,:), u_y(:,:), p(:,:), ct
        CHARACTER(LEN=4), INTENT(IN) :: solver
        REAL(DP), INTENT(OUT)        :: fh_x(:,:,:), fh_y(:,:,:)
        REAL(DP), ALLOCATABLE        :: qL_x(:,:,:), qR_x(:,:,:), qL_y(:,:,:),&
         qR_y(:,:,:), q(:,:,:)
        REAL(DP)                     :: beta(3), omega(3), cell_values(5), p_stencils(5)
        ALLOCATE(q(nx,ny,nv),qL_x(nx,ny,nv),qR_x(nx,ny,nv), qL_y(nx,ny,nv),qR_y(nx,ny,nv))
        q(:,:,:) = 0.0_dp
        qL_x(:,:,:) = 0.0_dp
        qR_x(:,:,:) = 0.0_dp
        qL_y(:,:,:) = 0.0_dp
        qR_y(:,:,:) = 0.0_dp

        !Build primitive variables array
        !$omp parallel do collapse(2) default(shared) private(i,j) schedule(static)
        DO i = 1, nx
            DO j = 1, ny
                q(i, j, 1) = rho(i, j)
                q(i, j, 2) = u_x(i, j)
                q(i, j, 3) = u_y(i, j)
                q(i, j, 4) = p(i, j)
            END DO
        END DO
        !$omp end parallel do
        
        !$omp parallel do collapse(2) default(shared) private(k,j,i,cell_values,p_stencils,beta,omega) schedule(static)
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
        !$omp end parallel do

        !$omp parallel do collapse(2) default(shared) private(k,i,j,cell_values,p_stencils,beta,omega) schedule(static)
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
        !$omp end parallel do

        ! Perform interpolation over ghost cells
        CALL ghostCellsInterpolation(qL_x, qL_y, qR_x, qR_y)

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
         qR_y(:,:,:), q(:,:,:)
        REAL(DP)                     :: beta(3), omega(3), cell_values(5), p_stencils(5)
        ALLOCATE(q(nx,ny,nv),qL_x(nx,ny,nv),qR_x(nx,ny,nv), qL_y(nx,ny,nv),qR_y(nx,ny,nv))
        q(:,:,:) = 0.0_dp
        qL_x(:,:,:) = 0.0_dp
        qR_x(:,:,:) = 0.0_dp
        qL_y(:,:,:) = 0.0_dp
        qR_y(:,:,:) = 0.0_dp

        !Build primitive variables array
        !$omp parallel do collapse(2) default(shared) private(i,j) schedule(static)
        DO i = 1, nx
            DO j = 1, ny
                q(i, j, 1) = rho(i, j)
                q(i, j, 2) = u_x(i, j)
                q(i, j, 3) = u_y(i, j)
                q(i, j, 4) = p(i, j)
            END DO
        END DO
        !$omp end parallel do
        
        !$omp parallel do collapse(2) default(shared) private(k,j,i,cell_values,p_stencils,beta,omega) schedule(static)
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
        !$omp end parallel do

        !$omp parallel do collapse(2) default(shared) private(k,i,j,cell_values,p_stencils,beta,omega) schedule(static)
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
        !$omp end parallel do

        ! Perform interpolation over ghost cells
        CALL ghostCellsInterpolation(qL_x, qL_y, qR_x, qR_y)

        ! Call approximate Riemann solver
        CALL riemann_solver(qL_x, qR_x, fh_x, qL_y, qR_y, fh_y, solver)
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
    SUBROUTINE ghostCellsInterpolation(qL_x, qL_y, qR_x, qR_y)
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

    END SUBROUTINE ghostCellsInterpolation
    
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
        q(:,:,:) = 0.0_dp
        qL_x(:,:,:) = 0.0_dp
        qR_x(:,:,:) = 0.0_dp
        qL_y(:,:,:) = 0.0_dp
        qR_y(:,:,:) = 0.0_dp
      
        eta = 1.d0/3.d0
        eps = 1e-20
      
        !Build primitive variables array
        !$omp parallel do collapse(2) default(shared) private(i,j) schedule(static)
        DO i = 1, nx
            DO j = 1, ny
                q(i, j, 1) = rho(i, j)
                q(i, j, 2) = u_x(i, j)
                q(i, j, 3) = u_y(i, j)
                q(i, j, 4) = p(i, j)
            END DO
        END DO
        !$omp end parallel do
      
        ! X reconstruction
        !$omp parallel do collapse(2) default(shared) &
        !$omp& private(k,j,i,qi,qip1,qim1,num,den,r,qLiphM,qRimhM,qLiphT,qRimhT,qLiphMT,qRimhMT) &
        !$omp& private(qmin,delq,theta,A,B,term,dpdx,drhodx,phip,phirho,phiprho,xsiph,arg1,arg2,xi) &
        !$omp& schedule(static)
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
        !$omp end parallel do

        ! Y reconstruction
        !$omp parallel do collapse(2) default(shared) &
        !$omp& private(k,i,j,qi,qip1,qim1,num,den,r,qLiphM,qRimhM,qLiphT,qRimhT,qLiphMT,qRimhMT) &
        !$omp& private(qmin,delq,theta,A,B,term,dpdx,drhodx,phip,phirho,phiprho,xsiph,arg1,arg2,xi) &
        !$omp& schedule(static)
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
        !$omp end parallel do
      
        qL_x(1,:,1) = qL_x(2,:,1)
        qL_x(1,:,2) = qL_x(2,:,2)
        qL_x(1,:,3) = qL_x(2,:,3)
        qL_x(1,:,4) = qL_x(2,:,4)
        qL_x(nx,:,1) = qL_x(nx-1,:,1)
        qL_x(nx,:,2) = qL_x(nx-1,:,2)
        qL_x(nx,:,3) = qL_x(nx-1,:,3)
        qL_x(nx,:,4) = qL_x(nx-1,:,4)
      
        qR_y(:,ny-1,1) = qR_y(:,ny-2,1)
        qR_y(:,ny-1,2) = qR_y(:,ny-2,2)
        qR_y(:,ny-1,3) = qR_y(:,ny-2,3)
        qR_y(:,ny-1,4) = qR_y(:,ny-2,4)
        qR_y(:,ny,1) = qR_y(:,ny-1,1)
        qR_y(:,ny,2) = qR_y(:,ny-1,2)
        qR_y(:,ny,3) = qR_y(:,ny-1,3)
        qR_y(:,ny,4) = qR_y(:,ny-1,4)
        ! Call approximate Riemann solver
        CALL riemann_solver(qL_x, qR_x, fh_x, qL_y, qR_y, fh_y, solver)
      
    END SUBROUTINE MUSCL_THINC

    REAL(DP) FUNCTION psi(arg)
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
        !$omp parallel do collapse(2) default(shared) &
        !$omp& private(i,j,pR,pL,cL,cR,hL,hR,uL1,uL2,uL3,uL4,uR1,uR2,uR3,uR4) &
        !$omp& private(fL1,fL2,fL3,fL4,fR1,fR2,fR3,fR4,SL,SR,sStar,delta_p,pStar,uL,uR,vL,vR) &
        !$omp& private(uL_cof,uR_cof,uL_s1,uL_s2,uL_s3,uL_s4,fL_s1,fL_s2,fL_s3,fL_s4) &
        !$omp& private(uR_s1,uR_s2,uR_s3,uR_s4,fR_s1,fR_s2,fR_s3,fR_s4) &
        !$omp& schedule(static)
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
        !$omp end parallel do
        
        !========
        !Y Fluxes
        !========
        !$omp parallel do collapse(2) default(shared) &
        !$omp& private(i,j,pR,pL,cL,cR,hL,hR,uL1,uL2,uL3,uL4,uR1,uR2,uR3,uR4) &
        !$omp& private(fL1,fL2,fL3,fL4,fR1,fR2,fR3,fR4,SL,SR,sStar,delta_p,pStar,uL,uR,vL,vR) &
        !$omp& private(uL_cof,uR_cof,uL_s1,uL_s2,uL_s3,uL_s4,fL_s1,fL_s2,fL_s3,fL_s4) &
        !$omp& private(uR_s1,uR_s2,uR_s3,uR_s4,fR_s1,fR_s2,fR_s3,fR_s4) &
        !$omp& schedule(static)
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
        !$omp end parallel do
        
        
    END SUBROUTINE riemann_HLLC
        
    !**************************************************************
    SUBROUTINE ApplyBoundaryConditions(rho, ux, uy, p)
            IMPLICIT NONE
            REAL(DP), INTENT(INOUT) :: rho(:,:), ux(:,:), uy(:,:), p(:,:)
            INTEGER :: i, j,n
        
            !=============================
            ! OUTFLOW BOUNDARY CONDITIONS
            !=============================
            ! Left boundary (x = 1)
            DO j = 1, ny
              DO n = 0, ng-1
                rho(1+n, j) = rho(ng+1, j)
                ux(1+n, j) = ux(ng+1, j)
                uy(1+n, j) = uy(ng+1, j)
                p(1+n, j) = p(ng+1, j)
              END DO
            END DO
        
            ! Right boundary (x = nx)
            DO j = 1, ny
              DO n = 0, ng-1
                rho(nx-n, j) = rho(nx-ng, j)
                ux(nx-n, j) = ux(nx-ng, j)
                uy(nx-n, j) = uy(nx-ng, j)
                p(nx-n, j) = p(nx-ng, j)
              END DO
            END DO
        
            ! Bottom boundary (y = 1)
            DO i = 1, nx
              DO n = 0, ng-1
                rho(i, 1+n) = rho(i, ng+1)
                ux(i, 1+n) = ux(i, ng+1)
                uy(i, 1+n) = uy(i, ng+1)
                p(i, 1+n) = p(i, ng+1)
              END DO
            END DO
        
            ! Top boundary (y = ny)
            DO i = 1, nx
              DO n = 0, ng-1
                rho(i, ny-n) = rho(i, ny-ng)
                ux(i, ny-n) = ux(i, ny-ng)
                uy(i, ny-n) = uy(i, ny-ng)
                p(i, ny-n) = p(i, ny-ng)
              END DO
            END DO
        
            ! Corners
            DO i = 1, ng
              DO j = 1, ng
                rho(i, j) = rho(ng, ng)
                ux(i, j) = ux(ng, ng)
                uy(i, j) = uy(ng, ng)
                p(i, j) = p(ng, ng)
    
                rho(nx+1-i, j) = rho(nx-ng, ng)
                ux(nx+1-i, j) = ux(nx-ng, ng)
                uy(nx+1-i, j) = uy(nx-ng, ng)
                p(nx+1-i, j) = p(nx-ng, ng)
    
                rho(i, ny+1-j) = rho(ng, ny-ng)
                ux(i, ny+1-j) = ux(ng, ny-ng)
                uy(i, ny+1-j) = uy(ng, ny-ng)
                p(i, ny+1-j) = p(ng, ny-ng)
    
                rho(nx+1-i, ny+1-j) = rho(nx-ng, ny-ng)
                ux(nx+1-i, ny+1-j) = ux(nx-ng, ny-ng)
                uy(nx+1-i, ny+1-j) = uy(nx-ng, ny-ng)
                p(nx+1-i, ny+1-j) = p(nx-ng, ny-ng)
              END DO
            END DO
            
    END SUBROUTINE ApplyBoundaryConditions
   
    !**************************************************************
    SUBROUTINE build_snapshot_paths(frame_dir, scheme, solver, nsnap, data_path, meta_path)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nsnap
        CHARACTER(LEN=*), INTENT(IN) :: frame_dir, scheme, solver
        CHARACTER(LEN=*), INTENT(OUT) :: data_path, meta_path
        CHARACTER(LEN=256) :: base

        WRITE(base, '("rho_",A,"_",A,"_grid_",I0,"_frames_",I0)') TRIM(scheme), TRIM(solver), nx, nsnap
        data_path = TRIM(frame_dir)//"/"//TRIM(base)//".bin"
        meta_path = TRIM(frame_dir)//"/"//TRIM(base)//".meta"
    END SUBROUTINE build_snapshot_paths

    !**************************************************************
    SUBROUTINE write_snapshot_meta(meta_path, scheme, solver, nsnap)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nsnap
        CHARACTER(LEN=*), INTENT(IN) :: meta_path, scheme, solver

        OPEN(UNIT=19, FILE=TRIM(meta_path), STATUS="REPLACE", ACTION="WRITE", FORM="FORMATTED")
        WRITE(19, '("scheme=",A)') TRIM(scheme)
        WRITE(19, '("solver=",A)') TRIM(solver)
        WRITE(19, '("nx_total=",I0)') nx
        WRITE(19, '("ny_total=",I0)') ny
        WRITE(19, '("nx_phys=",I0)') nx - 2*ng
        WRITE(19, '("ny_phys=",I0)') ny - 2*ng
        WRITE(19, '("ng=",I0)') ng
        WRITE(19, '("nsnap=",I0)') nsnap
        WRITE(19, '("tfinal=",ES24.16)') tfinal
        WRITE(19, '("xl=",ES24.16)') xl
        WRITE(19, '("xr=",ES24.16)') xr
        WRITE(19, '("yl=",ES24.16)') yl
        WRITE(19, '("yr=",ES24.16)') yr
        WRITE(19, '("dtype=float64")')
        WRITE(19, '("order=y-major")')
        CLOSE(19)
    END SUBROUTINE write_snapshot_meta

    !**************************************************************
    SUBROUTINE write_density_snapshot(snap_unit, rho)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: snap_unit
        REAL(DP), INTENT(IN) :: rho(:,:)
        INTEGER :: i, j

        WRITE(snap_unit) ((rho(i,j), i=ng+1,nx-ng), j=ng+1,ny-ng)
    END SUBROUTINE write_density_snapshot

    !**************************************************************
    SUBROUTINE RK3TVD(eu, cfl, ct, sharpness, rho, u_x, u_y, p, c, fh_x, fh_y, scheme, solver, nsnap, snap_unit)
            IMPLICIT NONE
            INTEGER i, j, next_snap
            INTEGER, INTENT(IN) :: nsnap, snap_unit
            REAL(DP)                  , INTENT(IN)    :: cfl, ct, sharpness
            CHARACTER(LEN=4)          , INTENT(IN)    :: scheme, solver
            REAL(DP), DIMENSION(:,:)  , INTENT(INOUT) :: rho, u_x, u_y, p, c
            REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: fh_x, fh_y, eu

            REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: eu0, eu1, eu2
            REAL(DP) :: snap_dt, snap_target
            
            ALLOCATE(eu0(nx,ny,nv), eu1(nx,ny,nv), eu2(nx,ny,nv))

            t = 0.d0
            next_snap = 1
            snap_dt = 0.0_dp
            IF (snap_unit > 0 .AND. nsnap > 0) THEN
                IF (nsnap > 1) snap_dt = tfinal / REAL(nsnap - 1, DP)
                CALL write_density_snapshot(snap_unit, rho)
                next_snap = 2
            END IF
            DO while (t<tfinal)
                print *, "Time: ",t
                wavespeed = maxval(SQRT(u_x**2 + u_y**2) + sqrt(gam*p/rho))
                print *, "wavespeed: ", wavespeed, "u:", maxval(u_x),"v:", maxval(u_y), "rho: ", maxval(rho) 
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
                eu1 = eu0

                CALL activateScheme(scheme, rho, u_x, u_y, p, fh_x, fh_y, ct, sharpness, solver)
                !$omp parallel do collapse(2) default(shared) private(i,j) schedule(static)
                DO i = 2, nx - 1
                    DO j = 2, ny - 1
                        eu1(i,j,:) = eu0(i,j,:) - dtodx*(fh_x(i,j,:) - fh_x(i-1,j,:)) - dtody*(fh_y(i,j,:) - fh_y(i,j-1,:))
                    END DO
                END DO
                !$omp end parallel do

                
                ! Stage 2
                CALL decomp(rho, u_x, u_y, p, eu1, c)
                CALL ApplyBoundaryConditions(rho, u_x, u_y, p)
                CALL solvec(rho, u_x, u_y, p, eu1)

                CALL activateScheme(scheme, rho, u_x, u_y, p, fh_x, fh_y, ct, sharpness, solver)
                eu2 = eu0
                !$omp parallel do collapse(2) default(shared) private(i,j) schedule(static)
                DO i = 2, nx - 1
                    DO j = 2, ny - 1
                        eu2(i,j,:) = (3.0_dp/4.0_dp)*eu0(i,j,:) + (1.0_dp/4.0_dp)*&
                        (eu1(i,j,:) - dtodx*(fh_x(i,j,:) - fh_x(i-1,j,:)) - dtody*(fh_y(i,j,:) - fh_y(i,j-1,:)))
                    END DO
                END DO
                !$omp end parallel do
                
                ! Stage 3
                CALL decomp(rho, u_x, u_y, p, eu2, c)
                CALL ApplyBoundaryConditions(rho, u_x, u_y, p)
                CALL solvec(rho, u_x, u_y, p, eu2)

                CALL activateScheme(scheme, rho, u_x, u_y, p, fh_x, fh_y, ct, sharpness, solver)
                !$omp parallel do collapse(2) default(shared) private(i,j) schedule(static)
                DO i = 2, nx - 1
                    DO j = 2, ny - 1
                        eu(i,j,:) = (1.0_dp/3.0_dp)*eu0(i,j,:) + &
                        (2.0_dp/3.0_dp)*(eu2(i,j,:) - dtodx*(fh_x(i,j,:) - fh_x(i-1,j,:)) - dtody*(fh_y(i,j,:) - fh_y(i,j-1,:)))
                    END DO
                END DO
                !$omp end parallel do
                
                CALL decomp(rho, u_x, u_y, p, eu, c)
                CALL ApplyBoundaryConditions(rho, u_x, u_y, p)
                CALL solvec(rho, u_x, u_y, p, eu)

                t = t + dt
                IF (snap_unit > 0 .AND. nsnap > 1) THEN
                    DO WHILE (next_snap <= nsnap)
                        snap_target = REAL(next_snap - 1, DP) * snap_dt
                        IF (t + 1.0e-12_dp < snap_target) EXIT
                        CALL write_density_snapshot(snap_unit, rho)
                        next_snap = next_snap + 1
                    END DO
                END IF
                if (t>=tfinal) then 
                    EXIT
                end if
            END DO
            IF (snap_unit > 0) THEN
                DO WHILE (next_snap <= nsnap)
                    CALL write_density_snapshot(snap_unit, rho)
                    next_snap = next_snap + 1
                END DO
            END IF
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
    SUBROUTINE output(x, y, rho, u_x, u_y, p, scheme, solver)
                    IMPLICIT NONE
                    INTEGER i, j
                    CHARACTER(LEN=50)            :: filename
                    CHARACTER(LEN=4), INTENT(IN) :: scheme, solver
                    REAL(DP)        , INTENT(IN) :: x(:), y(:), rho(:,:), u_x(:,:), u_y(:,:), p(:,:)
                    
                    WRITE(filename, '(A,A,A,A,I0,A)') scheme,"_",solver,"_grid_", nx, ".dat"
                    OPEN(UNIT=18, FILE=TRIM(filename), STATUS="REPLACE", ACTION="WRITE", FORM="FORMATTED")
                    DO i = 1, nx
                        DO j = 1, ny
                            WRITE(18,*) x(i), y(j), rho(i,j), u_x(i,j), u_y(i,j), p(i,j)
                        END DO
                    END DO
                    CLOSE(18)
                    
                END SUBROUTINE output
    END MODULE subroutines
            
    !**************************************************************
    PROGRAM myeuler1d 
        USE types_vars
        USE subroutines
        
        IMPLICIT NONE
        INTEGER                         :: grid_physical, nsnap, arg_count, ios
        REAL(DP)                        :: cfl, ct, sharpness, target_tfinal
        CHARACTER(LEN=4)                :: selected_scheme, solver
        CHARACTER(LEN=256)              :: arg, frame_dir
        
        !Number of ghost cells
        ng = 3
        grid_physical = 100
        nsnap = 10
        selected_scheme = "ALL "
        solver = "HLLC"
        frame_dir = "frames_preview"
        target_tfinal = 0.5d0

        arg_count = COMMAND_ARGUMENT_COUNT()
        IF (arg_count >= 1) THEN
            CALL GET_COMMAND_ARGUMENT(1, arg)
            selected_scheme = ADJUSTL(arg(1:MIN(4, LEN_TRIM(arg))))
        END IF
        IF (arg_count >= 2) THEN
            CALL GET_COMMAND_ARGUMENT(2, arg)
            READ(arg, *, IOSTAT=ios) grid_physical
            IF (ios /= 0) STOP "Invalid grid size argument"
        END IF
        IF (arg_count >= 3) THEN
            CALL GET_COMMAND_ARGUMENT(3, arg)
            READ(arg, *, IOSTAT=ios) nsnap
            IF (ios /= 0) STOP "Invalid snapshot count argument"
        END IF
        IF (arg_count >= 4) THEN
            CALL GET_COMMAND_ARGUMENT(4, frame_dir)
        END IF
        IF (arg_count >= 5) THEN
            CALL GET_COMMAND_ARGUMENT(5, arg)
            READ(arg, *, IOSTAT=ios) target_tfinal
            IF (ios /= 0) STOP "Invalid final time argument"
        END IF

        IF (selected_scheme == "ALL ") THEN
            CALL run_case("TENO")
            CALL run_case("WENO")
            CALL run_case("MUSC")
        ELSE
            CALL run_case(selected_scheme)
        END IF
                
    CONTAINS

        SUBROUTINE run_case(scheme)
            IMPLICIT NONE
            CHARACTER(LEN=4), INTENT(IN) :: scheme
            INTEGER :: snap_unit
            CHARACTER(LEN=512) :: data_path, meta_path
            REAL(DP), ALLOCATABLE :: x(:), y(:), rho(:,:),u_x(:,:), u_y(:,:), p(:,:), c(:,:)
            REAL(DP), ALLOCATABLE :: eu(:,:,:), fh_x(:,:,:), fh_y(:,:,:)

            nx = grid_physical + 2 * ng
            ny = grid_physical + 2 * ng
            
            ! Physical and domain/grid parameters
            gm1 = gam-1.d0
            gm1i = 1.d0/gm1
            nc_x = nx + 1
            nc_y = ny + 1
            nv = 4
            xl = -0.0d0
            xr = 1.0d0
            yl = -0.0d0
            yr = 1.0d0
            dx = (xr - xl) / REAL(grid_physical, DP)
            dy = (yr - yl) / REAL(grid_physical, DP)
            cfl = 0.5
            tfinal = target_tfinal
            ct = 10.d0**(-7.d0)
            sharpness = 2.d0

            ALLOCATE(x(nx), y(ny))
            ! Velocities array is defined as (x value | y values)
            ALLOCATE(rho(nx,ny),u_x(nx,ny), u_y(nx,ny), p(nx,ny), c(nx,ny))
            ! Converved values array is defined as (x co. | y co. | x values | y values)
            ALLOCATE(eu(nx,ny,nv), fh_x(nc_x,nc_y,nv), fh_y(nc_x,nc_y,nv))
            
            CALL grid(x, y)
            CALL init_riemann(x, y, rho, u_x, u_y, p)
            CALL ApplyBoundaryConditions(rho, u_x, u_y, p)
            CALL solvec(rho, u_x, u_y, p, eu)
            
            snap_unit = -1
            IF (nsnap > 0) THEN
                CALL build_snapshot_paths(frame_dir, scheme, solver, nsnap, data_path, meta_path)
                CALL write_snapshot_meta(meta_path, scheme, solver, nsnap)
                snap_unit = 40
                OPEN(UNIT=snap_unit, FILE=TRIM(data_path), STATUS="REPLACE", ACTION="WRITE", &
                     ACCESS="STREAM", FORM="UNFORMATTED")
            END IF

            PRINT *, "Running scheme: ", scheme, " grid: ", nx, " frames: ", nsnap
            CALL RK3TVD(eu, cfl, ct, sharpness, rho, u_x, u_y, p, c, fh_x, fh_y, scheme, solver, nsnap, snap_unit)
            IF (snap_unit > 0) CLOSE(snap_unit)
            CALL output(x, y, rho, u_x, u_y, p, scheme, solver)
            
            DEALLOCATE(x, y)
            DEALLOCATE(rho, u_x, u_y, p, c)
            DEALLOCATE(eu, fh_x, fh_y)
        END SUBROUTINE run_case
                 
    END PROGRAM myeuler1d
