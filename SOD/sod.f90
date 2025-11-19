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
    
    INTEGER :: nx, nc_x, ng, nv, nt, ntsteps, nf, kf
    REAL(DP) :: xl, xr, t, tfinal, dt, dx, dtodx, dxodt, gam=1.4d0, gm1, gm1i, wavespeed
    
    CONTAINS
    
    !**************************************************************
    SUBROUTINE grid(x)
        IMPLICIT NONE
        INTEGER :: i
        REAL(DP),INTENT(OUT) :: x(:)
        
        ! This subroutine generates a finite-volume grid in 1D
        ! Specifically, x(i) contains the coordinates of the 
        ! cell centers starting at h/2
        ! |---x---|---x---|.....|---x---|---x---|
        
        DO i = 1, nx
            x(i) = xl + REAL(i-1,DP)*dx + 0.5d0*dx
        END DO
        
    END SUBROUTINE grid
    
    !**************************************************************
    
    SUBROUTINE init_shocktube(x, rho, ux, p)
        IMPLICIT NONE
        REAL(DP), INTENT(in), DIMENSION(:) :: x
        REAL(DP), INTENT(out), DIMENSION(:) :: rho,ux,p
        INTEGER :: i
        ! Sod shock tube problem
        do i = 1,nx
            if (x(i) <= 0.5d0) then
                rho(i) = 1.0d0
                ux(i)  = 0.0d0
                p(i)   = 1.0d0
            else
                rho(i) = 0.125d0
                ux(i)  = 0.0d0
                p(i)   = 0.1d0
            end if
        end do
    END SUBROUTINE init_shocktube
    
    !**************************************************************
    
    SUBROUTINE solvec(rho, u_x, p, eu)
        IMPLICIT NONE
        REAL(DP), INTENT(IN)  :: rho(:), u_x(:), p(:)
        REAL(DP), INTENT(OUT) :: eu(:,:)
        INTEGER               :: i, j
        
        ! This subroutine constructs the conserved variable
        ! or solution vector from the primitive variables
        !
        ! U = [rho, rho*u, E]^T
        ! p = (gam-1)*(E-0.5*rho*u^2)
        ! E = p/(gam-1) + 0.5*rho*u^2
        
        DO i = 1, nx
            eu(i, 1) = rho(i)
            eu(i, 2) = rho(i) * u_x(i)
            eu(i, 3) = p(i) * gm1i + 0.5d0 * rho(i) * u_x(i)**2
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
    
    SUBROUTINE decomp(rho, u_x, p, eu, c)
        IMPLICIT NONE
        REAL(DP), INTENT(IN), DIMENSION(:,:) :: eu
        REAL(DP), INTENT(OUT), DIMENSION(:) :: rho, u_x, p, c
        INTEGER :: i, j
        
        ! This suborutine decomposes the solution vector U back
        ! to the primitive variables according to the defintion
        ! for 1D Euler equations
        !
        ! rho, rho*u, E
        ! p = (gam-1)*(E-0.5*rho*u^2)
        ! E = p/(gam-1) + 0.5*rho*u^2
        
        do i = 1, nx
            rho(i) = eu(i,1)
            u_x(i) = eu(i,2)/rho(i)
            p(i)   = gm1*(eu(i,3) - 0.5d0*rho(i)*u_x(i)**2)
            c(i)   = dsqrt(gam*p(i)/rho(i))
        end do
        
    END SUBROUTINE decomp
    
    !**************************************************************
    SUBROUTINE TENOFlux(fh_x, rho, u_x, p, ct, solver)
        IMPLICIT NONE
        INTEGER                      :: k, i, j
        REAL(DP), INTENT(IN)         :: rho(:), u_x(:), p(:), ct
        CHARACTER(LEN=4), INTENT(IN) :: solver
        REAL(DP), INTENT(OUT)        :: fh_x(:,:)
        REAL(DP), ALLOCATABLE        :: qL_x(:,:), qR_x(:,:), q(:,:), dqL(:,:), dqR(:,:), dq(:,:)
        REAL(DP)                     :: beta(3), omega(3), cell_values(5), p_stencils(5)
        ALLOCATE(q(nx,nv),qL_x(nx,nv),qR_x(nx,nv), dqL(nx,nv), dqR(nx,nv), dq(nx,nv))
        q(:,:) = 0.0_dp
        qL_x(:,:) = 0.0_dp
        qR_x(:,:) = 0.0_dp
        dqL(:,:) = 0.0_dp
        dqR(:,:) = 0.0_dp

        !Build primitive variables array
        DO i = 1, nx
            q(i, 1) = rho(i)
            q(i, 2) = u_x(i)
            q(i, 3) = p(i)
        END DO
        
        DO k = 1, nv                ! Loop over rho, u, v and p components
            DO i = 3, nx - 3
                ! Assemble left-biased stencil cell primitive values
                cell_values = (/q(i-2, k), q(i-1, k), q(i, k), q(i+1, k), q(i+2, k)/)
                
                ! Left-biased stencil primitive values before weights
                CALL stencilsValues(cell_values, p_stencils)
                
                ! Smoothness indicators for left-biased stencil
                CALL smoothnessIndicators(cell_values, beta)
                
                ! Compute non-linear weights
                CALL TENOWeights(beta, ct, omega)
                
                ! Compute the left-biased primitive values at the interface
                qL_x(i, k) = omega(1) * p_stencils(1) + omega(2) * p_stencils(2) + omega(3) * p_stencils(3)
                
                ! Assemble left-biased stencil cell primitive values
                cell_values = (/q(i+3, k), q(i+2, k), q(i+1, k), q(i, k), q(i-1, k)/)
                
                ! Right-biased stencil primitive values before weights
                CALL stencilsValues(cell_values, p_stencils)
                
                ! Smoothness indicators for right-biased stenci
                CALL smoothnessIndicators(cell_values, beta)
                
                ! Compute non-linear weights
                CALL TENOWeights(beta, ct, omega)
                
                ! Compute the right-biased primitive values at the interface
                qR_x(i, k) = omega(1) * p_stencils(1) + omega(2) * p_stencils(2) + omega(3) * p_stencils(3)  
            END DO    
        END DO

        ! Perform interpolation over ghost cells
        CALL WTENOghostCellsInterpolation(qL_x, qR_x)

        ! Call approximate Riemann solver
        CALL riemann_solver(qL_x, qR_x, fh_x, solver)
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
    SUBROUTINE WENOFlux(fh_x, rho, u_x, p, solver)
        IMPLICIT NONE
        INTEGER                      :: k, i
        REAL(DP), INTENT(IN)         :: rho(:), u_x(:), p(:)
        CHARACTER(LEN=4), INTENT(IN) :: solver
        REAL(DP), INTENT(OUT)        :: fh_x(:,:)
        REAL(DP), ALLOCATABLE        :: qL_x(:,:), qR_x(:,:), q(:,:), dqL(:,:), dqR(:,:), dq(:,:)
        REAL(DP)                     :: beta(3), omega(3), cell_values(5), p_stencils(5)
        ALLOCATE(q(nx,nv),qL_x(nx,nv),qR_x(nx,nv), dqL(nx,nv),dqR(nx,nv),dq(nx,nv))
        q(:,:) = 0.0_dp
        qL_x(:,:) = 0.0_dp
        qR_x(:,:) = 0.0_dp
        dqL(:,:) = 0.0_dp
        dqR(:,:) = 0.0_dp

        !Build primitive variables array
        DO i = 1, nx
            q(i, 1) = rho(i)
            q(i, 2) = u_x(i)
            q(i, 3) = p(i)
        END DO
        
        DO k = 1, nv                ! Loop over rho, u, v and p components
            DO i = 3, nx - 3
                ! Assemble left-biased stencil cell primitive values
                cell_values = (/q(i-2, k), q(i-1, k), q(i, k), q(i+1, k), q(i+2, k)/)
                
                ! Left-biased stencil primitive values before weights
                CALL stencilsValues(cell_values, p_stencils)
                
                ! Smoothness indicators for left-biased stencil
                CALL smoothnessIndicators(cell_values, beta)
                
                ! Compute non-linear weights
                CALL WENOWeights(beta, omega)
                
                ! Compute the left-biased primitive values at the interface
                qL_x(i, k) = omega(1) * p_stencils(1) + omega(2) * p_stencils(2) + omega(3) * p_stencils(3)
                
                ! Assemble left-biased stencil cell primitive values
                cell_values = (/q(i+3, k), q(i+2, k), q(i+1, k), q(i, k), q(i-1, k)/)
                
                ! Right-biased stencil primitive values before weights
                CALL stencilsValues(cell_values, p_stencils)
                
                ! Smoothness indicators for right-biased stenci
                CALL smoothnessIndicators(cell_values, beta)
                
                ! Compute non-linear weights
                CALL WENOWeights(beta, omega)
                
                ! Compute the right-biased primitive values at the interface
                qR_x(i, k) = omega(1) * p_stencils(1) + omega(2) * p_stencils(2) + omega(3) * p_stencils(3)  
            END DO    
        END DO

        ! Perform interpolation over ghost cells
        CALL WTENOghostCellsInterpolation(qL_x, qR_x)

        ! Call approximate Riemann solver
        CALL riemann_solver(qL_x, qR_x, fh_x, solver)
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
    SUBROUTINE WTENOghostCellsInterpolation(qL_x, qR_x)
        IMPLICIT NONE
        REAL(DP), DIMENSION(:,:), INTENT(INOUT)  :: qL_x, qR_x
        
        qL_x(2,:)     = qL_x(3,:)
        qL_x(1,:)     = qL_x(2,:)
        qR_x(2,:)     = qR_x(3,:)
        qR_x(1,:)     = qR_x(2,:)
        qL_x(nx-2,:)  = qL_x(nx-3,:)
        qL_x(nx-1,:)  = qL_x(nx-2,:)
        qL_x(nx,:)    = qL_x(nx-1,:)
        qR_x(nx-2,:)  = qR_x(nx-3,:)
        qR_x(nx-1,:)  = qR_x(nx-2,:)
        qR_x(nx,:)    = qR_x(nx-1,:)

    END SUBROUTINE WTENOghostCellsInterpolation
    
    !**************************************************************
    SUBROUTINE MUSCL_THINC(fh_x, rho, u_x, p, beta, solver)
      
        IMPLICIT NONE
        INTEGER :: k, i, j
        REAL(DP), INTENT(IN) :: rho(:), u_x(:), p(:), beta
        CHARACTER(4), INTENT(IN):: solver
        REAL(DP), INTENT(INOUT) :: fh_x(:,:)
        REAL(DP) :: qi,qip1,qim1
        REAL(DP) :: num,den,r
        REAL(DP) :: qLiphM,qRimhM
        REAL(DP) :: qLiphT,qRimhT
        REAL(DP) :: qLiphMT,qRimhMT
        REAL(DP) :: qmin,delq,theta,A,B,term
        REAL(DP) :: dpdx,drhodx,phip,phirho,phiprho,xsiph
        REAL(DP) :: eps, eta, arg1, arg2, xi
        REAL(DP), ALLOCATABLE :: q(:,:), qL_x(:,:), qR_x(:,:), qL_y(:,:), qR_y(:,:)
      
        ALLOCATE(q(nx,nv),qL_x(nx,nv),qR_x(nx,nv), qL_y(nx,nv),qR_y(nx,nv))
      
        eta = 1.d0/3.d0
        eps = 1e-20
      
        !Build primitive variables array
        DO i = 1, nx
            q(i, 1) = rho(i)
            q(i, 2) = u_x(i)
            q(i, 3) = p(i)
        END DO
      
        ! X reconstruction
        DO k = 1, nv
            ! Loop over each control volume
            DO i = 2, nx - 1
            
                qi   = q(i,k)
                qip1 = q(i+1,k)
                qim1 = q(i-1,k)
        
                num = qip1 - qi
                den = qi - qim1
                r = num/(den + 1e-08)
        
                ! MUSCL cell edge predictors
                qLiphM = qi + 0.25d0*((1.d0-eta)*psi(r)*(qi-qim1) &
                                +         (1.d0+eta)*psi(1.d0/(r+1e-08))*(qip1-qi))
                qRimhM = qi - 0.25d0*((1.d0+eta)*psi(r)*(qi-qim1) &
                                -         (1.d0-eta)*psi(1.d0/(r+1e-08))*(qip1-qi))
                qL_x(i,k)   = qLiphM
                qR_x(i-1,k) = qRimhM
        
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
                dpdx   = (p(i+1)-p(i))/dx
                drhodx = (rho(i+1)-rho(i))/dx
                phip     = max(p(i),p(i+1))/min(p(i),p(i+1))
                phirho   = max(rho(i),rho(i+1))/min(rho(i),rho(i+1))
                if (dpdx*drhodx >= 0) then
                    phiprho = phip/(phirho+1e-08)
                    xsiph = exp(-25.d0*(max(1.d0,phiprho)-1.d0))
                else
                    xsiph = 1.d0
                end if
        
            ! MUSCL-THINC hybrid cell edge predictor
                qLiphMT = (1.d0-xi*xsiph)*qLiphM + xi*xsiph*qLiphT
                qRimhMT = (1.d0-xi*xsiph)*qRimhM + xi*xsiph*qRimhT
        
                qL_x(i,k)   = qLiphMT
                qR_x(i-1,k) = qRimhMT
                end if
            END DO            
        END DO

        qL_x(1,:) = qL_x(2,:)
        qL_x(nx,:) = qL_x(nx-1,:)

        ! Call approximate Riemann solver
        CALL riemann_solver(qL_x, qR_x, fh_x, solver)
      
    END SUBROUTINE MUSCL_THINC

    REAL FUNCTION psi(arg)
        IMPLICIT NONE
        REAL(DP) :: arg
    
        ! This is needed for the MUSCL-THINC scheme
        psi = max(0.d0,min(1.d0,arg))
    
        RETURN
    END FUNCTION psi
    !**************************************************************
    SUBROUTINE riemann_HLLC(qL_x, qR_x, fh_x)
        IMPLICIT NONE
        INTEGER :: i
        REAL(DP), INTENT(IN), DIMENSION(:,:)  :: qL_x, qR_x
        REAL(DP), INTENT(INOUT), DIMENSION(:,:)  :: fh_x
        REAL(DP) :: pR, pL, cL, cR, eps, hL, hR
        REAL(DP) :: uL1,uL2,uL3
        REAL(DP) :: uR1,uR2,uR3
        REAL(DP) :: fL1,fL2,fL3
        REAL(DP) :: fR1,fR2,fR3
        REAL(DP) :: SL, SR, sStar, delta_p,pStar
        REAL(DP) :: uL, uR
        REAL(DP) :: uL_cof, uR_cof
        REAL(DP) :: uL_s1, uL_s2, uL_s3, fL_s1, fL_s2, fL_s3
        REAL(DP) :: uR_s1, uR_s2, uR_s3, fR_s1, fR_s2, fR_s3
        
        eps = 1.0e-6_dp
        
        DO i = 1, nx
            uL1 = max(qL_x(i,1), eps)
            uL2 = uL1 * qL_x(i,2) 
            uL3 = qL_x(i, 3) * gm1i + 0.5d0 * uL1 * qL_x(i,2)**2
            
            uR1 = max(qR_x(i,1), eps)
            uR2 = uR1 * qR_x(i,2)
            uR3 = qR_x(i,3) * gm1i + 0.5d0 * uR1 * qR_x(i,2)**2
            
            uL = qL_x(i,2)
            uR = qR_x(i,2)
            
            pL = MAX(gm1 * (uL3 - 0.5d0 * uL2**2 / (uL1 + eps)), eps)
            pR = MAX(gm1 * (uR3 - 0.5d0 * uR2**2 / (uR1 + eps)), eps)
            
            ! Compute left and right enthalpy
            hL = (uL3 + pL) / uL1
            hR = (uR3 + pR) / uR1
            
            cL = SQRT(gam * pL / uL1)
            cR = SQRT(gam * pR / uR1)
            
            ! Pressure-based entropy fix
            delta_p = 0.1_dp * ABS(pR - pL) / MAX(pL, pR, eps)
            SL = MIN(uL - cL, uR - cR) - delta_p
            SR = MAX(uL + cL, uR + cR) + delta_p
            
            fL1 = uL2                               ! rho_L * uL
            fL2 = uL2**2 / uL1 + pL                 ! rho_L * uL^2 + p_L
            fL3 = uL2 / uL1 * (uL3 + pL)            ! uL * (E_L + p_L)
            
            fR1 = uR2                               ! rho_R * uR
            fR2 = uR2**2 / uR1 + pR                 ! rho_R * uR^2 + p_R
            fR3 = uR2 / uR1 * (uR3 + pR)            ! uR * (E_R + p_R)
            
            IF (SL > 0.0_dp) THEN
                fh_x(i, 1) = fL1
                fh_x(i, 2) = fL2
                fh_x(i, 3) = fL3
            ELSE IF (SR < 0.0_dp) THEN
                fh_x(i, 1) = fR1
                fh_x(i, 2) = fR2
                fh_x(i, 3) = fR3
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
                uL_s3 = uL_cof * (uL3 / uL1 + (sStar - uL) * (sStar + pStar / (uL1 * (SL - uL) + eps)))
                
                fL_s1 = fL1 + SL * (uL_s1 - uL1)
                fL_s2 = fL2 + SL * (uL_s2 - uL2)
                fL_s3 = fL3 + SL * (uL_s3 - uL3)
                
                ! Right star region
                uR_cof = uR1 * (SR - uR) / (SR - sStar + eps)
                uR_s1 = uR_cof
                uR_s2 = uR_cof * sStar
                uR_s3 = uR_cof * (uR3 / uR1 + (sStar - uR) * (sStar + pStar / (uR1 * (SR - uR) + eps)))
                
                fR_s1 = fR1 + SR * (uR_s1 - uR1)
                fR_s2 = fR2 + SR * (uR_s2 - uR2)
                fR_s3 = fR3 + SR * (uR_s3 - uR3)
                
                ! If 0 <= sStar => use left star
                IF (sStar >= 0.0_dp) THEN
                    fh_x(i,1) = fL_s1
                    fh_x(i,2) = fL_s2
                    fh_x(i,3) = fL_s3
                ELSE
                    fh_x(i,1) = fR_s1
                    fh_x(i,2) = fR_s2
                    fh_x(i,3) = fR_s3
                END IF
            END IF
        END DO
    END SUBROUTINE riemann_HLLC
        
    !**************************************************************
    SUBROUTINE ApplyBoundaryConditions(rho, ux, p)
            IMPLICIT NONE
            REAL(DP), INTENT(INOUT) :: rho(:), ux(:), p(:)
            INTEGER :: n
        
            !=============================
            ! OUTFLOW BOUNDARY CONDITIONS
            !=============================
            ! Left boundary (x = 1)
            DO n = 0, ng-1
                rho(1+n) = rho(ng+1)
                ux(1+n) = ux(ng+1)
                p(1+n) = p(ng+1)
            END DO
        
            ! Right boundary (x = nx)
            DO n = 0, ng-1
                rho(nx-n) = rho(nx-ng)
                ux(nx-n) = ux(nx-ng)
                p(nx-n) = p(nx-ng)
            END DO            
    END SUBROUTINE ApplyBoundaryConditions
   
    !**************************************************************
    SUBROUTINE RK3TVD(eu, cfl, ct, sharpness, rho, u_x, p, c, fh_x, scheme, solver)
            IMPLICIT NONE
            INTEGER i
            REAL(DP)                  , INTENT(IN)    :: cfl, ct, sharpness
            CHARACTER(LEN=4)          , INTENT(IN)    :: scheme, solver
            REAL(DP), DIMENSION(:)  , INTENT(INOUT) :: rho, u_x, p, c
            REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: fh_x, eu

            REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: eu0, eu1, eu2
            
            ALLOCATE(eu0(nx,nv), eu1(nx,nv), eu2(nx,nv))

            t = 0.d0
            DO while (t<tfinal)
                print *, "Time: ",t
                wavespeed = maxval(abs(u_x) + sqrt(gam*p/rho))
                print *, "wavespeed: ", wavespeed, "u:", maxval(u_x), "rho: ", maxval(rho) 
                dt = cfl * dx / abs(wavespeed)
                
                if(t+dt > tfinal) then
                    dt = tfinal-t
                end if
                dtodx = dt/dx
                dxodt = dx/dt
                
                ! Stage 1
                eu0 = eu               ! Save old solution

                CALL activateScheme(scheme, rho, u_x, p, fh_x, ct, sharpness, solver)
                DO i = 2, nx - 1
                    eu1(i,:) = eu0(i,:) - dtodx*(fh_x(i,:) - fh_x(i-1,:))
                END DO

                
                ! Stage 2
                CALL decomp(rho, u_x, p, eu1, c)
                CALL ApplyBoundaryConditions(rho, u_x, p)
                CALL solvec(rho, u_x, p, eu1)

                CALL activateScheme(scheme, rho, u_x, p, fh_x, ct, sharpness, solver)
                DO i = 2, nx - 1
                    eu2(i,:) = (3.0_dp/4.0_dp)*eu0(i,:) + (1.0_dp/4.0_dp)*&
                        (eu1(i,:) - dtodx*(fh_x(i,:) - fh_x(i-1,:)))
                END DO
                
                ! Stage 3
                CALL decomp(rho, u_x, p, eu2, c)
                CALL ApplyBoundaryConditions(rho, u_x, p)
                CALL solvec(rho, u_x, p, eu2)

                CALL activateScheme(scheme, rho, u_x, p, fh_x, ct, sharpness, solver)
                DO i = 2, nx - 1
                    eu(i,:) = (1.0_dp/3.0_dp)*eu0(i,:) + &
                        (2.0_dp/3.0_dp)*(eu2(i,:) - dtodx*(fh_x(i,:) - fh_x(i-1,:)))
                END DO
                
                CALL decomp(rho, u_x, p, eu, c)
                CALL ApplyBoundaryConditions(rho, u_x, p)
                CALL solvec(rho, u_x, p, eu)

                t = t + dt
                if (t>=tfinal) then 
                    EXIT
                end if
            END DO
            DEALLOCATE(eu0, eu1, eu2)
    END SUBROUTINE RK3TVD
        
    !**************************************************************
    SUBROUTINE activateScheme(scheme, rho, u_x, p, fh_x, ct, sharpness, solver)
            IMPLICIT NONE
            REAL(DP)        , INTENT(IN) :: ct, sharpness, rho(:), u_x(:), p(:)
            CHARACTER(LEN=4), INTENT(IN) :: scheme, solver
            REAL(DP)        , INTENT(OUT):: fh_x(:,:)
            
            ! Determine which scheme to activate.
            IF (scheme == "TENO") THEN
                CALL TENOFlux(fh_x, rho, u_x, p, ct, solver)

             ELSEIF (scheme == "WENO") THEN
                CALL WENOFlux(fh_x, rho, u_x, p, solver)

             ELSEIF (scheme == "MUSC") THEN
                 CALL MUSCL_THINC(fh_x, rho, u_x, p, sharpness, solver)
            END IF    
    END SUBROUTINE activateScheme
                
    !**************************************************************
    SUBROUTINE riemann_solver(qL_x, qR_x, fh_x, solver)
        IMPLICIT NONE
        REAL(DP)        , INTENT(IN) :: qL_x(:,:), qR_x(:,:)
        CHARACTER(LEN=4), INTENT(IN) :: solver
        REAL(DP)        , INTENT(OUT):: fh_x(:,:)
        
        ! Determine which solver to activate.
        IF (solver == "HLLC") THEN
            CALL riemann_HLLC(qL_x, qR_x, fh_x)
        END IF  
    END SUBROUTINE riemann_solver
            
    !**************************************************************
    SUBROUTINE output(x, rho, u_x, p, scheme, solver)
                    IMPLICIT NONE
                    INTEGER i
                    CHARACTER(LEN=50)            :: filename
                    CHARACTER(LEN=4), INTENT(IN) :: scheme, solver
                    REAL(DP)        , INTENT(IN) :: x(:), rho(:), u_x(:), p(:)
                    
                    WRITE(filename, '(A,A,A,A,I0,A)') scheme,"_",solver,"_grid_", nx, ".dat"
                    OPEN(UNIT=18, FILE=TRIM(filename), STATUS="REPLACE", ACTION="WRITE", FORM="FORMATTED")
                    DO i = 1, nx
                        WRITE(18,*) x(i), rho(i), u_x(i), p(i)
                    END DO
                    CLOSE(18)
                    
                END SUBROUTINE output
    END MODULE subroutines
            
    !**************************************************************
    PROGRAM myeuler1d 
                USE types_vars
                USE subroutines
                
                IMPLICIT NONE
                INTEGER                                 :: scheme_index, solver_index, grid_index
                REAL(DP)                                :: cfl, ct, sharpness
                REAL(DP), DIMENSION(:), ALLOCATABLE     :: x, rho, u_x, p, c
                REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: eu, fh_x
                INTEGER, DIMENSION(1)                   :: grid_size = (/800/)
                CHARACTER(LEN=4)                        :: scheme(3) = (/"TENO", "MUSC", "WENO"/), solver(1) = (/"HLLC"/)
                
                !Number of ghost cells
                ng = 3

                DO grid_index = 1, SIZE(grid_size)
                    nx = grid_size(grid_index) + 2 * ng
                    
                    ! Physical and domain/grid parameters
                    gm1 = gam-1.d0
                    gm1i = 1.d0/gm1
                    nc_x = nx + 1
                    nv = 3
                    xl = -0.0d0
                    xr = 1.0d0
                    dx = (xr - xl) / (nx - 1)
                    cfl = 0.5
                    tfinal = 0.2
                    ct = 10.d0**(-5.d0)
                    sharpness = 2.d0

                    ! Run TENO, WENO and MUSCL schemes.
                    DO scheme_index = 1, SIZE(scheme)
                        DO solver_index = 1, SIZE(solver)
                            ALLOCATE(x(nx), rho(nx), u_x(nx), p(nx), c(nx))
                            ! Converved values array is defined as (x co. | y values)
                            ALLOCATE(eu(nx, nv), fh_x(nc_x, nv))
                            
                            CALL grid(x)
                            CALL init_shocktube(x, rho, u_x, p)
                            CALL solvec(rho, u_x, p, eu)
                            
                            
                            CALL RK3TVD(eu, cfl, ct, sharpness, rho, u_x, p, c, fh_x, &
                             scheme(scheme_index), solver(solver_index))
                            CALL output(x, rho, u_x, p, scheme(scheme_index), solver(solver_index))
                            
                            DEALLOCATE(x)
                            DEALLOCATE(rho, u_x, p, c)
                            DEALLOCATE(eu, fh_x)
                        END DO
                    END DO
                END DO
                
    END PROGRAM myeuler1d