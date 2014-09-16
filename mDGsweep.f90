! --------------------------------------------------------------------
! Modal DG sweep update for DG/FV hybrid 1D Advection
! By: Devin Light
! --------------------------------------------------------------------

SUBROUTINE mDGsweep(rhoq,u,uedge,dxel,nelem,N,wghts,DG_C,DG_LUC,DG_L,DG_DL,IPIV,dt, & 
			        doposlimit)
	USE mDGmod
	IMPLICIT NONE
	
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)

	! --- External functions
	REAL(KIND=DOUBLE), EXTERNAL :: B ! RHS function for evolution ODE for kth expansion coefficent

	! --- Inputs
	INTEGER, INTENT(IN) :: nelem,N ! Number of elements, highest degree of Legendre polynomial being used
	REAL(KIND=DOUBLE), INTENT(IN) :: dxel,dt ! Element width, Finite Volume sub-cell width, time step
	REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: wghts ! Quadrature weights and node locations
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:N), INTENT(IN) :: DG_C, DG_LUC,DG_L,DG_DL! C matrix, used to xfer between grids, LU decomp of C
	INTEGER, DIMENSION(0:N), INTENT(IN):: IPIV ! Pivot array for RHS when using DG_LUC
	REAL(KIND=DOUBLE), DIMENSION(1:3,1:(nelem*(N+1))), INTENT(IN) :: u ! Velocities at quadrature locations within each element
	REAL(KIND=DOUBLE), DIMENSION(1:3,1:nelem), INTENT(IN) :: uedge ! Edge velocities at RHS of each element
	LOGICAL, INTENT(IN) :: doposlimit

	! --- Outputs
	REAL(KIND=DOUBLE), DIMENSION(1:(nelem*(N+1))), INTENT(INOUT) :: rhoq ! Soln as sub-cell averages within each element at FV cell centers

	! --- Local Variables
	INTEGER :: i,j,k,stage
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: rqBAR! Reshaped cell averages
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: A,A1,A2 ! Reshaped DG coefficent matricies
	REAL(KIND=DOUBLE), DIMENSION(1:3,0:N,1:nelem) :: utild ! Reshaped velocities
	REAL(KIND=DOUBLE), DIMENSION(1:3,0:nelem) :: uedgetild ! Periodically extended edge velocities
	REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: flxrq ! Array of fluxes F(j,j+1) (flx(0) is left flux at left edge of domain)
	REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: fcfrq ! Flux correction factors for positivity limiting

	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: rqQuadVals
	REAL(KIND=DOUBLE), DIMENSION(0:1,0:nelem+1) :: rqEdgeVals
    REAL(KIND=DOUBLE), DIMENSION(0:N) :: currElemQuad,currElemU
    REAL(KIND=DOUBLE), DIMENSION(0:N) :: currLegendreDeriv
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: uTmpQuad
	REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: uTmpEdge

    REAL(KIND=DOUBLE) :: t0,tf,t1,t2

	! #####################################################
    ! A(k,j) gives kth coeff in the jth element for rhoq 
    ! #####################################################

    ! Reform incoming values to be more convienent
    DO j=1,nelem
        rqBAR(:,j) = rhoq(1+(N+1)*(j-1) : (N+1)*j)
        utild(1:3,:,j) = u(1:3,1+(N+1)*(j-1) : (N+1)*j)
    END DO

	uedgetild(1:3,1:nelem) = uedge(1:3,1:nelem)
	uedgetild(1:3,0) = uedge(1:3,nelem)

CALL CPU_TIME(t0)

	CALL projectAverages(A,DG_LUC,IPIV,rqBAR,N,nelem) ! Project incoming rhoq averages
	A1 = A

CALL CPU_TIME(t1)

	DO stage = 1,3
		uTmpQuad(:,:) = utild(stage,:,:)
		uTmpEdge(:) = uedgetild(stage,:)

		CALL evalExpansion(A1,DG_L,rqQuadVals,rqEdgeVals,N,nelem)
		CALL NUMFLUX(rqEdgeVals,uTmpEdge,nelem,flxrq)

		fcfrq = 1D0
		IF(doposlimit) THEN
			CALL FLUXCOR(A1,A,flxrq,DG_C,dxel,dt,nelem,N,stage,fcfrq)
		ENDIF

		! Forward step
		DO j=1,nelem
            currElemQuad = rqQuadVals(:,j)
            currElemU = uTmpQuad(:,j)
			DO k=0,N
                currLegendreDeriv = DG_DL(k,:)
				A2(k,j) = A1(k,j) + (dt/dxel)*B(currElemQuad,flxrq,currElemU,wghts,k,j,nelem,N,fcfrq,currLegendreDeriv)
			ENDDO
		ENDDO

		SELECT CASE(stage)
		CASE(1)
			A1 = A2
		CASE(2)
			A1 = 0.75D0*A + 0.25D0*A2
		CASE(3)
			A1 = A/3d0 + 2D0*A2/3D0
		END SELECT
	ENDDO ! stage

CALL CPU_TIME(t2)
    
    ! Average modal coefficients to get back to subcell averages
    DO j=1,nelem
    		rqBAR(:,j) = MATMUL(DG_C,A1(:,j))
    ENDDO

    IF(doposlimit) THEN
!        tmpMass = SUM(qBAR)
		CALL MFILL(rqBAR,N,nelem)
!        IF( abs(tmpMass - SUM(qBAR)) .gt. 1D-13) THEN
!            write(*,FMT='(A15,e12.4)') 'CHG After fill:',tmpMass - SUM(qBAR)
!        ENDIF
    ENDIF

CALL CPU_TIME(tf)
!write(*,FMT='(A11,e12.5,f7.2,A11)') 'Proj time: ',t1-t0,((t1-t0)/(tf-t0))*100D0,' % of total'
!write(*,FMT='(A11,e12.5,f7.2,A11)') 'Step time: ',t2-t1,((t2-t1)/(tf-t0))*100D0,' % of total'
!write(*,*) ''

	! Reform original shaped arrays
    DO j=1,nelem
        rhoq(1+(N+1)*(j-1) : (N+1)*j) = rqBAR(:,j)
    END DO

END SUBROUTINE mDGsweep

REAL(KIND=KIND(1D0)) FUNCTION B(quadVals,flx,uQuad,wghts,k,j,nelem,N,fluxcf,dLeg) 
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)

	! --- Inputs
	INTEGER, INTENT(IN) :: nelem, N
	INTEGER, INTENT(IN) :: k,j ! Mode number, element number
	REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: wghts,uQuad,dLeg
	REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: quadVals

	REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(IN) :: flx ! Array of element edge fluxes flx(j) = F(j,j+1)
	
	! Flux correction factor. fluxcf(j)=reduction in flux coming through face (j+1/2)
	! fluxcf(0)=reduction factor for flux coming through left boundary of domain.
	REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(IN) :: fluxcf
 
	B = SUM(wghts(:)*uQuad(:)*dLeg(:)*quadVals(:))

	IF(k .eq. 0) THEN
		B = B - fluxcf(j)*flx(j) + fluxcf(j-1)*flx(j-1)
	ELSE
		B = B - flx(j) + ((-1D0)**k)*flx(j-1)
	END IF

	B = (2D0*k+1D0)*B

END FUNCTION B

SUBROUTINE NUMFLUX(rqEdgeVals,uEdge,nelem,flxrq) 
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: nelem
	REAL(KIND=DOUBLE), DIMENSION(0:1,0:nelem+1), INTENT(IN) :: rqEdgeVals
	REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(IN) :: uEdge

	! -- Outputs	
	REAL(KIND=DOUBLE),DIMENSION(0:nelem), INTENT(OUT) :: flxrq

	! -- Local variables
	INTEGER :: j,which_el,which_sign

	DO j=0,nelem
        which_sign = 1-(1-(SIGN(1D0,uEdge(j))) )/2
        which_el = j + (1-(SIGN(1D0,uEdge(j))) )/2
!		flxrq(j) = 0.5D0*rqEdgeVals(1,j)*(uEdge(j)+DABS(uEdge(j)))+0.5D0*rqEdgeVals(0,j+1)*(uEdge(j)-DABS(uEdge(j)))
		flxrq(j) = uEdge(j)*rqEdgeVals(which_sign,which_el)
	ENDDO

END SUBROUTINE NUMFLUX

SUBROUTINE FLUXCOR(Acur,Apre,flx,DG_C,dxel,dt,nelem,N,substep,fluxcf)
	! Computes flux reductions factors to prevent total mass within each element from going negative
	! Outputs fluxcf. fluxcf(j) is the reduction factor for the right face of element j, 
	!  with fluxcf(0) being the factor for the left domain interface
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: N, nelem,substep
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem), INTENT(IN) :: Acur,Apre
	REAL(KIND=DOUBLE), INTENT(IN) :: dxel,dt
	REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(IN) :: flx
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:N), INTENT(IN) :: DG_C
	! -- Outputs
	REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(OUT) :: fluxcf
	! -- Local variables
	REAL(KIND=DOUBLE) :: Pj,Qj,eps
	REAL(KIND=DOUBLE), DIMENSION(0:nelem+1) :: R ! Reduction ratio for outward fluxes so that element j has non-negative values (1D0 indicates no limiting needed)
	INTEGER :: j,k

	eps = 1D-16 ! Small parameter used to ensure no division by 0

	DO j=1,nelem
		! Compute maximum allowable flux out of element j based on which substep of ssprk3 we are on
		SELECT CASE(substep)
			CASE(1) 
				Qj = (dxel/dt)*Acur(0,j)
!		IF(Qj .lt. 0D0) THEN
!			write(*,*) 'Stage 1 Qj < 0! j=',j,'Qj=',Qj
!			write(*,*) Acur(0,j)
!		END IF

			CASE(2) 
				Qj = (dxel/dt)*(3D0*Apre(0,j) + 1D0*Acur(0,j))
!		IF(Qj .lt. 0D0) THEN
!			write(*,*) 'Stage 2 Qj < 0! j=',j,'Qj=',Qj
!			write(*,*) Apre(0,j),Acur(0,j)
!		END IF

			CASE(3) 
				Qj = (dxel/(2D0*dt))*(1D0*Apre(0,j) + 2D0*Acur(0,j))
!		IF(Qj .lt. 0D0) THEN
!			write(*,*) 'Stage 3 Qj < 0! j=',j,'Qj=',Qj
!			write(*,*) Apre(0,j),Acur(0,j)
!		END IF

		END SELECT

		! Compute actual flux out of element j
		Pj = DMAX1(0D0,flx(j)) - DMIN1(0D0,flx(j-1)) + eps

		! Compute reduction ratio
		R(j) = DMIN1(1D0,Qj/Pj)
	END DO
	! Periodicity
	R(0) = R(nelem)
	R(nelem+1) = R(1)

	! Compute flux corection factors
!	fluxcf = R(0:nelem)
	DO j=0,nelem
		! If flux at right edge is negative, use limiting ratio in element to the right of current one
		! (since that is where we are pulling mass from)
		fluxcf(j) = R(j) - (1D0-INT(SIGN(1D0,flx(j))))/2D0*(R(j)-R(j+1))
	END DO

END SUBROUTINE FLUXCOR

SUBROUTINE MFILL(rhoq,N,nelem)
	! Subroutine for mass filling within an element to remove negative cell averaged values
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)

	! -- Inputs
	INTEGER, INTENT(IN) :: N,nelem
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem), INTENT(INOUT) :: rhoq
	! -- Local Variables
	INTEGER :: j,k
	REAL(KIND=DOUBLE) :: r,Mp,Mt,s

	DO j=1,nelem
		Mp = 0D0
		Mt = 0D0

		DO k=0,N
			Mt = Mt + rhoq(k,j)
			rhoq(k,j) = MAX(0D0,rhoq(k,j)) ! Zero out negative masses
			Mp = Mp + rhoq(k,j)
		ENDDO

		r = MAX(Mt,0D0)/MAX(Mp,TINY(1D0))

!        IF(Mt .lt. 0D0) THEN
!            write(*,*) 'WARNING: negative total mass when redistributing!',Mt
!        END IF

!		IF(r .gt. 1D0) THEN
!			write(*,*) 'WARNING REDUCTION RATIO > 1.0!!'
!		ENDIF
		rhoq(:,j) = r*rhoq(:,j) ! Reduce remaining positive masses by reduction factor

	ENDDO

END SUBROUTINE MFILL

SUBROUTINE evalExpansion(A,leg_quad,quadVals,edgeVals,N,nelem)
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: nelem,N
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:N), INTENT(IN) :: leg_quad
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem), INTENT(IN) :: A

	! -- Outputs
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem), INTENT(OUT) :: quadVals
	REAL(KIND=DOUBLE), DIMENSION(0:1,0:nelem+1), INTENT(OUT) :: edgeVals
	! -- Local variables
	INTEGER :: i,j

	DO j=1,nelem
		DO i=0,N
			quadVals(i,j) = SUM(A(:,j)*leg_quad(:,i))
		ENDDO
		edgeVals(0,j) = SUM(A(:,j)*(/ ((-1D0)**i , i=0,N) /)) ! Expansion value at left edge of element
		edgeVals(1,j) = SUM(A(:,j)) ! Expansion value at right edge of element
	ENDDO

	! Extend edgeVals periodically
	edgeVals(:,0) = edgeVals(:,nelem)
	edgeVals(:,nelem+1) = edgeVals(:,1)

END SUBROUTINE evalExpansion

SUBROUTINE projectAverages(A,DG_LUC,IPIV,avgs,N,nelem)
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: nelem,N
	REAL(KIND=DOUBLE), DIMENSION(0:N,0:N), INTENT(IN) :: DG_LUC
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem), INTENT(IN) :: avgs
	INTEGER, DIMENSION(0:N), INTENT(IN) :: IPIV
	! -- Outputs
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: A
	! -- Local variables
	INTEGER :: i,j,k
	REAL(KIND=DOUBLE), DIMENSION(1:nelem) :: hold
	REAL(KIND=DOUBLE), DIMENSION(0:N,1:nelem) :: fooBAR
	REAL(KIND=DOUBLE), DIMENSION(0:N) :: FOO_y

	fooBAR = avgs

	DO i=0,N ! Reorder RHS according to IPIV
        hold = fooBAR(i,1:nelem)
        	fooBAR(i,:) = fooBAR(IPIV(i)-1,:)
        	fooBAR(IPIV(i)-1,1:nelem) = hold
	ENDDO

	DO j=1,nelem
		FOO_y = 0D0
		! Solve Ly=RHS for y
		FOO_y(0) = fooBAR(0,j)
		DO k=1,N
			FOO_y(k) = fooBAR(k,j) - SUM(DG_LUC(k,0:k-1)*FOO_y(0:k-1))
		ENDDO
		! Solve Ux=y for x
		A(N,j) = (1D0/DG_LUC(N,N))*FOO_y(N)
		DO k=N-1,0,-1
			A(k,j) = (1D0/DG_LUC(k,k))*(FOO_y(k) - SUM(DG_LUC(k,k+1:N)*A(k+1:N,j)))
		ENDDO
	ENDDO

END SUBROUTINE projectAverages
