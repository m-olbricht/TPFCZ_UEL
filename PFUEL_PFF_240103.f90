!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ABAQUS user defined phase field element subroutine
!
! Stephan Roth, TU Bergakademie Freiberg, 30.07.2020
!
! 30.07.2020: Multi-phase multi-component
! 15.06.2021: concentrations independent on phases
! 01.11.2023: without NCP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE UEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                   mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                   ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                   pnewdt,jprops,njprop,period)

      USE ABQINTERFACE
      USE ABAModul

      IMPLICIT NONE

      ! UEL variables
      INTEGER(kind=AbqIK), INTENT(IN) :: nprops, njprop, mcrd, nnode, ndofel, &
                                         npredf, mlvarx, nrhs, nsvars, jtype, &
                                         kstep, kinc, jelem, ndload, mdload, &
                                         lflags(*), jdltyp(mdload,*), &
                                         jprops(njprop)
      REAL(kind=AbqRK), INTENT(IN) :: props(nprops), &
                                      coords(mcrd,nnode), u(ndofel), period, &
                                      v(ndofel), a(ndofel), time(2), dtime, &
                                      predef(2,npredf,nnode), du(mlvarx,*), &
                                      adlmag(mdload,*), ddlmag(mdload,*), params(*)
      REAL(kind=AbqRK), INTENT(INOUT) :: energy(8), pnewdt, svars(nsvars)
      REAL(kind=AbqRK), INTENT(OUT) :: amatrx(ndofel,ndofel), rhs(mlvarx,*)

      ! UMAT variables
      INTEGER(kind=AbqIK) :: ntens, layer, kspt, npt
      REAL(kind=AbqRK) :: predef_umat(npredf), dpred(npredf), rpl, drpldt, celent, F0(3,3), F1(3,3), GPcoords(3), drot(3,3)
      REAL(kind=AbqRK), ALLOCATABLE :: ddsddt(:), drplde(:)
      CHARACTER*45 :: cmname

      ! further variables
      REAL(kind=AbqRK) :: fint(ndofel), energy_gp(8), JacobiDet, prop_thickness, p, lambda
      REAL(kind=AbqRK), ALLOCATABLE :: Matrix_B(:,:), stran(:), dstran(:), Ct(:,:), stress(:)
      INTEGER(kind=AbqIK) :: i1, i2, i3
      INTEGER(kind=AbqIK) :: nNodalDOF, nPp, iT, iT_I
      LOGICAL :: isNlgeom, positivDetJ, k_RHS, k_K, k_M, positivDetJTemp
      

      ! initialisation
      positivDetJ=.TRUE.; isNlgeom=.FALSE.; k_RHS=.TRUE.; k_K=.TRUE.; k_M=.FALSE.; positivDetJTemp=.TRUE.
      ! umat dummies
      layer = 1; kspt = 1
      celent = one; F0 = zero; F1 = zero
      !
      
      IF (kinc .EQ. 1 .AND. jelem .EQ. 1) THEN
		  WRITE(*,*) ''
		  WRITE(*,*) 'u: ', u(:)
		  WRITE(*,*) ''
      END IF

      ! 'Wer bin ich, was kann ich ?'

      ! evaluation of lflags
      ! geometric nonlinear
      IF (lflags(2) .EQ. 1) isNlgeom = .TRUE.

      SELECT CASE(lflags(3))
      CASE(1)
        k_RHS = .TRUE.
        k_K   = .TRUE.
        k_M   = .FALSE.
      CASE(2)
        k_RHS = .FALSE.
        k_K   = .TRUE.
        k_M   = .FALSE.
      CASE(3)
        k_RHS = .FALSE.
        k_K   = .FALSE.
        k_M   = .FALSE.
      CASE(4)
        k_RHS = .FALSE.
        k_K   = .FALSE.
        k_M   = .TRUE.
      CASE(5)
        k_RHS = .TRUE.
        k_K   = .FALSE.
        k_M   = .FALSE.
      CASE(6)
        k_RHS = .FALSE.
        k_K   = .FALSE.
        k_M   = .TRUE.
      CASE DEFAULT
        k_RHS = .FALSE.
        k_K   = .FALSE.
        k_M   = .FALSE.
        write(7,*) 'no output required'
      END SELECT

      ! number of nodal DOF
      !nNodalDOF = ndofel/nnode
      nPp       = 1 ! number of order parameters, here: damage variable 

      ! number of tensor coordinates (generalised): strain, damage variable
      ntens=NDI+NSHR+nPp*(1+D) + (1+D)
	  
      ! allocation of UMAT-matrices
      ALLOCATE(ddsddt(ntens), drplde(ntens), Matrix_B(ntens,ndofel), Ct(ntens,ntens), stran(ntens), dstran(ntens), stress(ntens))
      ddsddt=zero; drplde=zero; Matrix_B = zero; stran = zero; dstran = zero; Ct = zero; stress = zero

      ! query at time 0: number of internal state variables per integration point and material parameters: answer from UMAT
      IF ( time(2) .EQ. zero ) THEN
        IF (nsvars .LT. NGP*numSDV) THEN
          WRITE(7,*) "A bigger number of solution-dependent state variables is required.", NGP*numSDV, nsvars
          WRITE(7,*) "Please note that this is a total number of SDVs, not per integration point."
          CALL XEXIT()
        END IF
        IF (nprops .LT. numMatPar) THEN
          WRITE(7,*) "A bigger number of real parameters is required.", numMatPar, nprops
          CALL XEXIT()
        END IF
        ! check of material parameters: answer from UMAT
        CALL CheckMaterialParameters(props)
      END IF


      ! thickness of 2D-Elements
      prop_thickness = one
      IF (D .EQ. 2) prop_thickness = props(thicknessIndex)

      ! Ende 'Wer bin ich, was kann ich ?'


      ! UMAT call only when neccessary
      IF (k_RHS .OR. k_K) THEN
        ! loop over integration points
        energy = zero
        amatrx = zero

        fint   = zero
        DO npt=1,NGP
          ! compute JACOBI determinante and B-matrix

          CALL BMatrixJacobian(coords(1:D,1:nnode),u,D,nnode,ndofel,NDI,NSHR,ntens,njprop,ShapeFunc(GPPos(npt,1:D)), &
                               ShapeFuncDeriv(GPPos(npt,1:D)),jprops,prop_thickness,isNlgeom,Matrix_B,JacobiDet,drot,GPcoords,positivDetJ)
		  

          IF (positivDetJ) THEN
          ! element OK
            ! compute local (transformed) separation
            stran  = matmul(Matrix_B, u)
            dstran = matmul(Matrix_B, du(1:ndofel,1))
            ! here extract stress from svars and compute F0,F1...
            stress = zero
            energy_gp = zero
            ! evaluation constitutive law
            
            CALL umat(stress,svars(numSDV*(npt-1)+1:numSDV*npt),Ct,energy_gp(2), &
                      energy_gp(4),energy_gp(3),rpl,ddsddt,drplde,drpldt,stran, &
                      dstran,time,dtime,predef_umat(1),dpred(1),predef_umat,dpred,cmname,NDI,NSHR,ntens, &
                      numSDV,props,nprops,GPcoords,drot,pnewdt,celent,F0,F1,jelem,npt,layer,kspt,kstep,kinc)

!            ! viscous dissipation = Helmholtz free energy bulk mix (SDV 31)
!            energy_gp(5) = svars(numSDV*(npt-1)+31)
!            ! electrostatic energy = gradient part of chemical energy (mix) (SDV 34)
!            energy_gp(7) = svars(numSDV*(npt-1)+34)
!            ! viscous dissipation = order parameter 1 (SDV 21)
!            energy_gp(5) = svars(numSDV*(npt-1)+21)
!            ! electrostatic energy = order parameter 2 (SDV 22)
!            energy_gp(7) = svars(numSDV*(npt-1)+22)
!            ! creep dissipation = order parameter 3 (SDV 23)
!            energy_gp(3) = svars(numSDV*(npt-1)+23)
!            ! electrostatic energy = sum of order parameters 2 (SDV 21 -- SDV 20+nOp)
!            energy_gp(7) = SUM(svars(numSDV*(npt-1)+21:numSDV*(npt-1)+20+nOp))
			
			
            ! integration over element
            ! energies
            energy(:) = energy(:) + GPWeight(npt)*JacobiDet * energy_gp(:)
            ! stiffness matrix
            amatrx = amatrx + GPWeight(npt)*JacobiDet * matmul(transpose(Matrix_B),matmul(Ct,Matrix_B))
            ! internal force vector
            fint = fint - GPWeight(npt)*JacobiDet * matmul(transpose(Matrix_B),stress)
		  ELSE
          ! distorted element
            ! stop loop
            EXIT
          END IF
        END DO
        
        
!~         WRITE(*,*) 'Perturbation UMAT call: kstep=', kstep, ' dtime=', dtime
!~ 		CALL PerturbationTest( &
!~ 			coords, u, du, D, nnode, ndofel, NDI, NSHR, ntens, njprop,  &
!~ 			jprops, prop_thickness, isNlgeom,                            &
!~ 			svars, numSDV, NGP, GPPos, GPWeight,                         &
!~ 			predef_umat, dpred, npredf, cmname,                          &
!~ 			props, nprops, pnewdt, celent,                               &
!~ 			F0, F1, jelem, layer, kspt, kstep, kinc,                     &
!~ 			time, dtime, amatrx)
        

        ! zero matrices for distorted elements
        IF (.NOT. positivDetJ) THEN
          WRITE(7,*) 'non-positive JACOBIAN det(J)<=0, distorted Element: ', jelem, ' matrices set to zero'
          amatrx = zero; fint = zero; energy = zero
        END IF

        rhs(1:ndofel,1) = fint

      END IF

      IF (k_M) THEN
      ! mass matrix
        amatrx = zero
      ELSE IF (.NOT. k_K) THEN
      ! no matrix required
        amatrx = zero
      END IF

      IF (.NOT. k_RHS) THEN
      ! no rhs required
        rhs(1:ndofel,1) = zero
      END IF
      
      IF (kinc .EQ. 1 .AND. jelem .EQ. 1) THEN
		  WRITE(*,*) ''
		  WRITE(*,'(A)') 'RHS: '
		  DO i1 = 1, ndofel
			WRITE(*,'(2X,I3,2X,ES14.6)') i1, rhs(i1,1)
	      END DO
		  !
		  WRITE(*,'(A)') 'amatrx:'
		  DO i1 = 1, ndofel
		  	WRITE(*,'(100(ES14.6,1X))') (amatrx(i1,i2), i2 = 1, ndofel)
		  END DO
		  WRITE(*,*) ''
      END IF
      
      DEALLOCATE(Matrix_B, ddsddt, drplde, stran, dstran, Ct, stress)

      RETURN

!      CONTAINS

!------------------------------------------------------------------------------------

    END SUBROUTINE UEL

!------------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!~ !------------------------------------------------------------------------------------!

!~ ! === SCHRITT 2: d am Gausspunkt aus U interpolieren ===========
!~ ! d_gp = Summe N_I(xgp) * U(d-DOF von Knoten I)
!~ ! d-DOF-Indices für 4 Knoten (u1 u2 d T): 3, 7, 11, 15

!~ subroutine get_d_gp(U, N, d_gp)
!~   implicit none
!~   double precision, intent(in)  :: U(16)   ! alle DOF des Elements
!~   double precision, intent(in)  :: N(4)    ! Formfunktionswerte am GP
!~   double precision, intent(out) :: d_gp

!~   d_gp = N(1)*U(3) + N(2)*U(7) + N(3)*U(11) + N(4)*U(15)
!~   d_gp = max(0.d0, min(1.d0, d_gp))        ! Klemmen [0,1]
!~ end subroutine
!~ ! =============================================================

!~ !------------------------------------------------------------------------------------!

!~ ! === SCHRITT 4: Phasenfeld-Residuum und Jacobi =================
!~ ! Lokaler Anteil (ohne Gradiententerm, der kommt aus ∫∇N·∇N dV):
!~ !
!~ !   R_d = Gc/l * d  -  2*(1-d)*H        (AT2-Modell)
!~ !   K_dd = Gc/l  +  2*H                 (immer > 0)
!~ !
!~ ! Gradientenanteil separat: Gc*l * B_d^T * B_d  (in deine B-Matrix)

!~ subroutine pff_local(d_gp, H_new, Gc, l_pf, R_d, K_dd)
!~   implicit none
!~   double precision, intent(in)  :: d_gp, H_new, Gc, l_pf
!~   double precision, intent(out) :: R_d, K_dd

!~   double precision :: gc_over_l
!~   gc_over_l = Gc / l_pf

!~   R_d  = gc_over_l * d_gp  -  2.d0*(1.d0 - d_gp)*H_new
!~   K_dd = gc_over_l          +  2.d0*H_new

!~   ! Sanity check
!~   if (.not.(R_d  > -1.d300 .and. R_d  < 1.d300)) then
!~     write(6,'(a,2es14.4)') '[INF] R_d! d_gp, H_new =', d_gp, H_new
!~     call xit
!~   end if
!~   if (K_dd <= 0.d0) then
!~     write(6,'(a,es14.4)') '[ERR] K_dd <= 0 ! K_dd =', K_dd
!~     call xit
!~   end if
!~ end subroutine
!~ ! =============================================================

!------------------------------------------------------------------------------------!

SUBROUTINE PerturbationTest( &
    coords, u, du, D, nnode, ndofel, NDI, NSHR, ntens, njprop,  &
    jprops, prop_thickness, isNlgeom,                            &
    svars, numSDV, NGP, GPPos, GPWeight,                         &
    predef_umat, dpred, npredf, cmname,                          &
    props, nprops, pnewdt, celent,                               &
    F0, F1, jelem, layer, kspt, kstep, kinc,                     &
    time, dtime, amatrx)

  USE ABQINTERFACE
  USE ABAModul
  IMPLICIT NONE

  !--- Eingaben ---
  INTEGER, INTENT(IN) :: D, nnode, ndofel, NDI, NSHR, ntens
  INTEGER, INTENT(IN) :: njprop, nprops, numSDV, NGP, npredf
  INTEGER, INTENT(IN) :: jelem, layer, kspt, kstep, kinc
  INTEGER, INTENT(IN) :: jprops(njprop)
  LOGICAL, INTENT(IN) :: isNlgeom
  CHARACTER*45, INTENT(IN) :: cmname
  REAL(8), INTENT(IN) :: coords(D, nnode)
  REAL(8), INTENT(IN) :: u(ndofel)
  REAL(8), INTENT(IN) :: du(ndofel, *)
  REAL(8), INTENT(IN) :: props(nprops)
  REAL(8), INTENT(IN) :: prop_thickness
  REAL(8), INTENT(IN) :: time(2), dtime, celent, pnewdt
  REAL(8), INTENT(IN) :: predef_umat(npredf), dpred(npredf)
  REAL(8), INTENT(IN) :: F0(3,3), F1(3,3)
  REAL(8), INTENT(IN) :: svars(*)
  REAL(8), INTENT(IN) :: GPPos(NGP, D)
  REAL(8), INTENT(IN) :: GPWeight(NGP)
  REAL(8), INTENT(IN) :: amatrx(ndofel, ndofel)

  !--- Lokale Variablen ---
  REAL(8), ALLOCATABLE :: fint_plus(:), fint_minus(:), U_pert(:)
  REAL(8), ALLOCATABLE :: K_num(:,:)
  REAL(8), ALLOCATABLE :: Matrix_B_pert(:,:)
  REAL(8), ALLOCATABLE :: stress_pert(:), Ct_pert(:,:)
  REAL(8), ALLOCATABLE :: stran_pert(:), dstran_pert(:)
  REAL(8) :: JacobiDet_pert, delta, err_rel, K_num_corr
  REAL(8) :: drot_pert(3,3), GPcoords_pert(3)
  REAL(8) :: energy_gp_pert(8)
  REAL(8) :: rpl_pert, drpldt_pert
  REAL(8) :: ddsddt_pert(ntens), drplde_pert(ntens)
  REAL(8) :: svars_pert(numSDV*NGP)
  LOGICAL :: positivDetJ_pert
  LOGICAL :: already_run
  INTEGER :: npt_pert, i, j, sign
  
  CHARACTER*45 :: cmname_pert
  REAL(8)      :: pnewdt_pert

  !--- Nur für definierten Zustand ---
  IF (kstep /= 2 .OR. kinc /= 20 .OR. jelem /= 1) RETURN
  IF (already_run) RETURN
  already_run = .TRUE.


  !--- Allokierung ---
  ALLOCATE(fint_plus(ndofel), fint_minus(ndofel), U_pert(ndofel))
  ALLOCATE(K_num(ndofel, ndofel))
  ALLOCATE(Matrix_B_pert(ntens, ndofel))
  ALLOCATE(stress_pert(ntens), Ct_pert(ntens, ntens))
  ALLOCATE(stran_pert(ntens), dstran_pert(ntens))
  
  cmname_pert  = cmname
  pnewdt_pert  = pnewdt


  delta  = 1.0D-6
  K_num  = 0.0D0

  WRITE(*,*) '====== PERTURBATIONSTEST: jelem=', jelem, ' kinc=', kinc, ' ======'

  ! ================================================================
  ! Schleife über alle DOFs: +delta und -delta
  ! ================================================================
  DO j = 1, ndofel

    fint_plus  = 0.0D0
    fint_minus = 0.0D0

    DO sign = -1, 1, 2   ! sign=-1 → minus, sign=+1 → plus

      U_pert    = u
      U_pert(j) = u(j) + sign * delta

      ! SVARs kopieren (kein Update)
      svars_pert = svars(1:numSDV*NGP)

      DO npt_pert = 1, NGP

        Matrix_B_pert  = 0.0D0
        stress_pert    = 0.0D0
        Ct_pert        = 0.0D0
        energy_gp_pert = 0.0D0
        rpl_pert       = 0.0D0
        drpldt_pert    = 0.0D0
        ddsddt_pert    = 0.0D0
        drplde_pert    = 0.0D0

        ! B-Matrix mit gestörtem U
        CALL BMatrixJacobian( &
          coords(1:D,1:nnode), U_pert, D, nnode, ndofel, NDI, NSHR, ntens, njprop, &
          ShapeFunc(GPPos(npt_pert,1:D)),                                            &
          ShapeFuncDeriv(GPPos(npt_pert,1:D)),                                       &
          jprops, prop_thickness, isNlgeom,                                          &
          Matrix_B_pert, JacobiDet_pert, drot_pert, GPcoords_pert, positivDetJ_pert)

        IF (.NOT. positivDetJ_pert) EXIT

        stran_pert  = matmul(Matrix_B_pert, U_pert)
        dstran_pert = matmul(Matrix_B_pert, du(1:ndofel,1))

		CALL umat(stress_pert,                                                     &
				  svars_pert(numSDV*(npt_pert-1)+1:numSDV*npt_pert),               &
				  Ct_pert, energy_gp_pert(2), energy_gp_pert(4),                   &
				  energy_gp_pert(3), rpl_pert, ddsddt_pert, drplde_pert,           &
				  drpldt_pert, stran_pert, dstran_pert,                             &
				  time, dtime, predef_umat(1), dpred(1),                            &
				  predef_umat, dpred, cmname_pert, NDI, NSHR, ntens,               &
				  numSDV, props, nprops, GPcoords_pert, drot_pert,                  &
				  pnewdt_pert, celent, F0, F1, jelem, npt_pert,                    &
				  layer, kspt, kstep, kinc)
		WRITE(*,*) 'nach UMAT: pnewdt_pert=', pnewdt_pert, ' j=', j, ' sign=', sign, ' npt_pert=', npt_pert
		IF (pnewdt_pert < 1.0D0) THEN
		  WRITE(*,*) 'WARNUNG: pnewdt_pert < 1 !'
		END IF
		!
        ! fint assemblieren
        IF (sign == 1) THEN
          fint_plus  = fint_plus  - GPWeight(npt_pert)*JacobiDet_pert &
                                  * matmul(transpose(Matrix_B_pert), stress_pert)
        ELSE
          fint_minus = fint_minus - GPWeight(npt_pert)*JacobiDet_pert &
                                  * matmul(transpose(Matrix_B_pert), stress_pert)
        END IF

      END DO  ! npt_pert

    END DO  ! sign

    ! Numerische Tangentenspalte j
    K_num(:,j) = (fint_plus - fint_minus) / (2.0D0 * delta)

  END DO  ! j
  
  
  ! Nach dem UMAT-Aufruf - pnewdt zurücksetzen:
  pnewdt_pert = 1.0D0

  ! ================================================================
  ! Ausgabe
  ! ================================================================
  WRITE(*,*) 'j-Loop fertig, starte Ausgabe'
  FLUSH(6)

  WRITE(*,*) '====== PERTURBATIONSTEST ERGEBNIS ======'
  WRITE(*,'(A6,A6,A14,A14,A14)') 'i', 'j', 'K_ana', 'K_num_corr', 'rel_err'
  FLUSH(6)

  DO j = 1, ndofel
    DO i = 1, ndofel
      K_num_corr = -K_num(i,j)
      IF (ABS(amatrx(i,j)) > 1.0D-6) THEN
        err_rel = ABS(amatrx(i,j) - K_num_corr) / ABS(amatrx(i,j))
      ELSE IF (ABS(K_num_corr) > 1.0D-6) THEN
        err_rel = 1.0D0
      ELSE
        CYCLE
      END IF
      IF (err_rel > 1.0D-3) THEN
        WRITE(*,'(I6,I6,3ES14.5,A)') i, j, amatrx(i,j), K_num_corr, err_rel, ' <-- FEHLER'
        FLUSH(6)
      END IF
    END DO
  END DO

  WRITE(*,*) '====== ENDE PERTURBATIONSTEST ======'
  FLUSH(6)


  DEALLOCATE(fint_plus, fint_minus, U_pert, K_num)
  DEALLOCATE(Matrix_B_pert, stress_pert, Ct_pert, stran_pert, dstran_pert)

END SUBROUTINE PerturbationTest
