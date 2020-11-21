!===============================================================================!
MODULE MOD_TimeDiscretization
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE TimeDiscretization
  MODULE PROCEDURE TimeDiscretization
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: TimeDiscretization
!-------------------------------------------------------------------------------!
!
!
!
!===============================================================================!
CONTAINS
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretization()
!-------------------------------------------------------------------------------!
USE MOD_Output,             ONLY: WriteMeshToDisk
USE MOD_Output,             ONLY: WriteSolutionToDisk
USE MOD_Equation,           ONLY: TimeStep
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: tEnd
USE MOD_FiniteVolume2D_vars,ONLY: dt_Analyze
USE MOD_FiniteVolume2D_vars,ONLY: nOutputFiles
USE MOD_FiniteVolume2D_vars,ONLY: WhichOutput
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL :: t, tAnalyze, dt_min
!-------------------------------------------------------------------------------!

t = 0.0
dt_Analyze = tEnd/REAL(nOutputFiles)
tAnalyze   = t + dt_Analyze
IF (WhichOutput .EQ. 1) THEN
  CALL WriteMeshToDisk()
END IF
CALL WriteSolutionToDisk(t)
IF (t .EQ. tEnd) THEN
  RETURN
END IF

DO
  dt_min = TimeStep()
  dt = MIN(MIN(dt_min,tAnalyze-t),MIN(dt_min,tEnd-t))
  CALL TimeDiscretizationBySSPRK4(t)
  t = t + dt
  IF (ABS(tAnalyze-t) .LT. 1.0E-10) THEN
    CALL WriteSolutionToDisk(t)
    tAnalyze = tAnalyze + dt_Analyze
    IF (tAnalyze .GT. tEnd) THEN
      tAnalyze = tEnd
    END IF
  END IF
  IF (ABS(t-tEnd) .LT. 1.0E-10) THEN
    EXIT
  END IF
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretization
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TimeDiscretizationBySSPRK4(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D,     ONLY: FVTimeDerivative
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: dt
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: K0(1:nVar,1:nElemsX,1:nElemsY)
REAL            :: K1(1:nVar,1:nElemsX,1:nElemsY)
REAL            :: K2(1:nVar,1:nElemsX,1:nElemsY)
REAL            :: K3(1:nVar,1:nElemsX,1:nElemsY)
REAL            :: K4(1:nVar,1:nElemsX,1:nElemsY)
REAL            :: K5(1:nVar,1:nElemsX,1:nElemsY)
REAL            :: tStage
INTEGER         :: ii, jj
!-------------------------------------------------------------------------------!

K0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! First Stage        !
!--------------------!
tStage = t + 0.0*dt
CALL FVTimeDerivative(tStage)
K1(1:nVar,1:nElemsX,1:nElemsY) = &
    1.00000000000000*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.39175222700392*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K1(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Second Stage       !
!--------------------!
tStage = t + 0.39175222700392*dt
CALL FVTimeDerivative(tStage)
K2(1:nVar,1:nElemsX,1:nElemsY) = &
    0.44437049406734*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.55562950593266*K1(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.36841059262959*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K2(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Third Stage        !
!--------------------!
tStage = t + 0.58607968896780*dt
CALL FVTimeDerivative(tStage)
K3(1:nVar,1:nElemsX,1:nElemsY) = &
    0.62010185138540*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.37989814861460*K2(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.25189177424738*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K3(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Fourth Stage       !
!--------------------!
tStage = t + 0.474542364687*dt
CALL FVTimeDerivative(tStage)
K4(1:nVar,1:nElemsX,1:nElemsY) = &
    0.17807995410773*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.82192004589227*K3(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.54497475021237*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K4(1:nVar,1:nElemsX,1:nElemsY)
K5(1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

!--------------------!
! Fifth Stage        !
!--------------------!
tStage  = t + 0.93501063100924*dt
CALL FVTimeDerivative(tStage)
U(1:nVar,1:nElemsX,1:nElemsY)  = &
    0.00683325884039*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.51723167208978*K2(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.12759831133288*K3(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.34833675773694*K4(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.08460416338212*K5(1:nVar,1:nElemsX,1:nElemsY)*dt &
  + 0.22600748319395*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt

!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationBySSPRK4
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_TimeDiscretization
!-------------------------------------------------------------------------------!
