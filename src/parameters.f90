!===============================================================================!
MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE InitializeParameters
  MODULE PROCEDURE InitializeParameters
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: InitializeParameters
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
SUBROUTINE InitializeParameters()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: PI
USE MOD_FiniteVolume2D_vars,ONLY: CFL
USE MOD_FiniteVolume2D_vars,ONLY: TEnd
USE MOD_FiniteVolume2D_vars,ONLY: Gravity
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
USE MOD_FiniteVolume2D_vars,ONLY: WhichOutput
USE MOD_FiniteVolume2D_vars,ONLY: nOutputFiles
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: VarNameVisu
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

InitialCondition = 211

SELECT CASE(InitialCondition)
  CASE(200) ! Constant State
    TEnd    = 1.0
    Gravity = 9.8
    nElemsX = 100
    nElemsY = 100
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/)
  CASE(211) ! Dam Break
    TEnd    = 5.0
    Gravity = 9.8
    nElemsX = 100
    nElemsY = 100
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/40.0,40.0/)
    BoundaryConditionsType = (/2,2,2,2/)
  CASE(212) ! Circular Dam Break
    TEnd    = 5.0
    Gravity = 9.8
    nElemsX = 100
    nElemsY = 100
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/40.0,40.0/)
    BoundaryConditionsType = (/2,2,2,2/)
  CASE DEFAULT
    ErrorMessage = "Initial condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

CFL      = 0.95

Reconstruction    = 4
ReconstructionFix = 3

WhichOutput  = 2
nOutputFiles = 100

VarNameVisu(1) = "Depth"
VarNameVisu(2) = "VelocityX"
VarNameVisu(3) = "VelocityY"

!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeParameters
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
