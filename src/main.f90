!===============================================================================!
PROGRAM FiniteVolume2D
!-------------------------------------------------------------------------------!
USE MOD_Mesh,              ONLY: BuildMesh
USE MOD_Parameters,        ONLY: InitializeParameters
USE MOD_FiniteVolume2D,    ONLY: FillInitialConditions
USE MOD_FiniteVolume2D,    ONLY: InitializeFiniteVolume
USE MOD_FiniteVolume2D,    ONLY: FinalizeFiniteVolume
USE MOD_TimeDiscretization,ONLY: TimeDiscretization
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!

CALL InitializeParameters()
CALL InitializeFiniteVolume()
CALL BuildMesh()
CALL FillInitialConditions()
CALL TimeDiscretization()
CALL FinalizeFiniteVolume()

!-------------------------------------------------------------------------------!
END PROGRAM FiniteVolume2D
!===============================================================================!
