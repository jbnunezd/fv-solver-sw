!===============================================================================!
MODULE MOD_Mesh
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE BuildMesh
  MODULE PROCEDURE BuildMesh
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: BuildMesh
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
SUBROUTINE BuildMesh()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: MESH_SX
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY
USE MOD_FiniteVolume2D_vars,ONLY: TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: MeshNodes
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP
!-------------------------------------------------------------------------------!

MeshNodes = 0.0
MeshBary  = 0.0

Mesh_SX    = ABS(Mesh_X1-Mesh_X0)
Mesh_DX(1) = ABS(Mesh_SX(1))/(REAL(nElemsX))
Mesh_DX(2) = ABS(Mesh_SX(2))/(REAL(nElemsY))

DO jj=0,nElemsY
  DO ii=0,nElemsX
    MeshNodes(1:nDims,ii,jj) = Mesh_X0(1:2) + (/REAL(ii),REAL(jj)/)*Mesh_DX(1:2)
  END DO
END DO

DO jj=1,nElemsY
  DO ii=1,nElemsX
    MeshBary(1:nDims,ii,jj) = MeshNodes(1:nDims,ii-1,jj-1) + 0.5*Mesh_DX(1:2)
  END DO
END DO

DO iGP=1,nGPs
  !------------------------------!
  ! Normal vectors: x-direction  !
  !------------------------------!
  DO jj=1,nElemsY
    DO ii=0,nElemsX
      NormVectX(1:nDims,iGP,ii,jj) = (/1.0,0.0/)
    END DO
  END DO

  !------------------------------!
  ! Normal vectors: y-direction  !
  !------------------------------!
  DO jj=0,nElemsY
    DO ii=1,nElemsX
      NormVectY(1:nDims,iGP,ii,jj) = (/0.0,1.0/)
    END DO
  END DO

  !------------------------------!
  ! Tangent vectors: x-direction !
  !------------------------------!
  DO jj=1,nElemsY
    DO ii=0,nElemsX
      TangVectX(1:nDims,iGP,ii,jj) = (/0.0,1.0/)
    END DO
  END DO

  !------------------------------!
  ! Tangent vectors: y-direction !
  !------------------------------!
  DO jj=0,nElemsY
    DO ii=1,nElemsX
      TangVectY(1:nDims,iGP,ii,jj) = (/-1.0,0.0/)
    END DO
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE BuildMesh
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Mesh
!-------------------------------------------------------------------------------!
