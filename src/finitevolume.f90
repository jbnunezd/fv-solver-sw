!===============================================================================!
MODULE MOD_FiniteVolume2D
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE InitializeFiniteVolume
  MODULE PROCEDURE InitializeFiniteVolume
END INTERFACE

INTERFACE FillInitialConditions
  MODULE PROCEDURE FillInitialConditions
END INTERFACE

INTERFACE FVTimeDerivative
  MODULE PROCEDURE FVTimeDerivative
END INTERFACE

INTERFACE FinalizeFiniteVolume
  MODULE PROCEDURE FinalizeFiniteVolume
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: InitializeFiniteVolume
PUBLIC :: FillInitialConditions
PUBLIC :: FVTimeDerivative
PUBLIC :: FinalizeFiniteVolume
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
SUBROUTINE InitializeFiniteVolume()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MeshNodes
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY
USE MOD_FiniteVolume2D_vars,ONLY: TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FluxX
USE MOD_FiniteVolume2D_vars,ONLY: FluxY
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: K1
USE MOD_FiniteVolume2D_vars,ONLY: K2
USE MOD_FiniteVolume2D_vars,ONLY: K3
USE MOD_FiniteVolume2D_vars,ONLY: K4
USE MOD_FiniteVolume2D_vars,ONLY: K5
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
INTEGER :: ii
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

SELECT CASE (Reconstruction)
  CASE(1) ! NONE
    nGhosts = 1
    nGPs    = 1
  CASE(2) ! MUSCL
    nGhosts = 1
    nGPs    = 1
  CASE(3) ! WENO3
    nGhosts = 1
    nGPs    = 2
  CASE(4) ! WENO5
    nGhosts = 2
    nGPs    = 2
  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

ALLOCATE(MeshNodes(1:nDims,       0:nElemsX,0:nElemsY))
ALLOCATE(MeshBary (1:nDims,       1:nElemsX,1:nElemsY))
ALLOCATE(NormVectX(1:nDims,1:nGPs,0:nElemsX,1:nElemsY))
ALLOCATE(TangVectX(1:nDims,1:nGPs,0:nElemsX,1:nElemsY))
ALLOCATE(NormVectY(1:nDims,1:nGPs,1:nElemsX,0:nElemsY))
ALLOCATE(TangVectY(1:nDims,1:nGPs,1:nElemsX,0:nElemsY))

ALLOCATE( U(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
ALLOCATE( V(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
ALLOCATE(Ut(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(WM(1:nVar,1:nGPs,0:nElemsX+1,0:nElemsY+1))
ALLOCATE(WP(1:nVar,1:nGPs,0:nElemsX+1,0:nElemsY+1))
ALLOCATE( S(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(FX(1:nVar,0:nElemsX,1:nElemsY))
ALLOCATE(FY(1:nVar,1:nElemsX,0:nElemsY))
ALLOCATE(FluxX(1:nVar,1:nGPs,0:nElemsX,1:nElemsY))
ALLOCATE(FluxY(1:nVar,1:nGPs,1:nElemsX,0:nElemsY))
ALLOCATE(Ind(1:2,0:nElemsX+1,0:nElemsY+1))

ALLOCATE(K0(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K1(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K2(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K3(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K4(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K5(1:nVar,1:nElemsX,1:nElemsY))

U  = 0.0
V  = 0.0
S  = 0.0
Ut = 0.0
FX = 0.0
FY = 0.0
WM = 0.0
WP = 0.0
FluxX = 0.0
FluxY = 0.0

Ind = .FALSE.

K0 = 0.0
K1 = 0.0
K2 = 0.0
K3 = 0.0
K4 = 0.0
K5 = 0.0

!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeFiniteVolume
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE FillInitialConditions()
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

DO jj=1,nElemsY
  DO ii=1,nElemsX
    CALL ExactFunction(&
      InitialCondition,0.0,MeshBary(1:nDims,ii,jj),U(1:nVar,ii,jj))
    CALL ConsToPrim(U(1:nVar,ii,jj),V(1:nVar,ii,jj))
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE FillInitialConditions
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE FVTimeDerivative(t)
!-------------------------------------------------------------------------------!
USE MOD_Equation,       ONLY: BoundaryConditions
USE MOD_Equation,       ONLY: SourceTerms
USE MOD_Reconstruction, ONLY: ReconstructionX
USE MOD_Reconstruction, ONLY: ReconstructionY
USE MOD_Reconstruction, ONLY: ReconstructionFixX
USE MOD_Reconstruction, ONLY: ReconstructionFixY
USE MOD_ShocksIndicator,ONLY: ShocksIndicatorX
USE MOD_ShocksIndicator,ONLY: ShocksIndicatorY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

CALL BoundaryConditions(t)

CALL ShocksIndicatorX()
CALL ReconstructionX()
CALL ReconstructionFixX()
CALL NumericalFluxFX()

CALL ShocksIndicatorY()
CALL ReconstructionY()
CALL ReconstructionFixY()
CALL NumericalFluxFY()

CALL SourceTerms(t)
CALL UpdateTimeDerivative()

!-------------------------------------------------------------------------------!
END SUBROUTINE FVTimeDerivative
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE UpdateTimeDerivative()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: Ut, S
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

DO jj=1,nElemsY
  DO ii=1,nElemsX
    Ut(1:nVar,ii,jj) = S(1:nVar,ii,jj) &
                     - (FX(1:nVar,ii+0,jj+0)-FX(1:nVar,ii-1,jj+0))/Mesh_DX(1) &
                     - (FY(1:nVar,ii+0,jj+0)-FY(1:nVar,ii+0,jj-1))/Mesh_DX(2)
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE UpdateTimeDerivative
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE NumericalFluxFX()
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: RiemannSolver
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FluxX
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX, TangVectX
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP
!-------------------------------------------------------------------------------!

FX    = 0.0
FluxX = 0.0

DO jj=1,nElemsY
  DO ii=0,nElemsX
    CALL RiemannSolver(&
                       WP(1:nVar,1:nGPs,ii+0,jj),&
                       WM(1:nVar,1:nGPs,ii+1,jj),&
                       NormVectX(1:nDims,1:nGPs,ii,jj),&
                       TangVectX(1:nDims,1:nGPs,ii,jj),&
                       FluxX(1:nVar,1:nGPs,ii,jj))
  END DO
END DO

SELECT CASE(Reconstruction)
  CASE(1,2)
    FX(:,:,:) = FX(:,:,:) + FluxX(:,nGPs,:,:)

  CASE(3,4)
    DO iGP = 1,nGPs
      FX(:,:,:) = FX(:,:,:) + 0.5*FluxX(:,iGP,:,:)
    END DO

END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE NumericalFluxFX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE NumericalFluxFY()
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: RiemannSolver
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FluxY
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY, TangVectY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP
!-------------------------------------------------------------------------------!

FY    = 0.0
FluxY = 0.0

DO jj=0,nElemsY
  DO ii=1,nElemsX
    CALL RiemannSolver(&
                       WP(1:nVar,1:nGPs,ii,jj+0),&
                       WM(1:nVar,1:nGPs,ii,jj+1),&
                       NormVectY(1:nDims,1:nGPs,ii,jj),&
                       TangVectY(1:nDims,1:nGPs,ii,jj),&
                       FluxY(1:nVar,1:nGPs,ii,jj))
  END DO
END DO

SELECT CASE(Reconstruction)
  CASE(1,2)
    FY(:,:,:) = FY(:,:,:) + FluxY(:,nGPs,:,:)

  CASE(3,4)
    DO iGP=1,nGPs
      FY(:,:,:) = FY(:,:,:) + 0.5*FluxY(:,iGP,:,:)
    END DO

END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE NumericalFluxFY
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE FinalizeFiniteVolume()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: MeshNodes
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY
USE MOD_FiniteVolume2D_vars,ONLY: TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FluxX
USE MOD_FiniteVolume2D_vars,ONLY: FluxY
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: K1
USE MOD_FiniteVolume2D_vars,ONLY: K2
USE MOD_FiniteVolume2D_vars,ONLY: K3
USE MOD_FiniteVolume2D_vars,ONLY: K4
USE MOD_FiniteVolume2D_vars,ONLY: K5
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!

DEALLOCATE(MeshNodes)
DEALLOCATE(MeshBary)
DEALLOCATE(NormVectX)
DEALLOCATE(TangVectX)
DEALLOCATE(NormVectY)
DEALLOCATE(TangVectY)
DEALLOCATE(U)
DEALLOCATE(V)
DEALLOCATE(Ut)
DEALLOCATE(WM)
DEALLOCATE(WP)
DEALLOCATE(S)
DEALLOCATE(FX)
DEALLOCATE(FY)
DEALLOCATE(FluxX)
DEALLOCATE(FluxY)
DEALLOCATE(Ind)
DEALLOCATE(K0)
DEALLOCATE(K1)
DEALLOCATE(K2)
DEALLOCATE(K3)
DEALLOCATE(K4)
DEALLOCATE(K5)

!-------------------------------------------------------------------------------!
END SUBROUTINE FinalizeFiniteVolume
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_FiniteVolume2D
!-------------------------------------------------------------------------------!
