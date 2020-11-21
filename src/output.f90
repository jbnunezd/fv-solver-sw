!===============================================================================!
MODULE MOD_Output
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE WriteMeshToDisk
  MODULE PROCEDURE WriteMeshToDisk
END INTERFACE

INTERFACE WriteSolutionToDisk
  MODULE PROCEDURE WriteSolutionToDisk
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: WriteMeshToDisk
PUBLIC :: WriteSolutionToDisk
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
SUBROUTINE WriteSolutionToDisk(OutputTime)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: VarNameVisu
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: WhichOutput
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: OutputTime
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
REAL,ALLOCATABLE   :: U_NVisu(:,:,:)
CHARACTER(LEN=255) :: FileName
!-------------------------------------------------------------------------------!

ALLOCATE(U_NVisu(1:nVar,1:nElemsX,1:nElemsY))

FileName = TIMESTAMP("SOLUTION",OutputTime)
FileName = TRIM(FileName)//".dat"

DO jj=1,nElemsY
  DO ii=1,nElemsX
    U_NVisu(1:nVar,ii,jj) = V(1:nVar,ii,jj)
  END DO
END DO

IF (WhichOutput .EQ. 1) THEN
  CALL WriteDataToDisk_OCTAVE(&
        nVar,                &
        (/nElemsX,nElemsY/), &
        'XY',                &
        MeshBary,            &
        VarNameVisu,         &
        U_NVisu,             &
        OutputTime,          &
        TRIM(FileName),      &
        "FV2D")
ELSEIF (WhichOutput .EQ. 2) THEN
  CALL WriteDataToDisk_TECPLOT(&
        nVar,                &
        (/nElemsX,nElemsY/), &
        'XY',                &
        MeshBary,            &
        VarNameVisu,         &
        U_NVisu,             &
        OutputTime,          &
        TRIM(FileName),      &
        "FV2D")
END IF

!-------------------------------------------------------------------------------!
END SUBROUTINE WriteSolutionToDisk
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WriteMeshToDisk()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: UNIT_FILE
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: FileName
!-------------------------------------------------------------------------------!

Filename = TRIM("MESH_CoordinateX.dat")
WRITE(*,'(A22,1X,A)') &
  "Writing DATA to Disk: ", "OCTAVE ASCII"//" => "//TRIM(FileName)
OPEN(UNIT_FILE, FILE=Filename)
WRITE(UNIT_FILE,'(A)') "CoordinateX"
DO ii=1,nElemsX
  WRITE(UNIT_FILE,'(SP,1(ES21.14E2))') MeshBary(1,ii,1)
END DO
CLOSE(UNIT_FILE)

Filename = TRIM("MESH_CoordinateY.dat")
WRITE(*,'(A22,1X,A)') &
  "Writing DATA to Disk: ", "OCTAVE ASCII"//" => "//TRIM(FileName)
OPEN(UNIT_FILE, FILE=Filename)
WRITE(UNIT_FILE,'(A)') "CoordinateY"
DO jj=1,nElemsY
  WRITE(UNIT_FILE,'(SP,1(ES21.14E2))') MeshBary(2,1,jj)
END DO
CLOSE(UNIT_FILE)

!-------------------------------------------------------------------------------!
END SUBROUTINE WriteMeshToDisk
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WriteDataToDisk_TECPLOT(nVar,nElems,CoordNames,Coords,VarNames,&
  OutputData,OutputTime,FileName,ZoneTitel)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: UNIT_FILE
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN)          :: nVar
INTEGER,INTENT(IN)          :: nElems(2)
REAL,INTENT(IN)             :: OutputTime
REAL,INTENT(IN)             :: Coords(2,1:nElems(1),1:nElems(2))
REAL,INTENT(IN)             :: OutputData(1:nVar,1:nElems(1),1:nElems(2))
CHARACTER(LEN=2),INTENT(IN) :: CoordNames
CHARACTER(LEN=*),INTENT(IN) :: VarNames(1:nVar)
CHARACTER(LEN=*),INTENT(IN) :: FileName
CHARACTER(LEN=*),INTENT(IN) :: ZoneTitel
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER                     :: ii, jj, iVar, Offset
INTEGER                     :: Width, Pos
CHARACTER(LEN=3)            :: space1
CHARACTER(LEN=7)            :: space2
CHARACTER(LEN=35)           :: VarString
CHARACTER(LEN=255)          :: FormatString
CHARACTER(LEN=512)          :: FormatTitle
CHARACTER(LEN=255)          :: temp1, temp2, temp3
!-------------------------------------------------------------------------------!

WRITE(*,'(A22,1X,A)') &
  "Writing DATA to Disk: ", "TECPLOT ASCII"//" => "//TRIM(FileName)

Offset = 38
FormatTitle = ""
FormatTitle(1:37) = &
  'VARIABLES="Coordinate'//CoordNames(1:1)//'","Coordinate'//CoordNames(2:2)//'"'

DO iVar=1,nVar
  WRITE(VarString,'(A2,A,A1)') ',"', TRIM(VarNames(iVar)), '"'
  FormatTitle(Offset:Offset+LEN(TRIM(VarString))) = TRIM(VarString)
  Offset = Offset + LEN(TRIM(VarString))
END DO

WRITE(temp1,'(F15.8)') OutputTime
WRITE(temp2,'(I8)') nElems(1)
WRITE(temp3,'(I8)') nElems(2)
OPEN(UNIT_FILE, FILE = TRIM(FileName), STATUS = "REPLACE")
WRITE(UNIT_FILE,'(A)') FormatTitle(1:Offset-1)
WRITE(UNIT_FILE,'(A,A,A,A,A,A,A,A)') &
  'ZONE T="', TRIM(ZoneTitel), '",STRANDID=1,SOLUTIONTIME=', &
  TRIM(ADJUSTL(TRIM(temp1))), ', DATAPACKING=POINT, ZONETYPE=ORDERED, I=', &
  TRIM(ADJUSTL(TRIM(temp2))), ', J=', TRIM(ADJUSTL(TRIM(temp3)))

WRITE(FormatString,'(A4,I2,A15)') "(SP,", nVar+2, "(ES21.14E2,2X))"

DO jj=1,nElems(2)
  DO ii=1,nElems(1)
    WRITE(UNIT_FILE,FormatString) Coords(1:2,ii,jj), OutputData(1:nVar,ii,jj)
  END DO
END DO

CLOSE(UNIT_FILE)

!-------------------------------------------------------------------------------!
END SUBROUTINE WriteDataToDisk_TECPLOT
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WriteDataToDisk_OCTAVE(nVar,nElems,CoordNames,Coords,VarNames,&
  OutputData,OutputTime,FileName,ZoneTitel)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: UNIT_FILE
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN)          :: nVar
INTEGER,INTENT(IN)          :: nElems(2)
REAL,INTENT(IN)             :: OutputTime
REAL,INTENT(IN)             :: Coords(2,1:nElems(1),1:nElems(2))
REAL,INTENT(IN)             :: OutputData(1:nVar,1:nElems(1),1:nElems(2))
CHARACTER(LEN=2),INTENT(IN) :: CoordNames
CHARACTER(LEN=*),INTENT(IN) :: VarNames(1:nVar)
CHARACTER(LEN=*),INTENT(IN) :: FileName
CHARACTER(LEN=*),INTENT(IN) :: ZoneTitel
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER                     :: ii, jj, iVar, Offset
INTEGER                     :: Width, Pos
CHARACTER(LEN=3)            :: space1
CHARACTER(LEN=7)            :: space2
CHARACTER(LEN=35)           :: VarString
CHARACTER(LEN=255)          :: FormatString
CHARACTER(LEN=512)          :: FormatTitle
CHARACTER(LEN=255)          :: temp1, temp2, temp3
!-------------------------------------------------------------------------------!

WRITE(*,'(A22,1X,A)') &
  "Writing DATA to Disk: ", "OCTAVE ASCII"//" => "//TRIM(FileName)

Pos   = 1
Width = 23
space1 = "   "
FormatTitle = ""

DO iVar=1,nVar
  WRITE(VarString,'(A)') TRIM(VarNames(iVar))
  FormatTitle(Pos:Pos+Width-1) = ADJUSTL(TRIM(VarString))
  Pos = Pos + Width
END DO

OPEN(UNIT_FILE, FILE = TRIM(FileName), STATUS = "REPLACE")
WRITE(UNIT_FILE,'(A)') FormatTitle(1:Pos-1)
WRITE(FormatString,'(A4,I2,A15)') "(SP,", nVar, "(ES21.14E2,2X))"

DO jj=1,nElems(2)
  DO ii=1,nElems(1)
    WRITE(UNIT_FILE,FormatString) OutputData(1:nVar,ii,jj)
  END DO
END DO

CLOSE(UNIT_FILE)

!-------------------------------------------------------------------------------!
END SUBROUTINE WriteDataToDisk_OCTAVE
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION TIMESTAMP(Filename,Time)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,            INTENT(IN) :: Time
CHARACTER(LEN=*),INTENT(IN) :: Filename
CHARACTER(LEN=255)          :: TimeStamp
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER                     :: i
!-------------------------------------------------------------------------------!

WRITE(TimeStamp,'(F15.7)') Time
DO i=1,LEN(TRIM(TimeStamp))
  IF(TimeStamp(i:i) .EQ. " ") TimeStamp(i:i) = "0"
END DO
TimeStamp=TRIM(Filename)//'_'//TRIM(TimeStamp)

!-------------------------------------------------------------------------------!
END FUNCTION TIMESTAMP
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Output
!-------------------------------------------------------------------------------!
