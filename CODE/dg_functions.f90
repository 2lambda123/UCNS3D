MODULE DG_FUNCTIONS

USE BASIS

IMPLICIT NONE

CONTAINS

FUNCTION DG_SOL(N, I_ELEM, X_IN, Y_IN, NUM_VARIABLES, IORDER, IDEGFREE, U_C_VALDG)
!> @brief
!> This function returns the DG solution at a given point (X_IN, Y_IN)\n
!> REQUIRES: X_IN, Y_IN: coordinates of the point where the solution is requested, NUM_VARIABLES: number of solution variables, IDEGFREE: number of basis terms
    INTEGER,INTENT(IN)::N,I_ELEM,IORDER,IDEGFREE,NUM_VARIABLES
    REAL,INTENT(IN)::X_IN,Y_IN ! Coordinates of the point where the solution is requested
    REAL,DIMENSION(:,:),INTENT(IN)::U_C_VALDG
    REAL,DIMENSION(IDEGFREE)::BASIS_TEMP
    INTEGER::I_DOF, I_VAR
    REAL,DIMENSION(NUM_VARIABLES)::DG_SOL

    IF(ALL(SHAPE(U_C_VALDG) /= (/ NUM_VARIABLES,IDEGFREE+1 /))) THEN
        WRITE(400+N,*) 'DG_SOL: U_C_VALDG WRONG DIMENSIONS:', SHAPE(U_C_VALDG)
        STOP
    END IF
        
    BASIS_TEMP = BASIS_REC2D(N,X_IN,Y_IN,IORDER,I_ELEM,IDEGFREE)

    DO I_VAR = 1, NUM_VARIABLES
        DG_SOL = U_C_VALDG(I_VAR,1) + DOT_PRODUCT(BASIS_TEMP(:), U_C_VALDG(I_VAR,2:))
    END DO

END FUNCTION DG_SOL

FUNCTION DG_RHS_INTEGRAL(N,I_ELEM,QP_X,QP_Y,QP_WEIGHT,NUM_VARS,ORDER,NUM_DOFS,CELL_VOL_OR_SURF,FLUX_TERM,VOL_OR_SURF)
!> @brief
!> Calculates the volume or surface integral term in the DG RHS for scalar linear advection with speed = 1
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N,I_ELEM,ORDER,NUM_VARS,NUM_DOFS,VOL_OR_SURF
    REAL,INTENT(IN)::QP_X,QP_Y,QP_WEIGHT,CELL_VOL_OR_SURF
    REAL,DIMENSION(:),INTENT(IN)::FLUX_TERM
    INTEGER::I_VAR
    REAL,DIMENSION(NUM_DOFS+1,NUM_VARS)::DG_RHS_INTEGRAL
    
    IF (SIZE(FLUX_TERM) /= NUM_VARS) THEN
        WRITE(400+N,*) 'DG_RHS_INTEGRAL: FLUX_TERM WRONG DIMENSIONS:', SHAPE(FLUX_TERM)
        STOP
    END IF
    
    IF (VOL_OR_SURF == 1) THEN ! VOLUME INTEGRAL
        DO I_VAR = 1, NUM_VARS
            DG_RHS_INTEGRAL(1,I_VAR) = 0
            DG_RHS_INTEGRAL(2:,I_VAR) = FLUX_TERM(I_VAR) * QP_WEIGHT * CELL_VOL_OR_SURF * (BASIS_REC2D_DERIVATIVE(N,QP_X,QP_Y,ORDER,I_ELEM,NUM_DOFS,1) + BASIS_REC2D_DERIVATIVE(N,QP_X,QP_Y,ORDER,I_ELEM,NUM_DOFS,2)) ! For linear advection with speed = 1
        END DO
    ELSE IF (VOL_OR_SURF == 2) THEN ! SURFACE INTEGRAL
        DO I_VAR = 1, NUM_VARS
            DG_RHS_INTEGRAL(1,I_VAR) = FLUX_TERM(I_VAR) * QP_WEIGHT * CELL_VOL_OR_SURF
            DG_RHS_INTEGRAL(2:,I_VAR) = FLUX_TERM(I_VAR) * QP_WEIGHT * CELL_VOL_OR_SURF * BASIS_REC2D(N,QP_X,QP_Y,ORDER,I_ELEM,NUM_DOFS)
        END DO
    END IF

END FUNCTION DG_RHS_INTEGRAL

FUNCTION CALC_DELTA_XYZ(NUM_NODES, NUM_DIMS, NODES_IN)
!> @brief
!> Calculates the "delta x/y/z" as in Luo 2012 eq 3.12
    IMPLICIT NONE
    INTEGER,INTENT(IN)::NUM_NODES, NUM_DIMS
    REAL,DIMENSION(:,:),INTENT(IN)::NODES_IN ! (NODE, DIMENSION)
    INTEGER::I_NODES, I_DIM
    REAL::XYZ_MAX,XYZ_MIN
    REAL,DIMENSION(NUM_DIMS)::CALC_DELTA_XYZ
    
    DO I_DIM = 1, NUM_DIMS
        XYZ_MAX = -1E12
        XYZ_MIN = 1E12
        
        DO I_NODES = 1, NUM_NODES
            IF (NODES_IN(I_NODES,I_DIM) > XYZ_MAX) XYZ_MAX = NODES_IN(I_NODES,I_DIM)
            IF (NODES_IN(I_NODES,I_DIM) < XYZ_MIN) XYZ_MIN = NODES_IN(I_NODES,I_DIM)
        END DO
        
        CALC_DELTA_XYZ(I_DIM) = 0.5 * ABS(XYZ_MAX - XYZ_MIN)
    END DO
    
END FUNCTION CALC_DELTA_XYZ

SUBROUTINE PRESTORE_DG
!> @brief
!> Prestores IELEM(N,I)%DELTA_XYZ, QP_ARRAY
    IMPLICIT NONE
    INTEGER::I,K,I_QP,N_QP
    
    DO I=1,XMPIELRANK(N)
        DO K=1,IELEM(N,I)%NONODES
            NODES_LIST(k,1:2)=INODER(IELEM(N,I)%NODES(K))%CORD(1:2)
            vext(k,1:2)=NODES_LIST(k,1:2)
        END DO
        
        IELEM(N,I)%DELTA_XYZ = CALC_DELTA_XYZ(IELEM(N,I)%NONODES, DIMENSIONA, NODES_LIST)
    
        SELECT CASE(ielem(n,i)%ishape)
        CASE(5)
            CALL QUADRATUREQUAD(N,IGQRULES)
            N_QP=QP_quad
        CASE(6)
            CALL QUADRATURETRIANGLE(N,IGQRULES)
            N_QP=QP_Triangle
        END SELECT
                
        DO I_QP=1,N_QP
            QP_ARRAY(I,I_QP)%X = QPOINTS(1,I_QP)
            QP_ARRAY(I,I_QP)%Y = QPOINTS(2,I_QP)
            QP_ARRAY(I,I_QP)%QP_WEIGHT = WEQUA3D(I_QP)
        END DO
    END DO
    
    CALL ASS_MASS_MATRIX(N)
    
END SUBROUTINE PRESTORE_DG


SUBROUTINE ASS_MASS_MATRIX(N)
!> @brief
!> Assembles the mass matrix
!> REQUIRES: Globals: IELEM, QP_QUAD, QP_TRIANGLE, MASS_MATRIX
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::I_ELEM, I_QP, N_QP, I_DOF, J_DOF, KMAXE
    REAL,DIMENSION(IDEGFREE)::BASIS_VECTOR
    
    KMAXE = XMPIELRANK(N)
    
    ALLOCATE(MASS_MATRIX_CENTERS(N:N,KMAXE,NUM_DG_DOFS,NUM_DG_DOFS)); MASS_MATRIX_CENTERS(N,:,:,:) = ZERO; MASS_MATRIX_CENTERS(N,:,1,1) = 1.0D0
    
    DO I_ELEM = 1, KMAXE
        SELECT CASE(IELEM(N,I_ELEM)%ISHAPE)
        CASE(5) ! Quadrilateral
            N_QP = QP_QUAD
        CASE(6) ! Triangle
            N_QP = QP_TRIANGLE
        END SELECT
        
        DO I_QP = 1, N_QP
            BASIS_VECTOR = BASIS_REC2D(N,QP_ARRAY(I_ELEM,I_QP)%X,QP_ARRAY(I_ELEM,I_QP)%Y,IORDER,I_ELEM,IDEGFREE)
            
            DO I_DOF = 1, IDEGFREE
                MASS_MATRIX_CENTERS(N, I_ELEM, 1, I_DOF+1) = MASS_MATRIX_CENTERS(N, I_ELEM, 1, I_DOF+1) + BASIS_VECTOR(I_DOF) * QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT
                MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF+1, 1) = MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF+1, 1) + BASIS_VECTOR(I_DOF) * QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT
                
                DO J_DOF = 1, IDEGFREE
                    MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF+1, J_DOF+1) = MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF, J_DOF) + BASIS_VECTOR(I_DOF) * BASIS_VECTOR(J_DOF) * QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT
                END DO
            END DO
        END DO
        
        MASS_MATRIX_CENTERS(N, I_ELEM, :, :) = MASS_MATRIX_CENTERS(N, I_ELEM, :, :) * IELEM(N,I_ELEM)%TOTVOLUME
    END DO
    
    ALLOCATE(INV_MASS_MATRIX(N:N,KMAXE,NUM_DG_DOFS,NUM_DG_DOFS)); INV_MASS_MATRIX(N,:,:,:) = ZERO
    
    CALL COMPMASSINV(MASS_MATRIX_CENTERS(N,:,:,:), INV_MASS_MATRIX(N,:,:,:), NUM_DG_DOFS, KMAXE)
    
    DO I_ELEM = 1, KMAXE
        WRITE(300+N,*) I_ELEM
!         DO I_QP = 1, N_QP
!             WRITE(300+N,*) 'QP_ARRAY', QP_ARRAY(I_ELEM,I_QP)%X, QP_ARRAY(I_ELEM,I_QP)%Y, QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT
!             WRITE(300+N,*) 'BASIS,', BASIS_REC2D(N,QP_ARRAY(I_ELEM,I_QP)%X,QP_ARRAY(I_ELEM,I_QP)%Y,IORDER,I_ELEM,IDEGFREE)
!         END DO
        WRITE(300+N,*) 'XYZ', IELEM(N,I_ELEM)%NODES
!         WRITE(300+N,*) 'DELTAXYZ', IELEM(N,I_ELEM)%DELTA_XYZ
        WRITE(300+N,*) 'MMC', MASS_MATRIX_CENTERS(N,I_ELEM,:,:)
        WRITE(300+N,*) 'Inverse,', INV_MASS_MATRIX(N,I_ELEM,:,:)
        WRITE(300+N,*) 'Identity', MATMUL(MASS_MATRIX_CENTERS(N,I_ELEM,:,:),INV_MASS_MATRIX(N,I_ELEM,:,:))
    END DO
    
END SUBROUTINE ASS_MASS_MATRIX

SUBROUTINE COMPMASSINV(totalMM,invMM,N_DOFS,kmaxe)
!Calculate the inverse of the input matrix with Gauss-Jordan Elimination
IMPLICIT NONE
 
integer :: i,j,k,l,m,irow,P
real:: big,dum
real,DIMENSION(N_DOFS,N_DOFS)::a,b
integer,INTENT(IN)::N_DOFS,kmaxe
REAL,DIMENSION(:,:,:),INTENT(IN)::totalMM
REAL,DIMENSION(:,:,:),INTENT(INOUT)::invMM

DO P=1,kmaxe

a(:,:)=totalMM(P,:,:)
b(:,:)=zero

do i = 1,N_DOFS
    do j = 1,N_DOFS
        b(i,j) = 0.0
    end do
    b(i,i) = 1.0
end do 

do i = 1,N_DOFS   
   big = a(i,i)
   do j = i,N_DOFS
     if (a(j,i).gt.big) then
       big = a(j,i)
       irow = j
     end if
   end do
   ! interchange lines i with irow for both a() and b() matrices
   if (big.gt.a(i,i)) then
     do k = 1,N_DOFS
       dum = a(i,k)                      ! matrix a()
       a(i,k) = a(irow,k)
       a(irow,k) = dum
       dum = b(i,k)                 ! matrix b()
       b(i,k) = b(irow,k)
       b(irow,k) = dum
     end do
   end if
   ! divide all entries in line i from a(i,j) by the value a(i,i); 
   ! same operation for the identity matrix
   dum = a(i,i)
   do j = 1,N_DOFS
     a(i,j) = a(i,j)/dum
     b(i,j) = b(i,j)/dum
   end do
   ! make zero all entries in the column a(j,i); same operation for indent()
   do j = i+1,N_DOFS
     dum = a(j,i)
     do k = 1,N_DOFS
       a(j,k) = a(j,k) - dum*a(i,k)
       b(j,k) = b(j,k) - dum*b(i,k)               
            
     end do
   end do
end do
  
 do i = 1,N_DOFS-1
   do j = i+1,N_DOFS
     dum = a(i,j)
     do l = 1,N_DOFS
       a(i,l) = a(i,l)-dum*a(j,l)
       b(i,l) = b(i,l)-dum*b(j,l)
     end do
   end do
 end do
 
 invMM(P,:,:)=b(:,:)
  
END DO
 
END SUBROUTINE COMPMASSINV 



! FUNCTION LUO_LSQ_RECONSTRUCT(N, N_DIM)
! !> @brief
! !> This function reconstructs an approximation of with NUM_DG_RECONSTRUCT_DOFS
! !> REQUIRES: IELEM, U_C as globals
!     INTEGER,INTENT(IN)::N, N_DIM
!     INTEGER::I_ELEM, I_FACE, NEIGHBOR_INDEX, I_DIM
!     REAL,DIMENSION(:,:)::LHS_MATRIX ! See eq. 3.19 of Luo 2012
!     REAL,DIMENSION(:)::RHS_DG_RECONSTRUCT ! See eq. 3.19 of Luo 2012
!     REAL,DIMENSION(:)::LUO_LSQ_RECONSTRUCT ! NUM_DG_RECONSTRUCT_DOFS
!     
!     KMAXE = XMPIELRANK(N)
!     
!     DO I_ELEM = 1, KMAXE
!         IF (IELEM(N, I_ELEM)%INTERIOR == 0)THEN ! Element is interior
!             DO I_FACE = 1, IELEM(N, I_ELEM)%IFCA
!                 NEIGHBOR_INDEX = IELEM(N,I)%INEIGH(I_FACE)
!             
!                 LHS_MATRIX = RESHAPE( (/ /), SHAPE(LHS_MATRIX) ) ! Inverse?
!                 RHS_DG_RECONSTRUCT(1) = DG_SOL(N, I_ELEM, IELEM(N, NEIGHBOR_INDEX)%XXC, IELEM(N,NEIGHBOR_INDEX)%YYC, NOF_VARIABLES, IELEM(N, I_ELEM)%IORDER, IELEM(N,I_ELEM)%IDEGREE, U_C(NEIGHBOR_INDEX)%VALDG(1,:,:) - DG_SOL(N, I_ELEM, IELEM(N, NEIGHBOR_INDEX)%XXC, IELEM(N,NEIGHBOR_INDEX)%YYC, NOF_VARIABLES, IELEM(N, I_ELEM)%IORDER, IELEM(N,I_ELEM)%IDEGREE, U_C(I_ELEM)%VALDG(1,:,:))
!                 
!                 DO I_DIM = 1, N_DIM
!                     RHS_DG_RECONSTRUCT(I_DIM+1) = IELEM(N,I_ELEM)%DELTA_XYZ(I_DIM) / IELEM(N,NEIGHBOR_INDEX)%DELTA_XYZ(I_DIM) * U_C(I_ELEM)%VALDG(1,:,I_DIM+1) - U_C(NEIGHBOR_INDEX)%VALDG(1,:,I_DIM+1)
!                 END DO
!             END DO
!         END IF
!     END DO
!     
! END FUNCTION LUO_LSQ_RECONSTRUCT

END MODULE DG_FUNCTIONS
