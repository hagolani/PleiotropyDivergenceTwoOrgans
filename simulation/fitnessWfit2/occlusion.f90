PROGRAM OCCLUSION

!USE OVERLAP

!calculate the optimal occlusion between two teeth
!it considers: overlapping area and distance between the cells.
!minimizes both 
!also option to maximize surface area of teeth (here of in selection?)

!assumes that the middle cusps of both teeth are always neighbors of 
!each other

!      /\   \    2   /       \    2   /  /\
!     /  \   \      /         \      /  /  \
!    /    \   \    /   or      \    /  /    \
!   /      \   \  /             \  /  /      \
!  /    1   \   \/               \/  /    1   \


!move the teeth in certain area until overlap happens
!stablish the limits determined by the overlap: 
!inside that space, find the minimal distance between cells


!1. find hard limits: area where we look for the optimal possitioing
!   area determined by: 
!       ----- cusp of 1. horizontal
!      ||\   \    2   /
!     / | \   \      / 
!    /  |  \   \    /  
!   /   |   \   \  /   
!  /    | 1  \   \/
!       |    ------ valley of 1.  
!       cusp of 1. vertical 
!
! next cusp of 1. 

!2. inside hard limit, find overlaping limits, due to teeth 1 and 2 
! overlaping
!3. inside overlaping limit, find minimal distance between cells.
!    divide in squares, keep square with minimal distance, repeat... 
!    TODO: how many times necessary to do this?

!to consider: although we are looking into this specific region, the
!overlaping and distance calculated is for the WHOLE teeth


IMPLICIT NONE
!read files
REAL*8, ALLOCATABLE :: POLY1(:,:), POLY2(:,:)  !to read the initial file, before number of cells is known
REAL*8, ALLOCATABLE :: profile1(:,:) , profile2(:,:) !to work with smaller matrices, transfer POLY 1 and 2 here
REAL*8, ALLOCATABLE :: z2(:) !safe the z-values from profile2

CHARACTER*40 ::  FNAM1, FNAM2, FNAM3, FNAM4     !name of files, 1 and 3 profiles, 2 and 4 toothmaker input files
INTEGER :: KTR1, KTR2 !number of cells in each profile
INTEGER :: cuspd  !cell in cusp of upper teeth

INTEGER :: gen !generation
CHARACTER(LEN=7) :: genc

INTEGER :: fittype !generation
CHARACTER(LEN=4) :: fittypecc

REAL*8 :: cusp2dv(2)  !coordinates of tallest cusp teeth 2 

REAL*8, DIMENSION(4) :: hlimits

REAL*8 :: fitt, dtf1 ,dtf2

REAL*8 :: xspam1 ,xspam2

ALLOCATE(POLY1(1000,2)) ; POLY1=0D0
ALLOCATE(POLY2(1000,2)) ; POLY2=0D0




!Reading input files of polygons

!FNAM1 and FNAM3 have the profiles to be "occluded". 
!FNAM2 FNAM4 has the name of the input file for toothmaker (needed to mutate later)

CALL GETARG(1,FNAM1) ; CALL GETARG(2,FNAM2) 
CALL GETARG(3,FNAM3) ; CALL GETARG(4,FNAM4)

!Generation
CALL getarg(5,genc)
READ (genc,*) gen

!Fitness type. 1: distance based. 2: surface based
CALL getarg(6,fittypecc)
READ (fittypecc,*) fittype

print*,"A ",FNAM1, FNAM2, FNAM3, FNAM4, gen, fittype

IF(fittype==4)THEN
   !development run out of time
   fitt=1984
   GOTO 687
ENDIF   


1000 FORMAT(A40)
     OPEN(1,FILE=FNAM1,STATUS='OLD')
     OPEN(2,FILE=FNAM3,STATUS='OLD')

!READ IN THE POLYGONS

     KTR1=1
     
1    CONTINUE
     READ(1,*,END=500)POLY1(KTR1,1),POLY1(KTR1,2)
     KTR1=KTR1+1 
     GOTO 1

500  CONTINUE
     IF(KTR1<5)THEN
        fitt=9999 ; GOTO 687
     ENDIF
     
     KTR2=1

2    CONTINUE
     READ(2,*,END=501)POLY2(KTR2,1),POLY2(KTR2,2)
     KTR2=KTR2+1
     GOTO 2

501 CONTINUE 
    IF(KTR2<5)THEN
        fitt=9999 ; GOTO 687
    ENDIF
    
    !number of cells in each teeth
    KTR1=KTR1-1 !; print*,KTR1
    KTR2=KTR2-1 !; print*,KTR2

    CLOSE(1)
    CLOSE(2)
    
     
    allocate(profile1(ktr1,2)) ; profile1 = POLY1(1:ktr1,1:2)
    DEALLOCATE(POLY1)

    allocate(profile2(ktr2,2)) ; profile2 = POLY2(1:ktr2,1:2)
    DEALLOCATE(POLY2)

  
      

    !Max x distance, make sure tooth do not get too small     
     xspam1 = profile1(ktr1,1) - profile1(1,1) 
     xspam2 = profile2(ktr2,1) - profile2(1,1) 
  !   PRINT*,"xspam1",xspam1,"xspam2",xspam2
     IF(xspam1<5.or.xspam2<5)THEN
       fitt=8888 ; GOTO 687
     ENDIF
      

   
  IF(fittype==3)THEN !conservative selection
    
  !  print*,"doing 140"

    CALL conservative(profile1,profile2,KTR1,KTR2,dtf1,dtf2,FNAM1, FNAM2)
    
    fitt = (dtf1 + dtf2) / 2 !simply the average of both tooth
    
   ! print*,"fit",fitt
    
  ELSE
    
    !Limits in which the "main" cusp of teeth 2 can be positioned
!    PRINT*,"Calling hard limits",ktr1,"***",profile1,"*-----*"    
    CALL HARD_LIMITS(profile1,ktr1,hlimits)
!    PRINT*,"Done with hard limits"
    ALLOCATE(Z2(ktr2))
    
    !Find best occluding position: 1. No overlap 2. Minimal distance between cells    
    z2=profile2(:,2) 
    cuspd=MINLOC(z2,DIM=1) !the cell in the lowest z-axis (in upper mandible)
    cusp2dv(1)=profile2(cuspd,1)
    cusp2dv(2)=profile2(cuspd,2)
    DEALLOCATE(z2)

    IF(fittype==1)THEN
!       print*,"type 1"
       CALL OCCLUSION_FITTING(hlimits,profile1,profile2,ktr1,ktr2,cusp2dv,FNAM1,FNAM3,fitt)!mindistance)
    ELSE
!       print*,"type 2"    
       CALL OCCLUSION_SURFACE(hlimits,profile1,profile2,ktr1,ktr2,cusp2dv,FNAM1,FNAM3,fitt)!maxocclusion)
    ENDIF 
    
  ENDIF  
    
  687 CONTINUE
    
  OPEN(2,FILE=TRIM(FNAM1)//TRIM(FNAM2)//"FITNESS.dat")
  WRITE(2,*)	FNAM1," ",FNAM2," ", FNAM3," ", FNAM4," ",fitt," ", gen
  CLOSE(2)
    

END PROGRAM OCCLUSION














SUBROUTINE conservative(profile1,profile2,ktr1,ktr2,dtf1,dtf2,FNAM1, FNAM2)

!Calculate how much a teeth has change comparing input 1 to 2

IMPLICIT NONE



REAL*8 , INTENT(IN):: profile1(ktr1,2) , profile2(ktr2,2) 
INTEGER, INTENT(IN) :: KTR1, KTR2 !number of cells in each profile
INTEGER :: KTRa, KTRb !number of cells in each profile

REAL*8, ALLOCATABLE :: profilea(:,:) , profileb(:,:)

REAL*8 :: d0, d1, dt, x1, y1 , dtf1 , dtf2 
INTEGER :: ii, jj ,i

REAL*8, ALLOCATABLE :: POLY1(:,:), POLY2(:,:)  !to read the initial file, before number of cells is known


CHARACTER*40,INTENT(IN) ::  FNAM1, FNAM2

ALLOCATE(POLY1(1000,2)) ; POLY1=0D0
ALLOCATE(POLY2(1000,2)) ; POLY2=0D0



1000 FORMAT(A40)
     OPEN(1,FILE="sealmaster1.txtprofile__up.txt",STATUS='OLD')

     OPEN(2,FILE="sealmaster2.txtprofile_down.txt",STATUS='OLD')

!READ IN THE POLYGONS

     KTRa=1
     
1    CONTINUE
     READ(1,*,END=500)POLY1(KTRa,1),POLY1(KTRa,2)
     KTRa=KTRa+1 
     GOTO 1

500  CONTINUE

     KTRb=1

2    CONTINUE
     READ(2,*,END=501)POLY2(KTRb,1),POLY2(KTRb,2)
     KTRb=KTRb+1
     GOTO 2

501 CONTINUE

    !number of cells in each teeth
    KTRa=KTRa-1
    KTRb=KTRb-1

    CLOSE(1)
    CLOSE(2)


    allocate(profilea(ktra,2)) ; profilea = POLY1(1:ktra,1:2)
    DEALLOCATE(POLY1)
    
    allocate(profileb(ktrb,2)) ; profileb = POLY2(1:ktrb,1:2)
    DEALLOCATE(POLY2)

!uppers

    !using EMD

     dt = 0
     DO ii=1,ktr1
         d1=1d5 ; x1 = profile1(ii,1) ; y1 = profile1(ii,2)
         DO jj=1,ktra
              d0=(profilea(jj,1)-x1)**2+(profilea(jj,2)-y1)**2
              IF(d0<d1)d1=d0
         ENDDO
!write(87,*)sqrt(d1)
         dt=dt+sqrt(d1)
      ENDDO
      !distance from teeth two to one
      DO ii=1,ktra
         d1=1d5 ; x1 = profilea(ii,1) ; y1 = profilea(ii,2)
         DO jj=1,ktr1
              d0=(profile1(jj,1)-x1)**2+(profile1(jj,2)-y1)**2
              IF(d0<d1)d1=d0
         ENDDO
!write(87,*)sqrt(d1)
         dt=dt+sqrt(d1)
      ENDDO
      dtf1=dt/(ktr1+ktra)



!lowers 

    !using EMD
     dt = 0
     DO ii=1,ktr2
         d1=1d5 ; x1 = profile2(ii,1) ; y1 = profile2(ii,2)
         DO jj=1,ktrb
              d0=(profileb(jj,1)-x1)**2+(profileb(jj,2)-y1)**2
              IF(d0<d1)d1=d0
         ENDDO
         dt=dt+sqrt(d1)
      ENDDO
      !distance from teeth two to one
      DO ii=1,ktrb
         d1=1d5 ; x1 = profileb(ii,1) ; y1 = profileb(ii,2)
         DO jj=1,ktr2
              d0=(profile2(jj,1)-x1)**2+(profile2(jj,2)-y1)**2
              IF(d0<d1)d1=d0
         ENDDO
         dt=dt+sqrt(d1)
      ENDDO
      dtf2=dt/(ktr2+ktrb)
  







END SUBROUTINE conservative





SUBROUTINE HARD_LIMITS(profile1,ktr1,hlimits)
!Find the area in which the teeth will move around to find the optimal occlusion
!Based only on teeth one, the "static" one, the second teeth will move around to 
!fit into the 1. one

IMPLICIT NONE

INTEGER*8 :: cusp1, cusp1b, cusp2 !cell in first cusp
INTEGER ,INTENT(IN) :: ktr1 
INTEGER :: i, down, celly1
REAL*8 , INTENT(IN):: profile1(ktr1,2) 
REAL*8 :: z1(ktr1)
REAL*8 :: cusp1v, cusp2v, z, x1, x2, y1, y2
REAL*8, DIMENSION(4), INTENT(OUT) :: hlimits


!find central and central adjacent cusps 
!print*,"l 190"
z1=profile1(:,2) 

!teeth one
 cusp1=MAXLOC(z1,DIM=1) !the cell in the highest z-axis (in lower mandible)
! print*,"l 195"
 cusp1v=z1(cusp1)
! print*,"cusp1",cusp1,"cusp1v",cusp1v
 !secondary cusps
!To consider: cells are in order. We look for the next cell that is at the peak of a cusp
!Add also that they need to be knot cells? TODO think ...
down=0 ; z=cusp1v
!DO i=cusp1,1,-1  !we start with the cell in the heighest cusp. We look for the next highest cusp, which should be between two valleys
		     !OJO look at order ***

!test skipping this part, occlude in any part from the middle towards end of the tooth                     
!if(1==2)then
!DO i=cusp1,ktr1!
!    !look for 1. valley
!    IF(z1(i)<z.AND.down==0)THEN
!	z=z1(i) 
!    ELSEIF(z1(i)>z)THEN
!        IF(down==0) celly1=i-1  ! lower part of valley between the two cusps, will be one of the limits
!    	z=z1(i) ; down=1 !; print*,"2A","celly",celly1 
!    ELSEIF(z1(i)<z.AND.down==1)THEN
!!        cusp1b=i+1 ; EXIT   !safe the previous value, its going downhill again
!        cusp1b=i-1 ; EXIT   !safe the previous value, its going downhill again
!    ENDIF
!    
!ENDDO
!endif

!IF Only one cusp or higher cusp after the lowest one
!IF(down==0)THEN
!     Print*,"no second cusp was found",1
     !celly1=1 ; cusp1b=1 !just the "end" of the teeth
 celly1=ktr1 ; cusp1b=ktr1 !just the "end" of the tooth
!ENDIF


! celly1=ktr1 ; cusp1b=ktr1

!limits are:
x1=profile1(cusp1,1) 
x2=profile1(cusp1b,1) !; print*,"x1",x1,"x2",x2

y1=cusp1v
y2=profile1(celly1,2) !; print*,"y1",y1,"y2",y2
! x1=lower x-axis ; x2=highest x-axis
! y1=highest z-axis ; y2=lowest z-axis
hlimits(1)=x1 ; hlimits(2)=x2 
hlimits(3)=y1 ; hlimits(4)=y2


 !open(182,file="limitssee.xyz")
 !write(182,*)x1,y1,0
 !write(182,*)x2,y2,0
 !close(182)


!PRINT*,"limits are:",hlimits

END SUBROUTINE HARD_LIMITS


 
SUBROUTINE OCCLUSION_FITTING(hlimits,profile1,profile2,ktr1,ktr2,cusp2dv,FNAM1,FNAM2,mindistance)

USE OVERLAP

IMPLICIT NONE

INTEGER ,INTENT(IN) :: ktr1 , ktr2
REAL*8, DIMENSION(4), INTENT(IN) :: hlimits
REAL*8 , INTENT(IN):: profile1(ktr1,2) , profile2(ktr2,2)
REAL*8 :: profile2_adjusted(ktr2,2)
REAL*8, INTENT(IN) :: cusp2dv(2)
CHARACTER*40, INTENT(IN) ::  FNAM1, FNAM2
REAL*8 :: trans(2)
REAL*8 :: overl

REAL*8 :: v1, v2, v1_100, v2_100, x0, y0

REAL*8 , DIMENSION(100,100,2) :: explore

INTEGER :: i, j, kk, ii, jj, k

CHARACTER(LEN=1024) :: filename

REAL*8 :: maxover, minover , distance, mindistance   

REAL*8, DIMENSION(10000,ktr2,2) :: all_profiles2

!EMD
REAL*8 :: d0, d1, dt, x1, y1

INTEGER :: theone !the chosen one with the best occlusion

!divide the area stablish by the limits into equal parts

v1=sqrt((hlimits(1)-hlimits(2))**2)
v2=sqrt((hlimits(3)-hlimits(4))**2)

v1_100=v1/100 ; v2_100=v2/100

mindistance=1d3

!1000 new coordinates
!btter to use fix steps?? so y and x equally explored OJO
x0=hlimits(2) ; y0=hlimits(4)
!print*,"v1_100",v1_100,"v2_100",v2_100,x0,y0

DO i=1,100 !optimize!! TODO vectorize?!
  DO j=1,100
         explore(i,j,1)=x0
	 y0=y0+v2_100 ; explore(i,j,2)=y0
  ENDDO
  x0=x0+v1_100	  
  y0=hlimits(4)
ENDDO

open(4,file="explore.xyz")
do i=1,100
  do j=1,100
      write(4,*) explore(i,j,:),0
  enddo
enddo
close(4)
!stop

!translocate the upper teeth and calculate fitness
maxover=0 ; minover=1d5
all_profiles2=0d0
kk=0 ; distance=0d0
DO i=1,100 !optimize!! TODO vectorize?!
  DO j=1,100
    kk=kk+1 !; print*,"kk",kk
    trans=cusp2dv-explore(i,j,:)  !vector necessary to move the cusp of teeth two to a position in explore matrix
    
    profile2_adjusted=profile2 - SPREAD(trans(1:2),1,ktr2) !apply this vector to all cells in teeth two
  !  print*,trans ; stop
!    WRITE (filename, "(A7,I1,A1,I1,A4)") "prof_ad",i,"_",j,".xyz"     
!    OPEN(4,FILE=filename)    
!    DO k=1,ktr2    
!       WRITE(4,*) profile2_adjusted(k,:),0    
!    ENDDO    
!    CLOSE(4)   
    
    all_profiles2(kk,:,:)=profile2_adjusted
   
   ! TODO : calculate overlap and distance ... keep the best ... see  
   
   CALL overlapping(profile2_adjusted,ktr2,profile1,ktr1,overl)
   
   IF(maxover<overl)maxover=overl
   IF(minover>overl)minover=overl
   !if no overlap, calculate the distance between cells (EMD?), keep the lowest one
   IF(overl<1d-20)THEN !how much overlap is allowed
   
   
      !safe files with adjusted coordinates, not necessary for program     
 !   WRITE (filename, "(A7,I2,A1,I2,A4)") "prof_ad",i,"_",j,".xyz"   
 !   OPEN(4,FILE=filename)    
 !   DO k=1,ktr2    
 !      WRITE(4,*) profile2_adjusted(k,:),0    
 !   ENDDO    
 !   CLOSE(4)    
   
   
     ! print*,"checking distance"
      dt=0d0 !total distance
      !distance from teeth one to two
      DO ii=1,ktr1
         d1=1d5 ; x1 = profile1(ii,1) ; y1 = profile1(ii,2)
         DO jj=1,ktr2
              d0=(profile2_adjusted(jj,1)-x1)**2+(profile2_adjusted(jj,2)-y1)**2
              IF(d0<d1)d1=d0
         ENDDO
         dt=dt+sqrt(d1)
      ENDDO
      !distance from teeth two to one
      DO ii=1,ktr2
         d1=1d5 ; x1 = profile2_adjusted(ii,1) ; y1 = profile2_adjusted(ii,2)
         DO jj=1,ktr1
              d0=(profile1(jj,1)-x1)**2+(profile1(jj,2)-y1)**2
              IF(d0<d1)d1=d0
         ENDDO
         dt=dt+sqrt(d1)
      ENDDO
      dt=dt/(ktr1+ktr2)

     IF(dt<mindistance)THEN      
       mindistance=dt
       theone=kk      
     ENDIF
   ENDIF
 
   
   !safe files with adjusted coordinates, not necessary for program     
!    WRITE (filename, "(A7,I1,A1,I1,A4)") "prof_ad",i,"_",j,".xyz"     
!    OPEN(4,FILE=filename)    
!    DO k=1,ktr2    
!       WRITE(4,*) profile2_adjusted(k,:),0    
!    ENDDO    
!    CLOSE(4)    
!    A = A - spread( B(1,:), 1, 250 )
  !  print*,1, trans
  !  print*,2, cusp2dv
  !  print*,3, explore(i,j,:)
  ! if(j==2)stop
  ENDDO
ENDDO

open(4,file=TRIM(FNAM2)//"adjusted.xyz")
do i=1,ktr2
  write(4,*)all_profiles2(theone,i,:),"0"
enddo
close(4)

open(4,file=TRIM(FNAM1)//"static.xyz")
do i=1,ktr1
  write(4,*)profile1(i,:),"0"
enddo
close(4)

END SUBROUTINE OCCLUSION_FITTING






SUBROUTINE OCCLUSION_SURFACE(hlimits,profile1,profile2,ktr1,ktr2,cusp2dv,FNAM1,FNAM2,maxocclusion)

USE OVERLAP

IMPLICIT NONE

INTEGER ,INTENT(IN) :: ktr1 , ktr2
REAL*8, DIMENSION(4), INTENT(IN) :: hlimits
REAL*8 , INTENT(IN):: profile1(ktr1,2) , profile2(ktr2,2)
REAL*8 :: profile2_adjusted(ktr2,2)
REAL*8, INTENT(IN) :: cusp2dv(2)
CHARACTER*40, INTENT(IN) ::  FNAM1, FNAM2
REAL*8 :: trans(2)
REAL*8 :: overl

REAL*8 :: v1, v2, v1_100, v2_100, x0, y0

REAL*8 , DIMENSION(100,100,2) :: explore

INTEGER :: i, j, kk, ii, jj, k, iii

CHARACTER(LEN=1024) :: filename

REAL*8 :: d !, mindistance   

REAL*8, DIMENSION(10000,ktr2,2) :: all_profiles2

!EMD
REAL*8 :: d0, d1, dt, x1, y1, contactDistance, x2, y2

INTEGER, DIMENSION(KTR1) :: cellsInContact1 !, cellsInContact1_SAFE !SAFE THE IDS OF CELLS IN MORPHOLOGY 1 THAT ARE CLOSE ENOUGH TO MORPHOLOGY 2
INTEGER, DIMENSION(KTR2) :: cellsInContact2 !, cellsInContact2_SAFE !SAFE THE IDS OF CELLS IN MORPHOLOGY 2 THAT ARE CLOSE ENOUGH TO MORPHOLOGY 1

INTEGER :: theone !the chosen one with the best occlusion

REAL*8 :: maxocclusion

!divide the area stablish by the limits into equal parts

v1=sqrt((hlimits(1)-hlimits(2))**2)
v2=sqrt((hlimits(3)-hlimits(4))**2)

v1_100=v1/100 ; v2_100=v2/100

!mindistance=1d3

contactDistance = 0.5

!1000 new coordinates
!btter to use fix steps?? so y and x equally explored OJO
x0=hlimits(1) ; y0=hlimits(4)
!print*,"v1_100",v1_100,"v2_100",v2_100,x0,y0

DO i=1,100 !optimize!! TODO vectorize?!
  DO j=1,100
         explore(i,j,1)=x0
	 y0=y0+v2_100 ; explore(i,j,2)=y0
  ENDDO
  x0=x0+v1_100	  
  y0=hlimits(4)
ENDDO

open(4,file="explore.xyz")
do i=1,100
  do j=1,100
      write(4,*) explore(i,j,:),0
  enddo
enddo
close(4)
!stop

!translocate the upper teeth and calculate fitness

maxocclusion = 0.0001
all_profiles2=0d0
kk=0 
DO i=1,100 !optimize!! TODO vectorize?!
  DO j=1,100
    kk=kk+1 !; print*,"kk",kk
    trans=cusp2dv-explore(i,j,:)  !vector necessary to move the cusp of teeth two to a position in explore matrix
    
    profile2_adjusted=profile2 - SPREAD(trans(1:2),1,ktr2) !apply this vector to all cells in teeth two
  !  print*,trans ; stop
!    WRITE (filename, "(A7,I1,A1,I1,A4)") "prof_ad",i,"_",j,".xyz"     
!    OPEN(4,FILE=filename)    
!    DO k=1,ktr2    
!       WRITE(4,*) profile2_adjusted(k,:),0    
!    ENDDO    
!    CLOSE(4)   
    
    all_profiles2(kk,:,:)=profile2_adjusted
   
   ! TODO : calculate overlap and distance ... keep the best ... see  
   
   CALL overlapping(profile2_adjusted,ktr2,profile1,ktr1,overl)
   
   !if no overlap, calculate the distance between cells (EMD?), keep the lowest one
   IF(overl<1d-20)THEN !how much overlap is allowed
   
   
      !safe files with adjusted coordinates, not necessary for program     
 !   WRITE (filename, "(A7,I2,A1,I2,A4)") "prof_ad",i,"_",j,".xyz"   
 !   OPEN(4,FILE=filename)    
 !   DO k=1,ktr2    
 !      WRITE(4,*) profile2_adjusted(k,:),0    
 !   ENDDO    
 !   CLOSE(4)    
   
   

      cellsInContact1 = 0 ; cellsInContact2 = 0
      !find cells which are close to other teeth
      DO ii=1,ktr1
         x1 = profile1(ii,1) ; y1 = profile1(ii,2)
         DO jj=1,ktr2
              d0=sqrt((profile2_adjusted(jj,1)-x1)**2+(profile2_adjusted(jj,2)-y1)**2) 
              IF(d0<contactDistance)THEN
		 cellsInContact1(ii)=1 ; EXIT ! IF CLOSE ENOUGH TO ONE CELL IN OTHER MORPHOLOGY NO NEED TO LOOK IF CLOSE TO OTHER CELLS
              ENDIF
         ENDDO
      ENDDO

      !surface of contact of those cells
      dt=0d0 !total distance
      DO ii=1,ktr1
         IF(cellsInContact1(ii)==0)CYCLE
	 x1 = profile1(ii,1) ; y1 = profile1(ii,2)
         iii=ii+1 !the next cell, they are in order
         IF(cellsInContact1(iii)==0)CYCLE !no cell after this one
	 x2 = profile1(iii,1) ; y2 = profile1(iii,2)
	 d = sqrt ( (x1-x2)**2 + (y1-y2)**2 ) !; print*,"*** ",d
         dt = dt + d

      ENDDO

 !find cells which are close to other teeth
      DO ii=1,ktr2
         x1 = profile2_adjusted(ii,1) ; y1 = profile2_adjusted(ii,2)
         DO jj=1,ktr1
              d0=sqrt((profile1(jj,1)-x1)**2+(profile1(jj,2)-y1)**2)
              IF(d0<contactDistance)THEN
		 cellsInContact2(ii)=1 ; EXIT! IF CLOSE ENOUGH TO ONE CELL IN OTHER MORPHOLOGY NO NEED TO LOOK IF CLOSE TO OTHER CELLS
              ENDIF
         ENDDO
      ENDDO
      !surface of contact of those cells

      DO ii=1,ktr2
         IF(cellsInContact2(ii)==0)CYCLE
	 x1 = profile2_adjusted(ii,1) ; y1 = profile2_adjusted(ii,2)
         iii=ii+1 !the next cell, they are in order
         IF(cellsInContact2(iii)==0)CYCLE !no cell after this one
	 x2 = profile2_adjusted(iii,1) ; y2 = profile2_adjusted(iii,2)
	 d = sqrt ( (x1-x2)**2 + (y1-y2)**2 )
         dt = dt + d

      ENDDO

   ENDIF
 

 
   IF(maxocclusion<dt)THEN
     maxocclusion=dt
     theone=kk
!     cellsInContact2_SAFE = cellsInContact2 !JUST FOR TESTING
!     cellsInContact1_SAFE = cellsInContact1 !JUST FOR TESTING
!     print*, "occlusion distance",dt,"sum1 ", sum(cellsInContact1),"sum2 ", sum(cellsInContact2)
   ENDIF
  
   !safe files with adjusted coordinates, not necessary for program     
!    WRITE (filename, "(A7,I1,A1,I1,A4)") "prof_ad",i,"_",j,".xyz"     
!    OPEN(4,FILE=filename)    
!    DO k=1,ktr2    
!       WRITE(4,*) profile2_adjusted(k,:),0    
!    ENDDO    
!    CLOSE(4)    
!    A = A - spread( B(1,:), 1, 250 )
  !  print*,1, trans
  !  print*,2, cusp2dv
  !  print*,3, explore(i,j,:)
  ! if(j==2)stop
  ENDDO
ENDDO

!print*, "MAX occlusion distance",maxocclusion

!open(4,file="POINTS1.xyz")
!do i=1,ktr1
!  IF(cellsInContact1_SAFE(i)==0)CYCLE
!  write(4,*)profile1(i,:),"0"
!enddo
!close(4)

!open(4,file="POINTS2.xyz")
!do i=1,ktr2
!  IF(cellsInContact2_SAFE(i)==0)CYCLE
!  write(4,*)all_profiles2(theone,i,:),"0"
!enddo
!close(4)


maxocclusion=1/maxocclusion ! only because evlution.sh script writen to take the 50% with lowest fitness

open(4,file=TRIM(FNAM2)//"adjusted.xyz")
do i=1,ktr2
  write(4,*)all_profiles2(theone,i,:),"0"
enddo
close(4)

open(4,file=TRIM(FNAM1)//"static.xyz")
do i=1,ktr1
  write(4,*)profile1(i,:),"0"
enddo
close(4)

END SUBROUTINE OCCLUSION_SURFACE






