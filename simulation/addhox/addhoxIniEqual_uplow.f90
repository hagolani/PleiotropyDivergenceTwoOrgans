PROGRAM addhox

IMPLICIT NONE

CHARACTER*40 :: FNAM, FNAM2

INTEGER :: i , j, KTR, k
REAL*8 :: a,b, hoxx

REAL*8 :: mutt !chances to mutate a parameter
CHARACTER(LEN=4) :: muttc

REAL*8 :: modi !mutation step hox


INTEGER, ALLOCATABLE :: seed(:)

INTEGER :: nseed

REAL*8, DIMENSION(27) :: step, hox1, hox2, oriparas1, oriparas2, diff
INTEGER*8, DIMENSION(27) :: epi_or_mesen


CHARACTER*2024 :: extras(27)
CHARACTER(LEN=2024) :: line1, line2


!hox gene upper will regulate epithelial related genes
!hox lower will regulate mesenchymal related gene
!rest genes random

!Random seed

CALL init_random_seed()

!Teeth to model, two independent teeth, read parameters of first teeth, adds random hox. Then uses the parameters without hox for the new tooth. Hox is then added to recover the original morphology. This way with no hox, both teeth have the same morphology
CALL getarg(1,FNAM) !upper

CALL getarg(2,FNAM2) !lower

!chances to be regulated by hox
!CALL getarg(3,muttc)
!READ (muttc,*) mutt

!mutation step hox
!CALL getarg(4,stephc)
!READ (muttc,*) steph


!Read Mutate both teeth ****

OPEN(1,FILE=FNAM,STATUS='OLD')
oriparas1 = 0
KTR=1
DO WHILE(KTR<28)	
     READ(1,*,END=500)oriparas1(KTR), extras(KTR)
     KTR=KTR+1
ENDDO
READ (1, '(A)') line1 
READ (1, '(A)') line2

 CLOSE(1)

OPEN(1,FILE=FNAM2,STATUS='OLD')
oriparas2 = 0
KTR=1
DO WHILE(KTR<28)	
     READ(1,*,END=500)oriparas2(KTR), extras(KTR)
     KTR=KTR+1
ENDDO
READ (1, '(A)') line1 
READ (1, '(A)') line2

 CLOSE(1)



!Introduce extraline that starts with "hox" at the end of the file which contains the hox vector. Each V[i] of the vector determines
 !how gene "i" is affected by the hox gene

hox1 = 0d0 ; hox2 = 0d0



!epithelial parameters
!1, 7,12 (ssh), 9,13 (fgf4), 15,16 (ssh treshold), diff (18)




!mesenchymal parameters
!2, 17, 19, 26
!6,11 (bmp4)

!both
!3, 5, 25 (degradation ... all?), 27

!all molecules (ssh, fgf and bmp in both?) see jernvall and thesleff 1999)

!nothing
!4,8,10,14

!biases
!20 (anterior), 21 (posterior), 22 (buccal), 23 (lingual)

!1= epi , 2=mesen , 0= both ; 3 = nothing


DATA epi_or_mesen /  1, 2, & ! 1 2
                  &  0, 3, & ! 3 4
                  &  0, 0, & ! 5 6 
                  &  0, 3, & ! 7 8
                  &  0, 3, & ! 9 10 
                  &  0, 0, & ! 11 12
                  &  0, 0, & ! 13 14
                  &  0, 0, & ! 15 16
                  &  2, 1, & ! 17 18
                  &  2, 0, & ! 19 20
                  &  0, 0, & ! 21 22
                  &  0, 3, & ! 23 24
                  &  0, 2, & ! 25 26
                  &  0/ ! 27

diff = oriparas2 - oriparas1


DO i=1,27
   IF(epi_or_mesen(i)==3)THEN !nothing
    
     hox1(i) = 0d0 ; hox2(i) = 0d0
     
   ELSEIF(epi_or_mesen(i)==0)THEN !both
     
     CALL random_number(a) !assign randomly to either upper or lower
     IF(a.le.0.5)THEN !to epi hox
       modi = diff(i)*-1
       hox1(i) = modi ; oriparas1(i)=oriparas1(i)-modi
       hox2(i) = 0d0
     
     ELSE             !to mesen hox
       
       hox2(i) = diff(i) ; oriparas2(i)=oriparas2(i)-diff(i)
       hox1(i) = 0d0 
     
     ENDIF 
  
   ELSEIF(epi_or_mesen(i)==1)THEN !epi
     modi = diff(i)*-1
     hox1(i) = modi ; oriparas1(i)=oriparas1(i)-modi
     hox2(i) = 0d0     
     
   ELSE ! == 2             ! mesen
      
     hox2(i) = diff(i) ; oriparas2(i)=oriparas2(i)-diff(i)
     hox1(i) = 0d0     
    
   ENDIF 

ENDDO


hoxx=1d0

!write files

 OPEN(4,FILE=trim(FNAM))
  DO k=1,27              
      WRITE(4,*) oriparas1(k), trim(extras(k))      
  ENDDO    
  WRITE(4,'(A)') trim(line1)
  WRITE(4,'(A)') trim(line2)          
  WRITE(4,*) "hox ", hoxx , hox1

 CLOSE(4)
 
 
 OPEN(4,FILE=trim(FNAM2))
  DO k=1,27              
      WRITE(4,*) oriparas2(k), trim(extras(k))      
  ENDDO    
  WRITE(4,'(A)') trim(line1)
  WRITE(4,'(A)') trim(line2)          
  WRITE(4,*) "hox ", hoxx , hox2

 CLOSE(4)





 STOP 
  
 500 PRINT*, "Somethig wrong with the file, not enough parameters"



END PROGRAM addhox


subroutine init_random_seed()
   !use iso_fortran_env, only: int64
   implicit none
   integer, allocatable :: seed(:)
   integer :: i, n, un, istat, dt(8), pid
   integer :: t

   call random_seed(size = n)
     allocate(seed(n))  
       !  XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365 * 24 * 60 * 60 * 1000 &
            + dt(2) * 31 * 24 * 60 * 60 * 1000 &
            + dt(3) * 24 * 60 * 60 * 1000 &
            + dt(5) * 60 * 60 * 1000 &
            + dt(6) * 60 * 1000 + dt(7) * 1000 &
            + dt(8)
       end if
       pid = getpid()
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
       call random_seed(put=seed)

       contains
      ! This simple PRNG might not be good enough for real work, but is
      ! sufficient for seeding a better PRNG.
        function lcg(s)
          integer :: lcg
          integer :: s
   
          if (s == 0) then
              s = 104729
          else
              s = mod(s, 2147483647)
          end if
          s = mod(s * 279470273, 2147483647)
          lcg = int(mod(s, int(huge(0))), kind(0))
        end function lcg
end subroutine init_random_seed












































