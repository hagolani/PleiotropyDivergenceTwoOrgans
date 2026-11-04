PROGRAM evol

!For teeth evolution. Introduces mutations or noise if it is the 1. generation. 

USE param_modifier !changes parameter values in a file and created a new file with the changes. Can introduce "noise" for initial population or mutate. 

IMPLICIT NONE

CHARACTER*40 :: FNAM, FNAM2, filename, filename2
CHARACTER*80 :: file2
CHARACTER*40, DIMENSION(2) ::  indis
CHARACTER*40, DIMENSION(2,2) ::new_indis
REAL*8, DIMENSION(27) :: paramen, paramen2
!CHARACTER(LEN=1024) :: extras(27) c
CHARACTER*2024 :: extras(27), extras2(27)
INTEGER :: ktr, KTR2
INTEGER :: psize !population size
CHARACTER(LEN=3) :: psizec 
CHARACTER(LEN=2024) :: line1, line2, lineA, lineB
INTEGER :: ini !indicates to param_mod if changes in parameters should be small (noise), when =0 or big (mutations), when =1.
CHARACTER(LEN=3) :: inic
REAL*8 :: mutt !chances to mutate a parameter
CHARACTER(LEN=10) :: muttc

REAL*8 :: stepmutation !chances to mutate a parameter
CHARACTER(LEN=10) :: stepmutationc

REAL*8 :: pleitropic !do teeth change independently? 1 independent ; else mutations the same for both
CHARACTER(LEN=4) :: pleitropicc

INTEGER :: gen !generation
CHARACTER(LEN=4) :: genc
INTEGER :: indi, i !individual in certain generation
CHARACTER(LEN=4) :: indic
CHARACTER(LEN=1) :: checkMut !A if file has mutated B if not

REAL*8 :: x_pfh,y_pfh,z_pfh

INTEGER, ALLOCATABLE :: seed(:)
INTEGER :: nseed

CHARACTER(LEN=4) :: qwerty !if lines starts with "hox", this is the hox line :P
REAL*8 :: hox0,hox02 !the parameter that multiplies the hox vector
REAL*8,DIMENSION(27) :: hoxv,hoxv2 !the hox vector, each position defines how a parameter is affected by hox gene

REAL*8,DIMENSION(27) :: step ! the size of the mutation for each parameter

REAL*8 :: hoxx !hox !mutae hox gene? 1 yes, else no. If trying to use with input file that does not have the line, it will simply run as if no hox gene
!if hox gene is present it will still be used, for now that line needs to be taken out if you do not want to use it in developing the teeth
CHARACTER(LEN=7) :: hoxc

REAL*8 :: hoxmutt !chance modifier to mutate hox vector
CHARACTER(LEN=10) :: hoxmuttc

REAL*8 :: hoxstep !steo multiplier for hox
CHARACTER(LEN=10) :: hoxstepc

REAL*8 :: hoxhnmut !chance modifier to mutate hox gene - hn -
CHARACTER(LEN=10) :: hoxhnmutc



CALL init_random_seed()

!Teeth to model, two independent teeth
CALL getarg(1,FNAM)
CALL getarg(2,FNAM2)

indis(1)=FNAM ; indis(2)=FNAM2 

!is it the initial generation? 0==yes ; 1==no
CALL getarg(3,inic)
READ (inic,*) ini

!chances to mutate a parameter
CALL getarg(4,muttc)
READ (muttc,*) mutt

!Generation
CALL getarg(5,genc)
READ (genc,*) gen

!individual
CALL getarg(6,indic)
READ (indic,*) indi

!step modifier
CALL getarg(7,stepmutationc)
READ (stepmutationc,*) stepmutation


!pleitropic?
CALL getarg(8,pleitropicc)
READ (pleitropicc,*) pleitropic

!hox?
CALL getarg(9,hoxc)
READ (hoxc,*) hoxx

!hox vector probability to mutate
CALL getarg(10,hoxmuttc)
READ (hoxmuttc,*) hoxmutt

!hox vector step size when mutates
CALL getarg(11,hoxstepc)
READ (hoxstepc,*) hoxstep

!hox gene (hn) probabibility to mutate
CALL getarg(12,hoxhnmutc)
READ (hoxhnmutc,*) hoxhnmut


!print*,FNAM,FNAM2,ini,mutt,gen,indi,stepmutation,pleitropic,hoxx,hoxmutt,hoxstep,hoxhnmut

IF(ini==0.and.pleitropic==1)THEN  !this is just to create duplicate population for 2 genome population based on the 1 genome popu
     stepmutation=0d0
ENDIF


!size of mutational steps for each parameter
!           1        2     3       4   5           6        7      8       9        10  11  
DATA step / 0.0001 , 1.1 , 0.006 , 0 , 0.0000080 , 0.007 ,  0.15 , 0.001 , 0.00018 , 0 , 0.0013 ,&
! 12      13      14  15       16       17       18           19     20     21 
& 0.0013 , 0.0013 ,0 ,0.0013 , 0.007  , 0.0008 , 0.0000045 , 0.007 , 0.10 , 0.13,&
! 22      23     24   25       26  27
& 0.0075 , 0.0075 , 0 , 0.0008 , 0 , 0.00000005 /


!Read and Mutate both teeth ****

hox0=1D4 ; hoxv=0D0
hox02=1D4 ; hoxv2=0D0

IF(pleitropic==1)THEN !TOOTH MUTATE INDEPENDENTLY 
        FNAM=indis(1) ; FNAM2=indis(2)
	OPEN(1,FILE=FNAM,STATUS='OLD')
	KTR=1
	DO WHILE(KTR<28)	
           READ(1,*,END=500)paramen(KTR), extras(KTR)
           KTR=KTR+1
	ENDDO
	read (1, '(A)') line1 
	read (1, '(A)') line2 
	
	IF(hoxx==1)THEN
          READ(1,*,END=500) qwerty
          IF(qwerty=="hox")THEN
             BACKSPACE(1)
             READ(1,*,END=500) qwerty, hox0, hoxv         
          ENDIF
        ENDIF
	
        CLOSE(1)    
       
        OPEN(1,FILE=FNAM2,STATUS='OLD')
	KTR2=1
	DO WHILE(KTR2<28)	
           READ(1,*,END=500)paramen2(KTR2), extras2(KTR2)
           KTR2=KTR2+1
	ENDDO
	read (1, '(A)') lineA 
	read (1, '(A)') lineB 
       
        IF(hoxx==1)THEN
          READ(1,*,END=500) qwerty
          IF(qwerty=="hox")THEN
             BACKSPACE(1)
             READ(1,*,END=500) qwerty, hox02, hoxv2         
          ENDIF
        ENDIF
       
        CLOSE(1)    
   
       !Generate population by introducing some random and very small changes in the parameters. 
       !Or mutate if ini==1

       CALL pchange(paramen,paramen2,extras,extras2,line1,line2,lineA,lineB,&
                &ini,mutt,gen,FNAM,FNAM2,indi,filename,filename2,checkMut,stepmutation,&
                &hox0,hoxv,hox02,hoxv2,hoxmutt,hoxstep,step,hoxhnmut)
                
       
       new_indis(1,1)=filename ; new_indis(1,2)=checkMut
       new_indis(2,1)=filename2 ; new_indis(2,2)=checkMut
       
ELSE  !MUTATIONS AFFECT UPPER AND LOWER EQUALLY

        FNAM=indis(1) ; FNAM2=indis(2)
	OPEN(1,FILE=FNAM,STATUS='OLD')
	KTR=1
	DO WHILE(KTR<28)	
           READ(1,*,END=500)paramen(KTR), extras(KTR)
           KTR=KTR+1
	ENDDO
	read (1, '(A)') line1 
	read (1, '(A)') line2 
	
	IF(hoxx==1)THEN
          READ(1,*,END=500) qwerty
          IF(qwerty=="hox")THEN
             BACKSPACE(1)
             READ(1,*,END=500) qwerty, hox0, hoxv         
          ENDIF
        ENDIF
	
        CLOSE(1)    
       
        OPEN(1,FILE=FNAM2,STATUS='OLD')
	KTR2=1
	DO WHILE(KTR2<28)	
           READ(1,*,END=500)paramen2(KTR2), extras2(KTR2)
           KTR2=KTR2+1
	ENDDO
	read (1, '(A)') lineA 
	read (1, '(A)') lineB 
       
        IF(hoxx==1)THEN
          READ(1,*,END=500) qwerty
          IF(qwerty=="hox")THEN
             BACKSPACE(1)
             READ(1,*,END=500) qwerty, hox02, hoxv2         
          ENDIF
        ENDIF
       
        CLOSE(1)    
   
       !Generate population by introducing some random and very small changes in the parameters. 
       !Or mutate if ini==1

       CALL pchange_dependent(paramen,paramen2,extras,extras2,line1,line2,lineA,lineB,&
                &ini,mutt,gen,FNAM,FNAM2,indi,filename,filename2,checkMut,stepmutation,&
                &hox0,hoxv,hox02,hoxv2,hoxmutt,hoxstep,step,hoxhnmut)
                
       
       new_indis(1,1)=filename ; new_indis(1,2)=checkMut
       new_indis(2,1)=filename2 ; new_indis(2,2)=checkMut

 


ENDIF



IF(new_indis(1,2)=="B".AND.new_indis(2,2)=="B")THEN !NO Mutations in any teeth
  OPEN(65,file=trim(new_indis(1,1))//trim(new_indis(2,1))//"genealoo.t00t")              
  WRITE(65,*)trim(indis(1))," ",trim(indis(2))," ",trim(new_indis(1,1))," ",trim(new_indis(2,1))!,&
!  &" ",trim(new_indis(2,1))," ",trim(new_indis(2,2))," "
  CLOSE(65)
ENDIF  
  
!ELSE !At least one teeth has mutated

  OPEN(65,file=trim(new_indis(1,1))//trim(new_indis(2,1))//"TORUN1")                
  WRITE(65,*)trim(indis(1))," ",trim(indis(2))," ",trim(new_indis(1,1))," ",trim(new_indis(1,2)),&
  &" ",trim(new_indis(2,1))," ",trim(new_indis(2,2))," "
  CLOSE(65)

!ENDIF


STOP

500 PRINT*, "Somethig wrong with the file, not enough parameters"

END PROGRAM evol






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

