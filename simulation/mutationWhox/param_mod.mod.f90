MODULE param_modifier

CONTAINS


function r8_normal_ab ( a, b )

!*****************************************************************************80
!
!! R8_NORMAL_AB returns a scaled pseudonormal R8.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, the mean of the PDF.
!
!    Input, real ( kind = rk ) B, the standard deviation of the PDF.
!
!    Output, real ( kind = rk ) R8_NORMAL_AB, a sample of the normal PDF.
!
  implicit none

  !integer, parameter :: rk = kind ( 1.0D+00 )

  real*8 :: a
  real*8 :: b
  real*8 :: r1
  real*8 :: r2
  real*8 :: r8_normal_ab
  real*8, parameter :: r8_pi = 3.141592653589793D+00
  real*8 :: x

  call random_number ( harvest = r1 )
  call random_number ( harvest = r2 )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  r8_normal_ab = a + b * x

  return
end
!*****************************************************************************80

!Modifies the input parameter values according to whether we want to introduce "noise" (for initial population) or mutations
!Generates the modified files and a list which can be used later to run development

!SUBROUTINE pchange(paramen,extras,line1,line2,ini,mutt,gen,FNAM,indi,filename,&
!                    &checkMut,updown,stepmutation,hox0i, hoxvi,hoxmutt,hoxstep,stepin,&
!                    &hoxhnmut) !FOR INDEPENDENT MUTATIONS

SUBROUTINE pchange(paramen,paramen2,extras,extras2,line1,line2,lineA,lineB&
              &,ini,mutt,gen,FNAM,FNAM2,indi,filename,filename2,checkMut,stepmutation&
              &,hox0i,hoxvi,hox0i2,hoxvi2,hoxmutt,hoxstep,stepin,hoxhnmut)

IMPLICIT NONE

INTEGER, INTENT(IN) :: ini, gen, indi !psize
REAL*8, DIMENSION(27), INTENT(IN) :: paramen, paramen2,stepin
CHARACTER(LEN=2024), INTENT(IN) :: line1, line2, lineA, lineB
CHARACTER(LEN=2024), INTENT(IN) :: extras(27), extras2(27)
REAL*8, INTENT(IN) :: mutt,stepmutation,hoxstep,hoxhnmut
CHARACTER*40, INTENT(IN) :: FNAM, FNAM2

REAL*8 :: hoxmutt !chance modifier to mutate hox vector
REAL*8 :: hox0i,hox0i2 !the parameter that multiplies the hox vector
REAL*8,DIMENSION(27) :: hoxvi,hoxvi2 !the hox vector, each position defines how a parameter is affected by hox gene
REAL*8, DIMENSION(28) :: mutaIDsHox, mutaIDsHox2 !safes if parameter mutated or not, 0 no, 1 yes, 28 is for hn
REAL*8, DIMENSION(27) :: mutaIDs, mutaIDs2 !safes if parameter mutated, in genome, or not, 0 no, 1 ye
REAL*8,DIMENSION(27) :: step

CHARACTER(LEN=40) :: filename, filename2
REAL*8, DIMENSION(27) :: randomvec, randomvec2, randomvechox, randomvechox2 !vector containing the modifiers for initial population
!REAL*8, DIMENSION(27) :: safemut !safes values of the mutattion on 1 and on 2 the new value of the mutated parameter
INTEGER :: i, ii, k, j, muttcheck!, mutcount
REAL*8 :: a, m, m2
INTEGER :: updown ! to specify in the name of the file if the teeth up or down in mandibule

INTEGER :: nhoxgene, nhoxvector, ngenome ! number of mutations of each type

CHARACTER(LEN=1) , INTENT(OUT) :: checkMut !1 if file has mutated 0 if not

!open(182,file="normal_test")
!do i=1,1000
    !m = r8_normal_ab ( 0d0, 1d0 )
!    print*,m
!    write(182,*)m 
!enddo
close(182)
!stop
nhoxgene=0 ; nhoxvector=0 ; ngenome=0
mutaIDs = 0 ; mutaIDs2 = 0 
mutaIDsHox = 0 ; mutaIDsHox2 = 0 

step = stepin*stepmutation!*hoxstep

IF(ini==0)then !make the initial population
         checkMut="A"        
         step=step*0.0001             
         CALL RANDOM_NUMBER(randomvec) !Generate vector of random numbers        
         CALL RANDOM_NUMBER(randomvec2) !Generate vector of random numbers
         DO j=1,27
            IF(j==4.OR.j==8.OR.j==10.OR.j==14.OR.j==24)THEN
            	    randomvec(j)=paramen(j) 			            
            	    randomvec2(j)=paramen2(j) 			            
            ELSE
	            m=randomvec(j)*step(j)!*0.01    !0.01 just to make the step smaller  
	            m2=randomvec2(j)*step(j)!*0.01  !0.01 just to make the step smaller 	            
	            CALL RANDOM_NUMBER(a)
	            IF(a>0.5)then
	                m=m*(-1)
	                m2=m2*(-1)
	            ENDIF    
	            randomvec(j)=paramen(j)+m          
	            randomvec2(j)=paramen2(j)+m2           
            ENDIF
         ENDDO    
         
         IF(hox0i<999)THEN !if higher than 999 no hox is actually present, see evol.f90 line 87 and 88
            randomvechox=0
            DO j=1,27
               IF(hoxvi(j)==0)CYCLE !in generation 0 we only introduce some noise in the vector, no changes in the hox network
               IF(j==4.OR.j==8.OR.j==10.OR.j==14.OR.j==24)THEN
            	    randomvechox(j)=hoxvi(j) 			            
               ELSE
                   CALL RANDOM_NUMBER(m)
	            !m=randomvechox(j)*step(j)*hoxstep!*0.01 !0.01 just to make the step smaller            
	            m=m*step(j)*hoxstep!*0.01 !0.01 just to make the step smaller            
	            CALL RANDOM_NUMBER(a)
	            IF(a>0.5)m=m*(-1)
	            randomvechox(j)=hoxvi(j)+m 
               ENDIF
            ENDDO  
         ENDIF
         IF(hox0i2<999)THEN !if higher than 999 is actually not present, see evol.f90 line 87 and 88
            randomvechox2=0
            DO j=1,27
               IF(hoxvi2(j)==0)CYCLE !in generation 0 we only introduce some noise in the vector
               IF(j==4.OR.j==8.OR.j==10.OR.j==14.OR.j==24)THEN
            	    randomvechox2(j)=hoxvi2(j) 			            
               ELSE
                   CALL RANDOM_NUMBER(m)
!	            m=randomvechox2(j)*step(j)*hoxstep!*0.01 !0.01 just to make the step smaller            
	            m=m*step(j)*hoxstep!*0.01 !0.01 just to make the step smaller            	            
	            CALL RANDOM_NUMBER(a)
	            IF(a>0.5)m=m*(-1)
	            randomvechox2(j)=hoxvi2(j)+m 
               ENDIF
            ENDDO  
         ENDIF

            
        !Write the file, both files
         
        !keep track of mutations         

         WRITE (filename, "(A5,A1,A2,I3,A2,I3,A6)") "seal_","3","__",indi,"__",gen,".mutas"
     	 filename=trim(adjustl(filename))
     	 DO ii=1,LEN(filename)
 	             IF(filename(ii:ii)==".")EXIT
 	             IF (filename(ii:ii)==" ") THEN ; filename(ii:ii)="_" ; ENDIF
 	        ENDDO
         OPEN(65,file=trim(adjustl(filename))      )      
	 WRITE(65,*)trim(adjustl(filename)), nhoxgene, nhoxvector, ngenome, gen
         CLOSE(65)

         
         !Upper

         WRITE (filename, "(A8,I3,A7)") "seal_1__",indi,"__0.txt"              
         filename=trim(filename)     
         DO ii=1,LEN(filename) !replace empty spaces with "_"
              IF(filename(ii:ii)==".")EXIT
              IF (filename(ii:ii)==" ") THEN ; filename(ii:ii)="_" ; ENDIF
         ENDDO
         !safe who the parent is and changes in parameters
         OPEN(65,file=trim(filename)//"genealoo")              
         WRITE(65,'(4a)',advance="no")trim(FNAM)," ",trim(filename)," "
         CLOSE(65)

         OPEN(4,FILE=trim(filename))
         DO k=1,27              
               WRITE(4,*) randomvec(k), trim(extras(k))        
         ENDDO    
         WRITE(4,'(A)') trim(line1)
         WRITE(4,'(A)') trim(line2)  
         IF(hox0i<999)WRITE(4,*) "hox ", hox0i , randomvechox 
         CLOSE(4)



	 !Lower

         WRITE (filename2, "(A8,I3,A7)") "seal_2__",indi,"__0.txt"              
         filename2=trim(filename2)     
         DO ii=1,LEN(filename2) !replace empty spaces with "_"
              IF(filename2(ii:ii)==".")EXIT
              IF (filename2(ii:ii)==" ") THEN ; filename2(ii:ii)="_" ; ENDIF
         ENDDO
         !safe who the parent is and changes in parameters
         OPEN(65,file=trim(filename2)//"genealoo")              
         WRITE(65,'(4a)',advance="no")trim(FNAM2)," ",trim(filename2)," "
         CLOSE(65)

         OPEN(4,FILE=trim(filename2))
         DO k=1,27              
               WRITE(4,*) randomvec2(k), trim(extras2(k))        
         ENDDO    
         WRITE(4,'(A)') trim(lineA)
         WRITE(4,'(A)') trim(lineB)  
         IF(hox0i2<999)WRITE(4,*) "hox ", hox0i2 , randomvechox2 
         CLOSE(4)


    
ELSE  ! Mutations on some parameters
     
        !Generate vector of random numbers        
         do i=1,27
            randomvec(i)=r8_normal_ab ( 0d0, 1d0 )
         enddo
         
         randomvechox=0
         randomvechox2=0

         do i=1,27
            randomvec2(i)=r8_normal_ab ( 0d0, 1d0 )
         enddo
         
         muttcheck=0 !checks if parameter has mutated, 0 if no mutation, 1 if it has mutated
         DO j=1,28
            IF(j==4.OR.j==8.OR.j==10.OR.j==14.OR.j==24)THEN
            	randomvec(j)=paramen(j)  !	these parameters do not change 			            
            	randomvec2(j)=paramen2(j)  !	these parameters do not change 			            
            ELSEIF(j==28)THEN !HOX GENE
                IF(hox0i<999)THEN
              	  CALL RANDOM_NUMBER(a) !; print*,"hox HN",mutt*hoxhnmut,mutt,hoxhnmut
                  IF(mutt*hoxhnmut>a)THEN !Chances for a parameter to mutate        
                     !CALL RANDOM_NUMBER(m)
                     !CALL RANDOM_NUMBER(a)
                     !IF(a>0.5) m=m*(-1)  !50% of decreasing value
                     m=r8_normal_ab ( 0d0, 1d0 )
                     hox0i=hox0i+m
                     muttcheck=1 ; nhoxgene = nhoxgene + 1 !; print*,"hn muta !!!"
                     mutaIDsHox(28)=1
                  ENDIF            	
                ENDIF
                IF(hox0i2<999)THEN
              	  CALL RANDOM_NUMBER(a) !; print*,"hox HN",mutt*hoxhnmut,mutt,hoxhnmut
                  IF(mutt*hoxhnmut>a)THEN !Chances for a parameter to mutate        
                     !CALL RANDOM_NUMBER(m)
                     !CALL RANDOM_NUMBER(a)
                     !IF(a>0.5) m=m*(-1)  !50% of decreasing value
                     m=r8_normal_ab ( 0d0, 1d0 )
                     hox0i2=hox0i2+m
                     muttcheck=1 ; nhoxgene = nhoxgene + 2 !; print*,"hn muta !!!"
                     mutaIDsHox2(28)=1
                  ENDIF            	
                ENDIF
            ELSE	
            	CALL RANDOM_NUMBER(a)
            	IF(mutt>a)THEN !Chances for a parameter to mutate     
                       muttcheck=1   ;  ngenome = ngenome + 1 !; print*,"mut genome1",ngenome
                       m=randomvec(j)*step(j)
               	      !CALL RANDOM_NUMBER(a)
               	      !IF(a>0.5)THEN
               	      !   m=m*(-1)  !50% of decreasing value
               	      !ENDIF     
               	      randomvec(j)=paramen(j)+m !mutation 
                             mutaIDs(j)=1  
                ELSE
               	      randomvec(j)=paramen(j)   !no mutation
                ENDIF 
                
                CALL RANDOM_NUMBER(a)
            	IF(mutt>a)THEN !Chances for a parameter to mutate     
                       muttcheck=1   ;  ngenome = ngenome + 1 !; print*,"mut genome2",ngenome
                       m2=randomvec2(j)*step(j)               	
               	 !     CALL RANDOM_NUMBER(a)
               	 !     IF(a>0.5)THEN
               	 !        m2=m2*(-1)  !50% of decreasing value
               	 !     ENDIF     
               	      randomvec2(j)=paramen2(j)+m2 !mutation
                      mutaIDs2(j)=1               	             
                ELSE
               	      randomvec2(j)=paramen2(j)   !no mutation
                ENDIF 
                
                
                !parameter in hox gene regulation
                IF(hox0i<999)THEN
                  CALL RANDOM_NUMBER(a) !; print*,"hox vector",mutt*hoxmutt,mutt,hoxmutt
                  IF(mutt*hoxmutt>a)THEN !Chances for hox vector to mutate
                     !CALL RANDOM_NUMBER(m)
                     m=r8_normal_ab ( 0d0, 1d0 )   
                     !m=randomvechox(j)*step(j)*hoxstep
                     m=m*step(j)*hoxstep                     
                     !CALL RANDOM_NUMBER(a)
                     !IF(a>0.5) m=m*(-1)  !50% of decreasing value
                     randomvechox(j)=hoxvi(j)+m !mutation
                     muttcheck=1 ; nhoxvector = nhoxvector +1 !; print*,"hox vector muta !!!",nhoxvector
                     mutaIDsHox(j)=1
                  ELSE
                     randomvechox(j)=hoxvi(j)   !no mutation
                  ENDIF
                ENDIF
               
               IF(hox0i2<999)THEN
                  CALL RANDOM_NUMBER(a)  !; print*,"hox vector",mutt*hoxmutt,mutt,hoxmutt
                  IF(mutt*hoxmutt>a)THEN !Chances for a parameter to mutate   !to modify hox gene structure more unlikely, we *0.1
                     !CALL RANDOM_NUMBER(m)   
                     !m=randomvechox2(j)*step(j)*hoxstep
                     m=r8_normal_ab ( 0d0, 1d0 )
                     m=m*step(j)*hoxstep
                     !CALL RANDOM_NUMBER(a)
                     !IF(a>0.5) m=m*(-1)  !50% of decreasing value
                     randomvechox2(j)=hoxvi2(j)+m !mutation
                     muttcheck=1 ; nhoxvector = nhoxvector +1 !; print*,"hox vector muta !!!",nhoxvector          
                     mutaIDsHox2(j)=1 
                  ELSE
                     randomvechox2(j)=hoxvi2(j)   !no mutation
                  ENDIF
               ENDIF                             
                 
            ENDIF  
         ENDDO       
         !keep track of mutations         

         WRITE (filename, "(A5,A1,A2,I3,A2,I3,A6)") "seal_","3","__",indi,"__",gen,".mutas"
     	 filename=trim(adjustl(filename))
     	 DO ii=1,LEN(filename)
 	             IF(filename(ii:ii)==".")EXIT
 	             IF (filename(ii:ii)==" ") THEN ; filename(ii:ii)="_" ; ENDIF
 	        ENDDO
         OPEN(65,file=trim(adjustl(filename))      )      
	 WRITE(65,*)trim(adjustl(filename)),nhoxgene,nhoxvector,ngenome,gen,mutaIDs,mutaIDs2,mutaIDsHox,mutaIDsHox2
         CLOSE(65)
     	 
         
         !Write the file
         IF(muttcheck==1)THEN !write new file and adds to list torun only if mutations have occured
                 checkMut="A"
	         
	         !Upper
	         WRITE (filename, "(A8,I3,A2,I3,A4)") "seal_1__",indi,"__",gen,".txt"
        	 filename=trim(adjustl(filename))
 	         DO ii=1,LEN(filename)
 	             IF(filename(ii:ii)==".")EXIT
 	             IF (filename(ii:ii)==" ") THEN ; filename(ii:ii)="_" ; ENDIF
 	        ENDDO
 	        !safe who the parent is
                OPEN(65,file=trim(filename)//"genealoo")              
	        WRITE(65,'(4a)',advance="no")trim(FNAM)," ",trim(filename)," "
               CLOSE(65)
               
 	        OPEN(4,FILE=filename)
 	        DO k=1,27              
 	              WRITE(4,*) randomvec(k), trim(extras(k))      
 	        ENDDO    
 	        WRITE(4,'(A)') trim(line1)
 	        WRITE(4,'(A)') trim(line2)  
 	        IF(hox0i<999)WRITE(4,*) "hox ", hox0i , randomvechox
 	        CLOSE(4)
 	        
	         !Lower 	        
                 WRITE (filename2, "(A8,I3,A2,I3,A4)") "seal_2__",indi,"__",gen,".txt"
        	 filename2=trim(adjustl(filename2))
 	         DO ii=1,LEN(filename2)
 	             IF(filename2(ii:ii)==".")EXIT
 	             IF (filename2(ii:ii)==" ") THEN ; filename2(ii:ii)="_" ; ENDIF
 	        ENDDO
 	        !safe who the parent is
                OPEN(65,file=trim(filename2)//"genealoo")              
	        WRITE(65,'(4a)',advance="no")trim(FNAM2)," ",trim(filename2)," "
               CLOSE(65)
               
 	        OPEN(4,FILE=filename2)
 	        DO k=1,27              
 	              WRITE(4,*) randomvec2(k), trim(extras2(k))      
 	        ENDDO    
 	        WRITE(4,'(A)') trim(lineA)
 	        WRITE(4,'(A)') trim(lineB)
 	        IF(hox0i2<999)WRITE(4,*) "hox ", hox0i2 , randomvechox2	          
 	        CLOSE(4) 	        
 	        
 	        
 	 ELSE !no muta, safe just the name for selection    
 	        checkMut="B"   
 	        
 	        !Upper
 	        WRITE (filename, "(A8,I3,A2,I3,A4)") "seal_1__",indi,"__",gen,".t0t" !.t0t: no mutation
        	filename=trim(filename)     	         	         
 	        DO ii=1,LEN(filename)
 	             IF(filename(ii:ii)==".")EXIT
 	             IF (filename(ii:ii)==" ") THEN ; filename(ii:ii)="_" ; ENDIF
 	        ENDDO
 	        !safe who the parent is
                OPEN(65,file=trim(filename)//"genealoo.t0t")              
	        WRITE(65,'(4a)',advance="no")trim(FNAM)," ",trim(filename)," "
                CLOSE(65)

               OPEN(4,FILE=filename)
 	        DO k=1,27              
 	              WRITE(4,*) paramen(k), trim(extras(k))      
 	        ENDDO    
 	        WRITE(4,'(A)') trim(line1)
 	        WRITE(4,'(A)') trim(line2)  
  	        IF(hox0i<999)WRITE(4,*) "hox ", hox0i , randomvechox
 	        CLOSE(4)
 	        
 	        !Lower
 	        WRITE (filename2, "(A8,I3,A2,I3,A4)") "seal_2__",indi,"__",gen,".t0t" !.t0t: no mutation
        	filename2=trim(filename2)     	         	         
 	        DO ii=1,LEN(filename2)
 	             IF(filename2(ii:ii)==".")EXIT
 	             IF (filename2(ii:ii)==" ") THEN ; filename2(ii:ii)="_" ; ENDIF
 	        ENDDO
 	        !safe who the parent is
                OPEN(65,file=trim(filename2)//"genealoo.t0t")              
	        WRITE(65,'(4a)',advance="no")trim(FNAM2)," ",trim(filename2)," "
                CLOSE(65)

               OPEN(4,FILE=filename2)
 	        DO k=1,27              
 	              WRITE(4,*) paramen2(k), trim(extras2(k))      
 	        ENDDO    
 	        WRITE(4,'(A)') trim(lineA)
 	        WRITE(4,'(A)') trim(lineB)  
 	        IF(hox0i2<999)WRITE(4,*) "hox ", hox0i2 , randomvechox2
 	        CLOSE(4)
 	        
 	        
	 ENDIF

ENDIF

END SUBROUTINE pchange




!CALL pchange(paramen,paramen2,extras,extras2,line1,line2,lineA,lineB,ini,mutt,gen,FNAM,FNAM2,indi,filename,filename2,checkMut,stepmutation)
SUBROUTINE pchange_dependent(paramen,paramen2,extras,extras2,line1,line2,lineA,lineB&
              &,ini,mutt,gen,FNAM,FNAM2,indi,filename,filename2,checkMut,stepmutation&
              &,hox0i,hoxvi,hox0i2,hoxvi2,hoxmutt,hoxstep,stepin,hoxhnmut) !FOR DEPENDENT MUTATIONS


IMPLICIT NONE

INTEGER, INTENT(IN) :: ini, gen, indi !psize
REAL*8, DIMENSION(27), INTENT(IN) :: paramen, paramen2,stepin
CHARACTER(LEN=2024), INTENT(IN) :: line1, line2, lineA, lineB
CHARACTER(LEN=2024), INTENT(IN) :: extras(27), extras2(27)
REAL*8, INTENT(IN) :: mutt,stepmutation,hoxstep,hoxhnmut
CHARACTER*40, INTENT(IN) :: FNAM, FNAM2

REAL*8 :: hoxmutt !chance modifier to mutate hox vector
REAL*8 :: hox0i,hox0i2 !the parameter that multiplies the hox vector
REAL*8,DIMENSION(27) :: hoxvi,hoxvi2 !the hox vector, each position defines how a parameter is affected by hox gene

REAL*8,DIMENSION(27) :: step

CHARACTER(LEN=40) :: filename, filename2
REAL*8, DIMENSION(27) :: randomvec, randomvec2, randomvechox, randomvechox2 !vector containing the modifiers for initial population
!REAL*8, DIMENSION(27) :: safemut !safes values of the mutattion on 1 and on 2 the new value of the mutated parameter
INTEGER :: i, ii, k, j, muttcheck!, mutcount
REAL*8 :: a, m, m2
INTEGER :: updown ! to specify in the name of the file if the teeth up or down in mandibule

REAL*8, DIMENSION(28) :: mutaIDsHox, mutaIDsHox2 !safes if parameter mutated or not, 0 no, 1 yes, 28 is for hn
REAL*8, DIMENSION(27) :: mutaIDs, mutaIDs2 !safes if parameter mutated, in genome, or not, 0 no, 1 ye


INTEGER :: nhoxgene, nhoxvector, ngenome ! number of mutations of each type

CHARACTER(LEN=1) , INTENT(OUT) :: checkMut !1 if file has mutated 0 if not

!size of mutational steps for each parameter
     !        1    , 2 ,   3  ,   4 ,    5   ,    6,        7,   8,  9  , 10, 11   
!DATA step / 0.0001 , 0.2 , 0.0001, 0 ,0.000001 , 0.000001,   0, 0, 0.00001, 0 , 0,&
!     !    12    ,   13  , 14, 15,   16   ,   17  ,   18   ,   19  ,    20	 , 21      ,
!     &  0.00002, 0.00002, 0 , 0 , 0.000095, 0.00001, 0.000002, 0.00001, 0.000015,  0.000018,&
!     !    22  ,   23   ,   24    ,     25   , 26,     27
!     & 0.00001, 0.00001 ,  0     ,   0.00001,  0, 0.000000001 /


!open(182,file="normal_test")
!do i=1,1000
!    m = r8_normal_ab ( 0d0, 1d0 )
!    print*, m
!    write(182,*)m 
!enddo
!close(182)
!stop

nhoxgene=0 ; nhoxvector=0 ; ngenome=0

mutaIDs = 0 ; mutaIDs2 = 0
mutaIDsHox = 0 ; mutaIDsHox2 = 0

step = stepin*stepmutation!*hoxstep

IF(ini==0)then !make the initial population
         step=step*0.0001 
         checkMut="A"                     
         CALL RANDOM_NUMBER(randomvec) !Generate vector of random numbers        
         randomvec2 = randomvec      
         DO j=1,27
            IF(j==4.OR.j==8.OR.j==10.OR.j==14.OR.j==24)THEN
            	    randomvec(j)=paramen(j) 			            
            	    randomvec2(j)=paramen2(j) 			            
            ELSE
	            m=randomvec(j)*step(j)!*0.01    !0.01 just to make the step smaller  
	            m2=randomvec2(j)*step(j)!*0.01  !0.01 just to make the step smaller 	            
	            CALL RANDOM_NUMBER(a)
	            IF(a>0.5)then
	                m=m*(-1)
	                m2=m2*(-1)
	            ENDIF    
	            randomvec(j)=paramen(j)+m          
	            randomvec2(j)=paramen2(j)+m2           
            ENDIF
         ENDDO    
         
         IF(hox0i<999)THEN !if higher than 999 no hox is actually present, see evol.f90 line 87 and 88
            randomvechox=0
            DO j=1,27
               IF(hoxvi(j)==0)CYCLE !in generation 0 we only introduce some noise in the vector, no changes in the hox network
               IF(j==4.OR.j==8.OR.j==10.OR.j==14.OR.j==24)THEN
            	    randomvechox(j)=hoxvi(j) 			            
               ELSE
                   CALL RANDOM_NUMBER(m)
	            !m=randomvechox(j)*step(j)*hoxstep!*0.01 !0.01 just to make the step smaller            
	            m=m*step(j)*hoxstep!*0.01 !0.01 just to make the step smaller            
	            CALL RANDOM_NUMBER(a)
	            IF(a>0.5)m=m*(-1)
	            randomvechox(j)=hoxvi(j)+m 
               ENDIF
            ENDDO  
         ENDIF
         IF(hox0i2<999)THEN !if higher than 999 is actually not present, see evol.f90 line 87 and 88
            randomvechox2=0
            DO j=1,27
               IF(hoxvi2(j)==0)CYCLE !in generation 0 we only introduce some noise in the vector
               IF(j==4.OR.j==8.OR.j==10.OR.j==14.OR.j==24)THEN
            	    randomvechox2(j)=hoxvi2(j) 			            
               ELSE
                   CALL RANDOM_NUMBER(m)
!	            m=randomvechox2(j)*step(j)*hoxstep!*0.01 !0.01 just to make the step smaller            
	            m=m*step(j)*hoxstep!*0.01 !0.01 just to make the step smaller            	            
	            CALL RANDOM_NUMBER(a)
	            IF(a>0.5)m=m*(-1)
	            randomvechox2(j)=hoxvi2(j)+m 
               ENDIF
            ENDDO  
         ENDIF

            
        !Write the file, both files
         
        !keep track of mutations         

         WRITE (filename, "(A5,A1,A2,I3,A2,I3,A6)") "seal_","3","__",indi,"__",gen,".mutas"
     	 filename=trim(adjustl(filename))
     	 DO ii=1,LEN(filename)
 	             IF(filename(ii:ii)==".")EXIT
 	             IF (filename(ii:ii)==" ") THEN ; filename(ii:ii)="_" ; ENDIF
 	        ENDDO
         OPEN(65,file=trim(adjustl(filename))      )      
	 WRITE(65,*)trim(adjustl(filename)), nhoxgene, nhoxvector, ngenome, gen
         CLOSE(65)

         
         !Upper

         WRITE (filename, "(A8,I3,A7)") "seal_1__",indi,"__0.txt"              
         filename=trim(filename)     
         DO ii=1,LEN(filename) !replace empty spaces with "_"
              IF(filename(ii:ii)==".")EXIT
              IF (filename(ii:ii)==" ") THEN ; filename(ii:ii)="_" ; ENDIF
         ENDDO
         !safe who the parent is and changes in parameters
         OPEN(65,file=trim(filename)//"genealoo")              
         WRITE(65,'(4a)',advance="no")trim(FNAM)," ",trim(filename)," "
         CLOSE(65)

         OPEN(4,FILE=trim(filename))
         DO k=1,27              
               WRITE(4,*) randomvec(k), trim(extras(k))        
         ENDDO    
         WRITE(4,'(A)') trim(line1)
         WRITE(4,'(A)') trim(line2)  
         IF(hox0i<999)WRITE(4,*) "hox ", hox0i , randomvechox 
         CLOSE(4)



	 !Lower

         WRITE (filename2, "(A8,I3,A7)") "seal_2__",indi,"__0.txt"              
         filename2=trim(filename2)     
         DO ii=1,LEN(filename2) !replace empty spaces with "_"
              IF(filename2(ii:ii)==".")EXIT
              IF (filename2(ii:ii)==" ") THEN ; filename2(ii:ii)="_" ; ENDIF
         ENDDO
         !safe who the parent is and changes in parameters
         OPEN(65,file=trim(filename2)//"genealoo")              
         WRITE(65,'(4a)',advance="no")trim(FNAM2)," ",trim(filename2)," "
         CLOSE(65)

         OPEN(4,FILE=trim(filename2))
         DO k=1,27              
               WRITE(4,*) randomvec2(k), trim(extras2(k))        
         ENDDO    
         WRITE(4,'(A)') trim(lineA)
         WRITE(4,'(A)') trim(lineB)  
         IF(hox0i2<999)WRITE(4,*) "hox ", hox0i2 , randomvechox2 
         CLOSE(4)


    
ELSE  ! Mutations on some parameters
     
         !CALL RANDOM_NUMBER(randomvec) !Generate vector of random numbers        
         do i=1,27
            randomvec(i)=r8_normal_ab ( 0d0, 1d0 )
         enddo
         !CALL RANDOM_NUMBER(randomvechox)
         randomvechox=0
         !CALL RANDOM_NUMBER(randomvechox2)
         randomvechox2=0
         randomvec2=randomvec
         muttcheck=0 !checks if parameter has mutated, 0 if no mutation, 1 if it has mutated
         DO j=1,28
            IF(j==4.OR.j==8.OR.j==10.OR.j==14.OR.j==24)THEN
            	randomvec(j)=paramen(j)  !	these parameters do not change 			            
            	randomvec2(j)=paramen2(j)  !	these parameters do not change 			            
            ELSEIF(j==28)THEN !HOX GENE
                IF(hox0i<999)THEN
              	  CALL RANDOM_NUMBER(a) !; print*,"hox HN",mutt*hoxhnmut,mutt,hoxhnmut
                  IF(mutt*hoxhnmut>a)THEN !Chances for a parameter to mutate        
                     !CALL RANDOM_NUMBER(m)
                     !CALL RANDOM_NUMBER(a)
                     !IF(a>0.5) m=m*(-1)  !50% of decreasing value
                     m=r8_normal_ab ( 0d0, 1d0 )
                     hox0i=hox0i+m
                     muttcheck=1 ; nhoxgene = nhoxgene + 1 !; print*,"hn muta !!!"
                     mutaIDsHox(28)=1 
                  ENDIF            	
                ENDIF
                IF(hox0i2<999)THEN
              	  CALL RANDOM_NUMBER(a) !; print*,"hox HN",mutt*hoxhnmut,mutt,hoxhnmut
                  IF(mutt*hoxhnmut>a)THEN !Chances for a parameter to mutate        
                     !CALL RANDOM_NUMBER(m)
                     !CALL RANDOM_NUMBER(a)
                     !IF(a>0.5) m=m*(-1)  !50% of decreasing value
                     m=r8_normal_ab ( 0d0, 1d0 )
                     hox0i2=hox0i2+m
                     muttcheck=1 ; nhoxgene = nhoxgene + 2 !; print*,"hn muta !!!"
                     mutaIDsHox2(28)=1
                  ENDIF            	
                ENDIF
            ELSE	
            	CALL RANDOM_NUMBER(a)
            	IF(mutt>a)THEN !Chances for a parameter to mutate     
                       muttcheck=1   ;  ngenome = ngenome + 1
                       m=randomvec(j)*step(j)
                       m2=randomvec2(j)*step(j)               	
               	      !CALL RANDOM_NUMBER(a)
               	      !IF(a>0.5)THEN
               	      !   m=m*(-1)  !50% of decreasing value
               	      !   m2=m2*(-1)  !50% of decreasing value
               	      !ENDIF     
               	      randomvec(j)=paramen(j)+m !mutation
               	      randomvec2(j)=paramen2(j)+m2 !mutation
                             mutaIDs(j)=1 ; mutaIDs2(1)=1               	             
                ELSE
               	      randomvec(j)=paramen(j)   !no mutation
               	      randomvec2(j)=paramen2(j)   !no mutation
                ENDIF 
                !parameter in hox gene regulation
                IF(hox0i<999)THEN
                  CALL RANDOM_NUMBER(a) !; print*,"hox vector",mutt*hoxmutt,mutt,hoxmutt
                  IF(mutt*hoxmutt>a)THEN !Chances for hox vector to mutate
                     !CALL RANDOM_NUMBER(m)   
                     !m=randomvechox(j)*step(j)*hoxstep
                     m=r8_normal_ab ( 0d0, 1d0 )
                     m=m*step(j)*hoxstep                     
                     !CALL RANDOM_NUMBER(a)
                     !IF(a>0.5) m=m*(-1)  !50% of decreasing value
                     randomvechox(j)=hoxvi(j)+m !mutation
                     muttcheck=1 ; nhoxvector = nhoxvector +1 !; print*,"hox vector muta !!!",nhoxvector
                     mutaIDsHox(j)=1
                  ELSE
                     randomvechox(j)=hoxvi(j)   !no mutation
                  ENDIF
                ENDIF
               
               IF(hox0i2<999)THEN
                  CALL RANDOM_NUMBER(a)  !; print*,"hox vector",mutt*hoxmutt,mutt,hoxmutt
                  IF(mutt*hoxmutt>a)THEN !Chances for a parameter to mutate   !to modify hox gene structure more unlikely, we *0.1
                     !CALL RANDOM_NUMBER(m)   
!                     m=randomvechox2(j)*step(j)*hoxstep
                     m=r8_normal_ab ( 0d0, 1d0 )
                     m=m*step(j)*hoxstep
                     !CALL RANDOM_NUMBER(a)
                     !IF(a>0.5) m=m*(-1)  !50% of decreasing value
                     randomvechox2(j)=hoxvi2(j)+m !mutation
                     muttcheck=1 ; nhoxvector = nhoxvector +1 !; print*,"hox vector muta !!!",nhoxvector      
                     mutaIDsHox2(j)=1     
                  ELSE
                     randomvechox2(j)=hoxvi2(j)   !no mutation
                  ENDIF
               ENDIF                             
                 
            ENDIF  
         ENDDO       
         !keep track of mutations         

         WRITE (filename, "(A5,A1,A2,I3,A2,I3,A6)") "seal_","3","__",indi,"__",gen,".mutas"
     	 filename=trim(adjustl(filename))
     	 DO ii=1,LEN(filename)
 	             IF(filename(ii:ii)==".")EXIT
 	             IF (filename(ii:ii)==" ") THEN ; filename(ii:ii)="_" ; ENDIF
 	        ENDDO
         OPEN(65,file=trim(adjustl(filename))      )      
	 WRITE(65,*)trim(adjustl(filename)),nhoxgene,nhoxvector,ngenome,gen,mutaIDs,mutaIDs2,mutaIDsHox,mutaIDsHox2
         CLOSE(65)
     	 
         
         !Write the file
         IF(muttcheck==1)THEN !write new file and adds to list torun only if mutations have occured
                 checkMut="A"
	         
	         !Upper
	         WRITE (filename, "(A8,I3,A2,I3,A4)") "seal_1__",indi,"__",gen,".txt"
        	 filename=trim(adjustl(filename))
 	         DO ii=1,LEN(filename)
 	             IF(filename(ii:ii)==".")EXIT
 	             IF (filename(ii:ii)==" ") THEN ; filename(ii:ii)="_" ; ENDIF
 	        ENDDO
 	        !safe who the parent is
                OPEN(65,file=trim(filename)//"genealoo")              
	        WRITE(65,'(4a)',advance="no")trim(FNAM)," ",trim(filename)," "
               CLOSE(65)
               
 	        OPEN(4,FILE=filename)
 	        DO k=1,27              
 	              WRITE(4,*) randomvec(k), trim(extras(k))      
 	        ENDDO    
 	        WRITE(4,'(A)') trim(line1)
 	        WRITE(4,'(A)') trim(line2)  
 	        IF(hox0i<999)WRITE(4,*) "hox ", hox0i , randomvechox
 	        CLOSE(4)
 	        
	         !Lower 	        
                 WRITE (filename2, "(A8,I3,A2,I3,A4)") "seal_2__",indi,"__",gen,".txt"
        	 filename2=trim(adjustl(filename2))
 	         DO ii=1,LEN(filename2)
 	             IF(filename2(ii:ii)==".")EXIT
 	             IF (filename2(ii:ii)==" ") THEN ; filename2(ii:ii)="_" ; ENDIF
 	        ENDDO
 	        !safe who the parent is
                OPEN(65,file=trim(filename2)//"genealoo")              
	        WRITE(65,'(4a)',advance="no")trim(FNAM2)," ",trim(filename2)," "
               CLOSE(65)
               
 	        OPEN(4,FILE=filename2)
 	        DO k=1,27              
 	              WRITE(4,*) randomvec2(k), trim(extras2(k))      
 	        ENDDO    
 	        WRITE(4,'(A)') trim(lineA)
 	        WRITE(4,'(A)') trim(lineB)
 	        IF(hox0i2<999)WRITE(4,*) "hox ", hox0i2 , randomvechox2	          
 	        CLOSE(4) 	        
 	        
 	        
 	 ELSE !no muta, safe just the name for selection    
 	        checkMut="B"   
 	        
 	        !Upper
 	        WRITE (filename, "(A8,I3,A2,I3,A4)") "seal_1__",indi,"__",gen,".t0t" !.t0t: no mutation
        	filename=trim(filename)     	         	         
 	        DO ii=1,LEN(filename)
 	             IF(filename(ii:ii)==".")EXIT
 	             IF (filename(ii:ii)==" ") THEN ; filename(ii:ii)="_" ; ENDIF
 	        ENDDO
 	        !safe who the parent is
                OPEN(65,file=trim(filename)//"genealoo.t0t")              
	        WRITE(65,'(4a)',advance="no")trim(FNAM)," ",trim(filename)," "
                CLOSE(65)

               OPEN(4,FILE=filename)
 	        DO k=1,27              
 	              WRITE(4,*) paramen(k), trim(extras(k))      
 	        ENDDO    
 	        WRITE(4,'(A)') trim(line1)
 	        WRITE(4,'(A)') trim(line2)  
  	        IF(hox0i<999)WRITE(4,*) "hox ", hox0i , randomvechox
 	        CLOSE(4)
 	        
 	        !Lower
 	        WRITE (filename2, "(A8,I3,A2,I3,A4)") "seal_2__",indi,"__",gen,".t0t" !.t0t: no mutation
        	filename2=trim(filename2)     	         	         
 	        DO ii=1,LEN(filename2)
 	             IF(filename2(ii:ii)==".")EXIT
 	             IF (filename2(ii:ii)==" ") THEN ; filename2(ii:ii)="_" ; ENDIF
 	        ENDDO
 	        !safe who the parent is
                OPEN(65,file=trim(filename2)//"genealoo.t0t")              
	        WRITE(65,'(4a)',advance="no")trim(FNAM2)," ",trim(filename2)," "
                CLOSE(65)

               OPEN(4,FILE=filename2)
 	        DO k=1,27              
 	              WRITE(4,*) paramen2(k), trim(extras2(k))      
 	        ENDDO    
 	        WRITE(4,'(A)') trim(lineA)
 	        WRITE(4,'(A)') trim(lineB)  
 	        IF(hox0i2<999)WRITE(4,*) "hox ", hox0i2 , randomvechox2
 	        CLOSE(4)
 	        
 	        
	 ENDIF

ENDIF


END SUBROUTINE pchange_dependent

















END MODULE param_modifier 
