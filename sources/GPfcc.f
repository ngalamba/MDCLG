      SUBROUTINE GPFCC
c
c     DEFINE INITIAL CONFIGURATION OF THE SYSTEM BY ASSIGNING ATOMIC POSITIONS ON A FCC NaCl TYPE OF LATTICE
c     ATTENTION!!! THE LENGHT OF THE SIDE OF THE MD BOX IS GIVEN BY DIST DEFINED IN THE MDINPUT.DAT FILE
c
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c
      PARAMETER(NCOMPMAX=5)
      PARAMETER(NMOLMAX=256)
      PARAMETER(NSPECMAX=10)
      PARAMETER(NATSPMAX=100) 
      PARAMETER(NPARTCL=500)
      PARAMETER(NLGAUSS=500)
c      
      PARAMETER(BANG=0.5291772083D0)   
      PARAMETER(PI=3.141592654D0)
c 
      CHARACTER*150 NLINEG
      LOGICAL LRAND
c                                                    
      COMMON/INTGRS/IN(14)
      COMMON/REALS /RL(18)     
c      
      COMMON/POSIT /RX(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              RY(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              RZ(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)
      COMMON/POSSD /XX(NPARTCL),YY(NPARTCL),ZZ(NPARTCL)
      COMMON/POST0 /RX0(NCOMPMAX,NSPECMAX,NATSPMAX),
     &              RY0(NCOMPMAX,NSPECMAX,NATSPMAX),
     &              RZ0(NCOMPMAX,NSPECMAX,NATSPMAX) 
      DIMENSION PFCCX(NPARTCL),PFCCY(NPARTCL),PFCCZ(NPARTCL)
      DIMENSION SPHPX(NCOMPMAX,NMOLMAX),SPHPY(NCOMPMAX,NMOLMAX),
     &          SPHPZ(NCOMPMAX,NMOLMAX) 
      COMMON/PARMOL/NMOL(NMOLMAX),NSPEC(NSPECMAX),
     &NATMSP(NCOMPMAX,NATSPMAX)
      COMMON/ATMASS/ZMASS(NCOMPMAX,NSPECMAX,NATSPMAX)
      COMMON/ATOMZ /NATZ(NCOMPMAX,NSPECMAX,NATSPMAX)   
      COMMON/GAUSSL/NLINEG(NLGAUSS)
c...  Number of atoms p/ molecule for each component      
      COMMON/ATCOMP/NATMOL(NPARTCL)
c
      EQUIVALENCE (IN(1),NCOMP),(IN(2),NLG),(IN(8),NMOLCS)
      EQUIVALENCE (IN(9),LTYPE),(IN(12),NGA),(IN(13),NGB)
      EQUIVALENCE (RL(11),FNATOMS),(RL(14),FNMOLCS),(RL(15),DIST)
c    
c     CALCULATE NUMBER OF UNIT CELLS p/ DIMENSION - NUNIT
c
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     !  NUMBER OF MOLECULES = (NUMBER OF MOLECULES p/ UNIT CELL)*(NUMBER OF UNIT CELLS p/ DIMENSION)**3 !
c     !  FCC NaCl TYPE OF LATTICE                                                                        !
c     !  NUMBER OF CATIONS p/ UNIT CELL = 12/4 + 1   = 4                                                 !
c     !  NUMBER OF ANIONS  p/ UNIT CELL =  8/8 + 6/2 = 4                                                 !
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
c                                    
      NUNIT = (FNMOLCS/8.0D0)**(1.0D0/3.0D0)+0.1D0                                
c     
c     CHECK WETHER USER ASKED FOR NON-FCC NUMBER OF MOLECULES
c             
    5 NCHECK = 8*(NUNIT**3)                                        
      IF (NCHECK.LT.NMOLCS) THEN                                    
         NUNIT = NUNIT + 1                                          
         GOTO 5                                                     
      ENDIF                                                        
c     
c     DEFINE LATTICE DISTANCE BASED ON CUBE OF SIDE DIST
c
      CUBE=DIST
      CUBEA=CUBE*BANG
c  
c     RE-DEFINE DIST TO HALF OF THE LENGHT OF THE SIDE OF THE UNIT CELL (FCC)
c             
      DIST=0.5D0*CUBE/DFLOAT(NUNIT)
      NMOLCS2=NMOLCS/2
c
c     CENTER OF THE MD BOX IS (DIST,DIST,DIST)_ALL COORDINATES ARE POSITIVE
c 
c     
c     ASSIGN POSITIONS OF FIRST FOUR ANIONS_COMPONENT_1                             
c     -------------------------------------------------      
      PFCCX(1)  = 0.0D0  
      PFCCY(1)  = 0.0D0 
      PFCCZ(1)  = 0.0D0
c     -----------------      
      PFCCX(2)  = 0.0D0  
      PFCCY(2)  =  DIST 
      PFCCZ(2)  =  DIST
c     -----------------      
      PFCCX(3)  =  DIST  
      PFCCY(3)  = 0.0D0 
      PFCCZ(3)  =  DIST
c     ----------------- 
      PFCCX(4)  =  DIST  
      PFCCY(4)  =  DIST 
      PFCCZ(4)  = 0.0D0
c     -----------------       
c     
c     ASSIGN POSITIONS OF FIRST FOUR CATIONS_COMPONENT_2                            
c     --------------------------------------------------
      PFCCX(NMOLCS2+1)  =  DIST  
      PFCCY(NMOLCS2+1)  =  DIST 
      PFCCZ(NMOLCS2+1)  =  DIST
c     -------------------------      
      PFCCX(NMOLCS2+2)  = 0.0D0  
      PFCCY(NMOLCS2+2)  = 0.0D0 
      PFCCZ(NMOLCS2+2)  =  DIST
c     -------------------------      
      PFCCX(NMOLCS2+3)  = 0.0D0  
      PFCCY(NMOLCS2+3)  =  DIST 
      PFCCZ(NMOLCS2+3)  = 0.0D0
c     ------------------------- 
      PFCCX(NMOLCS2+4)  =  DIST  
      PFCCY(NMOLCS2+4)  = 0.0D0 
      PFCCZ(NMOLCS2+4)  = 0.0D0
c     ------------------------- 
c     
c     REPLICATE FIRST FOUR POSITIONS OF ANIONS AND CATIONS OVER NUNITS  
c                     
      M=0                                                            
      MM = 0                                                        
      DO 15 IW=1,NUNIT                                                
         DO 25 JW=1,NUNIT                                                
            DO 35 KW=1,NUNIT                                                
               DO 45 IJW=1,4                                                 
c     
c     IF NUMBER OF POSITIONS ASSIGNED EXCEEDS NUMBER OF MOLECULES GET OUT
c    
*                  IF (MM.LT.NMOLCS2) THEN                                     
                     PFCCX(IJW+M)=PFCCX(IJW)+2.0D0*DIST*(KW-1)                         
                     PFCCY(IJW+M)=PFCCY(IJW)+2.0D0*DIST*(JW-1)                         
                     PFCCZ(IJW+M)=PFCCZ(IJW)+2.0D0*DIST*(IW-1)
                     PFCCX(IJW+M+NMOLCS2)=PFCCX(IJW+NMOLCS2)+
     &2.0D0*DIST*(KW-1)                         
                     PFCCY(IJW+M+NMOLCS2)=PFCCY(IJW+NMOLCS2)+
     &2.0D0*DIST*(JW-1)                         
                     PFCCZ(IJW+M+NMOLCS2)=PFCCZ(IJW+NMOLCS2)+
     &2.0D0*DIST*(IW-1)                         
*                  ENDIF                                                      
                  MM = MM + 1                                              
   45          CONTINUE                                                     
               M=M+4                                                        
   35       CONTINUE
   25    CONTINUE
   15 CONTINUE
c 
      WRITE(6,'(A9,F11.3)')'L (A) = ',CUBEA
      WRITE(6,*)'                                                     ' 
      WRITE(6,79)  
c       
      IF(LTYPE.EQ.2)THEN
      NN=1
      DO 55 I = 1,NCOMP
         JMAX = NMOL(I)
         DO 65 J=1,JMAX 
            SPHPX(I,J)=PFCCX(NN)
            SPHPY(I,J)=PFCCY(NN)
            SPHPZ(I,J)=PFCCZ(NN)
            NN=NN+1
   65    CONTINUE
   55 CONTINUE        
      ELSEIF(LTYPE.EQ.3)THEN
      NN=1
      DO 75 I = NCOMP,1,-1
         JMAX = NMOL(I)
         DO 85 J=1,JMAX 
            SPHPX(I,J)=PFCCX(NN)
            SPHPY(I,J)=PFCCY(NN)
            SPHPZ(I,J)=PFCCZ(NN)
            NN=NN+1
   85    CONTINUE
   75 CONTINUE        
      ENDIF 
c 
c     ASSIGN ATOMIC POSITIONS TO THE POINTS ON THE FCC NaCl TYPE OF LATTICE
c
      DO 95 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 105 J=1,JMAX 
            DO 115 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 125 L = 1,LMAX
                  RX(I,J,K,L)=RX0(I,K,L)+SPHPX(I,J)
                  RY(I,J,K,L)=RY0(I,K,L)+SPHPY(I,J)
                  RZ(I,J,K,L)=RZ0(I,K,L)+SPHPZ(I,J)
  125          CONTINUE
  115       CONTINUE      
  105    CONTINUE
   95 CONTINUE   
c
c     CALL MDROT TO RANDOMLY ROTATE MOLECULES 
c 
      LRAND=.TRUE.
      IF(LRAND)CALL MDROT
c      
*      IF(LRAND)THEN
*         DO 135 I = 1,NCOMP
*            IF((NMOL(I).GT.1).AND.(NATMOL(I).GT.1))CALL MDROT
*  135    CONTINUE       
*      ENDIF
c
c     ---------------------------------------------------------------------------------
c  
c     PRINT OUT THE TWO INITIAL CARTESIAN COORDINATES INPUT FILES 
c
c     gaussr.com - Positions at step 1/2..........[Gausian 03 input]
c     MDPOST.DAT - Positions at step 1/2..........[MDCLG input]
c
c     WRITE GAUSSIAN OPTIONS TO THE NEW POSITION FILE gaussr.dat 
c 
      DO 145 LL=1,NGA
         WRITE(5,'(A115)')NLINEG(LL)
  145 CONTINUE      
c
c     WRITE HEAD OF MDPOST.DAT
c
      WRITE(3,'((7X,A7),2(17X,A7))')'X(Bohr)','Y(Bohr)','Z(Bohr)'
c
      DO 155 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 165 J=1,JMAX 
            DO 175 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 185 L = 1,LMAX
                  WRITE(3,*)RX(I,J,K,L),RY(I,J,K,L),RZ(I,J,K,L)
                  WRITE(5,'(2X,I3,1X,3(1X,F19.14))')NATZ(I,K,L),
     &RX(I,J,K,L),RY(I,J,K,L),RZ(I,J,K,L)   
  185          CONTINUE
  175       CONTINUE      
  165    CONTINUE
  155 CONTINUE
c...  The last line of Gaussian 03 input file must be empty       
      WRITE(5,*)'                                                     '
c
c...  SPECIFY ATOMS'S BASIS SETS IF THE OPTION GEN IS USED IN GAUSSIAN
c
c     ATTENTION!!! NGB CAN BE > NLG IN WHICH CASE GEN WAS NOT USED IN THE GAUSSIAN INPUT FILE
c                  AND NO LINES ARE WRITTEN IN THE FOLLOWING LOOP
c  
      IF(NGB.LT.NLG)THEN
         DO 195 KK=NGB,NLG
            WRITE(5,'(A115)')NLINEG(KK)
  195    CONTINUE
      ENDIF          
c
c     ---------------------------------------------------------------------------------
c   
      CALL MDVISUAL
c
79    FORMAT(T1,107('_')) 
      RETURN
      END 
c
c
c      