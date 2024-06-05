c      
      SUBROUTINE GSTAR
c
c     DEFINE INITIAL CONFIGURATION OF THE SYSTEM BY ASSIGNING ATOMIC POSITIONS ON THE SURFACE OF CONCENTRIC SPHERES/CUBES
c     THE RADII/SIDES LENGHT OF THE SPHERES/CUBES ARE CONTROLED BY THE VARIABLE DIST DEFINED IN THE MDINPUT.DAT FILE
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                                                                 *
c     * ATTENTION!!! IF THE SOLVENT MOLECULES ARE TOO SMALL USE THE ROUTINE GPOST.f INSTEAD OF GSTAR.f  *
c     *                                                                                                 *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
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
      PARAMETER(NRADMAX=25)
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
      DIMENSION RADII(NRADMAX)     
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
c     DEFINE VECTORS FOR THE TWENTY SIX POSITIONS OF THREE CONCENTRIC SPHERES AND THE ADDITIONAL POINT AT (0.0,0.0,0.0)
c 
c     ATTENTION!!! SMSPHA, SMSPHB AND SMSPHC CAN BE CHOSEN SO THAT THE RADII OF THE 5TH, 6TH AND 7TH SPHERES ARE EQUAL                  
c     ......................................................... 
c     FIRST SPHERE RADIUS   = 0 >>>>>>>>>>>>>>>>>> 1 MOLECULE   
c     .........................................................
c     SECOND SPHERE RADIUS  = DIST >>>>>>>>>>>>>>> 6 MOLECULES
c     THIRD SPHERE RADIUS   = 2*DIST >>>>>>>>>>>>> 8 MOLECULES
c     FOURTH SPHERE RADIUS  = 3*DIST >>>>>>>>>>>>> 12 MOLECULES
c     .........................................................
c     FIFTH SPHERE RADIUS   = 4*DIST-SMSPHA >>>>>> 6 MOLECULES    
c     SIXTH  SPHERE RADIUS  = 5*DIST-SMSPHB >>>>>> 8 MOLECULES
c     SEVENTH SPHERE RADIUS = 6*DIST-SMSPHC >>>>>> 12 MOLECULES
c     etc,...                                                         
c     .........................................................                  
c
      DISTSQRT = DIST*2.0D0/DSQRT(3.0D0)
      DISTRTSQ = DIST*3.0D0/DSQRT(2.0D0)
c
c...  Point in the centre of coordinates - sphere of radius 0
      PFCCX(0)  = 0.0D0  
      PFCCY(0)  = 0.0D0 
      PFCCZ(0)  = 0.0D0
c...  Points on the six faces of an imaginary fcc unit cell of side lenght 2*DIST - sphere of radius DIST     
c     -----------------      
      PFCCX(1)  =  DIST  
      PFCCY(1)  = 0.0D0 
      PFCCZ(1)  = 0.0D0
c     -----------------      
      PFCCX(2)  = -DIST  
      PFCCY(2)  = 0.0D0 
      PFCCZ(2)  = 0.0D0
c     -----------------      
      PFCCX(3)  = 0.0D0  
      PFCCY(3)  =  DIST 
      PFCCZ(3)  = 0.0D0
c     ----------------- 
      PFCCX(4)  = 0.0D0  
      PFCCY(4)  = -DIST 
      PFCCZ(4)  = 0.0D0
c     ----------------- 
      PFCCX(5)  = 0.0D0  
      PFCCY(5)  = 0.0D0 
      PFCCZ(5)  =  DIST
c     -----------------  
      PFCCX(6)  = 0.0D0  
      PFCCY(6)  = 0.0D0 
      PFCCZ(6)  = -DIST
c     -----------------   
c...  Points on the eight corners of an imaginary fcc unit cell of side lenght (2/SQRT(3))*DIST - sphere of radius 2*DIST
c     ---------------------      
      PFCCX(7)  =  DISTSQRT  
      PFCCY(7)  =  DISTSQRT 
      PFCCZ(7)  =  DISTSQRT
c     ---------------------    
      PFCCX(8)  = -DISTSQRT  
      PFCCY(8)  =  DISTSQRT 
      PFCCZ(8)  =  DISTSQRT
c     ---------------------  
      PFCCX(9)  = -DISTSQRT  
      PFCCY(9)  =  DISTSQRT 
      PFCCZ(9)  = -DISTSQRT
c     ---------------------
      PFCCX(10) =  DISTSQRT  
      PFCCY(10) =  DISTSQRT 
      PFCCZ(10) = -DISTSQRT
c     ---------------------
      PFCCX(11) =  DISTSQRT  
      PFCCY(11) = -DISTSQRT 
      PFCCZ(11) =  DISTSQRT
c     --------------------- 
      PFCCX(12) = -DISTSQRT  
      PFCCY(12) = -DISTSQRT 
      PFCCZ(12) =  DISTSQRT
c     ---------------------   
      PFCCX(13) = -DISTSQRT  
      PFCCY(13) = -DISTSQRT 
      PFCCZ(13) = -DISTSQRT
c     --------------------- 
      PFCCX(14) =  DISTSQRT  
      PFCCY(14) = -DISTSQRT 
      PFCCZ(14) = -DISTSQRT
c     ---------------------
c...  Points on the twelve half sides of an imaginary cubic unit cell of side lenght (3/SQRT(2))*DIST - sphere of radius 3*DIST
c     ---------------------      
      PFCCX(15) =  DISTRTSQ  
      PFCCY(15) =  DISTRTSQ 
      PFCCZ(15) =  0.0D0
c     ---------------------    
      PFCCX(16) =  0.0D0  
      PFCCY(16) =  DISTRTSQ 
      PFCCZ(16) =  DISTRTSQ
c     ---------------------  
      PFCCX(17) = -DISTRTSQ  
      PFCCY(17) =  DISTRTSQ 
      PFCCZ(17) =  0.0D0
c     ---------------------
      PFCCX(18) =  0.0D0  
      PFCCY(18) =  DISTRTSQ 
      PFCCZ(18) = -DISTRTSQ
c     ---------------------
      PFCCX(19) =  DISTRTSQ  
      PFCCY(19) =  0.0D0 
      PFCCZ(19) =  DISTRTSQ
c     --------------------- 
      PFCCX(20) = -DISTRTSQ  
      PFCCY(20) =  0.0D0 
      PFCCZ(20) =  DISTRTSQ
c     ---------------------   
      PFCCX(21) = -DISTRTSQ  
      PFCCY(21) =  0.0D0 
      PFCCZ(21) = -DISTRTSQ
c     --------------------- 
      PFCCX(22) =  DISTRTSQ  
      PFCCY(22) =  0.0D0 
      PFCCZ(22) = -DISTRTSQ
c     ---------------------
      PFCCX(23) =  DISTRTSQ  
      PFCCY(23) = -DISTRTSQ 
      PFCCZ(23) =  0.0D0
c     ---------------------   
      PFCCX(24) =  0.0D0  
      PFCCY(24) = -DISTRTSQ 
      PFCCZ(24) =  DISTRTSQ
c     --------------------- 
      PFCCX(25) = -DISTRTSQ  
      PFCCY(25) = -DISTRTSQ 
      PFCCZ(25) =  0.0D0
c     ---------------------
      PFCCX(26) =  0.0D0  
      PFCCY(26) = -DISTRTSQ 
      PFCCZ(26) = -DISTRTSQ
c     ---------------------
c
      IF (NMOLCS.LE.26)THEN
         NSPHER  = 4
         RADII(1) = 0.0D0
         RADII(2) = DIST*BANG 
         RADII(3) = 2.0D0*DIST*BANG 
         RADII(4) = 3.0D0*DIST*BANG 
         WRITE(*,*)'                                                  '
         WRITE(*,'(A62)')'Four concentric spheres are enough to assign 
     &initial positions'
         WRITE(*,'(A23,4(F9.3))')'Spheres radii are (Ang)',RADII(1),
     &RADII(2),RADII(3),RADII(4)
         WRITE(*,*)'                                                  '
      ELSE  
         RADII(1) = 0.0D0
         RADII(2) = DIST*BANG 
         RADII(3) = 2.0D0*DIST*BANG 
         RADII(4) = 3.0D0*DIST*BANG  
         NLATTP   = IDNINT(DFLOAT(NMOLCS)/DFLOAT(26)+0.49D0)
         FNLATTP  = DFLOAT(NLATTP)
         NSPHER   = 1+3*NLATTP
c...Define additional positions on concentric spheres of larger radii
c   This loop creates in general a much larger number of positions than those needed, i.e., than the number of molecules
c   ATTENTION!!! Decrease the radius of additional spheres by subtracting SMPHA, SMPHB and SMPHC to condense the system
         KL=26
         KW=4
         KV=5
         KP=6
         KK=0
         KU=10
         KH=20
c        ............................................................
         DO 5 IL = 2,NLATTP
            DO 15 JL = 1,26
               IF(JL.LE.6)THEN
                  SMSPHA=DFLOAT(KK)/10.0D0
                  PFCCX(JL+KL) = PFCCX(JL)*(DFLOAT(KW)-SMSPHA)
                  PFCCY(JL+KL) = PFCCY(JL)*(DFLOAT(KW)-SMSPHA)
                  PFCCZ(JL+KL) = PFCCZ(JL)*(DFLOAT(KW)-SMSPHA)
                  RADII(IL*3-1)  = (DFLOAT(IL*3-2)-SMSPHA)*DIST*BANG
               ELSEIF((JL.GT.6).AND.(JL.LE.14))THEN
                  SMSPHB=(DFLOAT(KU)/2.0D0)/10.0D0
                  PFCCX(JL+KL) = PFCCX(JL)*(DFLOAT(KV)/2.0D0-SMSPHB)
                  PFCCY(JL+KL) = PFCCY(JL)*(DFLOAT(KV)/2.0D0-SMSPHB)
                  PFCCZ(JL+KL) = PFCCZ(JL)*(DFLOAT(KV)/2.0D0-SMSPHB)
                  RADII(IL*3)= (DFLOAT(IL*3-1)-2.0D0*SMSPHB)*DIST*BANG
               ELSEIF((JL.GT.14).AND.(JL.LE.26))THEN
                  SMSPHC=(DFLOAT(KH)/3.0D0)/10.0D0
                  PFCCX(JL+KL) = PFCCX(JL)*(DFLOAT(KP)/3.0D0-SMSPHC)
                  PFCCY(JL+KL) = PFCCY(JL)*(DFLOAT(KP)/3.0D0-SMSPHC)
                  PFCCZ(JL+KL) = PFCCZ(JL)*(DFLOAT(KP)/3.0D0-SMSPHC)
                  RADII(IL*3+1)= (DFLOAT(IL*3)-3.0D0*SMSPHC)*DIST*BANG
               ENDIF
   15       CONTINUE
            KL=KL+26
            KW=KW+3
            KV=KV+3
            KP=KP+3
            KK=KK+20
            KU=KU+20
            KH=KH+20
c   
    5    CONTINUE
c   
c...Aproximate volume of the smallest sphere that still contains the system 
         RADMAX  = DIST*BANG*(3.0D0*FNLATTP-3.0D0*SMSPHC)
         VSPHMAX = (4.0D0/3.0D0)*PI*(RADMAX)**3.0D0
c 
         WRITE(*,*)'                                                  '
         WRITE(*,'(A36,I3,A20)')'Positions assigned on the surface of',
     &NSPHER,'  concentric spheres'
         WRITE(*,'(A69,F15.3)')'Radii of the smallest sphere that 
     &still contains the system (A)  < = ',RADMAX
         WRITE(*,'(A69,E15.3)')'Volume of the smallest sphere that 
     &still contains the system (A³) < = ',VSPHMAX
         WRITE(*,*)'                                                  '
      ENDIF 
c...Print out to MDOUTPUT.DAT the radii of the spheres used to assign initial positions
         DO 25 N=1,NSPHER
            WRITE(6,'(A8,I2,A14,F11.3))')'Sphere ',N,' radius (A) = '
     &,RADII(N)
  25     CONTINUE
         WRITE(6,*)'                                                  ' 
         WRITE(6,69)     
c     
      IF(LTYPE.EQ.0)THEN
c...Pass NMOLCS positions on the surface of the concentric spheres to a two dimensional array SPHPX starting from component 1       
      NN=0
      DO 35 I = 1,NCOMP
         JMAX = NMOL(I)
         DO 45 J=1,JMAX 
            SPHPX(I,J)=PFCCX(NN)
            SPHPY(I,J)=PFCCY(NN)
            SPHPZ(I,J)=PFCCZ(NN)
            NN=NN+1
   45    CONTINUE
   35 CONTINUE  
c             
      ELSEIF(LTYPE.EQ.1)THEN
c...Pass NMOLCS positions on the surface of the concentric spheres to a two dimensional array SPHPX starting from component NCOMP      
      NN=0
      DO 55 I = NCOMP,1,-1
         JMAX = NMOL(I)
         DO 65 J=1,JMAX 
            SPHPX(I,J)=PFCCX(NN)
            SPHPY(I,J)=PFCCY(NN)
            SPHPZ(I,J)=PFCCZ(NN)
            NN=NN+1
   65    CONTINUE
   55 CONTINUE        
      ENDIF          
c
c     ASSIGN ATOMIC POSITIONS TO THE POINTS ON THE SURFACE OF THE CONCENTRIC SPHERES
c
      DO 75 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 85 J=1,JMAX 
            DO 95 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 105 L = 1,LMAX
                  RX(I,J,K,L)=RX0(I,K,L)+SPHPX(I,J)
                  RY(I,J,K,L)=RY0(I,K,L)+SPHPY(I,J)
                  RZ(I,J,K,L)=RZ0(I,K,L)+SPHPZ(I,J)
  105          CONTINUE
   95       CONTINUE      
   85    CONTINUE
   75 CONTINUE
c
c     CALL MDROT TO RANDOMLY ROTATE MOLECULES 
c 
      LRAND=.FALSE.
c      
      IF(LRAND)THEN
         DO 115 I = 1,NCOMP
            IF((NMOL(I).GT.3).AND.(NATMOL(I).GT.1))CALL MDROT
  115    CONTINUE       
      ENDIF
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
      DO 125 LL=1,NGA
         WRITE(5,'(A115)')NLINEG(LL)
  125 CONTINUE      
c
c     WRITE HEAD OF MDPOST.DAT
c
      WRITE(3,'((7X,A7),2(17X,A7))')'X(Bohr)','Y(Bohr)','Z(Bohr)'
c
      DO 135 I = 1,NCOMP
         JMAX = NMOL(I)
         KMAX = NSPEC(I)
         DO 145 J=1,JMAX 
            DO 155 K = 1,KMAX
               LMAX = NATMSP(I,K)
               DO 165 L = 1,LMAX
                  WRITE(3,*)RX(I,J,K,L),RY(I,J,K,L),RZ(I,J,K,L)
                  WRITE(5,'(2X,I3,1X,3(1X,F19.14))')NATZ(I,K,L),
     &RX(I,J,K,L),RY(I,J,K,L),RZ(I,J,K,L)   
  165          CONTINUE
  155       CONTINUE      
  145    CONTINUE
  135 CONTINUE
c...  The last line of Gaussian 03 input file must be empty       
      WRITE(5,*)'                                                     '
c
c...  SPECIFY ATOMS'S BASIS SETS IF THE OPTION GEN IS USED IN GAUSSIAN
c     ATTENTION!!! NGB CAN BE > NLG IN WHICH CASE GEN IS NOT USED IN THE GAUSSIAN INPUT FILE
c                  AND NO LINES ARE WRITTEN IN THE FOLLOWING LOOP
c  
      IF(NGB.LT.NLG)THEN
         DO 175 KK=NGB,NLG
            WRITE(5,'(A115)')NLINEG(KK)
  175    CONTINUE
      ENDIF          
c
c     ---------------------------------------------------------------------------------
c   
      CALL MDVISUAL
c
69    FORMAT(T1,107('_'))  
      RETURN
      END 
c
c
c      