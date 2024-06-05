c      
      SUBROUTINE GPOST
c
c     DEFINE INITIAL CONFIGURATION OF THE SYSTEM BY ASSIGNING ATOMIC POSITIONS ON THE SURFACE OF CONCENTRIC SPHERES
c     ATTENTION!!! THE RADIUS OF THE SPHERES IS CONTROLED BY THE VARIABLE DIST DEFINED IN THE MDINPUT.DAT FILE
c     ATTENTION!!! THIS ROUTINE IS OPTIMAL FOR SIXTY FOUR MOLECULES OF SOLVENT AND 1 MOLECULE OF SOLUTE
c     ATTENTION!!! FOR NMOLCS>64 A STAR TYPE OF CONFIGURATION IS OBTAINED BY ADDING LARGER SPHERES WITH 32 POSITIONS ASSIGNED
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                                                                 *
c     *  ATTENTION!!! IF THE SOLVENT MOLECULES ARE TOO BIG USE THE ROUTINE GSTAR.f INSTEAD OF GPOST.f   *
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
c     DEFINE VECTORS FOR SIX POINTS ON THE FACES OF A CUBE OF SIDE LENGHT DIST, TWENTY SIX POINTS ON THE 
c     SURFACE OF A SPHERE OF RADIUS 2*DIST AND THE ADDITIONAL POINT AT (0.0,0.0,0.0)               
c     ......................................................... 
c     FIRST SPHERE RADIUS        = 0 >>>>>>>>>>>>>> 1 MOLECULE   
c     .........................................................
c     SECOND SPHERE/CUBE RADIUS  = DIST >>>>>>>>>>> 6 MOLECULES
c     .........................................................
c     THIRD SPHERE RADIUS        = 2*DIST >>>>>>>> 26 MOLECULES
c     FOURTH SPHERE RADIUS       = 3*DIST >>>>>>>> 32 MOLECULES
c     FIFTH SPHERE RADIUS        = 4*DIST >>>>>>>> 32 MOLECULES
c     SIXTH SPHERE RADIUS        = 5*DIST >>>>>>>> 32 MOLECULES
c     etc,...                                                         
c     .........................................................                  
c
      SQRARG  = 3.0D0/2.0D0 
      DISTWO  = 2.0D0*DIST
      DISTSRA = 2.0D0*DIST/DSQRT(3.0D0)
      DISTSRB = 2.0D0*DIST/DSQRT(2.0D0)
      DISTSRC = 3.0D0*DIST/DSQRT(SQRARG)
      DISTSRD = 3.0D0*DIST/(2.0D0*DSQRT(SQRARG))
      DISTSRE = 3.0D0*DIST/DSQRT(3.0D0)
c     -------------------------------------------------------  
c     Point in the centre of coordinates - sphere of radius 0
c     -------------------------------------------------------
      PFCCX(0)  = 0.0D0  
      PFCCY(0)  = 0.0D0 
      PFCCZ(0)  = 0.0D0
c     -------------------------------------
c     Six points on a sphere of radius DIST     
c     -------------------------------------      
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
c     --------------------------------------- 
c     Six points on a sphere of radius 2*DIST     
c     ---------------------------------------      
      PFCCX(7)  =  DISTWO  
      PFCCY(7)  = 0.0D0 
      PFCCZ(7)  = 0.0D0
c     -------------------      
      PFCCX(8)  = -DISTWO  
      PFCCY(8)  = 0.0D0 
      PFCCZ(8)  = 0.0D0
c     -------------------     
      PFCCX(9)  = 0.0D0  
      PFCCY(9)  =  DISTWO 
      PFCCZ(9)  = 0.0D0
c     ------------------- 
      PFCCX(10) = 0.0D0  
      PFCCY(10) = -DISTWO 
      PFCCZ(10) = 0.0D0
c     ------------------- 
      PFCCX(11) = 0.0D0  
      PFCCY(11) = 0.0D0 
      PFCCZ(11) =  DISTWO
c     -------------------  
      PFCCX(12) = 0.0D0  
      PFCCY(12) = 0.0D0 
      PFCCZ(12) = -DISTWO
c     -----------------------------------------   
c     Eight points on a sphere of radius 2*DIST
c     -----------------------------------------      
      PFCCX(13) =  DISTSRA  
      PFCCY(13) =  DISTSRA 
      PFCCZ(13) =  DISTSRA
c     --------------------    
      PFCCX(14) = -DISTSRA  
      PFCCY(14) =  DISTSRA 
      PFCCZ(14) =  DISTSRA
c     --------------------  
      PFCCX(15) = -DISTSRA  
      PFCCY(15) =  DISTSRA 
      PFCCZ(15) = -DISTSRA
c     --------------------
      PFCCX(16) =  DISTSRA  
      PFCCY(16) =  DISTSRA 
      PFCCZ(16) = -DISTSRA
c     --------------------
      PFCCX(17) =  DISTSRA  
      PFCCY(17) = -DISTSRA 
      PFCCZ(17) =  DISTSRA
c     -------------------- 
      PFCCX(18) = -DISTSRA  
      PFCCY(18) = -DISTSRA 
      PFCCZ(18) =  DISTSRA
c     --------------------   
      PFCCX(19) = -DISTSRA  
      PFCCY(19) = -DISTSRA 
      PFCCZ(19) = -DISTSRA
c     -------------------- 
      PFCCX(20) =  DISTSRA  
      PFCCY(20) = -DISTSRA 
      PFCCZ(20) = -DISTSRA
c     ------------------------------------------
c     Twelve points on a sphere of radius 2*DIST
c     ------------------------------------------      
      PFCCX(21) =  DISTSRB  
      PFCCY(21) =  DISTSRB 
      PFCCZ(21) =  0.0D0
c     --------------------    
      PFCCX(22) =  0.0D0  
      PFCCY(22) =  DISTSRB 
      PFCCZ(22) =  DISTSRB
c     -------------------- 
      PFCCX(23) = -DISTSRB  
      PFCCY(23) =  DISTSRB 
      PFCCZ(23) =  0.0D0
c     --------------------
      PFCCX(24) =  0.0D0  
      PFCCY(24) =  DISTSRB 
      PFCCZ(24) = -DISTSRB
c     --------------------
      PFCCX(25) =  DISTSRB  
      PFCCY(25) =  0.0D0 
      PFCCZ(25) =  DISTSRB
c     -------------------- 
      PFCCX(26) = -DISTSRB  
      PFCCY(26) =  0.0D0 
      PFCCZ(26) =  DISTSRB
c     --------------------   
      PFCCX(27) = -DISTSRB  
      PFCCY(27) =  0.0D0 
      PFCCZ(27) = -DISTSRB
c     -------------------- 
      PFCCX(28) =  DISTSRB  
      PFCCY(28) =  0.0D0 
      PFCCZ(28) = -DISTSRB
c     --------------------
      PFCCX(29) =  DISTSRB  
      PFCCY(29) = -DISTSRB 
      PFCCZ(29) =  0.0D0
c     --------------------   
      PFCCX(30) =  0.0D0  
      PFCCY(30) = -DISTSRB 
      PFCCZ(30) =  DISTSRB
c     -------------------- 
      PFCCX(31) = -DISTSRB  
      PFCCY(31) = -DISTSRB 
      PFCCZ(31) =  0.0D0
c     --------------------
      PFCCX(32) =  0.0D0  
      PFCCY(32) = -DISTSRB 
      PFCCZ(32) = -DISTSRB
c     -----------------------------------------------
c     Twenty four points on a sphere of radius 3*DIST 
c     -----------------------------------------------     
      PFCCX(33) =  DISTSRC  
      PFCCY(33) =  DISTSRD 
      PFCCZ(33) = -DISTSRD
c     --------------------    
      PFCCX(34) =  DISTSRC  
      PFCCY(34) =  DISTSRD 
      PFCCZ(34) =  DISTSRD
c     --------------------  
      PFCCX(35) =  DISTSRC  
      PFCCY(35) = -DISTSRD 
      PFCCZ(35) = -DISTSRD
c     --------------------
      PFCCX(36) =  DISTSRC  
      PFCCY(36) = -DISTSRD 
      PFCCZ(36) =  DISTSRD
c     --------------------
      PFCCX(37) =  DISTSRD  
      PFCCY(37) =  DISTSRD 
      PFCCZ(37) =  DISTSRC
c     -------------------- 
      PFCCX(38) = -DISTSRD  
      PFCCY(38) =  DISTSRD 
      PFCCZ(38) =  DISTSRC
c     --------------------   
      PFCCX(39) =  DISTSRD  
      PFCCY(39) = -DISTSRD 
      PFCCZ(39) =  DISTSRC
c     -------------------- 
      PFCCX(40) = -DISTSRD  
      PFCCY(40) = -DISTSRD 
      PFCCZ(40) =  DISTSRC
c     --------------------     
      PFCCX(41) = -DISTSRC
      PFCCY(41) =  DISTSRD 
      PFCCZ(41) =  DISTSRD
c     --------------------    
      PFCCX(42) = -DISTSRC  
      PFCCY(42) =  DISTSRD 
      PFCCZ(42) = -DISTSRD
c     --------------------  
      PFCCX(43) = -DISTSRC  
      PFCCY(43) = -DISTSRD 
      PFCCZ(43) =  DISTSRD
c     --------------------
      PFCCX(44) = -DISTSRC  
      PFCCY(44) = -DISTSRD 
      PFCCZ(44) = -DISTSRD
c     --------------------
      PFCCX(45) = -DISTSRD  
      PFCCY(45) =  DISTSRD 
      PFCCZ(45) = -DISTSRC
c     -------------------- 
      PFCCX(46) =  DISTSRD  
      PFCCY(46) =  DISTSRD 
      PFCCZ(46) = -DISTSRC
c     --------------------   
      PFCCX(47) = -DISTSRD  
      PFCCY(47) = -DISTSRD 
      PFCCZ(47) = -DISTSRC
c     -------------------- 
      PFCCX(48) =  DISTSRD  
      PFCCY(48) = -DISTSRD 
      PFCCZ(48) = -DISTSRC
c     --------------------    
      PFCCX(49) =  DISTSRD  
      PFCCY(49) =  DISTSRC 
      PFCCZ(49) =  DISTSRD
c     --------------------    
      PFCCX(50) = -DISTSRD  
      PFCCY(50) =  DISTSRC 
      PFCCZ(50) =  DISTSRD
c     --------------------  
      PFCCX(51) =  DISTSRD  
      PFCCY(51) =  DISTSRC 
      PFCCZ(51) = -DISTSRD
c     --------------------
      PFCCX(52) = -DISTSRD  
      PFCCY(52) =  DISTSRC 
      PFCCZ(52) = -DISTSRD
c     --------------------
      PFCCX(53) =  DISTSRD  
      PFCCY(53) = -DISTSRC 
      PFCCZ(53) =  DISTSRD
c     -------------------- 
      PFCCX(54) = -DISTSRD  
      PFCCY(54) = -DISTSRC 
      PFCCZ(54) =  DISTSRD
c     --------------------   
      PFCCX(55) =  DISTSRD  
      PFCCY(55) = -DISTSRC 
      PFCCZ(55) = -DISTSRD
c     -------------------- 
      PFCCX(56) = -DISTSRD  
      PFCCY(56) = -DISTSRC 
      PFCCZ(56) = -DISTSRD
c     -----------------------------------------      
c     Eight points on a sphere of radius 3*DIST
c     -----------------------------------------      
      PFCCX(57) =  DISTSRE  
      PFCCY(57) =  DISTSRE 
      PFCCZ(57) =  DISTSRE
c     --------------------    
      PFCCX(58) = -DISTSRE  
      PFCCY(58) =  DISTSRE 
      PFCCZ(58) =  DISTSRE
c     --------------------  
      PFCCX(59) = -DISTSRE  
      PFCCY(59) =  DISTSRE 
      PFCCZ(59) = -DISTSRE
c     --------------------
      PFCCX(60) =  DISTSRE  
      PFCCY(60) =  DISTSRE 
      PFCCZ(60) = -DISTSRE
c     --------------------
      PFCCX(61) =  DISTSRE  
      PFCCY(61) = -DISTSRE 
      PFCCZ(61) =  DISTSRE
c     -------------------- 
      PFCCX(62) = -DISTSRE  
      PFCCY(62) = -DISTSRE 
      PFCCZ(62) =  DISTSRE
c     --------------------   
      PFCCX(63) = -DISTSRE  
      PFCCY(63) = -DISTSRE 
      PFCCZ(63) = -DISTSRE
c     -------------------- 
      PFCCX(64) =  DISTSRE  
      PFCCY(64) = -DISTSRE 
      PFCCZ(64) = -DISTSRE
c     --------------------
c
      IF (NMOLCS.LE.32)THEN
         NSPHER  = 3
         RADII(1) = 0.0D0
         RADII(2) = DIST*BANG
         RADII(3) = DISTWO*BANG
         WRITE(*,*)'                                                  '
         WRITE(*,'(A63)')'Three concentric spheres are enough to assign 
     &initial positions'
         WRITE(*,'(A23,3(F9.3))')'Spheres radii are (Ang)',RADII(1),
     &RADII(2),RADII(3)
         WRITE(*,*)'                                                  '
      ELSEIF (NMOLCS.LE.64)THEN
         NSPHER  = 4
         RADII(1) = 0.0D0
         RADII(2) = DIST*BANG
         RADII(3) = DISTWO*BANG
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
         RADII(3) = DISTWO*BANG
         RADII(4) = 3.0D0*DIST*BANG
         NLATTP   = IDNINT(DFLOAT(NMOLCS-32)/DFLOAT(32)+0.49D0)
         NSPHER   = 3+NLATTP
c...Define additional positions on concentric spheres of larger radii
c   This loop creates in general a much larger number of positions than those needed, i.e., than the number of molecules
c        ............................................................
         KL=64
         DO 5 IL = 2,NLATTP
            DO 15 JL = 1,32
               PFCCX(JL+KL) = PFCCX(JL+32)*DFLOAT(IL+2)/3.0D0
               PFCCY(JL+KL) = PFCCY(JL+32)*DFLOAT(IL+2)/3.0D0
               PFCCZ(JL+KL) = PFCCZ(JL+32)*DFLOAT(IL+2)/3.0D0
               RADII(IL+3)  = DFLOAT(IL+2)*DIST*BANG
   15       CONTINUE
            KL=KL+32   
    5    CONTINUE
c        ............................................................       
c   
c...Volume of the smallest sphere that still contains the system 
         RADMAX  = DIST*BANG*DFLOAT(NLATTP+2)
         VSPHMAX = (4.0D0/3.0D0)*PI*(RADMAX)**3.0D0
c 
         WRITE(*,*)'                                                  '
         WRITE(*,'(A36,I3,A20)')'Positions assigned on the surface of',
     &NSPHER,'  concentric spheres'
         WRITE(*,'(A67,F15.3)')'Radii of the smallest sphere that 
     &still contains the system (A)  = ',RADMAX
         WRITE(*,'(A67,E15.3)')'Volume of the smallest sphere that 
     &still contains the system (A³) = ',VSPHMAX
         WRITE(*,*)'                                                  '
      ENDIF 
c...Print out to MDOUTPUT.DAT the radii of the spheres used to assign initial positions
         DO 25 N=1,NSPHER
            WRITE(6,'(A8,I2,A14,F11.3)')'Sphere ',N,' radius (A) = '
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
      LRAND=.TRUE.
      IF(LRAND)CALL MDROT
c      
*      IF(LRAND)THEN
*         DO 115 I = 1,NCOMP
*            IF((NMOL(I).GT.1).AND.(NATMOL(I).GT.1))CALL MDROT
*  115    CONTINUE       
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