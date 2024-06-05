c     FORTRAN 77 FIXED FORMAT [7-72]
c
c     ***************************************************************************
c     **                                                                       **
c     **   PROGRAM: Molecular Dynamics of CLusters using Gaussian              **
c     **                                                                       **
c     **   FILENAME: MDCLG.f                                      Version 2.0  **
c     **                                                                       **
c     **   N. GALAMBA - November 2005                                          **
c     **                                                                       **
c     **   Modified November 2018                                              **
c     **                                                                       **
c     **   ALLOW RUN MD IN G09 -----> IN(14), NGAUSSV (03; 09)                 **
c     **                                                                       **
c     **   Modified February 2007                                              **
c     **                                                                       **
c     **   >>> ARRAY: NATMSP(NATMSPMAX) -----> NATMSP(NCOMPMAX,NATMSPMAX)      **
c     **                                                                       **
c     **   GRUPO DE FISICA-MATEMATICA DA UNIVERSIDADE DE LISBOA                **                                          
c     **                                                                       **                                           
c     ***************************************************************************
c     **                                                                       **
c     **   BORN-OPPENHEIMER MOLECULAR DYNAMICS SIMULATION OF CLUSTERS IN THE   **
c     **   MICROCANONICAL (N,V,E) ENSEMBLE FROM FIRST PRINCIPLES CALCULATION   **
c     **   OF THE ATOMIC FORCES USING LOCAL BASIS SETS WITH Gaussian 03        **
c     **                                                                       **
c     ***************************************************************************
c     **                                                                       **
c     **   PROGRAM DESCRIPTION                                                 **
c     **   -------------------                                                 **
c     **                                                                       **
c     **   A. SET UP INITIAL POSITIONS FOR ALL ATOMS                           **
c     **                                                                       ** 
c     **   B. DEFINE RANDOM VELOCITIES FOR ALL ATOMS AND SCALE VELOCITIES FOR  **
c     **      ZERO TOTAL LINEAR MOMENTUM AND A GIVEN TEMPERATURE               **
c     **                                                                       **
c     **   C. CALL GAUSSIAN TO CALCULATE Ab-Initio FORCES FOR SOME STARTING    **
c     **      CONFIGURATION                                                    **
c     **                                                                       **
c     **   D. READ FORCES FROM GAUSSIAN OUTPUT FILE AND INTEGRATE THE          **
c     **      EQUATIONS OF MOTION (Verlet "leap-frog" algorithm)               **
c     **                                                                       **
c     **   E. READ THE SCF ENERGY FROM GAUSSIAN OUTPUT TO CALCULATE THE TOTAL  **
c     **      ENERGY BY ADDING THE NUCLEAR KINETIC ENERGY                      **
c     **                                                                       **
c     **   F. SCALE VELOCITIES FOR DESIRED TEMPERATURE DURING EQUILIBRATION    **
c     **                                                                       **
c     **   G. WRITE TO FILE NEW POSITIONS (Gaussian FORMAT), NEW VELOCITIES    **
c     **      (FREE FORMAT) AND NEW POSITIONS (FREE FORMAT)                    **
c     **                                                                       **
c     **   H. CALCULATE NEW Ab-Initio FORCES BY CALLING GAUSSIAN AND READ OLD  **
c     **      POSITIONS AND VELOCITIES FROM FILE (FREE FORMAT) TO CONTINUE THE **
c     **      DYNAMICS                                                         **
c     **                                                                       **
c     ***************************************************************************
c     **                                                                       **
c     **   VARIABLES PASSED THROUGH THE ARRAY IN(...)                          **
c     **   ------------------------------------------                          **
c     **                                                                       **                                                             
c     **   IN(1) - NCOMP   - NUMBER OF CHEMICAL COMPONENTS                     **
c     **   IN(2) - NLG     - NUMBER OF LINES FOR GAUSSIAN OPTIONS              **
c     **   IN(3) - NATOMS  - TOTAL NUMBER OF ATOMS                             **  
c     **   IN(4) - NATOMS1 - NATOMS-1                                          **
c     **   IN(5) - MDSTEP  - EQUILIBRATION AND PRODUCTION MOLECULAR DYNAMICS   **
c     **                     STEP COUNTER - THIS COUNTER IS RESET TO ZERO AT   **
c     **                     END OF EQUILIBRATION                              **
c     **   IN(6) - MAXEQB  - NUMBER OF MD STEPS FOR EQUILIBRATION              ** 
c     **   IN(7) - NSTEP   - TOTAL no OF TIME-STEPS COUNTER - THIS COUNTER IS  **
c     **                     NOT RESET TO ZERO AT END OF EQUILIBRATION         **  
c     **   IN(8) - NMOLCS  - TOTAL NUMBER OF MOLECULES i.e., SUMMED OVER ALL   **  
c     **                     COMPONENTS                                        **
c     **   IN(9) - LTYPE   - FLAG TO DEFINE TYPE OF ORDERING OF COMPONENTS TO  **
c     **                     ASSIGN INITIAL POSITIONS IN GPOST AND TYPE OF     **
c     **                     GEOMETRICAL LATTICE i.e., SPHERICAL OR CUBIC FCC  **
c     **   IN(10)- NQUANTUM- FLAG TO DEFINE QUANTUM METHOD USED IN Gauss 03    **    
c     **   IN(11)- NGEN    - NUMBER OF OPTION LINES OF Gauss 03 TO SPECIFY     **
c     **                     ATOMS BASIS SETS IF GEN IS USED                   **
c     **   IN(12)- NGA     - NGA = NLG-NGEN                                    **
c     **   IN(13)- NGB     - NGB = NGA+1                                       **
c     **   IN(14)- NGAUSSV - GAUSSIAN VERION (03; 09)                          **
c     **                                                                       **
c     ***************************************************************************
c     **                                                                       **
c     **   VARIABLES PASSED THROUGH THE ARRAY LG(...)                          **
c     **   ------------------------------------------                          **
c     **                                                                       **                                                             
c     **   LG(1) - LPOST  - CONTROLS .TRUE./.FALSE. ON/OFF THE ROUTINE GPOST   **
c     **   LG(2) - LVEL   - CONTROLS .TRUE./.FALSE. ON/OFF THE ROUTINE MDVEL   **
c     **   LG(3) - LSCALE - CONTROLS .TRUE./.FALSE. ON/OFF VELOCITY SCALING    **
c     **                    FOR DESIRED TEMPERATURE - ON/OFF FOR EQUIL/PROD    **
c     **   LG(4) - LBOMD  - CONTROLS .TRUE./.FALSE. ON/OFF THE MD SIMULATION   **
c     **                    ALWAYS SET LBOMD .TRUE.                            **
c     **   LG(5) - LNBO   - CONTROLS .TRUE./.FALSE. ON/OFF NATURAL BOND        **
c     **                    ORBITAL ANALYSIS BY GAUSSIAN 03                    **
c     **                                                                       **
c     ***************************************************************************
c     **                                                                       **
c     **   VARIABLES PASSED THROUGH THE ARRAY RL(...)                          **
c     **   ------------------------------------------                          **
c     **                                                                       **                                                             
c     **   RL(1) - ENRNUC  - NUCLEAR REPULSION ENERGY CALCULATED BY GAUSSIAN   **
c     **   RL(2) - ENRSCF  - SCF ENERGY CALCULATED BY GAUSSIAN                 **
c     **   RL(3) - STEP    - TIME-STEP IN ATOMIC UNITS OF TIME                 **
c     **   RL(4) - STEPSQ  - SQUARE OF THE TIME STEP                           **
c     **   RL(5) - TKECONV - CONSTANT FACTOR TO CONVERT KINETIC ENERGY TO      **
c     **                     TEMPERATURE                                       ** 
c     **   RL(6) - ACTKE   - ACTUAL KINETIC ENERGY                             **
c     **   RL(7) - DTEMP   - DESIRED TEMPERATURE [K]                           **
c     **   RL(8) - ATEMP   - ACTUAL TEMPERATURE [K]                            **
c     **   RL(9) - ACTKEAV - AVERAGE ACTUAL KINETIC ENERGY ACCUMULATOR         **
c     **   RL(10)- ENSCFAV - AVERAGE SCF ENERGY ACCUMULATOR                    **
c     **   RL(11)- FNATOMS - NUMBER OF ATOMS CONVERTED TO DOUBLE PRECISION     **
c     **   RL(12)- FEKJM   - CONVERSION FACTOR FROM Ha TO kJ/mol atoms         **
c     **   RL(13)- COULV   - COULOMB PERMITTIVITY OF VACUUM TERM               **
c     **   RL(14)- FNMOLCS - NUMBER OF MOLECULES CONVERTED TO DOUBLE PRECISION ** 
c     **   RL(15)- DIST    - DISTANCE [Bohr] BETWEEN CONCENTRIC SPHERES WHERE  **
c     **                     TO INITIAL MOLECULAR POSITIONS ARE ASSIGNED       **
c     **   RL(16)- DTE     - TIME-STEP [fs] FOR EQUILIBRATION                  **
c     **   RL(17)- DTP     - TIME-STEP [fs] FOR PRODUCTION                     **
c     **   RL(18)- BOX     - LENGHT OF THE SIDE OF THE MD BOX IF LTYPE 2 OR 3  **
c     **                                                                       **
c     *************************************************************************** 
c     **                                                                       **
c     **   PARAMETERS                                                          **
c     **   ----------                                                          **
c     **                                                                       ** 
c     **   NCOMPMAX - MAXIMUM NUMBER OF COMPONENTS                             **
c     **   NMOLMAX  - MAXIMUM NUMBER OF MOLECULES/COMPONENT                    **
c     **   NSPECMAX - MAXIMUM NUMBER OF SPECIES/MOLECULE/COMPONENT             **
c     **   NATSPMAX - MAXIMUM NUMBER OF ATOMS/SPECIES/MOLECULE/COMPONENT       **
c     **   NLGAUSS  - MAXIMUM NUMBER OF GAUSSIAN OPTION LINES                  **
c     **   NPARTCL  - MAXIMUM NUMBER OF ATOMS                                  **
c     **   PTM      - MASS CONVERSION FACTOR FROM amu TO electron masses [me]  **
c     **   ATU      - TIME CONVERION FACTOR FROM atomic time units [atu] TO [s]**
c     **   THA      - TEMPERATURE CONVERSION FACTOR FROM Hartrees [Ha] TO [K]  **
c     **   BOLTZ    - BOLTZMANN CONSTANT IN [Ha/K]                             **
c     **   AVOGN    - AVOGADRO'S NUMBER                                        **
c     **   EHART    - ENERGY CONVERSION FACTOR FROM [Ha] TO [J]                **
c     **   PVACUM   - PERMITTIVITY OF VACUUM IN ATOMIC UNITS                   **
c     **   PI       - PI                                                       **
c     **   BANG     - VALUE OF ONE BOHR IN ANGSTROMS                           **
c     **   NRADMAX  - MAXIMUM NUMBER OF SPHERES RADII USED TO ASSIGN INITIAL   **
c     **              POSITIONS                                                **
c     **                                                                       **
c     ***************************************************************************
c     **                                                                       **
c     **   COMMON BLOCKS                                                       **
c     **   -------------                                                       **
c     **                                                                       ** 
c     **   POSIT  - POSITIONS AT STEP N                                        **
c     **   POSSD  - POSITIONS AT STEP N STORED IN A ONE-DIMENSION ARRAY        **
c     **   POSTD  - POSITIONS AT STEP N STORED IN A THREE-DIMENSION ARRAY      **
c     **   POST0  - ATOMIC POSITIONS INPUT FOR EACH MOLECULE OF EACH COMPONENT **
c     **   VELOC  - VELOCITIES AT STEP N-1/2                                   **   
c     **   VELOCN - VELOCITIES AT STEP N                                       **
c     **   FORCES - FORCES AT STEP N                                           **
c     **   PARMOL - MOLECULAR PARAMETERS                                       **
c     **   ATMASS - ATOMIC MASSES                                              **
c     **   ATNAME - ATOMIC CHEMICAL SYMBOLS                                    **
c     **   ATSYMB - ATOMIC CHEMICAL SYMBOLS                                    **
c     **   ATOMZ  - ATOMIC NUMBERS                                             **
c     **   ATNZSD - ATOMIC NUMBERS STORED IN A ONE-DIMENSION ARRAY             **
c     **   ATNUMB - ATOMS NUMBERING FROM 1 TO NATOMS                           **
c     **   ATCOMP - NUMBER OF ATOMS P/ MOLECULE FOR EACH COMPONENT OF THE      **
c     **            CLUSTER                                                    **
c     **   GAUSSL - GAUSSIAN 03 RUNNING OPTIONS                                **
c     **   MULLIK - MULLIKEN CHARGES                                           **
c     **   NBOCH  - NATURAL BOND ORBITAL ANALYSIS CHARGES                      **
c     **                                                                       ** 
c     ***************************************************************************
c     **                                                                       **
c     **   ATOMIC UNITS SYSTEM                                                 **
c     **   -------------------                                                 **
c     **                                                                       **
c     **   ALL CALCULATIONS INSIDE THIS CODE ARE CARRIED IN ATOMIC UNITS       **
c     **   THE INPUT OF THE CODE IS DEFINED IN ATOMIC UNITS EXCEPT FOR THE,    **
c     **                                                                       **
c     **   >>>>>>>> mass [atomic mass units]  e.g., m(Na) = 22.98977 [amu]     **
c     **   >>>>>>>> time-step [fs = 1.0E-15s] e.g., dt    = 5.0      [fs]      **
c     **   >>>>>>>> temperature [K]           e.g., T     = 300.0    [K]       **
c     **                                                                       **
c     **   THESE VARIABLES ARE CONVERTED TO ATOMIC UNITS INSIDE THE CODE       **
c     **                                                                       **
c     **   MASS   - ELECTRON MASSES                                            **
c     **            m[au] = m[amu]*(1.0E-3[Kg/amu]/Avog N)/me                  **
c     **            m[au] = m[amu]*(1.0E-3[Kg/amu]/6.022141990E+23)/me         **
c     **            m[au] = m[amu]*1.6605402E-27[Kg/amu]/9.10938188E-31[Kg]    **
c     **   >>>>>>>> m[au] = m[amu]*1822.88851 [electron masses/amu]            **
c     **                                                                       **
c     **   LENGHT - BOHRS                                                      **
c     **            Bohr = 0.5291772083E-10[m]                                 **
c     **                                                                       **
c     **   TIME   - ATOMIC TIME UNITS                                          **
c     **            atu = 2.418884E-17[s]                                      **
c     **   >>>>>>>> dt[atu] = dt[s]/2.418884E-17[s/atu]                        **
c     **                                                                       **
c     **   ENERGY - HARTREES                                                   **
c     *             Ha = 4.3597482E-18[J]                                      **
c     **                                                                       ** 
c     **   TEMP   - HARTREES                                                   **
c     **            Ha = 315773[K]                                             **
c     **            K  = Boltz k [J] = 1.38066E-23[J]                          **
c     **            K  = 1.38066E-23[J]/4.3597482E-18[J/Ha]                    **
c     **            K  = 3.166834268E-6[Ha]                                    **
c     **            Ha = 1.0[K]/3.166834268E-6 = 315773[K]                     **
c     **   >>>>>>>> T  = T[K]/315773[K/Ha]                                     **
c     **                                                                       **  
c     **   ATTENTION!!! This conversion is not used although the value of the  **  
c     **                desired T in [Ha] is printed out to MDOUTPUT.DAT file  **         
c     **                                                                       **
c     **   FORCE  - HARTREE/BOHR                                               **
c     **                                                                       **
c     **   VELOC  - BOHR/ATU                                                   **
c     **                                                                       **
c     ***************************************************************************
c          
      PROGRAM MDCLG
c          
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c          
      PARAMETER(NCOMPMAX=5)
      PARAMETER(NMOLMAX=256)
      PARAMETER(NSPECMAX=10)
      PARAMETER(NATSPMAX=100) 
      PARAMETER(NLGAUSS=500)
      PARAMETER(NPARTCL=500)     
c          
      CHARACTER*4  NAMEAT
      CHARACTER*150 NLINEG
c        
      LOGICAL LG(5),LPOST,LVEL,LSCALE,LBOMD,LNBO 
c                                               
      COMMON/INTGRS/IN(14)
      COMMON/LGCLS /LG      
      COMMON/REALS /RL(18)  
c         
      COMMON/POSIT /RX(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              RY(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              RZ(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)
      COMMON/POSSD /XX(NPARTCL),YY(NPARTCL),ZZ(NPARTCL) 
      COMMON/VELOC /VX(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              VY(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              VZ(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)      
      COMMON/FORCE /FX(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              FY(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX),
     &              FZ(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)         
      COMMON/PARMOL/NMOL(NMOLMAX),NSPEC(NSPECMAX),
     &NATMSP(NCOMPMAX,NATSPMAX)       
      COMMON/ATMASS/ZMASS(NCOMPMAX,NSPECMAX,NATSPMAX) 
      COMMON/ATNAME/NAMEAT(NCOMPMAX,NSPECMAX,NATSPMAX)
      COMMON/ATOMZ /NATZ(NCOMPMAX,NSPECMAX,NATSPMAX)
      COMMON/GAUSSL/NLINEG(NLGAUSS)
      COMMON/ATNUMB/NUMBAT(NCOMPMAX,NMOLMAX,NSPECMAX,NATSPMAX)    
c
      EQUIVALENCE (IN(1),NCOMP),(IN(2),NLG),(IN(3),NATOMS)
      EQUIVALENCE (IN(4),NATOMS1),(IN(5),MDSTEP),(IN(6),MAXEQB)
      EQUIVALENCE (IN(7),NSTEP),(IN(8),NMOLCS),(IN(9),LTYPE)
      EQUIVALENCE (IN(10),NQUANTUM),(IN(11),NGEN),(IN(12),NGA)
      EQUIVALENCE (IN(13),NGB),(IN(14),NGAUSSV)
      EQUIVALENCE (LG(1),LPOST),(LG(2),LVEL),(LG(3),LSCALE)
      EQUIVALENCE (LG(4),LBOMD),(LG(5),LNBO)
      EQUIVALENCE (RL(1),ENRNUC),(RL(2),ENRSCF),(RL(3),STEP)
      EQUIVALENCE (RL(4),STEPSQ),(RL(5),TKECONV),(RL(6),ACTKE)
      EQUIVALENCE (RL(7),DTEMP),(RL(8),ATEMP),(RL(9),ACTKEAV)
      EQUIVALENCE (RL(10),ENSCFAV),(RL(11),FNATOMS),(RL(12),FEKJM)
      EQUIVALENCE (RL(13),COULV),(RL(14),FNMOLCS),(RL(15),DIST)
      EQUIVALENCE (RL(16),DTE),(RL(17),DTP),(RL(18),BOX)         
c
c     OPEN I/O FILES FOR MD SIMULATION
c
c...  MD INPUT FILES
c
      OPEN(1,FILE='MDINPUT.DAT' ,ACCESS='SEQUENTIAL',STATUS='UNKNOWN')  
      OPEN(2,FILE='gaussf.com'  ,ACCESS='SEQUENTIAL',STATUS='UNKNOWN') 
      OPEN(3,FILE='MDPOST.DAT'  ,ACCESS='SEQUENTIAL',STATUS='UNKNOWN') 
      OPEN(4,FILE='MDVELT.DAT'  ,ACCESS='SEQUENTIAL',STATUS='UNKNOWN') 
c
c...  GAUSSIAN INPUT FILE/MD OUTPUT FILE
c       
      OPEN(5,FILE='gaussr.com'  ,ACCESS='SEQUENTIAL',STATUS='UNKNOWN')  
c
c...  MD OUTPUT FILES
c                
      OPEN(6,FILE='MDOUTPUT.DAT',ACCESS='APPEND'    ,STATUS='UNKNOWN')
      OPEN(7,FILE='mdcontrl.dat',ACCESS='SEQUENTIAL',STATUS='OLD')
      OPEN(8,FILE='MDVISUAL.PDB',ACCESS='APPEND'    ,STATUS='UNKNOWN') 
      OPEN(9,FILE='MDtrj.xyz',ACCESS='APPEND'    ,STATUS='UNKNOWN') 
c
c...  INITIALIZATION
c 
      CALL GPARM 
c    
      IF((LPOST).AND.(LTYPE.LE.1))CALL GPOST
      IF((LPOST).AND.(LTYPE.GE.2))CALL GPFCC
      IF(LVEL) CALL MDVEL
c
c     LPOST AND LVEL MUST BE .FALSE. DURING SIMULATION
c    
      LPOST = .FALSE.
      LVEL  = .FALSE.
c...  re-write part of the file mdcontrl.dat before the simulation starts     
      WRITE(7,*)LPOST
      WRITE(7,*)LVEL
      IF(.NOT.LBOMD)WRITE(7,*)0.0D0,0.0D0
c
      IF(LBOMD)THEN      
c
c     ATTENTION! RUN GAUSSIAN BEFORE EQUILIBRATION OR PRODUCTION STARTS
c    
*     g03<gaussr.com>gaussf.com 
c
c     ---------------------------------------
c     THIS CALL IS PERFORMED BY A BASH SCRIPT
c     ---------------------------------------
c  
c...  EQUILIBRATION/PRODUCTION
c          
      CALL GFORCE
*****      CALL GEWALD        
      CALL GVERLET
*****      CALL GPROPTY
      CALL MDVISUAL   
c      
      ELSE
      WRITE(*,'(A47)')'...............................................'
      WRITE(*,*)'                                                     '
      WRITE(*,'(A71)')'Finished setting up the initial positions and
     & velocities of the cluster'
      ENDIF  
c
      STOP
      END            
c
c
c      
      DOUBLE PRECISION FUNCTION RANF ( DUMMY )
c
c    *******************************************************************
c    **   RETURNS A UNIFORM RANDOM VARIATE IN THE RANGE 0 TO 1.       **
c    **                                                               **
c    **                 ***************                               **
c    **                 **  WARNING  **                               **
c    **                 ***************                               **
c    **                                                               **
c    **   GOOD RANDOM NUMBER GENERATORS ARE MACHINE SPECIFIC.         **
c    **   PLEASE USE THE ONE RECOMMENDED FOR YOUR MACHINE.            **
c    *******************************************************************
        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c        
      INTEGER     L, C, M
      PARAMETER ( L = 1029, C = 221591, M = 1048576 )
c
      INTEGER     SEED
      SAVE        SEED
      DATA        SEED / 0 /
      SEED = MOD ( SEED * L + C, M )
      RANF = REAL ( SEED ) / M
c
      RETURN
      END                                                                  
c
c
c