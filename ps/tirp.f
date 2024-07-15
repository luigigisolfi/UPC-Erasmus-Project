C*********************************************************************
C PROGRAM TIRP
C GIVEN A TABLE OF TIME, POS and VEL (IN ANY COORDINATES) STORED IN
C A FILE IT PERFORMS THE PARALLEL SHOOTING IN THE JPL SYSTEM USING 
C ROUTINE TPX (BUT CAN BE CHANGED TO ANY OTHER). 
C
C NOTE: SEE THE INPUT STATEMENTS FOR THE PRECISE FORMAT OF 
C       THE FILE (CHANNEL NCI).
C
C AT EACH ITERATION THE PROGRAM WRITES:
C          IN FILE ARXD (dibx.dat) TIME AND COORDINATES OF THE
C                       ORBIT FOR A PLOT. (IT IS WRITTEN BEFORE
C                       CALLING THE PARALLEL SHOOTING ROUTINE).
C          IN FILE ARXT (tpx.dat) EPOCHS AND POS+VEL IN JPL
C                       COORDINATES OF THE ACTUAL NODES OF THE 
C                       PARALLEL SHOOTING. IT IS WRITTEN BEFORE
C                       CALLING THE PARALLEL SHOOTING ROUTINE AND
C                       CAN BE USED AGAIN AS INPUT FILE OF THIS
C                       PROGRAM.
C          IN FILE tpx.itr THE PARALLEL SHOOTING ROUTINE WRITES
C                       INFORMATION ABOUT THE LAST ITERATION.
C
C NOTE: - THE VARIATIONAL MATRICES CAN BE WRITTEN VIA ROUTINE TPX.
C*********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C NP IS THE MAX NUMBER OF NODES OF THE PS.
      PARAMETER (NP=5000)
      CHARACTER*64 ARXI,ARXD,ARXT,ARXM
      DIMENSION VI(7*NP),FU(7*NP-7),PES(7*NP),CO(7*NP),IFI(7),XIN(9)
C
C VECTORFIELD CONTAINING ALSO THE VARIATIONAL EQUATIONS TO
C BE USED, SOME INPUT AND OUTPUT FILES AND JPL MODEL.
C
      EXTERNAL DERIVP
      ARXI='citp.dat'
      ARXD='dibx.dat'
      ARXT='tpx.dat'
      ARXM='model.dat'
      NCI=3
      NCR=4
      WRITE (*,*) 'JPL MODEL READ FROM FILE:',ARXM
      WRITE (*,*) 'IT IS IT SET UP (ICO, LI...) ?, YES=1'
      READ (*,*) ICON
      IF (ICON.NE.1) STOP
      OPEN(4,FILE=ARXM,STATUS='OLD')
      CALL MODJPL(4)
      CLOSE(4)
      CALL PODJE0(3,DJEZ)
      ORIG=0.D0
      WRITE (*,*) 'ORIGIN OF DAYS IN JULIAN EPOCHS: ',DJEZ
C
C THE INITIAL NODES ARE READ. THEN THEIR COORDINATES ARE 
C CHANGED INTO JPL ONES AND STORED IN VI READY FOR THE 
C PARALLEL SHOOTING ROUTINE.
C
      WRITE (*,*) 'GIVE THE FILE WHERE THE INITIAL NODES ARE STORED'
      WRITE (*,*) 'USUALLY: ',ARXI
      READ (*,*) ARXI
      OPEN(NCI,FILE=ARXI,STATUS='OLD')
C READING THE HEADER OF THE INPUT FILE.
C  * One header line for explanations
C  * T0 is the initial epoch (initial adim time or day)
C  * ITT, ITX the type of coordinates and NNO the number of nodes. 
      READ (NCI,*)
      READ (NCI,*) T0,XMU,ICO,LI,ITT,ITX,NNO ! molts no s'usen
      WRITE (*,*) 'THE FILE CONTAINS ',NNO,' NODES.'
      WRITE (*,*) 'GIVE INTIAL ONE, INPUT STEP FROM THE INITIAL'
      WRITE (*,*) 'ONE AND FINAL ONE '
      WRITE (*,*) '(i.e. 1 1 ',NNO,' reads all the file)'
      READ (*,*) INO,NSA,NFN
      IF (NFN.GT.NNO) NFN=NNO
      N=(NFN-INO)/NSA+1
      IF (N.GT.NP) THEN
      WRITE(*,*) 'MAIN. PARAMETER NP SMALL. NNO i NP: ',NNO,NP
      STOP
      ENDIF
      WRITE (*,*) 'IN THIS INPUT FILE ...'
      WRITE (*,*) 'TIME IS IN DAYS (1) or ADDIMENSIONAL (2) ?'
      WRITE (*,*) 'INDICATOR OF FILE SAYS: ',ITT
      READ (*,*) ITT
      WRITE (*,*) 'POS+VEL IS: ADIM-LI (1), ADIM-RTBP (2) or JPL (3) ?'
      WRITE (*,*) 'INDICATOR OF FILE SAYS: ',ITX
      READ (*,*) ITX
      WRITE (*,*) 'STORING INITIAL NODES IN SUITABLE COORDINATES'
C READING AND STORING THE DATA ABOUT THE NODES IN JPL COORDINATES
      N=-1
      DO 10 J=INO,NNO
      IF (J.GT.NFN) GOTO 10
      READ (NCI,*) NN,TJ,(XIN(K),K=1,6)
      IF (J.LT.INO.OR.JMOD(J-INO,NSA).NE.0) GOTO 10
      N=N+1
      TJ=T0+TJ
      VI(7*N+1)=TJ
      IF (ITT.EQ.2) CALL CDJTA(VI(7*N+1),TJ,ORIG,-1)
      CALL TRCABJ(VI(7*N+1),XIN,VI(7*N+2),ITX,3,2)
10    CONTINUE
      CLOSE(NCI)
      N=N+1
      WRITE (*,*) 'NUMBER OF NODES STORED: ',N
      WRITE (*,*) 'PARALLEL SHOOTING FROM, TO, AND NUMBER OF DAYS:'
      WRITE (*,*) VI(1),VI(7*(N-1)+1),VI(7*(N-1)+1)-VI(1)
C
C CONTROL PARAMETERS FOR THE PARALLEL SHOOTING AND OUTPUT.
C
      IPO=1
      IND=1 
      TOL=1.D-12
      WRITE (*,*) 'PARALLEL SHOOTING FORWARD (1) or BACKWARDS (-1)'
      WRITE (*,*) 'IN TIME (usually 1) ?'
      READ (*,*) IND
C      WRITE (*,*) ' IPO, IND, TOL ?'
C      READ (*,*) IPO,IND,TOL
      WRITE (*,*) 'MAXIMUM STEP IN DAYS FOR THE OUTPUT FILE WHICH'
      WRITE (*,*) 'CAN BE USED TO PLOT THE ORBIT ?'
      READ (*,*) HMAT
      WRITE (*,*) 'TYPE OF COORDINATES FOR THIS OUTPUT FILE ?'
      WRITE (*,*) '  ( 1=ADIM-LI; 2=BARIC-RTBP; 3=JPL )'
      READ (*,*) ITD
C PESOS I EXTREMS DEL TIR PARAL.LEL
      DO 1 I=0,N-1
      PES(7*I+1)=1.D0  ! Pes en component de t
      DO 2 K=2,4
      PES(7*I+K)=1.D0  ! Pes en pos 
      PES(7*I+K+3)=1.D0  ! Pes en vel
2     CONTINUE
1     CONTINUE
      IFI(1)=0
      IFI(2)=0
      IFI(3)=0
      IFI(4)=0
      IFI(5)=0
      IFI(6)=0
      IFI(7)=0
C
C ITERATIONS OF PARALLEL SHOOTING. WE DO A TABLE-PLOT
C BEFORE THE ITERATION USING TAUDPS, AND THEN WE CALL 
C THE MAIN ROUTINE OF PARALLEL SHOTING. THIS ROUTINE CAN
C BE SELECTED BETWEEN THE DIFFERENT POSSIBILITIES WE HAVE
C OF PARALLEL SHOOTING.
C
      ITER=0
100   ITER=ITER+1
      write (*,*) 'dibuixo trajectoria'
      CALL TAUDPS(N,VI,TOL,HMAT,ITD,NCR,ARXD,DERIVP)
      write (*,*) 'entro tpx'
      CALL TPX(N,VI,IPO,IND,IFI,PES,FU,CO,FNO,CNO,TOL,DERIVP)
      WRITE(*,*) ' ITERATION NUMBER ',ITER
C  WRITING THE NODES IN A FILE BEFORE DOING THE CORRECTION. 
      OPEN(NCR,FILE=ARXT)
      WRITE(NCR,*)'% OUTPUT NODES FROM SOME ITER OF TIRP OF PAR. SHOOT.'
      WRITE (NCR,200) T0,XMU,ICO,LI,1,3,N
200   FORMAT (2E24.16,4I3,I7)
      DO 101 I=1,N
      WRITE(NCR,201) I,VI(7*(I-1)+1)-T0,(VI(7*(I-1)+K),K=2,7)
201   FORMAT(I6,7E24.16)
101   CONTINUE
      CLOSE(NCR)
C  PERFORMING THE CORRECTION OF THE NODES
      DO 105 I=1,7*N
      VI(I)=VI(I)+CO(I)
105   CONTINUE
C  THERE IS NO STOP CONDITION IN THIS PROGRAM
      GOTO 100
      END