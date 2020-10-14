      SUBROUTINE DINTF2(N,N1,X,Y,X1,Y1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION  X(N),Y(N),X1(N1),Y1(N1)
      DIMENSION  IOP(2),W(400),A(400),B(400),C(400),TAB(3)
      IJ=1
      IOP(1)=5
      IOP(2)=5
      CALL SPLID1(N,X,Y,W,IOP,IJ,A,B,C)
      DO 1 I=1,N1
      X11=X1(I)
      CALL SPLID2(N,X,Y,W,IJ,X11,TAB)
      Y1(I)=TAB(1)
  1   CONTINUE
      RETURN
      END SUBROUTINE
C
      SUBROUTINE SPLID1(N,X,F,W,IOP,IJ,A,B,C)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     WHERE N=NUMBER OF POINTS IN THE INTERPOLATION
C           X=ORIGIN OF TABLE OF INDEPENDENT VARIABLE
C           F=ORIGIN OF TABLE OF DEPENDENT VARIABLE
C           W=AN ARRAY OF DIMENSION N WHICH CONTAINS THE CALCULATED
C             SECOND DERIVATIVES UPON RETURN
C           IOP=AN ARRAY OF DEMENSION 2 WHICH CONTAINS COMBINATIONS OF
C               THE INTEGERS 1 THRU 5 USED TO SPECIFY THE BOUNDARY
C               CONDITIONS
C           IJ=SPACING IN THE F AND W TABLES
C           A,B,C=ARRAYS OF DIMENSION N USED FOR TEMPORARY STORAGE
C
      LOGICAL MMMM
      DIMENSION IOP(2),X(N),F(N),W(N),A(N),B(N),C(N)
      K=N-1
      A(2)=-(X(2)-X(1))/6.
      B(2)=(X(3)-X(1))/3.
      W(IJ+1)=(F(2*IJ+1)-F(IJ+1))/(X(3)-X(2))-(F(IJ+1)-F(1))/(X(2)-X(1))
      IF (N-3)3,4,3
    3 DO 10 I=3,K
      M=(I-1)*IJ+1
      J1=M+IJ
      J2=M-IJ
      CON=(X(I+1)-X(I-1))/3.
      DON=(X(I)-X(I-1))/6.
      B(I)=CON-(DON**2)/B(I-1)
      E=(F(J1)-F(M))/(X(I+1)-X(I))-(F(M)-F(J2))/(X(I)-X(I-1))
      W(M)=E-(DON*W(J2))/B(I-1)
   10 A(I)=-(DON*A(I-1))/B(I-1)
    4 K1=(N-2)*IJ+1
      C(N-1)=-((X(N)-X(N-1))/6.)/B(N-1)
      W(K1)=W(K1)/B(N-1)
      A(N-1)=A(N-1)/B(N-1)
      K2=K-1
      IF (N-3)7,8,7
    7 DO 20 I=2,K2
      J=N-I
      CON=(X(J+1)-X(J))/6.
      A(J)=(A(J)-CON*A(J+1))/B(J)
      C(J)=-(CON*C(J+1))/B(J)
      K3=(J-1)*IJ+1
      M=K3+IJ
   20 W(K3)=(W(K3)-CON*W(M))/B(J)
    8 K4=(N-1)*IJ+1
      IF (IOP(1)-5) 201,200,201
  201 C1=W(1)
      IF (IOP(2)-5) 203,202,203
  203 C2=W(K4)
      GO TO 205
  200 IF (N-4)300,302,302
  302 A1=X(1)-X(2)
      A2=X(1)-X(3)
      A3=X(1)-X(4)
      A4=X(2)-X(3)
      A5=X(2)-X(4)
      A6=X(3)-X(4)
      W(1)=F(1)*(1./A1+1./A2+1./A3)-A2*A3*F(IJ+1)/(A1*A4*A5)+
     1A1*A3*F(2*IJ+1)/(A2*A4*A6)-A1*A2*F(3*IJ+1)/(A3*A5*A6)
      GO TO 201
  202 IF (N-4)300,303,303
  303 B1=X(N)-X(N-3)
      B2=X(N)-X(N-2)
      B3=X(N)-X(N-1)
      B4=X(N-1)-X(N-3)
      B5=X(N-1)-X(N-2)
      B6=X(N-2)-X(N-3)
      L1=K4-IJ
      L2=L1-IJ
      L3=L2-IJ
      W(K4)=-B2*B3*F(L3)/(B6*B4*B1)+B1*B3*F(L2)/(B6*B5*B2)
     1 -B1*B2*F(L1)/(B4*B5*B3)+F(K4)*(1./B1+1./B2+1./B3)
      GO TO 203
  205 DO 50 I=1,K
      M=(I-1)*IJ+1
      GO TO 60
   70 IF (I-1)80,50,80
   80 W(1)=W(1)-BOB*W(M)
      W(K4)=W(K4)-BILL*W(M)
      A(1)=A(1)-BOB*A(I)
      A(N)=A(N)-BILL*A(I)
      C(1)=C(1)-BOB*C(I)
      C(N)=C(N)-BILL*C(I)
      GO TO 50
   60 MK=IOP(1)
      GO TO (62,64,66,68,66),MK
   62 IF (I-1)71,63,71
   63 A(1)=-1.
      C(1)=0.
      GO TO 500
   71 BOB=0.
      GO TO 500
   64 IF (I-1)73,76,73
   76 A(1)=-1.
      C(1)=0.
      W(1)=0.
      GO TO 500
   73 IF (I-2)81,81,82
   81 BOB=-C1
      GO TO 500
   82 BOB=0.
      GO TO 500
   66 IF (I-1)83,84,83
   84 A(1)=-(X(2)-X(1))/3.
      C(1)=0.
      W(1)=-C1+(F(IJ+1)-F(1))/(X(2)-X(1))
      GO TO 500
   83 IF (I-2)85,85,86
   85 BOB=(X(2)-X(1))/6.
      GO TO 500
   86 BOB=0.
      GO TO 500
   68 IF (I-1)87,88,87
   88 A(1)=-1.
      C(1)=1.
      W(1)=0.
      GO TO 500
   87 BOB=0.
  500 ML=IOP(2)
      GO TO (120,130,140,150,140),ML
  120 IF (I-1)121,122,121
  122 A(N)=0.
      C(N)=-1.
      GO TO 70
  121 BILL=0.
      GO TO 70
  130 IF (I-1)131,132,131
  132 A(N)=0.
      C(N)=-1.
      W(K4)=0.
      GO TO 70
  131 IF (I-K)134,133,134
  133 BILL=-C2
      GO TO 70
  134 BILL=0.
      GO TO 70
  140 IF (I-1)141,142,141
  142 A(N)=0.
      C(N)=(X(N-1)-X(N))/3.
      W(K4)=C2-(F(K4)-F(K1))/(X(N)-X(N-1))
      GO TO 70
  141 IF (I-K)143,144,143
  144 BILL=(X(N)-X(N-1))/6.
      GO TO 70
  143 BILL=0.
      GO TO 70
  150 IF (I-1)151,152,151
  152 A(N)=0.
      C(N)=(X(N-1)+X(1)-X(N)-X(2))/3.
      W(K4)=(F(IJ+1)-F(1))/(X(2)-X(1))-(F(K4)-F(K1))/(X(N)-X(N-1))
      GO TO 70
  151 IF (I-2)153,154,153
  154 BILL=(X(2)-X(1))/6.
      GO TO 70
  153 IF (I-K)155,156,155
  156 BILL=(X(N)-X(N-1))/6.
      GO TO 70
  155 BILL=0.
      GO TO 70
   50 CONTINUE
      CON=A(1)*C(N)-C(1)*A(N)
      D1=-W(1)
      D2=-W(K4)
      W(1)=(D1*C(N)-C(1)*D2)/CON
      W(K4)=(A(1)*D2-D1*A(N))/CON
      DO 110 I=2,K
      M=(I-1)*IJ+1
  110 W(M)=W(M)+A(I)*W(1)+C(I)*W(K4)
      GO TO 305
  300 WRITE(6,333) MMMM
  333 FORMAT(1X,40HSPLID1 N LESS THAN 4. RESULTS INCORRECT.)
  305 RETURN
      END SUBROUTINE

      SUBROUTINE SPLID2(N,X,F,W,IJ,Y,TAB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     WHERE N=NUMBER OF POINTS IN THE INTERPOLATION
C           X=ORIGIN OF TABLE OF THE INDEPENDENT VARIABLE
C           F=ORIGIN OF TABLE OF THE DEPENDENT VARIABLE
C           W=ORIGIN OF TABLE OF SECOND DERIVATIVES AS CALCULATED
C             BY SPLID1
C           IJ=SPACING IN THE TABLES F AND W
C           Y=THE POINT AT WHICH INTERPOLATION IS DESIRED
C           TAB=AN ARRAY OF DIMENSION 3 WHICH CONTAINS THE FUNCTION
C               VALUE,FIRST DERIVATIVE,AND SECOND DERIVATIVE AT Y
C
      DIMENSION X(N),F(N),W(N),TAB(3)
C
C     LOCATE Y IN THE X TABLE
C
      IF (Y-X(1))10,10,20
   10 I=1
      GO TO 30
   20 IF (Y-X(N))15,40,40
   40 I=N-1
      GO TO 30
   15 CALL SEARCH(Y,X,N,I)
   30 MI=(I-1)*IJ+1
      KI=MI+IJ
      FLK=X(I+1)-X(I)
C
C     CALCULATE F(Y)
C
      A=(W(MI)*(X(I+1)-Y)**3+W(KI)*(Y-X(I))**3)/(6.*FLK)
      B=(F(KI)/FLK-W(KI)*FLK/6.)*(Y-X(I))
      C=(F(MI)/FLK-FLK*W(MI)/6.)*(X(I+1)-Y)
      TAB(1)=A+B+C
C
C     CALCULATE THE FIRST DERIVATIVE AT Y
C
      A=(W(KI)*(Y-X(I))**2-W(MI)*(X(I+1)-Y)**2)/(2.*FLK)
      B=(F(KI)-F(MI))/FLK
      C=FLK*(W(MI)-W(KI))/6.
      TAB(2)=A+B+C
C
C     CALCULATE THE SECOND DERIVATIVE AT Y
C
      TAB(3)=(W(MI)*(X(I+1)-Y)+W(KI)*(Y-X(I)))/FLK
      RETURN
      END SUBROUTINE

      SUBROUTINE SEARCH(Y,X,N,I)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(N)
      K=N-1
      DO 100 J=1,K
      IF (Y.GE.X(J).AND.Y.LT.X(J+1)) GO TO 101
      GO TO 100
  101 I=J
      GO TO 102
  100 CONTINUE
  102 RETURN
      END SUBROUTINE
