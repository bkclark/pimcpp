C     ALGORITHM 600, COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.9, NO. 2,
C     JUN., 1983, P. 258-259.
      SUBROUTINE QUINAT(N, X, Y, B, C, D, E, F)                         QUI   10
C
      INTEGER N
      REAL*8 X(N), Y(N), B(N), C(N), D(N), E(N), F(N)
C
C
C
C     QUINAT COMPUTES THE COEFFICIENTS OF A QUINTIC NATURAL QUINTIC SPLI
C     S(X) WITH KNOTS X(I) INTERPOLATING THERE TO GIVEN FUNCTION VALUES:
C               S(X(I)) = Y(I)  FOR I = 1,2, ..., N.
C     IN EACH INTERVAL (X(I),X(I+1)) THE SPLINE FUNCTION S(XX) IS A
C     POLYNOMIAL OF FIFTH DEGREE:
C     S(XX) = ((((F(I)*P+E(I))*P+D(I))*P+C(I))*P+B(I))*P+Y(I)    (*)
C           = ((((-F(I)*Q+E(I+1))*Q-D(I+1))*Q+C(I+1))*Q-B(I+1))*Q+Y(I+1)
C     WHERE  P = XX - X(I)  AND  Q = X(I+1) - XX.
C     (NOTE THE FIRST SUBSCRIPT IN THE SECOND EXPRESSION.)
C     THE DIFFERENT POLYNOMIALS ARE PIECED TOGETHER SO THAT S(X) AND
C     ITS DERIVATIVES UP TO S"" ARE CONTINUOUS.
C
C        INPUT:
C
C     N          NUMBER OF DATA POINTS, (AT LEAST THREE, I.E. N > 2)
C     X(1:N)     THE STRICTLY INCREASING OR DECREASING SEQUENCE OF
C                KNOTS.  THE SPACING MUST BE SUCH THAT THE FIFTH POWER
C                OF X(I+1) - X(I) CAN BE FORMED WITHOUT OVERFLOW OR
C                UNDERFLOW OF EXPONENTS.
C     Y(1:N)     THE PRESCRIBED FUNCTION VALUES AT THE KNOTS.
C
C        OUTPUT:
C
C     B,C,D,E,F  THE COMPUTED SPLINE COEFFICIENTS AS IN (*).
C         (1:N)  SPECIFICALLY
C                B(I) = S'(X(I)), C(I) = S"(X(I))/2, D(I) = S"'(X(I))/6,
C                E(I) = S""(X(I))/24,  F(I) = S""'(X(I))/120.
C                F(N) IS NEITHER USED NOR ALTERED.  THE FIVE ARRAYS
C                B,C,D,E,F MUST ALWAYS BE DISTINCT.
C
C        OPTION:
C
C     IT IS POSSIBLE TO SPECIFY VALUES FOR THE FIRST AND SECOND
C     DERIVATIVES OF THE SPLINE FUNCTION AT ARBITRARILY MANY KNOTS.
C     THIS IS DONE BY RELAXING THE REQUIREMENT THAT THE SEQUENCE OF
C     KNOTS BE STRICTLY INCREASING OR DECREASING.  SPECIFICALLY:
C
C     IF X(J) = X(J+1) THEN S(X(J)) = Y(J) AND S'(X(J)) = Y(J+1),
C     IF X(J) = X(J+1) = X(J+2) THEN IN ADDITION S"(X(J)) = Y(J+2).
C
C     NOTE THAT S""(X) IS DISCONTINUOUS AT A DOUBLE KNOT AND, IN
C     ADDITION, S"'(X) IS DISCONTINUOUS AT A TRIPLE KNOT.  THE
C     SUBROUTINE ASSIGNS Y(I) TO Y(I+1) IN THESE CASES AND ALSO TO
C     Y(I+2) AT A TRIPLE KNOT.  THE REPRESENTATION (*) REMAINS
C     VALID IN EACH OPEN INTERVAL (X(I),X(I+1)).  AT A DOUBLE KNOT,
C     X(J) = X(J+1), THE OUTPUT COEFFICIENTS HAVE THE FOLLOWING VALUES:
C       Y(J) = S(X(J))          = Y(J+1)
C       B(J) = S'(X(J))         = B(J+1)
C       C(J) = S"(X(J))/2       = C(J+1)
C       D(J) = S"'(X(J))/6      = D(J+1)
C       E(J) = S""(X(J)-0)/24     E(J+1) = S""(X(J)+0)/24
C       F(J) = S""'(X(J)-0)/120   F(J+1) = S""'(X(J)+0)/120
C     AT A TRIPLE KNOT, X(J) = X(J+1) = X(J+2), THE OUTPUT
C     COEFFICIENTS HAVE THE FOLLOWING VALUES:
C       Y(J) = S(X(J))         = Y(J+1)    = Y(J+2)
C       B(J) = S'(X(J))        = B(J+1)    = B(J+2)
C       C(J) = S"(X(J))/2      = C(J+1)    = C(J+2)
C       D(J) = S"'((X(J)-0)/6    D(J+1) = 0  D(J+2) = S"'(X(J)+0)/6
C       E(J) = S""(X(J)-0)/24    E(J+1) = 0  E(J+2) = S""(X(J)+0)/24
C       F(J) = S""'(X(J)-0)/120  F(J+1) = 0  F(J+2) = S""'(X(J)+0)/120
C
      INTEGER I, M
      REAL*8 B1, P, PQ, PQQR, PR, P2, P3, Q, QR, Q2, Q3, R, R2, S,T,U,V
C
      IF (N.LE.2) GO TO 190
C
C     COEFFICIENTS OF A POSITIVE DEFINITE, PENTADIAGONAL MATRIX,
C     STORED IN D,E,F FROM 2 TO N-2.
C
      M = N - 2
      Q = X(2) - X(1)
      R = X(3) - X(2)
      Q2 = Q*Q
      R2 = R*R
      QR = Q + R
      D(1) = 0.
      E(1) = 0.
      D(2) = 0.
      IF (Q.NE.0.) D(2) = 6.*Q*Q2/(QR*QR)
C
      IF (M.LT.2) GO TO 40
      DO 30 I=2,M
        P = Q
        Q = R
        R = X(I+2) - X(I+1)
        P2 = Q2
        Q2 = R2
        R2 = R*R
        PQ = QR
        QR = Q + R
        IF (Q) 20, 10, 20
   10   D(I+1) = 0.
        E(I) = 0.
        F(I-1) = 0.
        GO TO 30
   20   Q3 = Q2*Q
        PR = P*R
        PQQR = PQ*QR
        D(I+1) = 6.*Q3/(QR*QR)
        D(I) = D(I) + (Q+Q)*(15.*PR*PR+(P+R)*Q*(20.*PR+7.*Q2)+Q2*(8.*
     *   (P2+R2)+21.*PR+Q2+Q2))/(PQQR*PQQR)
        D(I-1) = D(I-1) + 6.*Q3/(PQ*PQ)
        E(I) = Q2*(P*QR+3.*PQ*(QR+R+R))/(PQQR*QR)
        E(I-1) = E(I-1) + Q2*(R*PQ+3.*QR*(PQ+P+P))/(PQQR*PQ)
        F(I-1) = Q3/PQQR
   30 CONTINUE
C
   40 IF (R.NE.0.) D(M) = D(M) + 6.*R*R2/(QR*QR)
C
C     FIRST AND SECOND ORDER DIVIDED DIFFERENCES OF THE GIVEN FUNCTION
C     VALUES, STORED IN B FROM 2 TO N AND IN C FROM 3 TO N
C     RESPECTIVELY. CARE IS TAKEN OF DOUBLE AND TRIPLE KNOTS.
C
      DO 60 I=2,N
        IF (X(I).NE.X(I-1)) GO TO 50
        B(I) = Y(I)
        Y(I) = Y(I-1)
        GO TO 60
   50   B(I) = (Y(I)-Y(I-1))/(X(I)-X(I-1))
   60 CONTINUE
      DO 80 I=3,N
        IF (X(I).NE.X(I-2)) GO TO 70
        C(I) = B(I)*0.5
        B(I) = B(I-1)
        GO TO 80
   70   C(I) = (B(I)-B(I-1))/(X(I)-X(I-2))
   80 CONTINUE
C
C     SOLVE THE LINEAR SYSTEM WITH C(I+2) - C(I+1) AS RIGHT-HAND SIDE.
C
      IF (M.LT.2) GO TO 100
      P = 0.
      C(1) = 0.
      E(M) = 0.
      F(1) = 0.
      F(M-1) = 0.
      F(M) = 0.
      C(2) = C(4) - C(3)
      D(2) = 1./D(2)
C
      IF (M.LT.3) GO TO 100
      DO 90 I=3,M
        Q = D(I-1)*E(I-1)
        D(I) = 1./(D(I)-P*F(I-2)-Q*E(I-1))
        E(I) = E(I) - Q*F(I-1)
        C(I) = C(I+2) - C(I+1) - P*C(I-2) - Q*C(I-1)
        P = D(I-1)*F(I-1)
   90 CONTINUE
C
  100 I = N - 1
      C(N-1) = 0.
      C(N) = 0.
      IF (N.LT.4) GO TO 120
      DO 110 M=4,N
C        I = N-2, ..., 2
        I = I - 1
        C(I) = (C(I)-E(I)*C(I+1)-F(I)*C(I+2))*D(I)
  110 CONTINUE
C
C     INTEGRATE THE THIRD DERIVATIVE OF S(X).
C
  120 M = N - 1
      Q = X(2) - X(1)
      R = X(3) - X(2)
      B1 = B(2)
      Q3 = Q*Q*Q
      QR = Q + R
      IF (QR) 140, 130, 140
  130 V = 0.
      T = 0.
      GO TO 150
  140 V = C(2)/QR
      T = V
  150 F(1) = 0.
      IF (Q.NE.0.) F(1) = V/Q
      DO 180 I=2,M
        P = Q
        Q = R
        R = 0.
        IF (I.NE.M) R = X(I+2) - X(I+1)
        P3 = Q3
        Q3 = Q*Q*Q
        PQ = QR
        QR = Q + R
        S = T
        T = 0.
        IF (QR.NE.0.) T = (C(I+1)-C(I))/QR
        U = V
        V = T - S
        IF (PQ) 170, 160, 170
  160   C(I) = C(I-1)
        D(I) = 0.
        E(I) = 0.
        F(I) = 0.
        GO TO 180
  170   F(I) = F(I-1)
        IF (Q.NE.0.) F(I) = V/Q
        E(I) = 5.*S
        D(I) = 10.*(C(I)-Q*S)
        C(I) = D(I)*(P-Q) + (B(I+1)-B(I)+(U-E(I))*P3-(V+E(I))*Q3)/PQ
        B(I) = (P*(B(I+1)-V*Q3)+Q*(B(I)-U*P3))/PQ -
     *   P*Q*(D(I)+E(I)*(Q-P))
  180 CONTINUE
C
C     END POINTS X(1) AND X(N).
C
      P = X(2) - X(1)
      S = F(1)*P*P*P
      E(1) = 0.
      D(1) = 0.
      C(1) = C(2) - 10.*S
      B(1) = B1 - (C(1)+S)*P
C
      Q = X(N) - X(N-1)
      T = F(N-1)*Q*Q*Q
      E(N) = 0.
      D(N) = 0.
      C(N) = C(N-1) + 10.*T
      B(N) = B(N) + (C(N)-T)*Q
  190 RETURN
      END
C                                                                       QUI   10
C                                                                       QUI   20
      SUBROUTINE QUINEQ(N, Y, B, C, D, E, F)                            QUI   30
C
      INTEGER N
      REAL*8 Y(N), B(N), C(N), D(N), E(N), F(N)
C
C
C
C     QUINEQ COMPUTES THE COEFFICIENTS OF QUINTIC NATURAL QUINTIC SPLINE
C     S(X) WITH EQUIDISTANT KNOTS X(I) INTERPOLATING THERE TO GIVEN
C     FUNCTION VALUES:
C               S(X(I)) = Y(I)  FOR I = 1,2, ..., N.
C     IN EACH INTERVAL (X(I),X(I+1)) THE SPLINE FUNCTION S(XX) IS
C     A POLYNOMIAL OF FIFTH DEGREE:
C     S(XX)=((((F(I)*P+E(I))*P+D(I))*P+C(I))*P+B(I))*P+Y(I)    (*)
C          =((((-F(I)*Q+E(I+1))*Q-D(I+1))*Q+C(I+1))*Q-B(I+1))*Q+Y(I+1)
C     WHERE  P = (XX - X(I))/X(I+1) - X(I))
C     AND    Q = (X(I+1) - XX)/(X(I+1) - X(I)).
C     (NOTE THE FIRST SUBSCRIPT IN THE SECOND EXPRESSION.)
C     THE DIFFERENT POLYNOMIALS ARE PIECED TOGETHER SO THAT S(X) AND
C     ITS DERIVATIVES UP TO S"" ARE CONTINUOUS.
C
C        INPUT:
C
C     N          NUMBER OF DATA POINTS, (AT LEAST THREE, I.E. N > 2)
C     Y(1:N)     THE PRESCRIBED FUNCTION VALUES AT THE KNOTS
C
C        OUTPUT:
C
C     B,C,D,E,F  THE COMPUTED SPLINE COEFFICIENTS AS IN (*).
C         (1:N)  IF X(I+1) - X(I) = 1., THEN SPECIFICALLY:
C                B(I) = S'X(I)), C(I) = S"(X(I))/2, D(I) = S"'(X(I))/6,
C                E(I) = S""(X(I))/24,  F(I) = S""'(X(I)+0)/120.
C                F(N) IS NEITHER USED NOR ALTERED.  THE ARRAYS
C                Y,B,C,D MUST ALWAYS BE DISTINCT.  IF E AND F ARE
C                NOT WANTED, THE CALL QUINEQ(N,Y,B,C,D,D,D) MAY
C                BE USED TO SAVE STORAGE LOCATIONS.
C
      INTEGER I, M
      REAL*8 P, Q, R, S, T, U, V
C
      IF (N.LE.2) GO TO 50
C
      M = N - 3
      P = 0.
      Q = 0.
      R = 0.
      S = 0.
      T = 0.
      D(M+1) = 0.
      D(M+2) = 0.
      IF (M.LE.0) GO TO 30
      DO 10 I=1,M
        U = P*R
        B(I) = 1./(66.-U*R-Q)
        R = 26. - U
        C(I) = R
        D(I) = Y(I+3) - 3.*(Y(I+2)-Y(I+1)) - Y(I) - U*S - Q*T
        Q = P
        P = B(I)
        T = S
        S = D(I)
   10 CONTINUE
C
      I = N - 2
      DO 20 M=4,N
C        I    = N-3, ..., 1
        I = I - 1
        D(I) = (D(I)-C(I)*D(I+1)-D(I+2))*B(I)
   20 CONTINUE
C
   30 M = N - 1
      Q = 0.
      R = D(1)
      T = R
      V = R
      DO 40 I=2,M
        P = Q
        Q = R
        R = D(I)
        S = T
        T = P - Q - Q + R
        F(I) = T
        U = 5.*(-P+Q)
        E(I) = U
        D(I) = 10.*(P+Q)
        C(I) = 0.5*(Y(I+1)+Y(I-1)+S-T) - Y(I) - U
        B(I) = 0.5*(Y(I+1)-Y(I-1)-S-T) - D(I)
   40 CONTINUE
C
      F(1) = V
      E(1) = 0.
      D(1) = 0.
      C(1) = C(2) - 10.*V
      B(1) = Y(2) - Y(1) - C(1) - V
      E(N) = 0.
      D(N) = 0.
      C(N) = C(N-1) + 10.*T
      B(N) = Y(N) - Y(N-1) + C(N) - T
   50 RETURN
      END
C                                                                       QUI   10
C                                                                       QUI   20
C                                                                       QUI   30
C                                                                       QUI   40
      SUBROUTINE QUINDF(N, X, Y, B, C, D, E, F)                         QUI   50
C
      INTEGER N
      REAL*8 X(N), Y(N), B(N), C(N), D(N), E(N), F(N)
C
C
C
C     QUINDF COMPUTES THE COEFFICIENTS OF A QUINTIC NATURAL QUINTIC
C     SPLINE S(X) WITH KNOTS X(I) FOR WHICH THE FUNCTION VALUES Y(I) AND
C     THE FIRST DERIVATIVES B(I) ARE SPECIFIED AT X(I), I = 1,2, ...,N.
C     IN EACH INTERVAL (X(I),X(I+1)) THE SPLINE FUNCTION S(XX) IS
C     A POLYNOMIAL OF FIFTH DEGREE:
C     S(XX)=((((F(I)*P+E(I))*P+D(I))*P+C(I))*P+B(I))*P+Y(I)    (*)
C     WHERE  P = XX - X(I).
C
C        INPUT:
C
C     N          NUMBER OF DATA POINTS, (AT LEAST TWO, I.E. N > 1)
C     X(1:N)     THE STRICTLY INCREASING OR DECREASING SEQUENCE OF
C                KNOTS.  THE SPACING MUST BE SUCH THAT THE FIFTH
C                POWER OF X(I+1) - X(I) CAN BE FORMED WITHOUT
C                OVERFLOW OR UNDERFLOW OF EXPONENTS
C     Y(1:N)     THE PRESCRIBED FUNCTION VALUES AT THE KNOTS
C     B(1:N)     THE PRESCRIBED DERIVATIVE VALUES AT THE KNOTS
C
C        OUTPUT:
C
C     C,D,E,F    THE COMPUTED SPLINE COEFFICIENTS AS IN (*).
C         (1:N)  E(N) AND F(N) ARE NEITHER USED NOR ALTERED.
C                THE ARRAYS C,D,E,F MUST ALWAYS BE DISTINCT.
C
      INTEGER I, M, N1
      REAL*8 CC, G, H, HH, H2, P, PP, Q, QQ, R, RR
C
      IF (N.LE.1) GO TO 40
      N1 = N - 1
      CC = 0.
      HH = 0.
      PP = 0.
      QQ = 0.
      RR = 0.
      G = 0.
      DO 10 I=1,N1
        H = 1./(X(I+1)-X(I))
        H2 = H*H
        D(I) = 3.*(HH+H) - G*HH
        P = (Y(I+1)-Y(I))*H2*H
        Q = (B(I+1)+B(I))*H2
        R = (B(I+1)-B(I))*H2
        CC = 10.*(P-PP) - 5.*(Q-QQ) + R + RR + G*CC
        C(I) = CC
        G = H/D(I)
        HH = H
        PP = P
        QQ = Q
        RR = R
   10 CONTINUE
C
      C(N) = (-10.*PP+5.*QQ+RR+G*CC)/(3.*HH-G*HH)
      I = N
      DO 20 M=1,N1
C        I      = N-1, ..., 1
        I = I - 1
        D(I+1) = 1./(X(I+1)-X(I))
        C(I) = (C(I)+C(I+1)*D(I+1))/D(I)
   20 CONTINUE
C
      DO 30 I=1,N1
        H = D(I+1)
        P = (((Y(I+1)-Y(I))*H-B(I))*H-C(I))*H
        Q = ((B(I+1)-B(I))*H-C(I)-C(I))*H
        R = (C(I+1)-C(I))*H
        G = Q - 3.*P
        RR = R - 3.*(P+G)
        QQ = -RR - RR + G
        F(I) = RR*H*H
        E(I) = QQ*H
        D(I) = -RR - QQ + P
   30 CONTINUE
C
      D(N) = 0.
      E(N) = 0.
      F(N) = 0.
   40 RETURN
      END
c$$$C    DRIVER PROGRAM FOR TEST OF QUINAT, QUINEQ AND QUINDF               MAN   10
c$$$C       FOLLOWS.                                                        MAN   20
c$$$C                                                                       MAN   30
c$$$      INTEGER N,NM1,M,MM,MM1,I,K,J,JJ                                   MAN   40
c$$$      REAL*8 Z                                                          MAN   50
c$$$      REAL*8 X(200),Y(200),B(200),BB(200),CC(200),DD(200),EE(200),      MAN   60
c$$$     *                 FF(200),A(200,6),C(6),DIFF(5),COM(5)             MAN   70
c$$$C                                                                       MAN   80
c$$$C     N          NUMBER OF DATA POINTS.                                 MAN   90
c$$$C     M          2*M-1 IS ORDER OF SPLINE.                              MAN  100
c$$$C                   M = 3 ALWAYS FOR QUINTIC SPLINE.                    MAN  110
c$$$C     NN,NM1,MM,                                                        MAN  120
c$$$C     MM1,I,K,                                                          MAN  130
c$$$C     J,JJ       TEMPORARY INTEGER VARIABLES.                           MAN  140
c$$$C     Z,P        TEMPORARY REAL*8 VARIABLES.                            MAN  150
c$$$C     X(1:N)     THE SEQUENCE OF KNOTS.                                 MAN  160
c$$$C     Y(1:N)     THE PRESCRIBED FUNCTION VALUES AT THE KNOTS.           MAN  170
c$$$C     B(1:N)     THE PRESCRIBED DERIVATIVE VALUES AT THE KNOTS.         MAN  180
c$$$C     BB,CC,DD,                                                         MAN  190
c$$$C     EE,FF(1:N) THE COMPUTED SPLINE COEFFICIENTS                       MAN  200
c$$$C     A(1:N,1:6) TWO DIMENSIONAL ARRAY WHOSE COLUMNS ARE                MAN  210
c$$$C                   Y, BB, CC, DD, EE, FF.                              MAN  220
c$$$C     DIFF(1:5)  MAXIMUM VALUES OF DIFFERENCES OF VALUES AND            MAN  230
c$$$C                   DERIVATIVES TO RIGHT AND LEFT OF KNOTS.             MAN  240
c$$$C     COM(1:5)   MAXIMUM VALUES OF COEFFICIENTS.                        MAN  250
c$$$C                                                                       MAN  260
c$$$C                                                                       MAN  270
c$$$C     TEST OF QUINAT WITH NONEQUIDISTANT KNOTS AND                      MAN  280
c$$$C        EQUIDISTANT KNOTS FOLLOWS.                                     MAN  290
c$$$C                                                                       MAN  300
c$$$      NOUT=6                                                            MAN  310
c$$$      WRITE(NOUT,10)                                                    MAN  320
c$$$   10 FORMAT(50H1         TEST OF QUINAT WITH NONEQUIDISTANT KNOTS)     MAN  330
c$$$      N = 5                                                             MAN  340
c$$$      X(1) = -3.0                                                       MAN  350
c$$$      X(2) = -1.0                                                       MAN  360
c$$$      X(3) =  0.0                                                       MAN  370
c$$$      X(4) =  3.0                                                       MAN  380
c$$$      X(5) =  4.0                                                       MAN  390
c$$$      Y(1) =  7.0                                                       MAN  400
c$$$      Y(2) = 11.0                                                       MAN  410
c$$$      Y(3) = 26.0                                                       MAN  420
c$$$      Y(4) = 56.0                                                       MAN  430
c$$$      Y(5) = 29.0                                                       MAN  440
c$$$      M = 3                                                             MAN  450
c$$$      MM = 2*M                                                          MAN  460
c$$$      MM1 = MM - 1                                                      MAN  470
c$$$      WRITE(NOUT,15) N,M                                                MAN  480
c$$$   15 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN  490
c$$$      CALL QUINAT(N,X,Y,BB,CC,DD,EE,FF)                                 MAN  500
c$$$      DO 20 I = 1,N                                                     MAN  510
c$$$         A(I,1) = Y(I)                                                  MAN  520
c$$$         A(I,2) = BB(I)                                                 MAN  530
c$$$         A(I,3) = CC(I)                                                 MAN  540
c$$$         A(I,4) = DD(I)                                                 MAN  550
c$$$         A(I,5) = EE(I)                                                 MAN  560
c$$$         A(I,6) = FF(I)                                                 MAN  570
c$$$   20 CONTINUE                                                          MAN  580
c$$$      DO 30 I = 1,MM1                                                   MAN  590
c$$$         DIFF(I) = 0.0                                                  MAN  600
c$$$         COM(I) = 0.0                                                   MAN  610
c$$$   30 CONTINUE                                                          MAN  620
c$$$      DO 70 K = 1,N                                                     MAN  630
c$$$         DO 35 I = 1,MM                                                 MAN  640
c$$$            C(I) = A(K,I)                                               MAN  650
c$$$   35 CONTINUE                                                          MAN  660
c$$$      WRITE(NOUT,40) K                                                  MAN  670
c$$$   40 FORMAT(40H ---------------------------------------,I3,            MAN  680
c$$$     *  45H --------------------------------------------)               MAN  690
c$$$      WRITE(NOUT,45) X(K)                                               MAN  700
c$$$   45 FORMAT(F12.8)                                                     MAN  710
c$$$      IF (K .EQ. N) WRITE(NOUT,55) C(1)                                 MAN  720
c$$$      IF (K .EQ. N) GO TO 75                                            MAN  730
c$$$      WRITE(NOUT,55) (C(I),I=1,MM)                                      MAN  740
c$$$      DO 50 I = 1,MM1                                                   MAN  750
c$$$         IF (ABS(A(K,I)) .GT. COM(I)) COM(I) = ABS(A(K,I))              MAN  760
c$$$   50 CONTINUE                                                          MAN  770
c$$$   55 FORMAT(6F16.8)                                                    MAN  780
c$$$      Z = X(K+1) - X(K)                                                 MAN  790
c$$$      DO 60 I = 2,MM                                                    MAN  800
c$$$         DO 60 JJ = I,MM                                                MAN  810
c$$$            J = MM + I - JJ                                             MAN  820
c$$$            C(J-1) = C(J)*Z + C(J-1)                                    MAN  830
c$$$   60 CONTINUE                                                          MAN  840
c$$$      WRITE(NOUT,55) (C(I),I=1,MM)                                      MAN  850
c$$$      DO 65 I = 1,MM1                                                   MAN  860
c$$$         IF ((K .GE. N-1) .AND. (I .NE. 1)) GO TO 65                    MAN  870
c$$$         Z = ABS(C(I) - A(K+1,I))                                       MAN  880
c$$$         IF (Z .GT. DIFF(I)) DIFF(I) = Z                                MAN  890
c$$$   65 CONTINUE                                                          MAN  900
c$$$   70 CONTINUE                                                          MAN  910
c$$$   75 WRITE(NOUT,80)                                                    MAN  920
c$$$   80 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN  930
c$$$      WRITE(NOUT,85) (DIFF(I),I=1,MM1)                                  MAN  940
c$$$   85 FORMAT(5E18.9)                                                    MAN  950
c$$$      WRITE(NOUT,90)                                                    MAN  960
c$$$   90 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN  970
c$$$      IF (ABS(C(1)) .GT. COM(1)) COM(1) =ABS(C(1))                      MAN  980
c$$$      WRITE(NOUT,95) (COM(I),I=1,MM1)                                   MAN  990
c$$$   95 FORMAT(5F16.8)                                                    MAN 1000
c$$$      M = 3                                                             MAN 1010
c$$$      DO 200 N = 10,100,10                                              MAN 1020
c$$$      MM = 2*M                                                          MAN 1030
c$$$      MM1 = MM - 1                                                      MAN 1040
c$$$      NM1 = N -1                                                        MAN 1050
c$$$      DO 100 I = 1,NM1,2                                                MAN 1060
c$$$         X(I)   = I                                                     MAN 1070
c$$$         X(I+1) = I + 1                                                 MAN 1080
c$$$         Y(I)   = 1.                                                    MAN 1090
c$$$         Y(I+1) = 0.                                                    MAN 1100
c$$$  100 CONTINUE                                                          MAN 1110
c$$$      IF (MOD(N,2) .EQ. 0) GOTO 105                                     MAN 1120
c$$$      X(N) = N                                                          MAN 1130
c$$$      Y(N) = 1.                                                         MAN 1140
c$$$  105 CONTINUE                                                          MAN 1150
c$$$      WRITE(NOUT,110) N,M                                               MAN 1160
c$$$  110 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 1170
c$$$      CALL QUINAT(N,X,Y,BB,CC,DD,EE,FF)                                 MAN 1180
c$$$      DO 115 I = 1,N                                                    MAN 1190
c$$$         A(I,1) = Y(I)                                                  MAN 1200
c$$$         A(I,2) = BB(I)                                                 MAN 1210
c$$$         A(I,3) = CC(I)                                                 MAN 1220
c$$$         A(I,4) = DD(I)                                                 MAN 1230
c$$$         A(I,5) = EE(I)                                                 MAN 1240
c$$$         A(I,6) = FF(I)                                                 MAN 1250
c$$$  115 CONTINUE                                                          MAN 1260
c$$$      DO 120 I = 1, MM1                                                 MAN 1270
c$$$         DIFF(I) = 0.0                                                  MAN 1280
c$$$         COM(I) = 0.0                                                   MAN 1290
c$$$  120 CONTINUE                                                          MAN 1300
c$$$      DO 165 K = 1,N                                                    MAN 1310
c$$$         DO 125 I = 1,MM                                                MAN 1320
c$$$            C(I) = A(K,I)                                               MAN 1330
c$$$  125    CONTINUE                                                       MAN 1340
c$$$         IF (N .GT. 10) GOTO 140                                        MAN 1350
c$$$         WRITE(NOUT,130) K                                              MAN 1360
c$$$  130 FORMAT(40H ---------------------------------------,I3,            MAN 1370
c$$$     *       45H --------------------------------------------)          MAN 1380
c$$$         WRITE(NOUT,135) X(K)                                           MAN 1390
c$$$  135 FORMAT(F12.8)                                                     MAN 1400
c$$$         IF (K .EQ. N) WRITE(NOUT,150) C(1)                             MAN 1410
c$$$  140    CONTINUE                                                       MAN 1420
c$$$         IF (K .EQ. N) GO TO 170                                        MAN 1430
c$$$         IF (N .LE. 10) WRITE(NOUT,150) (C(I), I=1,MM)                  MAN 1440
c$$$         DO 145 I = 1,MM1                                               MAN 1450
c$$$            IF (ABS(A(K,I)) .GT.  COM(I)) COM(I) = ABS(A(K,I))          MAN 1460
c$$$  145    CONTINUE                                                       MAN 1470
c$$$  150 FORMAT(6F16.8)                                                    MAN 1480
c$$$         Z = X(K+1) - X(K)                                              MAN 1490
c$$$         DO 155 I = 2,MM                                                MAN 1500
c$$$            DO 155 JJ = I,MM                                            MAN 1510
c$$$               J = MM + I - JJ                                          MAN 1520
c$$$               C(J-1) = C(J)*Z + C(J-1)                                 MAN 1530
c$$$  155    CONTINUE                                                       MAN 1540
c$$$         IF (N .LE. 10) WRITE(NOUT,150) (C(I), I=1,MM)                  MAN 1550
c$$$         DO 160 I = 1,MM1                                               MAN 1560
c$$$            IF ((K .GE. N-1) .AND. (I .NE. 1)) GO TO 160                MAN 1570
c$$$            Z = ABS(C(I) - A(K+1,I))                                    MAN 1580
c$$$            IF (Z .GT. DIFF(I)) DIFF(I) = Z                             MAN 1590
c$$$  160    CONTINUE                                                       MAN 1600
c$$$  165 CONTINUE                                                          MAN 1610
c$$$  170 WRITE(NOUT,175)                                                   MAN 1620
c$$$  175 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 1630
c$$$      WRITE(NOUT,180) (DIFF(I),I=1,MM1)                                 MAN 1640
c$$$  180 FORMAT(5E18.9)                                                    MAN 1650
c$$$      WRITE(NOUT,185)                                                   MAN 1660
c$$$  185 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 1670
c$$$      IF (ABS(C(1)) .GT. COM(1)) COM(1) = ABS(C(1))                     MAN 1680
c$$$      WRITE(NOUT,190) (COM(I),I=1,MM1)                                  MAN 1690
c$$$  190 FORMAT(5E18.9)                                                    MAN 1700
c$$$  200 CONTINUE                                                          MAN 1710
c$$$C                                                                       MAN 1720
c$$$C                                                                       MAN 1730
c$$$C     TEST OF QUINEQ FOLLOWS.                                           MAN 1740
c$$$C                                                                       MAN 1750
c$$$      WRITE(NOUT,210)                                                   MAN 1760
c$$$  210 FORMAT(18H1   TEST OF QUINEQ)                                     MAN 1770
c$$$      M = 3                                                             MAN 1780
c$$$      DO 400 N = 10,100,10                                              MAN 1790
c$$$      MM = 2*M                                                          MAN 1800
c$$$      MM1 = MM - 1                                                      MAN 1810
c$$$      NM1 = N -1                                                        MAN 1820
c$$$      DO 300 I = 1,NM1,2                                                MAN 1830
c$$$         X(I)   = I                                                     MAN 1840
c$$$         X(I+1) = I + 1                                                 MAN 1850
c$$$         Y(I)   = 1.                                                    MAN 1860
c$$$         Y(I+1) = 0.                                                    MAN 1870
c$$$  300 CONTINUE                                                          MAN 1880
c$$$      IF (MOD(N,2) .EQ. 0) GOTO 305                                     MAN 1890
c$$$      X(N) = FLOAT(N)                                                   MAN 1900
c$$$      Y(N) = 1.                                                         MAN 1910
c$$$  305 CONTINUE                                                          MAN 1920
c$$$      WRITE(NOUT,310) N,M                                               MAN 1930
c$$$  310 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 1940
c$$$      CALL QUINEQ(N,Y,BB,CC,DD,EE,FF)                                   MAN 1950
c$$$      DO 315 I = 1,N                                                    MAN 1960
c$$$         A(I,1) = Y(I)                                                  MAN 1970
c$$$         A(I,2) = BB(I)                                                 MAN 1980
c$$$         A(I,3) = CC(I)                                                 MAN 1990
c$$$         A(I,4) = DD(I)                                                 MAN 2000
c$$$         A(I,5) = EE(I)                                                 MAN 2010
c$$$         A(I,6) = FF(I)                                                 MAN 2020
c$$$  315 CONTINUE                                                          MAN 2030
c$$$      DO 320 I = 1, MM1                                                 MAN 2040
c$$$         DIFF(I) = 0.0                                                  MAN 2050
c$$$         COM(I) = 0.0                                                   MAN 2060
c$$$  320 CONTINUE                                                          MAN 2070
c$$$      DO 365 K = 1,N                                                    MAN 2080
c$$$         DO 325 I = 1,MM                                                MAN 2090
c$$$            C(I) = A(K,I)                                               MAN 2100
c$$$  325    CONTINUE                                                       MAN 2110
c$$$         IF (N .GT. 10) GOTO 340                                        MAN 2120
c$$$         WRITE(NOUT,330) K                                              MAN 2130
c$$$  330 FORMAT(40H ---------------------------------------,I3,            MAN 2140
c$$$     *       45H --------------------------------------------)          MAN 2150
c$$$         WRITE(NOUT,335) X(K)                                           MAN 2160
c$$$  335 FORMAT(F12.8)                                                     MAN 2170
c$$$         IF (K .EQ. N) WRITE(NOUT,350) C(1)                             MAN 2180
c$$$  340    CONTINUE                                                       MAN 2190
c$$$         IF (K .EQ. N) GO TO 370                                        MAN 2200
c$$$         IF (N .LE. 10) WRITE(NOUT,350) (C(I), I=1,MM)                  MAN 2210
c$$$         DO 345 I = 1,MM1                                               MAN 2220
c$$$            IF (ABS(A(K,I)) .GT.  COM(I)) COM(I) = ABS(A(K,I))          MAN 2230
c$$$  345    CONTINUE                                                       MAN 2240
c$$$  350 FORMAT(6F16.8)                                                    MAN 2250
c$$$         Z = 1.                                                         MAN 2260
c$$$         DO 355 I = 2,MM                                                MAN 2270
c$$$            DO 355 JJ = I,MM                                            MAN 2280
c$$$               J = MM + I - JJ                                          MAN 2290
c$$$               C(J-1) = C(J)*Z + C(J-1)                                 MAN 2300
c$$$  355    CONTINUE                                                       MAN 2310
c$$$         IF (N .LE. 10) WRITE(NOUT,350) (C(I), I=1,MM)                  MAN 2320
c$$$         DO 360 I = 1,MM1                                               MAN 2330
c$$$            IF ((K .GE. N-1) .AND. (I .NE. 1)) GO TO 360                MAN 2340
c$$$            Z = ABS(C(I) - A(K+1,I))                                    MAN 2350
c$$$            IF (Z .GT. DIFF(I)) DIFF(I) = Z                             MAN 2360
c$$$  360    CONTINUE                                                       MAN 2370
c$$$  365 CONTINUE                                                          MAN 2380
c$$$  370 WRITE(NOUT,375)                                                   MAN 2390
c$$$  375 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 2400
c$$$      WRITE(NOUT,380) (DIFF(I),I=1,MM1)                                 MAN 2410
c$$$  380 FORMAT(5E18.9)                                                    MAN 2420
c$$$      WRITE(NOUT,385)                                                   MAN 2430
c$$$  385 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 2440
c$$$      IF (ABS(C(1)) .GT. COM(1)) COM(1) = ABS(C(1))                     MAN 2450
c$$$      WRITE(NOUT,390) (COM(I),I=1,MM1)                                  MAN 2460
c$$$  390 FORMAT(5E18.9)                                                    MAN 2470
c$$$  400 CONTINUE                                                          MAN 2480
c$$$C                                                                       MAN 2490
c$$$C                                                                       MAN 2500
c$$$C     TEST OF QUINDF WITH NONEQUIDISTANT KNOTS FOLLOWS.                 MAN 2510
c$$$C                                                                       MAN 2520
c$$$      WRITE(NOUT,410)                                                   MAN 2530
c$$$  410 FORMAT(50H1         TEST OF QUINDF WITH NONEQUIDISTANT KNOTS)     MAN 2540
c$$$      N = 5                                                             MAN 2550
c$$$      X(1) = -3.0                                                       MAN 2560
c$$$      X(2) = -1.0                                                       MAN 2570
c$$$      X(3) =  0.0                                                       MAN 2580
c$$$      X(4) =  3.0                                                       MAN 2590
c$$$      X(5) =  4.0                                                       MAN 2600
c$$$      Y(1) =  7.0                                                       MAN 2610
c$$$      Y(2) = 11.0                                                       MAN 2620
c$$$      Y(3) = 26.0                                                       MAN 2630
c$$$      Y(4) = 56.0                                                       MAN 2640
c$$$      Y(5) = 29.0                                                       MAN 2650
c$$$      B(1) =  2.0                                                       MAN 2660
c$$$      B(2) = 15.0                                                       MAN 2670
c$$$      B(3) = 10.0                                                       MAN 2680
c$$$      B(4) =-27.0                                                       MAN 2690
c$$$      B(5) =-30.0                                                       MAN 2700
c$$$      M = 3                                                             MAN 2710
c$$$      MM = 2*M                                                          MAN 2720
c$$$      MM1 = MM - 1                                                      MAN 2730
c$$$      WRITE(NOUT,415) N,M                                               MAN 2740
c$$$  415 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 2750
c$$$      CALL QUINDF(N,X,Y,B,CC,DD,EE,FF)                                  MAN 2760
c$$$      DO 420 I = 1,N                                                    MAN 2770
c$$$         A(I,1) = Y(I)                                                  MAN 2780
c$$$         A(I,2) = B(I)                                                  MAN 2790
c$$$         A(I,3) = CC(I)                                                 MAN 2800
c$$$         A(I,4) = DD(I)                                                 MAN 2810
c$$$         A(I,5) = EE(I)                                                 MAN 2820
c$$$         A(I,6) = FF(I)                                                 MAN 2830
c$$$  420 CONTINUE                                                          MAN 2840
c$$$      DO 430 I = 1,MM1                                                  MAN 2850
c$$$         DIFF(I) = 0.0                                                  MAN 2860
c$$$         COM(I) = 0.0                                                   MAN 2870
c$$$  430 CONTINUE                                                          MAN 2880
c$$$      DO 470 K = 1,N                                                    MAN 2890
c$$$         DO 435 I = 1,MM                                                MAN 2900
c$$$            C(I) = A(K,I)                                               MAN 2910
c$$$  435 CONTINUE                                                          MAN 2920
c$$$      WRITE(NOUT,440) K                                                 MAN 2930
c$$$  440 FORMAT(40H ---------------------------------------,I3,            MAN 2940
c$$$     *  45H --------------------------------------------)               MAN 2950
c$$$      WRITE(NOUT,445) X(K)                                              MAN 2960
c$$$  445 FORMAT(F12.8)                                                     MAN 2970
c$$$      IF (K .EQ. N) WRITE(NOUT,455) C(1)                                MAN 2980
c$$$      IF (K .EQ. N) GO TO 475                                           MAN 2990
c$$$      WRITE(NOUT,455) (C(I),I=1,MM)                                     MAN 3000
c$$$      DO 450 I = 1,MM1                                                  MAN 3010
c$$$         IF (ABS(A(K,I)) .GT. COM(I)) COM(I) = ABS(A(K,I))              MAN 3020
c$$$  450 CONTINUE                                                          MAN 3030
c$$$  455 FORMAT(6F16.8)                                                    MAN 3040
c$$$      Z = X(K+1) - X(K)                                                 MAN 3050
c$$$      DO 460 I = 2,MM                                                   MAN 3060
c$$$         DO 460 JJ = I,MM                                               MAN 3070
c$$$            J = MM + I - JJ                                             MAN 3080
c$$$            C(J-1) = C(J)*Z + C(J-1)                                    MAN 3090
c$$$  460 CONTINUE                                                          MAN 3100
c$$$      WRITE(NOUT,455) (C(I),I=1,MM)                                     MAN 3110
c$$$      DO 465 I = 1,MM1                                                  MAN 3120
c$$$         IF ((K .GE. N-1) .AND. (I .NE. 1)) GO TO 465                   MAN 3130
c$$$         Z = ABS(C(I) - A(K+1,I))                                       MAN 3140
c$$$         IF (Z .GT. DIFF(I)) DIFF(I) = Z                                MAN 3150
c$$$  465 CONTINUE                                                          MAN 3160
c$$$  470 CONTINUE                                                          MAN 3170
c$$$  475 WRITE(NOUT,480)                                                   MAN 3180
c$$$  480 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 3190
c$$$      WRITE(NOUT,485) (DIFF(I),I=1,MM1)                                 MAN 3200
c$$$  485 FORMAT(5E18.9)                                                    MAN 3210
c$$$      WRITE(NOUT,490)                                                   MAN 3220
c$$$  490 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 3230
c$$$      IF (ABS(C(1)) .GT. COM(1)) COM(1) =ABS(C(1))                      MAN 3240
c$$$      WRITE(NOUT,495) (COM(I),I=1,MM1)                                  MAN 3250
c$$$  495 FORMAT(5F16.8)                                                    MAN 3260
c$$$      M = 3                                                             MAN 3270
c$$$      DO 600 N = 10,100,10                                              MAN 3280
c$$$      MM = 2*M                                                          MAN 3290
c$$$      MM1 = MM - 1                                                      MAN 3300
c$$$      P = 0.0                                                           MAN 3310
c$$$      DO 500 I = 1,N                                                    MAN 3320
c$$$         P = P + ABS(SIN((FLOAT(I)))) + 0.001*FLOAT(I)                  MAN 3330
c$$$         X(I) = P                                                       MAN 3340
c$$$         Y(I) = COS((FLOAT(I))) - 0.5                                   MAN 3350
c$$$         B(I) = COS((FLOAT(2*I))) - 0.5                                 MAN 3360
c$$$  500 CONTINUE                                                          MAN 3370
c$$$      WRITE(NOUT,510) N,M                                               MAN 3380
c$$$  510 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 3390
c$$$      CALL QUINDF(N,X,Y,B,CC,DD,EE,FF)                                  MAN 3400
c$$$      DO 515 I = 1,N                                                    MAN 3410
c$$$         A(I,1) = Y(I)                                                  MAN 3420
c$$$         A(I,2) = B(I)                                                  MAN 3430
c$$$         A(I,3) = CC(I)                                                 MAN 3440
c$$$         A(I,4) = DD(I)                                                 MAN 3450
c$$$         A(I,5) = EE(I)                                                 MAN 3460
c$$$         A(I,6) = FF(I)                                                 MAN 3470
c$$$  515 CONTINUE                                                          MAN 3480
c$$$      DO 520 I = 1, MM1                                                 MAN 3490
c$$$         DIFF(I) = 0.0                                                  MAN 3500
c$$$         COM(I) = 0.0                                                   MAN 3510
c$$$  520 CONTINUE                                                          MAN 3520
c$$$      DO 565 K = 1,N                                                    MAN 3530
c$$$         DO 525 I = 1,MM                                                MAN 3540
c$$$            C(I) = A(K,I)                                               MAN 3550
c$$$  525    CONTINUE                                                       MAN 3560
c$$$         IF (N .GT. 10) GOTO 540                                        MAN 3570
c$$$         WRITE(NOUT,530) K                                              MAN 3580
c$$$  530 FORMAT(40H ---------------------------------------,I3,            MAN 3590
c$$$     *       45H --------------------------------------------)          MAN 3600
c$$$         WRITE(NOUT,535) X(K)                                           MAN 3610
c$$$  535 FORMAT(F12.8)                                                     MAN 3620
c$$$         IF (K .EQ. N) WRITE(NOUT,550) C(1)                             MAN 3630
c$$$  540    CONTINUE                                                       MAN 3640
c$$$         IF (K .EQ. N) GO TO 570                                        MAN 3650
c$$$         IF (N .LE. 10) WRITE(NOUT,550) (C(I), I=1,MM)                  MAN 3660
c$$$         DO 545 I = 1,MM1                                               MAN 3670
c$$$            IF (ABS(A(K,I)) .GT.  COM(I)) COM(I) = ABS(A(K,I))          MAN 3680
c$$$  545    CONTINUE                                                       MAN 3690
c$$$  550 FORMAT(6F16.8)                                                    MAN 3700
c$$$         Z = X(K+1) - X(K)                                              MAN 3710
c$$$         DO 555 I = 2,MM                                                MAN 3720
c$$$            DO 555 JJ = I,MM                                            MAN 3730
c$$$               J = MM + I - JJ                                          MAN 3740
c$$$               C(J-1) = C(J)*Z + C(J-1)                                 MAN 3750
c$$$  555    CONTINUE                                                       MAN 3760
c$$$         IF (N .LE. 10) WRITE(NOUT,550) (C(I), I=1,MM)                  MAN 3770
c$$$         DO 560 I = 1,MM1                                               MAN 3780
c$$$            IF ((K .GE. N-1) .AND. (I .NE. 1)) GO TO 560                MAN 3790
c$$$            Z = ABS(C(I) - A(K+1,I))                                    MAN 3800
c$$$            IF (Z .GT. DIFF(I)) DIFF(I) = Z                             MAN 3810
c$$$  560    CONTINUE                                                       MAN 3820
c$$$  565 CONTINUE                                                          MAN 3830
c$$$  570 WRITE(NOUT,575)                                                   MAN 3840
c$$$  575 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 3850
c$$$      WRITE(NOUT,580) (DIFF(I),I=1,MM1)                                 MAN 3860
c$$$  580 FORMAT(5E18.9)                                                    MAN 3870
c$$$      WRITE(NOUT,585)                                                   MAN 3880
c$$$  585 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 3890
c$$$      IF (ABS(C(1)) .GT. COM(1)) COM(1) = ABS(C(1))                     MAN 3900
c$$$      WRITE(NOUT,590) (COM(I),I=1,MM1)                                  MAN 3910
c$$$  590 FORMAT(5E18.9)                                                    MAN 3920
c$$$  600 CONTINUE                                                          MAN 3930
c$$$C                                                                       MAN 3940
c$$$C                                                                       MAN 3950
c$$$C     TEST OF QUINAT WITH NONEQUIDISTANT DOUBLE KNOTS FOLLOWS.          MAN 3960
c$$$C                                                                       MAN 3970
c$$$      WRITE(NOUT,610)                                                   MAN 3980
c$$$  610 FORMAT(50H1  TEST OF QUINAT WITH NONEQUIDISTANT DOUBLE KNOTS)     MAN 3990
c$$$      N = 5                                                             MAN 4000
c$$$      X(1) = -3.                                                        MAN 4010
c$$$      X(2) = -3.                                                        MAN 4020
c$$$      X(3) = -1.                                                        MAN 4030
c$$$      X(4) = -1.                                                        MAN 4040
c$$$      X(5) = 0.                                                         MAN 4050
c$$$      X(6) = 0.                                                         MAN 4060
c$$$      X(7) = 3.                                                         MAN 4070
c$$$      X(8) = 3.                                                         MAN 4080
c$$$      X(9) = 4.                                                         MAN 4090
c$$$      X(10) = 4.                                                        MAN 4100
c$$$      Y(1) = 7.                                                         MAN 4110
c$$$      Y(2) = 2.                                                         MAN 4120
c$$$      Y(3) = 11.                                                        MAN 4130
c$$$      Y(4) = 15.                                                        MAN 4140
c$$$      Y(5) = 26.                                                        MAN 4150
c$$$      Y(6) = 10.                                                        MAN 4160
c$$$      Y(7) = 56.                                                        MAN 4170
c$$$      Y(8) = -27.                                                       MAN 4180
c$$$      Y(9) = 29.                                                        MAN 4190
c$$$      Y(10) = -30.                                                      MAN 4200
c$$$      M = 3                                                             MAN 4210
c$$$      NN = 2*N                                                          MAN 4220
c$$$      MM = 2*M                                                          MAN 4230
c$$$      MM1 = MM - 1                                                      MAN 4240
c$$$      WRITE(NOUT,615) N,M                                               MAN 4250
c$$$  615 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 4260
c$$$      CALL QUINAT(NN,X,Y,BB,CC,DD,EE,FF)                                MAN 4270
c$$$      DO 620 I = 1,NN                                                   MAN 4280
c$$$         A(I,1) = Y(I)                                                  MAN 4290
c$$$         A(I,2) = BB(I)                                                 MAN 4300
c$$$         A(I,3) = CC(I)                                                 MAN 4310
c$$$         A(I,4) = DD(I)                                                 MAN 4320
c$$$         A(I,5) = EE(I)                                                 MAN 4330
c$$$         A(I,6) = FF(I)                                                 MAN 4340
c$$$  620 CONTINUE                                                          MAN 4350
c$$$      DO 630 I = 1,MM1                                                  MAN 4360
c$$$         DIFF(I) = 0.0                                                  MAN 4370
c$$$         COM(I) = 0.0                                                   MAN 4380
c$$$  630 CONTINUE                                                          MAN 4390
c$$$      DO 670 K = 1,NN                                                   MAN 4400
c$$$         DO 635 I = 1,MM                                                MAN 4410
c$$$            C(I) = A(K,I)                                               MAN 4420
c$$$  635 CONTINUE                                                          MAN 4430
c$$$      WRITE(NOUT,640) K                                                 MAN 4440
c$$$  640 FORMAT(40H ---------------------------------------,I3,            MAN 4450
c$$$     *  45H --------------------------------------------)               MAN 4460
c$$$      WRITE(NOUT,645) X(K)                                              MAN 4470
c$$$  645 FORMAT(F12.8)                                                     MAN 4480
c$$$      IF (K .EQ. NN) WRITE(NOUT,655) C(1)                               MAN 4490
c$$$      IF (K .EQ. NN) GO TO 675                                          MAN 4500
c$$$      WRITE(NOUT,655) (C(I),I=1,MM)                                     MAN 4510
c$$$      DO 650 I = 1,MM1                                                  MAN 4520
c$$$         IF (ABS(A(K,I)) .GT. COM(I)) COM(I) = ABS(A(K,I))              MAN 4530
c$$$  650 CONTINUE                                                          MAN 4540
c$$$  655 FORMAT(6F16.8)                                                    MAN 4550
c$$$      Z = X(K+1) - X(K)                                                 MAN 4560
c$$$      DO 660 I = 2,MM                                                   MAN 4570
c$$$         DO 660 JJ = I,MM                                               MAN 4580
c$$$            J = MM + I - JJ                                             MAN 4590
c$$$            C(J-1) = C(J)*Z + C(J-1)                                    MAN 4600
c$$$  660 CONTINUE                                                          MAN 4610
c$$$      WRITE(NOUT,655) (C(I),I=1,MM)                                     MAN 4620
c$$$      DO 665 I = 1,MM1                                                  MAN 4630
c$$$         IF ((K .GE. NN-1) .AND. (I .NE. 1)) GO TO 665                  MAN 4640
c$$$         Z = ABS(C(I) - A(K+1,I))                                       MAN 4650
c$$$         IF (Z .GT. DIFF(I)) DIFF(I) = Z                                MAN 4660
c$$$  665 CONTINUE                                                          MAN 4670
c$$$  670 CONTINUE                                                          MAN 4680
c$$$  675 WRITE(NOUT,680)                                                   MAN 4690
c$$$  680 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 4700
c$$$      WRITE(NOUT,685) (DIFF(I),I=1,MM1)                                 MAN 4710
c$$$  685 FORMAT(5E18.9)                                                    MAN 4720
c$$$      WRITE(NOUT,690)                                                   MAN 4730
c$$$  690 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 4740
c$$$      IF (ABS(C(1)) .GT. COM(1)) COM(1) =ABS(C(1))                      MAN 4750
c$$$      WRITE(NOUT,695) (COM(I),I=1,MM1)                                  MAN 4760
c$$$  695 FORMAT(5F16.8)                                                    MAN 4770
c$$$      M = 3                                                             MAN 4780
c$$$      DO 800 N = 10,100,10                                              MAN 4790
c$$$      NN = 2*N                                                          MAN 4800
c$$$      MM = 2*M                                                          MAN 4810
c$$$      MM1 = MM - 1                                                      MAN 4820
c$$$      P = 0.0                                                           MAN 4830
c$$$      DO 700 I = 1,N                                                    MAN 4840
c$$$         P = P + ABS(SIN(FLOAT(I)))                                     MAN 4850
c$$$         X(2*I-1) = P                                                   MAN 4860
c$$$         X(2*I)   = P                                                   MAN 4870
c$$$         Y(2*I-1) = COS(FLOAT(I)) - 0.5                                 MAN 4880
c$$$         Y(2*I)   = COS(FLOAT(2*I)) - 0.5                               MAN 4890
c$$$  700 CONTINUE                                                          MAN 4900
c$$$      WRITE(NOUT,710) N,M                                               MAN 4910
c$$$  710 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 4920
c$$$      CALL QUINAT(NN,X,Y,BB,CC,DD,EE,FF)                                MAN 4930
c$$$      DO 715 I = 1,NN                                                   MAN 4940
c$$$         A(I,1) = Y(I)                                                  MAN 4950
c$$$         A(I,2) = BB(I)                                                 MAN 4960
c$$$         A(I,3) = CC(I)                                                 MAN 4970
c$$$         A(I,4) = DD(I)                                                 MAN 4980
c$$$         A(I,5) = EE(I)                                                 MAN 4990
c$$$         A(I,6) = FF(I)                                                 MAN 5000
c$$$  715 CONTINUE                                                          MAN 5010
c$$$      DO 720 I = 1, MM1                                                 MAN 5020
c$$$         DIFF(I) = 0.0                                                  MAN 5030
c$$$         COM(I) = 0.0                                                   MAN 5040
c$$$  720 CONTINUE                                                          MAN 5050
c$$$      DO 765 K = 1,NN                                                   MAN 5060
c$$$         DO 725 I = 1,MM                                                MAN 5070
c$$$            C(I) = A(K,I)                                               MAN 5080
c$$$  725    CONTINUE                                                       MAN 5090
c$$$         IF (N .GT. 10)  GOTO 740                                       MAN 5100
c$$$         WRITE(NOUT,730) K                                              MAN 5110
c$$$  730 FORMAT(40H ---------------------------------------,I3,            MAN 5120
c$$$     *       45H --------------------------------------------)          MAN 5130
c$$$         WRITE(NOUT,735) X(K)                                           MAN 5140
c$$$  735 FORMAT(F12.8)                                                     MAN 5150
c$$$         IF (K .EQ. NN) WRITE(NOUT,750) C(1)                            MAN 5160
c$$$  740    CONTINUE                                                       MAN 5170
c$$$         IF (K .EQ. NN) GO TO 770                                       MAN 5180
c$$$         IF (N .LE. 10) WRITE(NOUT,750) (C(I), I=1,MM)                  MAN 5190
c$$$         DO 745 I = 1,MM1                                               MAN 5200
c$$$            IF (ABS(A(K,I)) .GT.  COM(I)) COM(I) = ABS(A(K,I))          MAN 5210
c$$$  745    CONTINUE                                                       MAN 5220
c$$$  750 FORMAT(6F16.8)                                                    MAN 5230
c$$$         Z = X(K+1) - X(K)                                              MAN 5240
c$$$         DO 755 I = 2,MM                                                MAN 5250
c$$$            DO 755 JJ = I,MM                                            MAN 5260
c$$$               J = MM + I - JJ                                          MAN 5270
c$$$               C(J-1) = C(J)*Z + C(J-1)                                 MAN 5280
c$$$  755    CONTINUE                                                       MAN 5290
c$$$         IF (N .LE. 10) WRITE(NOUT,750) (C(I), I=1,MM)                  MAN 5300
c$$$         DO 760 I = 1,MM1                                               MAN 5310
c$$$            IF ((K .GE. NN-1) .AND. (I .NE. 1)) GO TO 760               MAN 5320
c$$$            Z = ABS(C(I) - A(K+1,I))                                    MAN 5330
c$$$            IF (Z .GT. DIFF(I)) DIFF(I) = Z                             MAN 5340
c$$$  760    CONTINUE                                                       MAN 5350
c$$$  765 CONTINUE                                                          MAN 5360
c$$$  770 WRITE(NOUT,775)                                                   MAN 5370
c$$$  775 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 5380
c$$$      WRITE(NOUT,780) (DIFF(I),I=1,MM1)                                 MAN 5390
c$$$  780 FORMAT(5E18.9)                                                    MAN 5400
c$$$      WRITE(NOUT,785)                                                   MAN 5410
c$$$  785 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 5420
c$$$      IF (ABS(C(1)) .GT. COM(1)) COM(1) = ABS(C(1))                     MAN 5430
c$$$      WRITE(NOUT,790) (COM(I),I=1,MM1)                                  MAN 5440
c$$$  790 FORMAT(5E18.9)                                                    MAN 5450
c$$$  800 CONTINUE                                                          MAN 5460
c$$$C                                                                       MAN 5470
c$$$C                                                                       MAN 5480
c$$$C     TEST OF QUINAT WITH NONEQUIDISTANT KNOTS, ONE DOUBLE KNOT,        MAN 5490
c$$$C        ONE TRIPLE KNOT, FOLLOWS.                                      MAN 5500
c$$$C                                                                       MAN 5510
c$$$      WRITE(NOUT,805)                                                   MAN 5520
c$$$      WRITE(NOUT,810)                                                   MAN 5530
c$$$  805 FORMAT(51H1         TEST OF QUINAT WITH NONEQUIDISTANT KNOTS,)    MAN 5540
c$$$  810 FORMAT(40H             ONE DOUBLE, ONE TRIPLE KNOT)               MAN 5550
c$$$      N = 8                                                             MAN 5560
c$$$      X(1) = -3.                                                        MAN 5570
c$$$      X(2) = -1.                                                        MAN 5580
c$$$      X(3) =  -1.                                                       MAN 5590
c$$$      X(4) =  0.                                                        MAN 5600
c$$$      X(5) =  3.                                                        MAN 5610
c$$$      X(6) = 3.                                                         MAN 5620
c$$$      X(7) = 3.                                                         MAN 5630
c$$$      X(8) = 4.                                                         MAN 5640
c$$$      Y(1) =  7.                                                        MAN 5650
c$$$      Y(2) = 11.                                                        MAN 5660
c$$$      Y(3) = 15.                                                        MAN 5670
c$$$      Y(4) = 26.                                                        MAN 5680
c$$$      Y(5) = 56.                                                        MAN 5690
c$$$      Y(6) = -30.                                                       MAN 5700
c$$$      Y(7) =  -7.                                                       MAN 5710
c$$$      Y(8) =  29.                                                       MAN 5720
c$$$      M = 3                                                             MAN 5730
c$$$      MM = 2*M                                                          MAN 5740
c$$$      MM1 = MM - 1                                                      MAN 5750
c$$$      WRITE(NOUT,815) N,M                                               MAN 5760
c$$$  815 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 5770
c$$$      CALL QUINAT(N,X,Y,BB,CC,DD,EE,FF)                                 MAN 5780
c$$$      DO 820 I = 1,N                                                    MAN 5790
c$$$         A(I,1) = Y(I)                                                  MAN 5800
c$$$         A(I,2) = BB(I)                                                 MAN 5810
c$$$         A(I,3) = CC(I)                                                 MAN 5820
c$$$         A(I,4) = DD(I)                                                 MAN 5830
c$$$         A(I,5) = EE(I)                                                 MAN 5840
c$$$         A(I,6) = FF(I)                                                 MAN 5850
c$$$  820 CONTINUE                                                          MAN 5860
c$$$      DO 830 I = 1,MM1                                                  MAN 5870
c$$$         DIFF(I) = 0.0                                                  MAN 5880
c$$$         COM(I) = 0.0                                                   MAN 5890
c$$$  830 CONTINUE                                                          MAN 5900
c$$$      DO 870 K = 1,N                                                    MAN 5910
c$$$         DO 835 I = 1,MM                                                MAN 5920
c$$$            C(I) = A(K,I)                                               MAN 5930
c$$$  835 CONTINUE                                                          MAN 5940
c$$$      WRITE(NOUT,840) K                                                 MAN 5950
c$$$  840 FORMAT(40H ---------------------------------------,I3,            MAN 5960
c$$$     *  45H --------------------------------------------)               MAN 5970
c$$$      WRITE(NOUT,845) X(K)                                              MAN 5980
c$$$  845 FORMAT(F12.8)                                                     MAN 5990
c$$$      IF (K .EQ. N) WRITE(NOUT,855) C(1)                                MAN 6000
c$$$      IF (K .EQ. N) GO TO 875                                           MAN 6010
c$$$      WRITE(NOUT,855) (C(I),I=1,MM)                                     MAN 6020
c$$$      DO 850 I = 1,MM1                                                  MAN 6030
c$$$         IF (ABS(A(K,I)) .GT. COM(I)) COM(I) = ABS(A(K,I))              MAN 6040
c$$$  850 CONTINUE                                                          MAN 6050
c$$$  855 FORMAT(6F16.8)                                                    MAN 6060
c$$$      Z = X(K+1) - X(K)                                                 MAN 6070
c$$$      DO 860 I = 2,MM                                                   MAN 6080
c$$$         DO 860 JJ = I,MM                                               MAN 6090
c$$$            J = MM + I - JJ                                             MAN 6100
c$$$            C(J-1) = C(J)*Z + C(J-1)                                    MAN 6110
c$$$  860 CONTINUE                                                          MAN 6120
c$$$      WRITE(NOUT,855) (C(I),I=1,MM)                                     MAN 6130
c$$$      DO 865 I = 1,MM1                                                  MAN 6140
c$$$         IF ((K .GE. N-1) .AND. (I .NE. 1)) GO TO 865                   MAN 6150
c$$$         Z = ABS(C(I) - A(K+1,I))                                       MAN 6160
c$$$         IF (Z .GT. DIFF(I)) DIFF(I) = Z                                MAN 6170
c$$$  865 CONTINUE                                                          MAN 6180
c$$$  870 CONTINUE                                                          MAN 6190
c$$$  875 WRITE(NOUT,880)                                                   MAN 6200
c$$$  880 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 6210
c$$$      WRITE(NOUT,885) (DIFF(I),I=1,MM1)                                 MAN 6220
c$$$  885 FORMAT(5E18.9)                                                    MAN 6230
c$$$      WRITE(NOUT,890)                                                   MAN 6240
c$$$  890 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 6250
c$$$      IF (ABS(C(1)) .GT. COM(1)) COM(1) =ABS(C(1))                      MAN 6260
c$$$      WRITE(NOUT,895) (COM(I),I=1,MM1)                                  MAN 6270
c$$$  895 FORMAT(5F16.8)                                                    MAN 6280
c$$$C                                                                       MAN 6290
c$$$C                                                                       MAN 6300
c$$$C     TEST OF QUINAT WITH NONEQUIDISTANT KNOTS, TWO DOUBLE KNOTS,       MAN 6310
c$$$C        ONE TRIPLE KNOT,FOLLOWS.                                       MAN 6320
c$$$C                                                                       MAN 6330
c$$$      WRITE(NOUT,905)                                                   MAN 6340
c$$$      WRITE(NOUT,910)                                                   MAN 6350
c$$$  905 FORMAT(51H1         TEST OF QUINAT WITH NONEQUIDISTANT KNOTS,)    MAN 6360
c$$$  910 FORMAT(40H             TWO DOUBLE, ONE TRIPLE KNOT)               MAN 6370
c$$$      N = 10                                                            MAN 6380
c$$$      X(1) = 0.                                                         MAN 6390
c$$$      X(2) = 2.                                                         MAN 6400
c$$$      X(3) = 2.                                                         MAN 6410
c$$$      X(4) = 3.                                                         MAN 6420
c$$$      X(5) = 3.                                                         MAN 6430
c$$$      X(6) = 3.                                                         MAN 6440
c$$$      X(7) = 5.                                                         MAN 6450
c$$$      X(8) = 8.                                                         MAN 6460
c$$$      X(9) = 9.                                                         MAN 6470
c$$$      X(10)= 9.                                                         MAN 6480
c$$$      Y(1) = 163.                                                       MAN 6490
c$$$      Y(2) = 237.                                                       MAN 6500
c$$$      Y(3) = -127.                                                      MAN 6510
c$$$      Y(4) = 119.                                                       MAN 6520
c$$$      Y(5) = -65.                                                       MAN 6530
c$$$      Y(6) = 192.                                                       MAN 6540
c$$$      Y(7) = 293.                                                       MAN 6550
c$$$      Y(8) =  326.                                                      MAN 6560
c$$$      Y(9) = 0.                                                         MAN 6570
c$$$      Y(10)= -414.0                                                     MAN 6580
c$$$      M = 3                                                             MAN 6590
c$$$      MM = 2*M                                                          MAN 6600
c$$$      MM1 = MM - 1                                                      MAN 6610
c$$$      WRITE(NOUT,915) N,M                                               MAN 6620
c$$$  915 FORMAT(5H-N = ,I3,7H    M =,I2)                                   MAN 6630
c$$$      CALL QUINAT(N,X,Y,BB,CC,DD,EE,FF)                                 MAN 6640
c$$$      DO 920 I = 1,N                                                    MAN 6650
c$$$         A(I,1) = Y(I)                                                  MAN 6660
c$$$         A(I,2) = BB(I)                                                 MAN 6670
c$$$         A(I,3) = CC(I)                                                 MAN 6680
c$$$         A(I,4) = DD(I)                                                 MAN 6690
c$$$         A(I,5) = EE(I)                                                 MAN 6700
c$$$         A(I,6) = FF(I)                                                 MAN 6710
c$$$  920 CONTINUE                                                          MAN 6720
c$$$      DO 930 I = 1,MM1                                                  MAN 6730
c$$$         DIFF(I) = 0.0                                                  MAN 6740
c$$$         COM(I) = 0.0                                                   MAN 6750
c$$$  930 CONTINUE                                                          MAN 6760
c$$$      DO 970 K = 1,N                                                    MAN 6770
c$$$         DO 935 I = 1,MM                                                MAN 6780
c$$$            C(I) = A(K,I)                                               MAN 6790
c$$$  935 CONTINUE                                                          MAN 6800
c$$$      WRITE(NOUT,940) K                                                 MAN 6810
c$$$  940 FORMAT(40H ---------------------------------------,I3,            MAN 6820
c$$$     *  45H --------------------------------------------)               MAN 6830
c$$$      WRITE(NOUT,945) X(K)                                              MAN 6840
c$$$  945 FORMAT(F12.8)                                                     MAN 6850
c$$$      IF (K .EQ. N) WRITE(NOUT,955) C(1)                                MAN 6860
c$$$      IF (K .EQ. N) GO TO 975                                           MAN 6870
c$$$      WRITE(NOUT,955) (C(I),I=1,MM)                                     MAN 6880
c$$$      DO 950 I = 1,MM1                                                  MAN 6890
c$$$         IF (ABS(A(K,I)) .GT. COM(I)) COM(I) = ABS(A(K,I))              MAN 6900
c$$$  950 CONTINUE                                                          MAN 6910
c$$$  955 FORMAT(6F16.8)                                                    MAN 6920
c$$$      Z = X(K+1) - X(K)                                                 MAN 6930
c$$$      DO 960 I = 2,MM                                                   MAN 6940
c$$$         DO 960 JJ = I,MM                                               MAN 6950
c$$$            J = MM + I - JJ                                             MAN 6960
c$$$            C(J-1) = C(J)*Z + C(J-1)                                    MAN 6970
c$$$  960 CONTINUE                                                          MAN 6980
c$$$      WRITE(NOUT,955) (C(I),I=1,MM)                                     MAN 6990
c$$$      DO 965 I = 1,MM1                                                  MAN 7000
c$$$         IF ((K .GE. N-1) .AND. (I .NE. 1)) GO TO 965                   MAN 7010
c$$$         Z = ABS(C(I) - A(K+1,I))                                       MAN 7020
c$$$         IF (Z .GT. DIFF(I)) DIFF(I) = Z                                MAN 7030
c$$$  965 CONTINUE                                                          MAN 7040
c$$$  970 CONTINUE                                                          MAN 7050
c$$$  975 WRITE(NOUT,980)                                                   MAN 7060
c$$$  980 FORMAT(41H  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES )              MAN 7070
c$$$      WRITE(NOUT,985) (DIFF(I),I=1,MM1)                                 MAN 7080
c$$$  985 FORMAT(5E18.9)                                                    MAN 7090
c$$$      WRITE(NOUT,990)                                                   MAN 7100
c$$$  990 FORMAT(42H  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS )             MAN 7110
c$$$      IF (ABS(C(1)) .GT. COM(1)) COM(1) =ABS(C(1))                      MAN 7120
c$$$      WRITE(NOUT,995) (COM(I),I=1,MM1)                                  MAN 7130
c$$$  995 FORMAT(5F16.8)                                                    MAN 7140
c$$$      STOP                                                              MAN 7150
c$$$      END                                                               MAN 7160
