    subroutine matrix_inv ( A, n, NP )
!-------------------------------------------------------------------------
!
!	      Taken from "Numeric recipes".  The original program was
!       GAUSSJ which solves linear equations by the Gauss_Jordon
!       elimination method.  Only the parts required to invert
!	      matrices have been retained.
!
!	      J.P. Griffith  6/88
!
!-------------------------------------------------------------------------
    
    PARAMETER (NMAX=5200)
    implicit double precision(a-h,o-z)

    DIMENSION A(NP,NP), IPIV(NMAX), INDXR(NMAX), INDXC(NMAX)

    double precision :: time_begin_0
    double precision :: time_elasped_0
    double precision :: time_end_0
    
    call cpu_time(time_begin_0)

    DO 11 J=1,N
        IPIV(J)=0
        CONTINUE
    11 END DO
    DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
            IF(IPIV(J) /= 1)THEN
                DO 12 K=1,N
                    IF (IPIV(K) == 0) THEN
                        IF (dABS(A(J,K)) >= BIG)THEN
                            BIG=dABS(A(J,K))
                            IROW=J
                            ICOL=K
                        ENDIF
                    ELSE IF (IPIV(K) > 1) THEN
                        print *, 'Singular matrix'
                    ENDIF
                    CONTINUE
                12 END DO
            ENDIF
            CONTINUE
        13 END DO
        
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW /= ICOL) THEN
            DO 14 L=1,N
                DUM=A(IROW,L)
                A(IROW,L)=A(ICOL,L)
                A(ICOL,L)=DUM
                CONTINUE
            14 END DO
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL) == 0.) print *, 'Singular matrix.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
            A(ICOL,L)=A(ICOL,L)*PIVINV
            CONTINUE
        16 END DO
        DO 21 LL=1,N
            IF(LL /= ICOL)THEN
                DUM=A(LL,ICOL)
                A(LL,ICOL)=0.
                DO 18 L=1,N
                    A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
                    CONTINUE
                18 END DO
            ENDIF
            CONTINUE
        21 END DO
        CONTINUE
    22 END DO
    DO 24 L=N,1,-1
        IF(INDXR(L) /= INDXC(L))THEN
            DO 23 K=1,N
                DUM=A(K,INDXR(L))
                A(K,INDXR(L))=A(K,INDXC(L))
                A(K,INDXC(L))=DUM
                CONTINUE
            23 END DO
        ENDIF
        CONTINUE
    24 END DO


    call cpu_time(time_end_0)
    
    time_elasped_0 = time_end_0 - time_begin_0
    print *,'==========================================='
    print *,'Inversion takes :',time_elasped_0,' seconds'
    print *,'==========================================='


    RETURN
    END
