subroutine tau_lines3(w,s,n,nlines,tauk,f1,stepj,nstep,nstep0,vgt,nmax)
implicit none
integer :: i
integer :: l
integer :: ind
double precision :: fr
double precision :: sr
integer :: nr
integer :: n1
integer :: n2
integer, intent(in) :: nstep
integer, intent(in) :: nstep0
integer, intent(in) :: nlines
double precision, intent(in) :: stepj
integer, intent(in) :: nmax
double precision, intent(in) :: f1
double precision, dimension(nlines) :: w
double precision, dimension(nlines) :: s
integer,dimension(nlines) :: n
double precision, dimension(5,nmax)    :: vgt
double precision, dimension(nstep0), intent(out) :: tauk

do 1 l=1,nlines
      fr=w(l)
      sr=s(l)
      nr=n(l)
      ind=idnint((fr-f1)/stepj)+1
!~    Label 16: fr > f1 and fr <= f1+cutoff
!~    Label 17: fr > f1+cutoff and <= f2-cutoff
!~    Label 18: fr > f2-cutoff and <= f2
!~    Label 19: fr > f2
      if(ind > nstep) goto 19
      if(ind > (nstep-nmax+1)) goto 18
      if(ind > nmax) goto 17
      if(ind > 1) goto 16
        n1=1
        n2=ind+nmax-1
        if(n1 <= n2) THEN
            do 30 i=n1,n2
                tauk(i)=tauk(i)+sr*vgt(nr,i+1-ind)
            30 END DO
        ENDIF
        goto 1
        16 n1=ind
        n2=ind+nmax-1
        do 31 i=n1,n2
            tauk(i)=tauk(i)+sr*vgt(nr,i+1-ind)
        31 END DO
        n2=1
        do 32 i=n1-1,n2,-1
            tauk(i)=tauk(i)+sr*vgt(nr,ind+1-i)
        32 END DO
        goto 1
        17 n1=ind
        n2=ind+nmax-1
        do 33 i=n1,n2
            tauk(i)=tauk(i)+sr*vgt(nr,i+1-ind)
        33 END DO
        n2=ind-nmax+1
        do 34 i=n1-1,n2,-1
            tauk(i)=tauk(i)+sr*vgt(nr,ind+1-i)
        34 END DO
        goto 1
        18 n1=ind
        n2=nstep
        do 35 i=n1,n2
            tauk(i)=tauk(i)+sr*vgt(nr,i+1-ind)
        35 END DO
        n2=ind-nmax+1
        do 36 i=n1-1,n2,-1
            tauk(i)=tauk(i)+sr*vgt(nr,ind+1-i)
        36 END DO
        goto 1
        19 n1=nstep
        n2=ind-nmax+1
        if(n1 >= n2) THEN
            do 37 i=n1,n2,-1
                tauk(i)=tauk(i)+sr*vgt(nr,ind+1-i)
            37 END DO
        ELSE
            return
        ENDIF
    1 END DO
return
end subroutine tau_lines3
