subroutine tau_lines1(w,s,n,nlines,tauk,f1,stepj,nstep,nstep0,vgt,nmax)
implicit none
integer :: i = 0
integer :: l = 0
integer :: ind = 0
double precision :: fr = 0.
double precision :: sr = 0.
integer :: nr = 0
integer :: n1 = 0
integer :: n2 = 0
integer, intent(in) :: nstep
integer, intent(in) :: nstep0
integer, intent(in) :: nlines
double precision, intent(in) :: stepj
integer, intent(in) :: nmax
double precision, intent(in) :: f1
double precision, dimension(nlines) :: w
double precision, dimension(nlines) :: s
integer, dimension(nlines) :: n
double precision, dimension(5,nmax)    :: vgt
double precision, dimension(nstep0), intent(out) :: tauk

do 1 l=1,nlines
      fr=w(l)
      sr=s(l)
      nr=n(l)
      ind=idnint((fr-f1)/stepj)+1
!~       Label 15: fr <= f2-cutoff
!~       Label 16: fr > f2-cutoff and fr <= f1
!~       Label 17: fr > f1 and <= f2
!~       Label 18: fr > f2 and <= f1+cutoff
!~       Label 19: fr > f1+cutoff
        if(ind.gt.nmax) goto 19
        if(ind.gt.nstep) goto 18
        if(ind.gt.1) goto 17
        if(ind.gt.nstep-nmax+1) goto 16
      n1=1
      n2=ind+nmax-1
      if(n1.le.n2) THEN
    do 30 i=n1,n2
30  tauk(i)=tauk(i)+sr*vgt(nr,i+1-ind)
        ENDIF
            goto 1
16    n1=1
      n2=nstep
      do 31 i=n1,n2
31    tauk(i)=tauk(i)+sr*vgt(nr,i+1-ind)
            goto 1
17    n1=ind
      n2=nstep
      do 32 i=n1,n2
32    tauk(i)=tauk(i)+sr*vgt(nr,i+1-ind)
      n2=1
      do 33 i=n1-1,n2,-1
33    tauk(i)=tauk(i)+sr*vgt(nr,ind+1-i)
            goto 1
18    n1=nstep
      n2=1
      do 34 i=n1,n2,-1
34    tauk(i)=tauk(i)+sr*vgt(nr,ind+1-i)
        goto 1
19    n1=nstep
      n2=ind-nmax+1
      if(n1.ge.n2) THEN
    do 35 i=n1,n2,-1
35  tauk(i)=tauk(i)+sr*vgt(nr,ind+1-i)
        ELSE
    return
        ENDIF
1     continue
return
end
