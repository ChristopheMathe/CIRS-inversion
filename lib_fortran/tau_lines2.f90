subroutine tau_lines2(w,s,n,nlines,tauk,f1,stepj,nstep,nstep0,vgt,nmax)
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

do l = 1, nlines ! 1
  fr=w(l)
  sr=s(l)
  nr=n(l)
  ind=idnint((fr-f1)/stepj)+1
!~       Label 15: fr <= f1
!~       Label 16: fr > f1 and fr < f2-cutoff
!~       Label 17: fr >= f2-cutoff and <= f1+cutoff
!~       Label 18: fr > f1+cutoff and <= f2
!~       Label 19: fr > f2
  if(ind>nstep) then ! 19
      n1 = nstep
      n2 = ind - nmax + 1
      if(n1>=n2) then
      do i = n1, n2, -1
        tauk(i) = tauk(i) + sr*vgt(nr,ind+1-i)
      end do
    else
      return
    end if
    cycle
  end if
  if(ind>nmax) then ! 18
    n1 = ind
    n2 = nstep
    do i = n1, n2
      tauk(i) = tauk(i) + sr*vgt(nr,i+1-ind)
    end do
    n2 = ind - nmax + 1
    do i = n1-1, n2, -1
      tauk(i) = tauk(i) + sr*vgt(nr,ind+1-i)
    end do
    cycle
  end if
  if(ind>=(nstep-nmax+1)) then ! 17
    n1 = ind
    n2 = nstep
    do i = n1, n2
      tauk(i) = tauk(i) + sr * vgt(nr,i+1-ind)
    end do
    n2=1
    do i = n1-1, n2, -1
      tauk(i)=tauk(i)+sr*vgt(nr,ind+1-i)
    end do
    cycle
  end if
  if(ind>1) then! 16
    n1 = ind
    n2 = ind + nmax - 1
    do i = n1, n2
      tauk(i) = tauk(i) + sr*vgt(nr,i+1-ind)
    end do
    n2=1
    do i = n1-1, n2, -1
      tauk(i)=tauk(i)+sr*vgt(nr,ind+1-i)
    end do
    cycle
  end if
  n1=1
  n2=ind+nmax-1
  if(n1.le.n2) THEN
    do i = n1, n2
      tauk(i)=tauk(i)+sr*vgt(nr,i+1-ind)
    end do
  end if
end do
return
end

