subroutine tau_lines1(w,s,n,nlines,tauk,f1,stepj,nstep,nstep0,vgt,nmax)
  implicit none
  integer, intent(in) :: nstep
  integer, intent(in) :: nstep0
  integer, intent(in) :: nlines
  double precision, intent(in) :: stepj
  integer, intent(in) :: nmax
  double precision, intent(in) :: f1
  double precision, dimension(nstep0), intent(out) :: tauk

  double precision, dimension(nlines) :: w
  double precision, dimension(nlines) :: s
  integer :: i
  integer :: l
  integer :: ind
  double precision :: fr
  double precision :: sr
  integer :: nr
  integer :: n1
  integer :: n2
  integer, dimension(nlines) :: n
  double precision, dimension(5,nmax)    :: vgt

  do l = 1, nlines
    fr = w(l)
    sr = s(l)
    nr = n(l)
    ind = nint((fr-f1)/stepj)+1
!~       Label 15: fr <= f2-cutoff
!~       Label 16: fr > f2-cutoff and fr <= f1
!~       Label 17: fr > f1 and <= f2
!~       Label 18: fr > f2 and <= f1+cutoff
!~       Label 19: fr > f1+cutoff
    if(ind>nmax) then
      n1=nstep! 19
      n2=ind-nmax+1
      if(n1>=n2) THEN
        do i=n1,n2,-1
          tauk(i)=tauk(i)+sr*vgt(nr,ind+1-i)! 35
        end do
      else
        return
      end if
      cycle
    end if
    if(ind>nstep) then
      n1=nstep ! 18
      n2=1
      do i = n1, n2, -1
        tauk(i) = tauk(i)+sr*vgt(nr,ind+1-i)! 34
      end do
      cycle
    end if
    if(ind>1) then
      n1=ind ! 17
      n2=nstep
      do i = n1, n2 ! 32
        tauk(i)=tauk(i)+sr*vgt(nr,i+1-ind)
      end do
      n2=1
      do i = n1-1, n2, -1 !33
        tauk(i)=tauk(i)+sr*vgt(nr,ind+1-i)
      end do
      cycle
    end if
    if(ind>nstep-nmax+1) then
      n1=1 ! 16
      n2=nstep
      do i = n1, n2 ! 31
        tauk(i) = tauk(i)+sr*vgt(nr,i+1-ind)
      end do
      cycle
    end if
    n1=1
    n2=ind+nmax-1
    if(n1<=n2) then
      do i = n1, n2
        tauk(i)=tauk(i)+sr*vgt(nr,i+1-ind)
      end do
    end if
  end do ! 1
  return
end
