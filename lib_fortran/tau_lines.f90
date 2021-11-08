module tau_lines
  implicit none
  contains
  subroutine get_line_absorption_cross_sections(w_out, s_out, e_out, n_out, nlines, vgt, nvoigt, nfreq_cut_off, &
                                               nbl_sp_max, nstep, nstep0, stepj, f1, f2, hckt, hck296, fdop, dtauk)
      ! """
      ! Get the absorption cross section of a line with Voigt profile.
      ! There must be one Voigt profile per line broadening, and the Voigt profile must be calculated from the
      ! center of the line to the cutoff.
      ! :param line_wavenumber: (cm-1) wavenumber of the line
      ! :param line_intensity: (cm-1/(molecules.cm-2)) intensity of the line
      ! :param i_broadening_line: line broadening index
      ! :param f_voigt: Voigt function of shape (n_broadening, n_wavenumbers_cutoff)
      ! :param wavenumber_min: (cm-1) minimum wavenumber of the spectral interval
      ! :param wavenumber_max: (cm-1) maximum wavenumber of the spectral interval
      ! :param wavenumber_step: (cm-1) step between to wavenumbers of the wavenumber interval
      ! :param n_wavenumbers: number of wavenumbers in the wavenumber interval
      ! :param cutoff: (cm-1) distance from the line center where the line is no more calculated
      ! :param n_wavenumbers_cutoff: number of wavenumbers from the center of the line profile to the cutoff
      ! :return absorption_cross_sections: (cm2) absorption cross section of the line
      ! """
      implicit none

  !    integer, intent(in) :: n_wavenumbers, n_wavenumbers_cutoff, i_broadening_line
  !    integer, dimension(:), intent(in) :: nlines
  !    double precision, intent(in) :: line_wavenumber, line_intensity, wavenumber_min, wavenumber_max, wavenumber_step, &
  !                                    cutoff
  !    double precision, dimension(:, :), intent(in) :: f_voigt
  !    doubleprecision, intent(out) :: absorption_cross_sections(n_wavenumbers)
      integer, intent(in) :: nstep, nstep0, nvoigt, nfreq_cut_off, nlines, nbl_sp_max
      integer, intent(in) :: n_out(nbl_sp_max)
      double precision, intent(in) :: vgt(5,nvoigt), w_out(nbl_sp_max), s_out(nbl_sp_max), e_out(nbl_sp_max),  stepj, f1, f2
      double precision, intent(in) :: hckt, hck296, fdop
      double precision, intent(out) :: dtauk(nstep0)
      double precision :: cut_off, s_over, dt
      integer :: &
          l, &
          nr, &
          wr, &
          er, &
          sr, &
          i, &                ! index
          i_line_center, &    ! index of the center of the line profile in the wavenumber interval
          i_line_min, &       ! index of the left cutoff of the line profile in the wavenumber interval
          i_line_max, &       ! index of the right cutoff of the line profile in the wavenumber interval
          i_min, &            ! minimum index of the line in the wavenumber interval where to calculate the abs.
          i_mid, &            ! index separating the left and right sides of the line profile in the wvn. int.
          i_max               ! maximum index of the line in the wavenumber interval where to calculate the abs.

      cut_off = fdop * stepj * (nfreq_cut_off-1)
      dt = hckt - hck296
      do l = 1, nlines
        wr = w_out(l)
        nr = n_out(l)
        er = e_out(l)
        sr = s_out(l)
        ! Exclude lines too far from the wavenumber interval
        if(wr > f2 + cut_off) then
          cycle
        else if(wr < f1 - cut_off) then
          cycle
        end if

        s_over = sr * dexp(-er*dt) * (1.d0-dexp(-wr*hckt)) / &
                 (1.d0 - dexp(-wr*hck296))

        ! Find the index of the wavenumber interval where the line center is located
        i_line_center = nint((wr - f1) / stepj) + 1

        ! Find the index of the borders of the line profile
        i_line_min = i_line_center - nfreq_cut_off + 1
        i_line_max = i_line_center + nfreq_cut_off - 1

        ! Find the absorption cross section calculation indices
        ! 0 <= i_mid <= n + 1 (and not 1 <= i_mid <= n) because we need to calculate at i = 1 and i = n
        i_min = max(i_line_min, 1)
        i_mid = min(max(i_line_center, 0), nstep + 1)
        i_max = min(i_line_max, nfreq_cut_off)

        ! Left side of the line profile (center of the line not included)
        if(i_line_center > 1) then  ! line_wavenumber > wavenumber_min
          do i = i_min, i_mid - 1
            dtauk(i) = dtauk(i) + s_over * vgt(nr, i_line_center - i + 1)
          end do
        end if

        ! Right side of the line profile (center of the line not included)
        if(i_line_center < nfreq_cut_off) then  ! line_wavenumber < wavenumber_max
          do i = i_mid + 1, i_max
            dtauk(i) = dtauk(i) + s_over * vgt(nr, i - i_line_center + 1)
          end do
        end if

        ! Center of the line
        if(i_line_center >= 1 .and. i_line_center <= nfreq_cut_off) then
          dtauk(i_mid) = dtauk(i_mid) + s_over * vgt(nr, 1)
        end if
      end do

      return
  end subroutine get_line_absorption_cross_sections

end module tau_lines
