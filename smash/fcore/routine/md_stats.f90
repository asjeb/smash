!%      (MD) Module Differentiated.
!%
!%      Subroutine
!%      ----------
!%
!%      - heap_sort
!%
!%      Function
!%      --------
!%
!%      - quantile1d_r

module md_stats

    use md_constant !% only: sp

    implicit none

contains

    subroutine heap_sort(n, arr)

        !% Notes
        !% -----
        !%
        !% Implement heap sort algorithm
        !%
        !% Computational complexity is O(n log n)

        implicit none

        integer, intent(in) :: n
        real(sp), dimension(n), intent(inout) :: arr

        integer :: l, ir, i, j
        real(sp) :: arr_l

        l = n/2 + 1

        ir = n

10      continue

        if (l .gt. 1) then

            l = l - 1

            arr_l = arr(l)

        else

            arr_l = arr(ir)

            arr(ir) = arr(1)

            ir = ir - 1

            if (ir .eq. 1) then

                arr(1) = arr_l

                return

            end if

        end if

        i = l

        j = l + l

20      if (j .le. ir) then

            if (j .lt. ir) then

                if (arr(j) .lt. arr(j + 1)) j = j + 1

            end if

            if (arr_l .lt. arr(j)) then

                arr(i) = arr(j)

                i = j; j = j + j

            else

                j = ir + 1

            end if

            goto 20

        end if

        arr(i) = arr_l

        goto 10

    end subroutine heap_sort

    function quantile1d_r(dat, p) result(res)

        !% Notes
        !% -----
        !%
        !% Quantile function for real 1d array using linear interpolation
        !%
        !% Similar to numpy.quantile

        implicit none

        real(sp), intent(in) :: p
        real(sp), dimension(:), intent(in) :: dat
        real(sp), dimension(size(dat)) :: sorted_dat
        integer :: n
        real(sp) :: res, q1, q2, frac

        res = dat(1)

        n = size(dat)

        if (n .gt. 1) then

            sorted_dat = dat

            call heap_sort(n, sorted_dat)

            frac = (n - 1)*p + 1

            if (frac .le. 1) then

                res = sorted_dat(1)

            else if (frac .ge. n) then

                res = sorted_dat(n)

            else
                q1 = sorted_dat(int(frac))

                q2 = sorted_dat(int(frac) + 1)

                res = q1 + (q2 - q1)*(frac - int(frac)) ! linear interpolation

            end if

        end if

    end function quantile1d_r

end module md_stats
