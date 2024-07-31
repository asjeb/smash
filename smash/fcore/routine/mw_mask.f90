!%      (MW) Module Wrapped.
!%
!%      Subroutine
!%      ----------
!%
!%      - mask_upstream_cells

module mw_mask

    use md_constant, only: sp
    use mwd_mesh, only: MeshDT

contains

    recursive subroutine mask_upstream_cells(mesh, row, col, mask)

        implicit none

        type(MeshDT), intent(in) :: mesh
        integer, intent(in) :: row, col
        logical, dimension(mesh%nrow, mesh%ncol), intent(inout) &
        & :: mask

        integer :: i, row_imd, col_imd
        integer, dimension(8) :: drow = (/1, 1, 0, -1, -1, -1, 0, 1/)
        integer, dimension(8) :: dcol = (/0, -1, -1, -1, 0, 1, 1, 1/)
        integer, dimension(8) :: opposite_dir = (/5, 6, 7, 8, 1, 2, 3, 4/)
        
        mask(row, col) = .true.

        do i = 1, 8

            row_imd = row + drow(i)
            col_imd = col + dcol(i)
            
            if (row_imd .lt. 1 .or. row_imd .gt. mesh%nrow .or. col_imd .lt. 1 .or. col_imd .gt. mesh%ncol) cycle
            
!~             if ( mesh%flwdir(row_imd, col_imd) .eq. i .and. mesh%flwdir(row, col) .eq. opposite_dir(i) ) then
                
!~                 write(*,*) 'A well is detected at row/col',row,col,&
!~                 &' Xcoord=', mesh%xmin+col*mesh%xres, 'Ycoord=', mesh%ymax-row*mesh%yres
                
!~                 cycle
                
!~             end if
            
            if (mesh%flwdir(row_imd, col_imd) .eq. i) call mask_upstream_cells(mesh, row_imd, col_imd, mask)

        end do

    end subroutine mask_upstream_cells

end module mw_mask
