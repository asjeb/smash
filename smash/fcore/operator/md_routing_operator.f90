!%      (MD) Module Differentiated.
!%
!%      Subroutine
!%      ----------
!%
!%      - upstream_discharge
!%      - linear_routing
!%      - kinematic_wave1d
!%      - lag0_time_step
!%      - lr_time_step
!%      - kw_time_step

module md_routing_operator

    use md_constant !% only : sp
    use mwd_setup !% only: SetupDT
    use mwd_mesh !% only: MeshDT
    use mwd_input_data !% only: Input_DataDT !% lie au prcp
    use mwd_options !% only: OptionsDT
    use mwd_returns !% only: ReturnsDT
    use mwd_atmos_manipulation !% get_ac_atmos_data_time_step

    implicit none

contains

    subroutine upstream_discharge(mesh, row, col, ac_q, qup)

        implicit none

        type(MeshDT), intent(in) :: mesh
        integer, intent(in) :: row, col
        real(sp), dimension(mesh%nac), intent(in) :: ac_q
        real(sp), intent(out) :: qup

        integer :: i, row_imd, col_imd, k
        integer, dimension(8) :: drow = (/1, 1, 0, -1, -1, -1, 0, 1/)
        integer, dimension(8) :: dcol = (/0, -1, -1, -1, 0, 1, 1, 1/)

        qup = 0._sp

        do i = 1, 8

            row_imd = row + drow(i)
            col_imd = col + dcol(i)

            if (row_imd .lt. 1 .or. row_imd .gt. mesh%nrow .or. col_imd .lt. 1 .or. col_imd .gt. mesh%ncol) cycle
            k = mesh%rowcol_to_ind_ac(row_imd, col_imd)

            if (mesh%flwdir(row_imd, col_imd) .eq. i) qup = qup + ac_q(k)

        end do

    end subroutine upstream_discharge

    subroutine linear_routing(dx, dy, dt, flwacc, llr, hlr, qup, q)

        implicit none

        real(sp), intent(in) :: dx, dy, dt, flwacc
        real(sp), intent(in) :: llr
        real(sp), intent(inout) :: hlr, qup, q

        real(sp) :: hlr_imd

        qup = (qup*dt)/(1e-3_sp*(flwacc - dx*dy))

        hlr_imd = hlr + qup

        hlr = hlr_imd*exp(-dt/(llr*60._sp))

        q = q + (hlr_imd - hlr)*1e-3_sp*(flwacc - dx*dy)/dt

    end subroutine linear_routing

    subroutine kinematic_wave1d(dx, dy, dt, akw, bkw, qlijm1, qlij, qim1j, qijm1, qij)

        implicit none

        real(sp), intent(in) :: dx, dy, dt
        real(sp), intent(in) :: akw, bkw
        real(sp), intent(in) :: qlijm1, qlij, qim1j, qijm1
        real(sp), intent(inout) :: qij

        real(sp) :: wqlijm1, wqlij, wqim1j, wqijm1
        real(sp) :: dtddx, n1, n2, n3, d1, d2, rhs, rsd, rsd_d
        integer :: iter, maxiter

        !% Avoid numerical issues
        wqlijm1 = max(1e-6_sp, qlijm1)
        wqlij = max(1e-6_sp, qlij)
        wqim1j = max(1e-6_sp, qim1j)
        wqijm1 = max(1e-6_sp, qijm1)

        dtddx = dt/dx

        d1 = dtddx
        d2 = akw*bkw*((wqijm1 + wqim1j)/2._sp)**(bkw - 1._sp)

        n1 = dtddx*wqim1j
        n2 = wqijm1*d2
        n3 = dtddx*(wqlijm1 + wqlij)/2._sp

        !% Linearized solution
        qij = (n1 + n2 + n3)/(d1 + d2)

        !% Non-Linear solution solved with Newton-Raphson
        !% Commented while testing Linearized solution

!~         rhs = n1 + akw*wqijm1**bkw + n3

!~         iter = 0
!~         maxiter = 2
!~         rsd = 1._sp

!~         do while (abs(rsd) > 1e-6 .and. iter < maxiter)

!~             rsd = dtddx*qij + akw*qij**bkw - rhs
!~             rsd_d = dtddx + akw*bkw*qij**(bkw - 1._sp)

!~             qij = qij - rsd/rsd_d

!~             qij = max(qij, 0._sp)

!~             iter = iter + 1

!~         end do

    end subroutine kinematic_wave1d

    subroutine lag0_time_step(setup, mesh, options, returns, time_step, ac_qtz, ac_qz)

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(OptionsDT), intent(in) :: options
        type(ReturnsDT), intent(inout) :: returns
        integer, intent(in) :: time_step
        real(sp), dimension(mesh%nac, setup%nqz), intent(in) :: ac_qtz
        real(sp), dimension(mesh%nac, setup%nqz), intent(inout) :: ac_qz

        integer :: i, j, row, col, k, time_step_returns
        real(sp) :: qup

        ac_qz(:, setup%nqz) = ac_qtz(:, setup%nqz)

        ! Skip the first partition because boundary cells are not routed
        do i = 2, mesh%npar

            ! Tapenade does not accept 'IF' condition within OMP directive. Therefore, the routing loop
            ! is duplicated ... Maybe there is another way to do it.
            if (mesh%ncpar(i) .ge. options%comm%ncpu) then
#ifdef _OPENMP
                !$OMP parallel do schedule(static) num_threads(options%comm%ncpu) &
                !$OMP& shared(setup, mesh, ac_qtz, ac_qz, i) &
                !$OMP& private(j, row, col, k, qup)
#endif
                do j = 1, mesh%ncpar(i)

                    row = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 1)
                    col = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 2)
                    k = mesh%rowcol_to_ind_ac(row, col)

                    if (mesh%active_cell(row, col) .eq. 0 .or. mesh%local_active_cell(row, col) .eq. 0) cycle

                    call upstream_discharge(mesh, row, col, ac_qz(:, setup%nqz), qup)

                    ac_qz(k, setup%nqz) = ac_qz(k, setup%nqz) + qup

                    !$AD start-exclude
                    !internal fluxes
                    if (returns%internal_fluxes_flag) then
                        if (allocated(returns%mask_time_step)) then
                            if (returns%mask_time_step(time_step)) then
                                time_step_returns = returns%time_step_to_returns_time_step(time_step)
                                returns%internal_fluxes( &
                                    row, &
                                    col, &
                                    time_step_returns, &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + 1: &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + setup%n_routing_fluxes &
                                    ) = (/qup/)
                            end if
                        end if
                    end if
                    !$AD end-exclude

                end do
#ifdef _OPENMP
                !$OMP end parallel do
#endif
            else

                do j = 1, mesh%ncpar(i)

                    row = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 1)
                    col = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 2)
                    k = mesh%rowcol_to_ind_ac(row, col)

                    if (mesh%active_cell(row, col) .eq. 0 .or. mesh%local_active_cell(row, col) .eq. 0) cycle

                    call upstream_discharge(mesh, row, col, ac_qz(:, setup%nqz), qup)

                    ac_qz(k, setup%nqz) = ac_qz(k, setup%nqz) + qup

                    !$AD start-exclude
                    !internal fluxes
                    if (returns%internal_fluxes_flag) then
                        if (allocated(returns%mask_time_step)) then
                            if (returns%mask_time_step(time_step)) then
                                time_step_returns = returns%time_step_to_returns_time_step(time_step)
                                returns%internal_fluxes( &
                                    row, &
                                    col, &
                                    time_step_returns, &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + 1: &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + setup%n_routing_fluxes &
                                    ) = (/qup/)
                            end if
                        end if
                    end if
                    !$AD end-exclude

                end do

            end if

        end do

    end subroutine lag0_time_step

    subroutine lr_time_step(setup, mesh, options, returns, time_step, ac_qtz, ac_llr, ac_hlr, ac_qz)

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(OptionsDT), intent(in) :: options
        type(ReturnsDT), intent(inout) :: returns
        integer, intent(in) :: time_step
        real(sp), dimension(mesh%nac, setup%nqz), intent(in) :: ac_qtz
        real(sp), dimension(mesh%nac), intent(in) :: ac_llr
        real(sp), dimension(mesh%nac), intent(inout) :: ac_hlr
        real(sp), dimension(mesh%nac, setup%nqz), intent(inout) :: ac_qz

        integer :: i, j, row, col, k, time_step_returns
        real(sp) :: qup

        ac_qz(:, setup%nqz) = ac_qtz(:, setup%nqz)

        ! Skip the first partition because boundary cells are not routed
        do i = 2, mesh%npar

            ! Tapenade does not accept 'IF' condition within OMP directive. Therefore, the routing loop
            ! is duplicated ... Maybe there is another way to do it.
            if (mesh%ncpar(i) .ge. options%comm%ncpu) then
#ifdef _OPENMP
                !$OMP parallel do schedule(static) num_threads(options%comm%ncpu) &
                !$OMP& shared(setup, mesh, ac_qtz, ac_llr, ac_hlr, ac_qz, i) &
                !$OMP& private(j, row, col, k, qup)
#endif
                do j = 1, mesh%ncpar(i)

                    row = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 1)
                    col = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 2)
                    k = mesh%rowcol_to_ind_ac(row, col)

                    if (mesh%active_cell(row, col) .eq. 0 .or. mesh%local_active_cell(row, col) .eq. 0) cycle

                    call upstream_discharge(mesh, row, col, ac_qz(:, setup%nqz), qup)

                    call linear_routing(mesh%dx(row, col), mesh%dy(row, col), setup%dt, mesh%flwacc(row, col), &
                    & ac_llr(k), ac_hlr(k), qup, ac_qz(k, setup%nqz))

                    !$AD start-exclude
                    !internal fluxes
                    if (returns%internal_fluxes_flag) then
                        if (allocated(returns%mask_time_step)) then
                            if (returns%mask_time_step(time_step)) then
                                time_step_returns = returns%time_step_to_returns_time_step(time_step)
                                returns%internal_fluxes( &
                                    row, &
                                    col, &
                                    time_step_returns, &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + 1: &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + setup%n_routing_fluxes &
                                    ) = (/qup/)
                            end if
                        end if
                    end if
                    !$AD end-exclude

                end do
#ifdef _OPENMP
                !$OMP end parallel do
#endif
            else

                do j = 1, mesh%ncpar(i)

                    row = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 1)
                    col = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 2)
                    k = mesh%rowcol_to_ind_ac(row, col)

                    if (mesh%active_cell(row, col) .eq. 0 .or. mesh%local_active_cell(row, col) .eq. 0) cycle

                    call upstream_discharge(mesh, row, col, ac_qz(:, setup%nqz), qup)

                    call linear_routing(mesh%dx(row, col), mesh%dy(row, col), setup%dt, mesh%flwacc(row, col), &
                    & ac_llr(k), ac_hlr(k), qup, ac_qz(k, setup%nqz))

                    !$AD start-exclude
                    !internal fluxes
                    if (returns%internal_fluxes_flag) then
                        if (allocated(returns%mask_time_step)) then
                            if (returns%mask_time_step(time_step)) then
                                time_step_returns = returns%time_step_to_returns_time_step(time_step)
                                returns%internal_fluxes( &
                                    row, &
                                    col, &
                                    time_step_returns, &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + 1: &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + setup%n_routing_fluxes &
                                    ) = (/qup/)
                            end if
                        end if
                    end if
                    !$AD end-exclude

                end do

            end if

        end do

    end subroutine lr_time_step

    subroutine kw_time_step(setup, mesh, options, returns, t, ac_qtz, ac_akw, ac_bkw, ac_qz)
        ! Bug fixed on the parallel adjoint linked to variable names by changing of variable names :
        ! Added local variable t instead of time_step as input to kw_time_step -> no change
        ! Added local variable t_returns instead of time_step_returns to kw_time_step
        ! -> modified openMP, unmodified non-openMP

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(OptionsDT), intent(in) :: options
        type(ReturnsDT), intent(inout) :: returns
        integer, intent(in) :: t
        real(sp), dimension(mesh%nac, setup%nqz), intent(in) :: ac_qtz
        real(sp), dimension(mesh%nac), intent(in) :: ac_akw, ac_bkw
        real(sp), dimension(mesh%nac, setup%nqz), intent(inout) :: ac_qz

        integer :: i, j, row, col, k, t_returns
        real(sp) :: qlijm1, qlij, qim1j, qijm1

        ac_qz(:, setup%nqz) = ac_qtz(:, setup%nqz)

        ! Skip the first partition because boundary cells are not routed
        do i = 2, mesh%npar

            ! Tapenade does not accept 'IF' condition within OMP directive. Therefore, the routing loop
            ! is duplicated ... Maybe there is another way to do it.
            if (mesh%ncpar(i) .ge. options%comm%ncpu) then
#ifdef _OPENMP
                !$OMP parallel do schedule(static) num_threads(options%comm%ncpu) &
                !$OMP& shared(setup, mesh, ac_qtz, ac_akw, ac_bkw, ac_qz, i) &
                !$OMP& private(j, row, col, k, qlijm1, qlij, qim1j, qijm1)
#endif
                do j = 1, mesh%ncpar(i)

                    row = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 1)
                    col = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 2)
                    k = mesh%rowcol_to_ind_ac(row, col)

                    if (mesh%active_cell(row, col) .eq. 0 .or. mesh%local_active_cell(row, col) .eq. 0) cycle

                    qlijm1 = ac_qtz(k, setup%nqz - 1)
                    qlij = ac_qtz(k, setup%nqz)
                    qijm1 = ac_qz(k, setup%nqz - 1)

                    call upstream_discharge(mesh, row, col, ac_qz(:, setup%nqz), qim1j)

                    call kinematic_wave1d(mesh%dx(row, col), mesh%dy(row, col), setup%dt, &
                    & ac_akw(k), ac_bkw(k), qlijm1, qlij, qim1j, qijm1, ac_qz(k, setup%nqz))

                    !$AD start-exclude
                    !internal fluxes
                    if (returns%internal_fluxes_flag) then
                        if (allocated(returns%mask_time_step)) then
                            if (returns%mask_time_step(t)) then
                                t_returns = returns%time_step_to_returns_time_step(t)
                                returns%internal_fluxes( &
                                    row, &
                                    col, &
                                    t_returns, &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + 1: &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + setup%n_routing_fluxes &
                                    ) = (/qim1j/)
                            end if
                        end if
                    end if
                    !$AD end-exclude
                end do
#ifdef _OPENMP
                !$OMP end parallel do
#endif
            else

                do j = 1, mesh%ncpar(i)

                    row = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 1)
                    col = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 2)
                    k = mesh%rowcol_to_ind_ac(row, col)

                    if (mesh%active_cell(row, col) .eq. 0 .or. mesh%local_active_cell(row, col) .eq. 0) cycle

                    qlijm1 = ac_qtz(k, setup%nqz - 1)
                    qlij = ac_qtz(k, setup%nqz)
                    qijm1 = ac_qz(k, setup%nqz - 1)

                    call upstream_discharge(mesh, row, col, ac_qz(:, setup%nqz), qim1j)

                    call kinematic_wave1d(mesh%dx(row, col), mesh%dy(row, col), setup%dt, &
                    & ac_akw(k), ac_bkw(k), qlijm1, qlij, qim1j, qijm1, ac_qz(k, setup%nqz))

                    !$AD start-exclude
                    !internal fluxes
                    if (returns%internal_fluxes_flag) then
                        if (allocated(returns%mask_time_step)) then
                            if (returns%mask_time_step(t)) then
                                t_returns = returns%time_step_to_returns_time_step(t)
                                returns%internal_fluxes( &
                                    row, &
                                    col, &
                                    t_returns, &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + 1: &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + setup%n_routing_fluxes &
                                    ) = (/qim1j/)
                            end if
                        end if
                    end if
                    !$AD end-exclude
                end do

            end if

        end do

    end subroutine kw_time_step

    subroutine apply_simple_canal(nrow, ncol, h, qx, qy, c)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        integer, intent(in) :: c

        !                            2
        !  |------------------------------------------------------|
        ! 1|                                                      |3
        !  |------------------------------------------------------|
        !                            4
        
        ! 1
        if (c .lt. 10) then
            h(:, 1) = 2._sp
        else
            h(:, 1) = 1._sp
        end if
        ! 2
        qy(1, :) = 0._sp

        ! 3
        qx(:, ncol+1) = 0._sp

        ! 4
        qy(nrow+1, :) = 0._sp

    end subroutine apply_simple_canal

    subroutine initial_macdonal(nrow, ncol, h, qx, qy, zb)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        real(sp), dimension(nrow, ncol), intent(in) :: zb
        
        h = 0._sp

        qx = 0._sp
        qy = 0._sp
    end subroutine initial_macdonal

    subroutine bc_mac_donald(nrow, ncol, h, qx, qy, zb, c)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        real(sp), dimension(nrow, ncol), intent(in) :: zb
        integer, intent(in) :: c

        !                            2
        !  |------------------------------------------------------|
        ! 1|                                                      |3
        !  |------------------------------------------------------|
        !                            4
        
        ! 1
        ! h(:, 1) = (4._sp / gravity) ** (1._sp / 3._sp) * (1 + 0.5 * exp(-4._sp))
        qx(:, 1) = 2._sp

        ! 2
        qy(1, :) = 0._sp

        ! 3
        h(:, ncol) = (4._sp / gravity) ** (1._sp / 3._sp) * (1 + 0.5 * exp(-4._sp))

        ! 4
        qy(nrow+1, :) = 0._sp

    end subroutine bc_mac_donald
    
    subroutine initial_bump(nrow, ncol, h, qx, qy, zb)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        real(sp), dimension(nrow, ncol), intent(in) :: zb
        
        h = 2._sp - zb

        qx = 4.42_sp
        qy = 0._sp
    end subroutine initial_bump
    
    subroutine apply_bump_bc(nrow, ncol, h, qx, qy, zb, c)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        real(sp), dimension(nrow, ncol), intent(in) :: zb
        integer, intent(in) :: c

        !                            2
        !  |------------------------------------------------------|
        ! 1|                                                      |3
        !  |------------------------------------------------------|
        !                            4
        
        ! 1
        qx(:, 1) = 4.42_sp

        ! 2
        qy(1, :) = 0._sp

        ! 3
        h(:, ncol) = 2._sp

        ! 4
        qy(nrow+1, :) = 0._sp

    end subroutine apply_bump_bc



    subroutine initial_lake_at_rest_immersive_bump(nrow, ncol, h, qx, qy, zb)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        real(sp), dimension(nrow, ncol), intent(in) :: zb
        
        h = 2._sp - zb
        qx = 0._sp
        qy = 0._sp
    end subroutine initial_lake_at_rest_immersive_bump

    subroutine initial_lake_at_rest_emersive_bump(nrow, ncol, h, qx, qy, zb)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        real(sp), dimension(nrow, ncol), intent(in) :: zb
       
        h(:, :) = max(0.1_sp, zb(:, :)) - zb(:, :)
        qx = 0._sp
        qy = 0._sp
    end subroutine initial_lake_at_rest_emersive_bump


    subroutine bc_wall(nrow, ncol, h, qx, qy, c)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        integer, intent(in) :: c

        !                            2
        !  |------------------------------------------------------|
        ! 1|                                                      |3
        !  |------------------------------------------------------|
        !                            4
        
        ! 1
        qx(:, 1) = 0._sp

        ! 2
        qy(1, :) = 0._sp

        ! 3
        qx(:, ncol+1) = 0._sp

        ! 4
        qy(nrow+1, :) = 0._sp

    end subroutine bc_wall

    ! subroutine free_outflow(h, qx, qy, nrow, ncol)
        ! implicit none
        ! real(sp), dimension(nrow, ncol) :: h
        ! real(sp), dimension(nrow, ncol+1) :: qx
        ! real(sp), dimension(nrow+1, ncol) :: qy
        ! do col = 1, mesh%ncol
                
        !     ! up
        !     ! dzb = zb(1, col) - zb(2, col)
        !     ! qy(1, col) = hsw(1, col) ** (5._sp / 3._sp) * &
        !     !     sqrt(abs(dzb) / mesh%dy(1, col)) / manning(1, col)
        !     qy(1, col) = 0._sp ! wall

        !     ! down
        !     ! dzb = zb(mesh%nrow-1, col) - zb(mesh%nrow, col)
        !     ! qy(mesh%nrow+1, col) = hsw(mesh%nrow, col) ** (5._sp / 3._sp) * &
        !     !     sqrt(abs(dzb) / mesh%dy(mesh%nrow, col)) / manning(mesh%nrow, col)
        !     qy(mesh%nrow-1, col) = 0._sp ! wall

        ! end do

        ! do row = 1, mesh%nrow
                            
        !     ! ! left
        !     ! dzb = zb(row, 1) - zb(row, 2)
        !     ! qx(row, 1) = hsw(row, 1) ** (5._sp / 3._sp) * &
        !     !     sqrt(abs(dzb)) / mesh%dx(row, 1) / manning(row, 1)
        !     qx(row, 1) = 2._sp ! qin 
        !     qy(row, 1) = 0._sp ! qin

            ! right
            ! dzb = zb(row, mesh%ncol-1) - zb(row, mesh%ncol)
            ! qx(row, mesh%ncol+1) = hsw(row, mesh%ncol) ** (5._sp / 3._sp) * &
            !     sqrt(abs(dzb)) / mesh%dx(row, mesh%ncol) / manning(row, mesh%ncol)
            ! qx(row, mesh%ncol+1) = 0._sp ! wall
        ! end do
    ! do row=1, mesh%nrow
    !     hsw(row, mesh%ncol) = (4._sp / gravity) ** (1._sp / 3._sp) * (1 + 0.5 * exp(-4._sp))
    ! end do
    ! end subroutine free_outflow


    subroutine shallow_water_2d_time_step(setup, mesh, input_data, options, returns, &
            time_step, ac_qtz, zb, manning, ac_qz)
        !hsw, qx, qy a mettre dans les flux internes
        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(OptionsDT), intent(in) :: options
        type(ReturnsDT), intent(inout) :: returns
        integer, intent(in) :: time_step
        real(sp), dimension(mesh%nac, setup%nqz), intent(in) :: ac_qtz
        real(sp), dimension(mesh%nrow, mesh%ncol), intent(in) :: zb, manning
        real(sp), dimension(mesh%nac, setup%nqz), intent(inout) :: ac_qz
        ! real(sp), dimension(:, :, :), allocatable :: hsw_t, eta_t, qx_t, qy_t
        real(sp), dimension(mesh%nrow, mesh%ncol) :: hsw, eta
        real(sp), dimension(mesh%nrow, mesh%ncol+1) :: qx
        real(sp), dimension(mesh%nrow+1, mesh%ncol) :: qy
        real(sp), dimension(mesh%nac) :: ac_prcp

        integer :: i, j, row, col, k, time_returns, nt_sw, nt_sw_old, ctt
        real(sp) :: t, dt
        real(sp) :: heps, hfx, hfy, maxhsw, hxm, hxp, hym, hyp, dzb

        !$AD start-exclude

        ! initialisation
        t = (time_step - 1) * setup%dt
        
        call initial_macdonal(mesh%nrow, mesh%ncol, hsw, qx, qy, zb)
        
        
        heps = 1.e-6
        ctt = 1

        ! print *, manning
        
        do while (t .lt. time_step * setup%dt .and. ctt .le. returns%nt_sw) 
            eta = zb + hsw
            ! print *, eta 
            ! write(*,*) "qx = "
            ! do i = 1, size(qx, 1)
            !     write(*,"(*(f8.4))") qx(i,:)
            ! end do
            ! write(*,*) "h = "
            ! do i = 1, size(hsw, 1)
            !     write(*,"(*(f8.4))") hsw(i,:)
            ! end do
            returns%sw2d(:, :, ctt, 1) = hsw(:, :)
            returns%sw2d(:, :, ctt, 2) = eta(:, :)
            returns%sw2d(:, :, ctt, 3) = qx(:, :)
            returns%sw2d(:, :, ctt, 4) = qy(:, :)

            maxhsw = max(heps, maxval(hsw))
            ! print *, maxhsw
            ! print *, mesh%dx(0,0), mesh%dy(0,0)
            dt = 0.5 * min(minval(mesh%dx), minval(mesh%dy)) &
            / sqrt(gravity * maxhsw)

            ! print *, eta(1, 1)
            ! print *, hsw(1, 1)
            ! print *, manning(1, 1)
            ! update fluxes 
            do row = 1, mesh%nrow
                do col = 2, mesh%ncol

                    hfx = max(heps, max(eta(row, col-1), eta(row, col)) - max(zb(row, col-1), zb(row, col)))
                    ! print *, row, col, eta(row, col), eta(row, col-1), eta(row, col) - eta(row, col-1)
                    qx(row, col) = (qx(row, col) - dt * gravity * hfx * &
                        (eta(row, col) - eta(row, col-1)) / mesh%dx(row, col)) / &
                        (1 + dt * gravity * manning(row, col) ** 2 * abs(qx(row, col)) &
                        / hfx ** (7._sp / 3._sp))
                end do
            end do
            ! print *, "POUET"
            do row = 2, mesh%nrow
                do col = 1, mesh%ncol

                    hfy = max(heps, max(eta(row-1, col), eta(row, col)) - max(zb(row-1, col), zb(row, col)))
                    qy(row, col) = (qy(row, col) - dt * gravity * hfy * &
                            (eta(row, col) - eta(row-1, col)) / mesh%dy(row, col)) / &
                            (1 + dt * gravity * manning(row, col) ** 2 * abs(qy(row, col)) &
                            / hfy ** (7._sp / 3._sp))

                end do
            end do

            
            ! update water height
            do row = 1, mesh%nrow
                do col = 1, mesh%ncol
                    k = mesh%rowcol_to_ind_ac(row, col)
                            
                    hsw(row, col) = hsw(row, col) + dt / mesh%dx(row, col) * (qx(row, col) - qx(row, col+1)) & 
                    + dt / mesh%dy(row, col) * (qy(row, col) - qy(row+1, col)) 

                    ! hsw(row, col) = max(0._sp, hsw(row, col))
                    ! if (t .eq. ((time_step - 1) * setup%dt)) then
                    !     hsw(row, col) = hsw(row, col) + dt * ac_qtz(k, setup%nqz) &
                    !         / mesh%dx(row, col) / mesh%dy(row, col)
                    ! end if
                    
                end do
            end do

            call bc_mac_donald(mesh%nrow, mesh%ncol, hsw, qx, qy, zb, ctt)

            ! call apply_simple_canal(mesh%nrow, mesh%ncol, hsw, qx, qy, ctt)
            
            ! call apply_mac_donald(mesh%nrow, mesh%ncol, hsw, qx, qy, ctt)

            ! call bc_wall(mesh%nrow, mesh%ncol, hsw, qx, qy, ctt)
            ! print *, dt
            ! ! write(*,*) "qx_t = "
            ! ! do i = 1, size(returns%qx_t, 1)
            ! !     write(*,"(*(f8.4))") returns%qx_t(i, :, ctt)
            ! ! end do 
            
            ! ! write(*,*) "qy_t = "
            ! ! do i = 1, size(returns%qy_t, 1)
            ! !     write(*,"(*(f8.4))") returns%qy_t(i, :, ctt)
            ! ! end do 
            ! write(*,*) "qy = "
            ! do i = 1, size(qy, 1)
            !     write(*,"(*(f8.4))") qy(i,:)
            ! end do 
            ! ! write(*,*) "hsw_t = "
            ! ! do i = 1, size(returns%hsw_t, 1)
            ! !     write(*,"(*(f8.4))") returns%hsw_t(i, :, ctt)
            ! ! end do 
            ! ! if (ctt .lt. 38 .and. ctt .gt. 28) then
            ! !     print *, ctt 
            ! write(*,*) "hsw = "
            ! do i = 1, size(hsw, 1)
            !     write(*,"(*(f8.4))") hsw(i,:)
            ! end do 
            ! ! end if
            ctt = ctt + 1
            t = t + dt
        end do 

        !update volume discharge
        do row = 1, mesh%nrow
            do col = 1, mesh%ncol
                k = mesh%rowcol_to_ind_ac(row, col)
                ac_qz(k, setup%nqz) = hsw(row, col) * mesh%dx(row, col) &
                    * mesh%dy(row, col) /setup%dt ! voir les unites
            end do
        end do
        
        !$AD end-exclude
    end subroutine shallow_water_2d_time_step

end module md_routing_operator
