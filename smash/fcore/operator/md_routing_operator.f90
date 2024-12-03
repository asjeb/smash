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

    subroutine initial_zero(nrow, ncol, h, qx, qy, zb)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        real(sp), dimension(nrow, ncol), intent(in) :: zb
        
        h = 0._sp

        qx = 0._sp
        qy = 0._sp
    end subroutine initial_zero

    subroutine bc_simple_wave(nrow, ncol, h, qx, qy, c)
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
            h(:, 1) = 0.3_sp
        else
            h(:, 1) = 0._sp
        end if
        ! 2
        qy(1, :) = 0._sp

        ! 3
        qx(:, ncol+1) = 0._sp

        ! 4
        qy(nrow+1, :) = 0._sp

    end subroutine bc_simple_wave



    subroutine bc_gaussian_wave(nrow, ncol, h, qx, qy, time)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        real(sp), intent(in) :: time
        real(sp) :: sigma, pi

        !                            2
        !  |------------------------------------------------------|
        ! 1|                                                      |3
        !  |------------------------------------------------------|
        !                            4

        sigma = sqrt(0.2_sp)
        pi = 3.14_sp

        ! 1
        h(:, 1) = 1._sp / sigma / sqrt(2 * pi) * exp(- time  ** 2 / 2._sp / sigma ** 2)

        ! 2
        qy(1, :) = 0._sp

        ! 3
        qx(:, ncol+1) = 0._sp

        ! 4
        qy(nrow+1, :) = 0._sp

    end subroutine bc_gaussian_wave



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
        qx(:, 1) = 2._sp

        ! 2
        qy(1, :) = 0._sp

        ! 3
        h(:, ncol) = (4._sp / gravity) ** (1._sp / 3._sp) * (1 + 0.5 * exp(-4._sp))

        ! 4
        qy(nrow+1, :) = 0._sp

    end subroutine bc_mac_donald


    subroutine bc_mac_donald_rain(nrow, ncol, h, qx, qy, c)
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
        qx(:, 1) = 1._sp

        ! 2
        qy(1, :) = 0._sp

        ! 3
        h(:, ncol) = (4._sp / gravity) ** (1._sp / 3._sp) * (1 + 0.5 * exp(-4._sp))

        ! 4
        qy(nrow+1, :) = 0._sp

    end subroutine bc_mac_donald_rain
    
    
    
    subroutine macdonald_rainfall(nrow, ncol, h, dt)
        implicit none

        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), intent(in) :: dt
        
        h = h + dt * 0.001_sp

    end subroutine macdonald_rainfall
    
    
    
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


    subroutine initial_bump_drain(nrow, ncol, h, qx, qy, zb)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        real(sp), dimension(nrow, ncol), intent(in) :: zb
       
        h = 0.5_sp - zb
        qx = 0._sp
        qy = 0._sp
    end subroutine initial_bump_drain


    subroutine initial_wet_dambreak(nrow, ncol, h, qx, qy, zb)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        real(sp), dimension(nrow, ncol), intent(in) :: zb
        real(sp) :: hl, hr
        integer :: col

        hl = 0.005_sp
        hr = 0.001_sp
        
        h = hl
        do col = ncol/2, ncol  
            h(:, col) = hr
        end do

        qx = 0._sp
        qy = 0._sp
    end subroutine initial_wet_dambreak


    subroutine initial_dry_dambreak(nrow, ncol, h, qx, qy, zb)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        real(sp), dimension(nrow, ncol), intent(in) :: zb
        real(sp) :: hl, hr
        integer :: col

        hl = 0.005_sp
        hr = 0._sp
        
        h = hl
        do col = ncol/2, ncol  
            h(:, col) = hr
        end do

        qx = 0._sp
        qy = 0._sp
    end subroutine initial_dry_dambreak
    

    subroutine bc_bump_drain(nrow, ncol, h, qx, qy, c)
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
        qx(:, ncol+1) = 0.1

        ! 4
        qy(nrow+1, :) = 0._sp

    end subroutine bc_bump_drain


    subroutine bc_wet_dambreak(nrow, ncol, h, qx, qy)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        real(sp) :: hl, hr

        !                            2
        !  |------------------------------------------------------|
        ! 1|                                                      |3
        !  |------------------------------------------------------|
        !                            4
        
        hl = 0.005_sp
        hr = 0.001_sp

        ! 1
        h(:, 1) = hl

        ! 2
        qy(1, :) = 0._sp

        ! 3
        h(:, ncol) = hr

        ! 4
        qy(nrow+1, :) = 0._sp

    end subroutine bc_wet_dambreak


    subroutine bc_dry_dambreak(nrow, ncol, h, qx, qy)
        implicit none
        integer, intent(in) :: nrow, ncol
        real(sp), dimension(nrow, ncol), intent(inout) :: h
        real(sp), dimension(nrow, ncol+1), intent(inout) :: qx
        real(sp), dimension(nrow+1, ncol), intent(inout) :: qy
        real(sp) :: hl, hr

        !                            2
        !  |------------------------------------------------------|
        ! 1|                                                      |3
        !  |------------------------------------------------------|
        !                            4
        
        hl = 0.005_sp
        hr = 0._sp

        ! 1
        h(:, 1) = hl

        ! 2
        qy(1, :) = 0._sp

        ! 3
        h(:, ncol) = hr

        ! 4
        qy(nrow+1, :) = 0._sp

    end subroutine bc_dry_dambreak


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


    subroutine free_outflow(mesh, hsw, zb, qx, qy, manning, dt)
        implicit none

        type(MeshDT), intent(in) :: mesh
        real(sp), dimension(mesh%nrow, mesh%ncol), intent(in) :: hsw, zb
        real(sp), dimension(mesh%nrow, mesh%ncol+1), intent(inout) :: qx
        real(sp), dimension(mesh%nrow+1, mesh%ncol), intent(inout) :: qy
        real(sp), dimension(mesh%nrow, mesh%ncol), intent(in) :: manning
        real(sp), intent(in) :: dt

        integer :: row, col
        real (sp) :: sgn
        real (sp) :: slope
        real (sp) :: qout, volume_out


        do col = 1, mesh%ncol
                
            ! up
            slope = (zb(1, col) - zb(2, col)) /  mesh%dy(1, col)
            
            if (slope .gt. 0) then
                sgn = 1._sp
            else
                sgn = -1._sp
            end if

            qy(1, col) = sgn * hsw(1, col) ** (5._sp / 3._sp) * &
                sqrt(abs(slope)) / manning(1, col)
            
                
            qout = qout + sgn * qy(1, col) 

            volume_out = volume_out + qout * mesh%dx(1, col) * dt

            ! down
            slope = (zb(mesh%nrow, col) - zb(mesh%nrow-1, col)) / mesh%dy(mesh%nrow, col)

            if (slope .gt. 0) then
                sgn = 1._sp
            else
                sgn = -1._sp
            end if

            qy(mesh%nrow+1, col) = sgn * hsw(mesh%nrow, col) ** (5._sp / 3._sp) * &
                sqrt(abs(slope)) / manning(mesh%nrow, col)
            

            qout = qout + sgn * qy(mesh%nrow+1, col)

            volume_out = volume_out + qout * mesh%dx(mesh%nrow, col) * dt

        end do



        do row = 1, mesh%nrow
                            
            ! left
            slope = (zb(row, 2) - zb(row, 1)) / mesh%dx(row, 1)

            if (slope .gt. 0) then
                sgn = 1
            else
                sgn = -1
            end if
            
            qx(row, 1) = sgn * hsw(row, 1) ** (5._sp / 3._sp) * &
                sqrt(abs(slope)) / manning(row, 1)


            qout = qout + sgn * qx(row, 1)

            volume_out = volume_out + qout * mesh%dy(row, 1) * dt

            ! right
            slope = (zb(row, mesh%ncol-1) - zb(row, mesh%ncol)) / mesh%dx(row, mesh%ncol)

            if (slope .gt. 0) then
                sgn = 1
            else
                sgn = -1
            end if

            qx(row, mesh%ncol+1) = hsw(row, mesh%ncol) ** (5._sp / 3._sp) * &
                sqrt(abs(slope)) / manning(row, mesh%ncol)


            qout = qout + sgn * qx(row, mesh%ncol+1)

            volume_out = volume_out + qout * mesh%dy(row, mesh%ncol) * dt
            
        end do
    end subroutine free_outflow



    subroutine bc_thacker2d(nrow, ncol, h, qx, qy, c)
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
        h(:, 1) = 0._sp

        ! 2
        h(1, :) = 0._sp

        ! 3
        h(:, ncol) = 0._sp

        ! 4
        h(nrow, :) = 0._sp

    end subroutine bc_thacker2d



    subroutine initial_thacker2d(mesh, h, qx, qy, zb)
        implicit none
        type(MeshDT), intent(in) :: mesh
        real(sp), dimension(mesh%nrow, mesh%ncol), intent(inout) :: h
        real(sp), dimension(mesh%nrow, mesh%ncol+1), intent(inout) :: qx
        real(sp), dimension(mesh%nrow+1, mesh%ncol), intent(inout) :: qy
        real(sp), dimension(mesh%nrow, mesh%ncol), intent(in) :: zb
        real(sp) :: a, b, hstar, v, L
        integer :: col

        a = 1._sp
        b = 0.5_sp
        hstar = 0.1_sp
        gravity = 9.81_sp
        L = 4._sp

        do col=1, mesh%ncol
            h(:, col) = max(0._sp, &
            (b * hstar) / (a * a) * (2 * ((col-1) * mesh%dx(:, col) - L / 2.) - b) - zb(:, col))
        end do
        
        v = sqrt(2 * gravity * hstar) / a
        
        qx(:, :) = 0._sp
        qy(:, :) = h(:, :) * v 

    end subroutine initial_thacker2d



    subroutine shallow_water_2d_time_step(setup, mesh, input_data, options, returns, &
            time_step, ac_qtz, zb, manning, ac_qz)
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

        integer :: i, j, row, col, k, time_returns, nt_sw, nt_sw_old, c_routing_time
        real(sp) :: t, dt
        real(sp) :: heps, qeps, hfx, hfy, maxhsw, hxm, hxp, hym, hyp, dzb, h, hh, z
        real(sp) :: qxc, qyc, theta
        real(sp) :: bb, hh1, hh2, ee1, ee2
        !$AD start-exclude

        ! initialisation
        t = (time_step - 1) * setup%dt
        
        call initial_zero(mesh%nrow, mesh%ncol, hsw, qx, qy, zb)
        ! call initial_thacker2d(mesh, hsw, qx, qy, zb)

        heps = 1.e-6
        qeps = 1.e-8
        c_routing_time = 1
        theta = 1._sp
        

        do while (t .lt. time_step * setup%dt .and. c_routing_time .le. returns%nt_sw) ! .and. time_step .le. 10
            ! print *, t, dt
            ! print *, c_routing_time, t, time_step * setup%dt, setup%dt, time_step
            eta = zb + hsw
            
            returns%sw2d(:, :, c_routing_time, 1) = hsw(:, :)
            returns%sw2d(:, :, c_routing_time, 2) = eta(:, :)
            returns%sw2d(:, :, c_routing_time, 3) = qx(:, :)
            returns%sw2d(:, :, c_routing_time, 4) = qy(:, :)
            returns%sw2d_times(c_routing_time) = t


            maxhsw = max(heps, maxval(hsw))
           
            dt = 0.7 * min(minval(mesh%dx), minval(mesh%dy)) &
            / sqrt(gravity * maxhsw)
            
            ! print *, "c est le tab de dt qui plante !!!"
            ! returns%sw2d_dt(c_routing_time) = dt
            ! print *, "ne s affiche pas "
            print *, dt, maxhsw

            do row = 1, mesh%nrow
                do col = 2, mesh%ncol

                    hfx = max(eta(row, col-1), eta(row, col)) - max(zb(row, col-1), zb(row, col))
                    
                    if (hfx .ge. heps) then

                        !theta = dt / mesh%dx(row, col) * min(abs(qx(row, col) / hfx), sqrt(gravity * hfx))
                        
                        qxc = theta * qx(row, col) + 0.5_sp * (1._sp - theta) * (qx(row, col-1) + qx(row, col+1))

                        qx(row, col) = (qxc - dt * gravity * hfx * &
                            (eta(row, col) - eta(row, col-1)) / mesh%dx(row, col)) / &
                            (1 + dt * gravity * manning(row, col) ** 2 * abs(qx(row, col)) &
                            / hfx ** (7._sp / 3._sp))
                    else 
                        qx(row, col) = 0._sp
                    end if

                    if (abs(qx(row, col)) .lt. qeps) then
                        qx(row, col) = 0._sp
                    end if

                end do
            end do

            do row = 2, mesh%nrow
                do col = 1, mesh%ncol

                    hfy = max(eta(row-1, col), eta(row, col)) - max(zb(row-1, col), zb(row, col))

                    if (hfy .ge. heps) then

                        qyc = theta * qy(row, col) + 0.5_sp * (1._sp - theta) * (qy(row-1, col) + qy(row+1, col))

                        qy(row, col) = (qyc - dt * gravity * hfy * &
                                (eta(row, col) - eta(row-1, col)) / mesh%dy(row, col)) / &
                                (1 + dt * gravity * manning(row, col) ** 2 * abs(qy(row, col)) &
                                / hfy ** (7._sp / 3._sp))
                    else 
                        qy(row, col) = 0._sp
                    end if

                    if (abs(qy(row, col)) .lt. qeps) then
                        qy(row, col) = 0._sp
                    end if

                end do
            end do
            
            ! update water height
            do row = 1, mesh%nrow
                do col = 1, mesh%ncol
                    k = mesh%rowcol_to_ind_ac(row, col)
                            
                    hsw(row, col) = hsw(row, col) + dt / mesh%dx(row, col) * (qx(row, col) - qx(row, col+1)) & 
                    + dt / mesh%dy(row, col) * (qy(row, col) - qy(row+1, col))

                    hsw(row, col) = max(0._sp, hsw(row, col))

                    if (t .eq. ((time_step - 1) * setup%dt)) then
                        hsw(row, col) = hsw(row, col) + dt * ac_qtz(k, setup%nqz) &
                            / mesh%dx(row, col) / mesh%dy(row, col)
                    end if
                    
                end do
            end do

            call free_outflow(mesh, hsw, zb, qx, qy, manning, dt)
            
            !call bc_wall(mesh%nrow, mesh%ncol, hsw, qx, qy, c_routing_time)

            c_routing_time = c_routing_time + 1
            t = t + dt
        end do

        !update volume discharge
        ! avoir les fluxes aux interfaces h au centre
        ! smash c'est un flux a l interieur
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
