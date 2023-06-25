!%      (MWD) Module Wrapped and Differentiated
!%
!%      Subroutine
!%      ----------
!%
!%      - map_control_to_parameters

module mwd_parameters_manipulation

    use md_constant !% only: sp, nopr_parameters, nopr_states
    use mwd_setup !% only: SetupDT
    use mwd_mesh !% only: MeshDT
    use mwd_input_data !% only: Input_DataDT
    use mwd_parameters !% only: ParametersDT
    use mwd_options !% only: OptionsDT
    use mwd_control !% only: ControlDT_initialise, ControlDT_finalise

    implicit none
    
    public :: parameters_to_control

contains
    
    subroutine uniform_get_control_size(setup, options, n)
        
        implicit none
        
        type(SetupDT), intent(in) :: setup
        type(OptionsDT), intent(in) :: options
        integer, intent(inout) :: n
        
        n = sum(options%optimize%opr_parameters) + sum(options%optimize%opr_initial_states)
    
    end subroutine uniform_get_control_size
    
    subroutine distributed_get_control_size(setup, mesh, options, n)
            
        implicit none
        
        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(OptionsDT), intent(in) :: options
        integer, intent(inout) :: n
        
        n = (sum(options%optimize%opr_parameters) + &
        & sum(options%optimize%opr_initial_states)) * mesh%nac
    
    end subroutine distributed_get_control_size
    
    subroutine multi_linear_get_control_size(setup, options, n)
            
        implicit none
        
        type(SetupDT), intent(in) :: setup
        type(OptionsDT), intent(in) :: options
        integer, intent(inout) :: n
        
        integer :: i
        
        n = 0
        
        do i = 1, nopr_parameters
        
            if (options%optimize%opr_parameters(i) .ne. 1) cycle
            
            n = n + 1 + sum(options%optimize%opr_parameters_descriptor(:, i))
        
        end do
        
        do i = 1, nopr_states
        
            if (options%optimize%opr_initial_states(i) .ne. 1) cycle
            
            n = n + 1 + sum(options%optimize%opr_initial_states_descriptor(:, i))
        
        end do
    
    end subroutine multi_linear_get_control_size
    
    
    subroutine sigmoide(x, res)
    
        implicit none
        
        real(sp), intent(in) :: x
        real(sp), intent(inout) :: res
        
        res = 1._sp / (1._sp + exp(-x))
    
    end subroutine sigmoide
    
    subroutine inv_sigmoide(x, res)
        
        implicit none
        
        real(sp), intent(in) :: x
        real(sp), intent(inout) :: res
        
        res = log(x / (1._sp - x))
    
    end subroutine inv_sigmoide
    
    
    subroutine scaled_sigmoide(x, l, u, res)
        
        implicit none
        
        real(sp), intent(in) :: x, l, u
        real(sp), intent(inout) :: res
        
        call sigmoide(x, res)
        
        res = res * (u - l) + l 
    
    end subroutine scaled_sigmoide
    
    subroutine inv_scaled_sigmoid(x, l, u, res)
    
        implicit none
        
        real(sp), intent(in) :: x, l, u
        real(sp), intent(inout) :: res
        
        real(sp) :: xw, eps = 1e-3_sp
        
        xw = max(x, l + eps)
        xw = min(x, u - eps)
        xw = (xw - l) / (u - l)
        
        call inv_sigmoide(xw, res)
    
    end subroutine inv_scaled_sigmoid
    
    
    subroutine sigmoide2d(x, res)
    
        implicit none
        
        real(sp), dimension(:, :), intent(in) :: x
        real(sp), dimension(:, :), intent(inout) :: res
        
        res = 1._sp / (1._sp + exp(-x))
    
    end subroutine sigmoide2d
    
    
    subroutine scaled_sigmoide2d(x, l, u, res)
    
        implicit none
        
        real(sp), dimension(:, :), intent(in) :: x
        real(sp), intent(in) :: l, u
        real(sp), dimension(:, :), intent(inout) :: res
        
        call sigmoide2d(x, res)
        
        res = res * (u - l) + l
    
    end subroutine scaled_sigmoide2d
    
    
    subroutine opr_parameters_to_matrix(setup, mesh, parameters, matrix)

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(ParametersDT), intent(in) :: parameters
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_parameters), intent(inout) :: matrix

        matrix(:, :, 1) = parameters%opr_parameters%ci
        matrix(:, :, 2) = parameters%opr_parameters%cp
        matrix(:, :, 3) = parameters%opr_parameters%cft
        matrix(:, :, 4) = parameters%opr_parameters%cst
        matrix(:, :, 5) = parameters%opr_parameters%kexc
        matrix(:, :, 6) = parameters%opr_parameters%llr

    end subroutine opr_parameters_to_matrix

    subroutine opr_initial_states_to_matrix(setup, mesh, parameters, matrix)

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(ParametersDT), intent(in) :: parameters
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_states), intent(inout) :: matrix

        matrix(:, :, 1) = parameters%opr_initial_states%hi
        matrix(:, :, 2) = parameters%opr_initial_states%hp
        matrix(:, :, 3) = parameters%opr_initial_states%hft
        matrix(:, :, 4) = parameters%opr_initial_states%hst
        matrix(:, :, 5) = parameters%opr_initial_states%hlr

    end subroutine opr_initial_states_to_matrix
    
    
    subroutine matrix_to_opr_parameters(setup, mesh, matrix, parameters)

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_parameters), intent(in) :: matrix
        type(ParametersDT), intent(inout) :: parameters

        parameters%opr_parameters%ci = matrix(:, :, 1)
        parameters%opr_parameters%cp = matrix(:, :, 2)
        parameters%opr_parameters%cft = matrix(:, :, 3)
        parameters%opr_parameters%cst = matrix(:, :, 4)
        parameters%opr_parameters%kexc = matrix(:, :, 5)
        parameters%opr_parameters%llr = matrix(:, :, 6)

    end subroutine matrix_to_opr_parameters

    subroutine matrix_to_opr_initial_states(setup, mesh, matrix, parameters)

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_states), intent(in) :: matrix
        type(ParametersDT), intent(inout) :: parameters

        parameters%opr_initial_states%hi = matrix(:, :, 1)
        parameters%opr_initial_states%hp = matrix(:, :, 2)
        parameters%opr_initial_states%hft = matrix(:, :, 3)
        parameters%opr_initial_states%hst = matrix(:, :, 4)
        parameters%opr_initial_states%hlr = matrix(:, :, 5)

    end subroutine matrix_to_opr_initial_states
    
    
    subroutine sbs_control_tfm(parameters)
    
        implicit none
        
        type(ParametersDT), intent(inout) :: parameters
        
        integer :: i
        logical, dimension(size(parameters%control%x)) :: nbd_mask
        
        !% Need lower and upper bound to sbs tfm
        nbd_mask = (parameters%control%nbd(:) .eq. 2)
        
        do i = 1, size(parameters%control%x)
        
            if (.not. nbd_mask(i)) cycle
        
            if (parameters%control%l_bkg(i) .lt. 0._sp) then
                
                parameters%control%x(i) = asinh(parameters%control%x(i))
                parameters%control%l(i) = asinh(parameters%control%l_bkg(i))
                parameters%control%u(i) = asinh(parameters%control%u_bkg(i))
            
            else if (parameters%control%l_bkg(i) .ge. 0._sp .and. parameters%control%u_bkg(i) .le. 1._sp) then
            
                parameters%control%x(i) = log(parameters%control%x(i)/(1._sp - parameters%control%x(i)))
                parameters%control%l(i) = log(parameters%control%l_bkg(i)/(1._sp - parameters%control%l_bkg(i)))
                parameters%control%u(i) = log(parameters%control%u_bkg(i)/(1._sp - parameters%control%u_bkg(i)))
                
            else
            
                parameters%control%x(i) = log(parameters%control%x(i))
                parameters%control%l(i) = log(parameters%control%l_bkg(i))
                parameters%control%u(i) = log(parameters%control%u_bkg(i))
            
            end if
        
        end do
    
    end subroutine sbs_control_tfm
    
    subroutine sbs_inv_control_tfm(parameters)
    
        implicit none
        
        type(ParametersDT), intent(inout) :: parameters
        
        integer :: i
        logical, dimension(size(parameters%control%x)) :: nbd_mask
        
        !% Need lower and upper bound to sbs tfm
        nbd_mask = (parameters%control%nbd(:) .eq. 2)
        
        do i = 1, size(parameters%control%x)
        
            if (.not. nbd_mask(i)) cycle
        
            if (parameters%control%l_bkg(i) .lt. 0._sp) then
                
                parameters%control%x(i) = sinh(parameters%control%x(i))
            
            else if (parameters%control%l_bkg(i) .ge. 0._sp .and. parameters%control%u_bkg(i) .le. 1._sp) then
            
                parameters%control%x(i) = exp(parameters%control%x(i))/(1._sp + exp(parameters%control%x(i))) 
                
            else
            
                parameters%control%x(i) = exp(parameters%control%x(i))
            
            end if
        
        end do
        
        parameters%control%l = parameters%control%l_bkg
        parameters%control%u = parameters%control%u_bkg

    end subroutine sbs_inv_control_tfm
    

    subroutine normalize_control_tfm(parameters)
    
        implicit none

        type(ParametersDT), intent(inout) :: parameters
 
        logical, dimension(size(parameters%control%x)) :: nbd_mask
        
        !% Need lower and upper bound to normalize
        nbd_mask = (parameters%control%nbd(:) .eq. 2)
        
        where (nbd_mask)
        
            parameters%control%x = (parameters%control%x - parameters%control%l_bkg) / &
            (parameters%control%u_bkg - parameters%control%l_bkg)
            parameters%control%l = 0._sp
            parameters%control%u = 1._sp
        
        end where
    
    end subroutine normalize_control_tfm
    
    subroutine normalize_inv_control_tfm(parameters)
        
        implicit none

        type(ParametersDT), intent(inout) :: parameters
 
        logical, dimension(size(parameters%control%x)) :: nbd_mask
        
        !% Need lower and upper bound to denormalize
        nbd_mask = (parameters%control%nbd(:) .eq. 2)
        
        where (nbd_mask)
        
            parameters%control%x = parameters%control%x * &
            & (parameters%control%u_bkg - parameters%control%l_bkg) + parameters%control%l_bkg
            parameters%control%l = parameters%control%l_bkg
            parameters%control%u = parameters%control%u_bkg
        
        end where
    
    end subroutine normalize_inv_control_tfm
    
    
    subroutine uniform_parameters_to_control(setup, mesh, parameters, options)
    
        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(ParametersDT), intent(inout) :: parameters
        type(OptionsDT), intent(in) :: options
        
        integer :: n, i, j
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_parameters) :: opr_parameters_matrix
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_states) :: opr_initial_states_matrix
        logical, dimension(mesh%nrow, mesh%ncol) :: ac_mask
        
        call opr_parameters_to_matrix(setup, mesh, parameters, opr_parameters_matrix)
        call opr_initial_states_to_matrix(setup, mesh, parameters, opr_initial_states_matrix)
        
        call uniform_get_control_size(setup, options, n)
        
        call ControlDT_initialise(parameters%control, n)
        
        ac_mask = (mesh%active_cell(:,:) .eq. 1)
        
        j = 0
        
        do i = 1, nopr_parameters
        
            if (options%optimize%opr_parameters(i) .ne. 1) cycle
            
            j = j + 1

            parameters%control%x(j) = sum(opr_parameters_matrix(:, :, i), mask=ac_mask) / mesh%nac
            parameters%control%l(j) = options%optimize%l_opr_parameters(i)
            parameters%control%u(j) = options%optimize%u_opr_parameters(i)

        end do
        
        do i = 1, nopr_states
        
            if (options%optimize%opr_initial_states(i) .ne. 1) cycle
            
            j = j + 1

            parameters%control%x(j) = sum(opr_initial_states_matrix(:, :, i), mask=ac_mask) / mesh%nac
            parameters%control%l(j) = options%optimize%l_opr_initial_states(i)
            parameters%control%u(j) = options%optimize%u_opr_initial_states(i)

        end do
        
!~         parameters%control%x_bkg = parameters%control%x
        parameters%control%l_bkg = parameters%control%l
        parameters%control%u_bkg = parameters%control%u
        parameters%control%nbd = 2
        
    end subroutine uniform_parameters_to_control
        
    subroutine uniform_control_to_parameters(setup, mesh, parameters, options)

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(ParametersDT), intent(inout) :: parameters
        type(OptionsDT), intent(in) :: options
        
        integer :: i, j
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_parameters) :: opr_parameters_matrix
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_states) :: opr_initial_states_matrix
        logical, dimension(mesh%nrow, mesh%ncol) :: ac_mask
        
        
        call opr_parameters_to_matrix(setup, mesh, parameters, opr_parameters_matrix)
        call opr_initial_states_to_matrix(setup, mesh, parameters, opr_initial_states_matrix)
        
        ac_mask = (mesh%active_cell(:,:) .eq. 1)
        
        j = 0
        
        do i = 1, nopr_parameters
        
            if (options%optimize%opr_parameters(i) .ne. 1) cycle
            
            j = j + 1
            
            where (ac_mask)
            
                opr_parameters_matrix(:, :, i) = parameters%control%x(j)
            
            end where

        end do
        
        do i = 1, nopr_states
        
            if (options%optimize%opr_initial_states(i) .ne. 1) cycle
            
            j = j + 1

            where (ac_mask)
            
                opr_initial_states_matrix(:, :, i) = parameters%control%x(j)
            
            end where

        end do
        
        call matrix_to_opr_parameters(setup, mesh, opr_parameters_matrix, parameters)
        call matrix_to_opr_initial_states(setup, mesh, opr_initial_states_matrix, parameters)
    
    end subroutine uniform_control_to_parameters
    
    
    subroutine distributed_parameters_to_control(setup, mesh, parameters, options)
    
        implicit none
        
        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(ParametersDT), intent(inout) :: parameters
        type(OptionsDT), intent(in) :: options
        
        integer :: n, i, j, row, col
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_parameters) :: opr_parameters_matrix
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_states) :: opr_initial_states_matrix
        
        call opr_parameters_to_matrix(setup, mesh, parameters, opr_parameters_matrix)
        call opr_initial_states_to_matrix(setup, mesh, parameters, opr_initial_states_matrix)
        
        call distributed_get_control_size(setup, mesh, options, n)
        
        call ControlDT_initialise(parameters%control, n)
        
        j = 0
        
        do i = 1, nopr_parameters
        
            if (options%optimize%opr_parameters(i) .ne. 1) cycle
            
            do col = 1, mesh%ncol
            
                do row = 1, mesh%nrow
                
                    if (mesh%active_cell(row, col) .eq. 0) cycle
                    
                    j = j + 1
                    
                    parameters%control%x(j) = opr_parameters_matrix(row, col, i)
                    parameters%control%l(j) = options%optimize%l_opr_parameters(i)
                    parameters%control%u(j) = options%optimize%u_opr_parameters(i)
                
                end do
            
            end do
        
        end do
        
        do i = 1, nopr_states
        
            if (options%optimize%opr_initial_states(i) .ne. 1) cycle
            
            do col = 1, mesh%ncol
            
                do row = 1, mesh%nrow
                
                    if (mesh%active_cell(row, col) .eq. 0) cycle
                    
                    j = j + 1
                    
                    parameters%control%x(j) = opr_initial_states_matrix(row, col, i)
                    parameters%control%l(j) = options%optimize%l_opr_initial_states(i)
                    parameters%control%u(j) = options%optimize%u_opr_initial_states(i)
                
                end do
            
            end do
        
        end do
        
        !~         parameters%control%x_bkg = parameters%control%x
        parameters%control%l_bkg = parameters%control%l
        parameters%control%u_bkg = parameters%control%u
        parameters%control%nbd = 2
    
    end subroutine distributed_parameters_to_control

    subroutine distributed_control_to_parameters(setup, mesh, parameters, options)
    
        implicit none
        
        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(ParametersDT), intent(inout) :: parameters
        type(OptionsDT), intent(in) :: options
        
        integer :: i, j, row, col
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_parameters) :: opr_parameters_matrix
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_states) :: opr_initial_states_matrix
        
        call opr_parameters_to_matrix(setup, mesh, parameters, opr_parameters_matrix)
        call opr_initial_states_to_matrix(setup, mesh, parameters, opr_initial_states_matrix)
        
        j = 0
        
        do i = 1, nopr_parameters
        
            if (options%optimize%opr_parameters(i) .ne. 1) cycle
            
            do col = 1, mesh%ncol
            
                do row = 1, mesh%nrow
                
                    if (mesh%active_cell(row, col) .eq. 0) cycle
                    
                    j = j + 1
                    
                    opr_parameters_matrix(row, col, i) = parameters%control%x(j)
                
                end do
            
            end do
        
        end do
        
        do i = 1, nopr_states
        
            if (options%optimize%opr_initial_states(i) .ne. 1) cycle
            
            do col = 1, mesh%ncol
            
                do row = 1, mesh%nrow
                
                    if (mesh%active_cell(row, col) .eq. 0) cycle
                    
                    j = j + 1
                    
                    opr_initial_states_matrix(row, col, i) = parameters%control%x(j)
                                    
                end do
            
            end do
        
        end do
        
        call matrix_to_opr_parameters(setup, mesh, opr_parameters_matrix, parameters)
        call matrix_to_opr_initial_states(setup, mesh, opr_initial_states_matrix, parameters)
    
    end subroutine distributed_control_to_parameters
    
    
    subroutine multi_linear_parameters_to_control(setup, mesh, input_data, parameters, options)
    
        implicit none
        
        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(inout) :: parameters
        type(OptionsDT), intent(in) :: options
        
        integer :: n, i, j, k
        real(sp) :: y, l, u
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_parameters) :: opr_parameters_matrix
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_states) :: opr_initial_states_matrix
        logical, dimension(mesh%nrow, mesh%ncol) :: ac_mask
        
        call opr_parameters_to_matrix(setup, mesh, parameters, opr_parameters_matrix)
        call opr_initial_states_to_matrix(setup, mesh, parameters, opr_initial_states_matrix)
        
        call multi_linear_get_control_size(setup, options, n)
        call ControlDT_initialise(parameters%control, n)
        
        ac_mask = (mesh%active_cell(:,:) .eq. 1)
        
        j = 0
        
        do i = 1, nopr_parameters
        
            if (options%optimize%opr_parameters(i) .ne. 1) cycle
            
            j = j + 1
            
            y = sum(opr_parameters_matrix(:, :, i), mask=ac_mask) / mesh%nac
            l = options%optimize%l_opr_parameters(i)
            u = options%optimize%u_opr_parameters(i)
            
            call inv_scaled_sigmoid(y, l, u, parameters%control%x(j))
            
            do k = 1, setup%nd
            
                if (options%optimize%opr_parameters_descriptor(k, i) .ne. 1) cycle
                
                j = j + 1
                
                parameters%control%x(j) = 0._sp
            
            end do

        end do
        
        do i = 1, nopr_states
        
            if (options%optimize%opr_initial_states(i) .ne. 1) cycle
            
            j = j + 1
            
            y = sum(opr_initial_states_matrix(:, :, i), mask=ac_mask) / mesh%nac
            l = options%optimize%l_opr_initial_states(i)
            u = options%optimize%u_opr_initial_states(i)
            
            call inv_scaled_sigmoid(y, l, u, parameters%control%x(j))
            
            do k = 1, setup%nd
            
                if (options%optimize%opr_initial_states_descriptor(k, i) .ne. 1) cycle
                
                j = j + 1
                
                parameters%control%x(j) = 0._sp
            
            end do

        end do
        
        parameters%control%nbd = 0
    
    end subroutine multi_linear_parameters_to_control
    
    subroutine multi_linear_control_to_parameters(setup, mesh, input_data, parameters, options)
    
        implicit none
        
        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(inout) :: parameters
        type(OptionsDT), intent(in) :: options
        
        integer :: i, j, k
        real(sp) :: l, u
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_parameters) :: opr_parameters_matrix
        real(sp), dimension(mesh%nrow, mesh%ncol, nopr_states) :: opr_initial_states_matrix
        real(sp), dimension(mesh%nrow, mesh%ncol) :: wa2d, norm_desc
        logical, dimension(mesh%nrow, mesh%ncol) :: ac_mask
        
        call opr_parameters_to_matrix(setup, mesh, parameters, opr_parameters_matrix)
        call opr_initial_states_to_matrix(setup, mesh, parameters, opr_initial_states_matrix)
        
        j = 0
        
        do i = 1, nopr_parameters
        
            if (options%optimize%opr_parameters(i) .ne. 1) cycle
            
            j = j + 1
            
            wa2d = parameters%control%x(j)
            
            do k = 1, setup%nd
            
                if (options%optimize%opr_parameters_descriptor(k, i) .ne. 1) cycle
                
                j = j + 1
                
                norm_desc = (input_data%physio_data%descriptor(:, :, k) - input_data%physio_data%l_descriptor(k)) / &
                & (input_data%physio_data%u_descriptor(k) - input_data%physio_data%l_descriptor(k))
                
                wa2d = wa2d + parameters%control%x(j) * norm_desc
            
            end do

            l = options%optimize%l_opr_parameters(i)
            u = options%optimize%u_opr_parameters(i)
        
            call scaled_sigmoide2d(wa2d, l, u, opr_parameters_matrix(:, :, i))
            
        end do
        
        do i = 1, nopr_states
        
            if (options%optimize%opr_initial_states(i) .ne. 1) cycle
            
            j = j + 1
            
            wa2d = parameters%control%x(j)
            
            do k = 1, setup%nd
            
                if (options%optimize%opr_initial_states_descriptor(k, i) .ne. 1) cycle
                
                j = j + 1
                
                norm_desc = (input_data%physio_data%descriptor(:, :, k) - input_data%physio_data%l_descriptor(k)) / &
                & (input_data%physio_data%u_descriptor(k) - input_data%physio_data%l_descriptor(k))
                
                wa2d = wa2d + parameters%control%x(j) * norm_desc
            
            end do
            
            l = options%optimize%l_opr_initial_states(i)
            u = options%optimize%u_opr_initial_states(i)
        
            call scaled_sigmoide2d(wa2d, l, u, opr_initial_states_matrix(:, :, i))
        
        end do
        
        call matrix_to_opr_parameters(setup, mesh, opr_parameters_matrix, parameters)
        call matrix_to_opr_initial_states(setup, mesh, opr_initial_states_matrix, parameters)
    
    end subroutine multi_linear_control_to_parameters
    
    
    subroutine control_tfm(parameters, options)
    
        implicit none
        
        type(ParametersDT), intent(inout) :: parameters
        type(OptionsDT), intent(in) :: options
        
        select case (options%optimize%control_tfm)
        
        case ("sbs")
        
            call sbs_control_tfm(parameters)
        
        case ("normalize")
        
            call normalize_control_tfm(parameters)
        
        end select
    
    end subroutine control_tfm
    
    subroutine inv_control_tfm(parameters, options)
    
        implicit none
        
        type(ParametersDT), intent(inout) :: parameters
        type(OptionsDT), intent(in) :: options
        
        select case (options%optimize%control_tfm)
        
        case ("sbs")
        
            call sbs_inv_control_tfm(parameters)
        
        case ("normalize")
        
            call normalize_inv_control_tfm(parameters)
        
        end select
    
    end subroutine inv_control_tfm
    
    
    subroutine parameters_to_control(setup, mesh, input_data, parameters, options)
    
        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(inout) :: parameters
        type(OptionsDT), intent(in) :: options
        
        select case (options%optimize%mapping)
        
        case ("uniform")
        
            call uniform_parameters_to_control(setup, mesh, parameters, options)
            
        case ("distributed")
        
            call distributed_parameters_to_control(setup, mesh, parameters, options)
            
        case ("multi-linear")
        
            call multi_linear_parameters_to_control(setup, mesh, input_data, parameters, options)
        
        case ("multi-polynomial")
        
        end select
        
        call control_tfm(parameters, options)
    
    end subroutine parameters_to_control

    subroutine control_to_parameters(setup, mesh, input_data, parameters, options)

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(Input_DataDT), intent(in) :: input_data
        type(ParametersDT), intent(inout) :: parameters
        type(OptionsDT), intent(in) :: options
        
        if (.not. allocated(parameters%control%x)) return
        
        call inv_control_tfm(parameters, options)
        
        select case (options%optimize%mapping)
        
        case ("uniform")
        
            call uniform_control_to_parameters(setup, mesh, parameters, options)
            
        case ("distributed")
        
            call distributed_control_to_parameters(setup, mesh, parameters, options)
            
        case ("multi-linear")
        
            call multi_linear_control_to_parameters(setup, mesh, input_data, parameters, options)
        
        end select
        
    end subroutine control_to_parameters

end module mwd_parameters_manipulation
