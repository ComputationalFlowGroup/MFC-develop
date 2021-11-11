!>
!! @file m_time_steppers.f90
!! @brief Contains module m_time_steppers
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The following module features a variety of time-stepping schemes.
!!              Currently, it includes the following Runge-Kutta (RK) algorithms:
!!                   1) 1st Order TVD RK
!!                   2) 2nd Order TVD RK
!!                   3) 3rd Order TVD RK
!!              where TVD designates a total-variation-diminishing time-stepper.
module m_time_steppers

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_rhs                  !< Right-hand-side (RHS) evaluation procedures

    use m_data_output          !< Run-time info & solution data output procedures

    use m_bubbles              !< Bubble dynamics routines

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy
    ! ==========================================================================

    implicit none

    type(vector_field), allocatable, dimension(:) :: q_cons_ts !<
    !! Cell-average conservative variables at each time-stage (TS)

    type(scalar_field), private, allocatable, dimension(:) :: q_prim_vf !<
    !! Cell-average primitive variables at the current time-stage

    type(scalar_field), allocatable, dimension(:) :: rhs_vf !<
    !! Cell-average RHS variables at the current time-stage

    type(vector_field), allocatable, dimension(:) :: q_prim_ts !<
    !! Cell-average primitive variables at consecutive TIMESTEPS

    integer, private :: num_ts !<
    !! Number of time stages in the time-stepping scheme

!$acc declare create(q_cons_ts,q_prim_vf,rhs_vf,q_prim_ts)

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_time_steppers_module() ! -----------------------

        type(bounds_info) :: ix, iy, iz !<
            !! Indical bounds in the x-, y- and z-directions

        integer :: i, j !< Generic loop iterators

        ! Setting number of time-stages for selected time-stepping scheme
        if (time_stepper == 1) then
            num_ts = 1
        elseif (any(time_stepper == (/2, 3/))) then
            num_ts = 2
        end if

        ! Setting the indical bounds in the x-, y- and z-directions
        ix%beg = -buff_size; ix%end = m + buff_size

        if (n > 0) then

            iy%beg = -buff_size; iy%end = n + buff_size

            if (p > 0) then
                iz%beg = -buff_size; iz%end = p + buff_size
            else
                iz%beg = 0; iz%end = 0
            end if

        else

            iy%beg = 0; iy%end = 0
            iz%beg = 0; iz%end = 0

        end if

        ! Allocating the cell-average conservative variables
        allocate (q_cons_ts(1:num_ts))

        do i = 1, num_ts
            allocate (q_cons_ts(i)%vf(1:sys_size))
!$acc enter data create(q_cons_ts(i)%vf(1:sys_size))
        end do

        do i = 1, num_ts
            do j = 1, sys_size
                allocate (q_cons_ts(i)%vf(j)%sf(ix%beg:ix%end, &
                                                iy%beg:iy%end, &
                                                iz%beg:iz%end))
!$acc enter data create(q_cons_ts(i)%vf(j)%sf(ix%beg:ix%end,iy%beg:iy%end,iz%beg:iz%end))
            end do
        end do

        ! Allocating the cell-average primitive ts variables
        if (any(com_wrt) .or. any(cb_wrt) .or. probe_wrt) then
            allocate (q_prim_ts(0:3))

            do i = 0, 3
                allocate (q_prim_ts(i)%vf(1:sys_size))
!$acc enter data create(q_prim_ts(i)%vf(1:sys_size))
            end do

            do i = 0, 3
                do j = 1, sys_size
                    allocate (q_prim_ts(i)%vf(j)%sf(ix%beg:ix%end, &
                                                    iy%beg:iy%end, &
                                                    iz%beg:iz%end))
!$acc enter data create(q_prim_ts(i)%vf(j)%sf(ix%beg:ix%end,iy%beg:iy%end,iz%beg:iz%end))
                end do
            end do
        end if

        ! Allocating the cell-average primitive variables
        allocate (q_prim_vf(1:sys_size))

        do i = mom_idx%beg, E_idx
            allocate (q_prim_vf(i)%sf(ix%beg:ix%end, &
                                      iy%beg:iy%end, &
                                      iz%beg:iz%end))
!$acc enter data create(q_prim_vf(i)%sf(ix%beg:ix%end,iy%beg:iy%end,iz%beg:iz%end))
        end do

        if (bubbles) then
            do i = bub_idx%beg, sys_size
                allocate (q_prim_vf(i)%sf(ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
!$acc enter data create(q_prim_vf(i)%sf(ix%beg:ix%end,iy%beg:iy%end,iz%beg:iz%end))
            end do
        end if


        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                allocate (q_prim_vf(i)%sf(ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end))
!$acc enter data create(q_prim_vf(i)%sf(ix%beg:ix%end,iy%beg:iy%end,iz%beg:iz%end))
            end do
        end if

        ! Allocating the cell-average RHS variables
        allocate (rhs_vf(1:sys_size))

        do i = 1, sys_size
            allocate (rhs_vf(i)%sf(0:m, 0:n, 0:p))
!$acc enter data create(rhs_vf(i)%sf(0:m, 0:n, 0:p))
        end do

        ! Opening and writing the header of the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_open_run_time_information_file()
        end if

    end subroutine s_initialize_time_steppers_module ! ---------------------

    !> 1st order TVD RK time-stepping algorithm
        !! @param t_step Current time step
    subroutine s_1st_order_tvd_rk(t_step) ! --------------------------------

        integer, intent(IN) :: t_step

        integer :: i !< Generic loop iterator

        ! Stage 1 of 1 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
!$acc enter data attach(q_prim_vf(i)%sf)
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
!$acc enter data attach(q_prim_vf(i)%sf)
        end do


        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)
        if (DEBUG) print *, 'got rhs'

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if
        if (DEBUG) print *, 'wrote runtime info'

        if (any(com_wrt) .or. any(cb_wrt) .or. probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        do i = 1, sys_size
            q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + dt*rhs_vf(i)%sf
!$acc update device(q_cons_ts(1)%vf(i)%sf)
        end do

        !print *, q_cons_ts(1)%vf(cont_idx%beg)%sf(50,30,0)
        !print *, q_cons_ts(1)%vf(E_idx)%sf(50,30,0)
        !print *, q_cons_ts(1)%vf(adv_idx%end)%sf(50,30,0)
        !print *, q_cons_ts(1)%vf(mom_idx%beg)%sf(50,30,0)

        
        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        do i = 1, cont_idx%end
!$acc exit data detach(q_prim_vf(i)%sf)
            q_prim_vf(i)%sf => null()
        end do


        do i = adv_idx%beg, adv_idx%end
!$acc exit data detach(q_prim_vf(i)%sf)
            q_prim_vf(i)%sf => null()
        end do

        ! ==================================================================

    end subroutine s_1st_order_tvd_rk ! ------------------------------------

    !> 2nd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_2nd_order_tvd_rk(t_step) ! --------------------------------

        integer, intent(IN) :: t_step

        integer :: i !< Generic loop iterator

        ! Stage 1 of 2 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
!$acc enter data attach(q_prim_vf(i)%sf)
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
!$acc enter data attach(q_prim_vf(i)%sf)
        end do

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (any(com_wrt) .or. any(cb_wrt) .or. probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        do i = 1, sys_size
            q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + dt*rhs_vf(i)%sf
        end do


        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        ! ==================================================================

        ! Stage 2 of 2 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
        end do

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) = &
                (q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + dt*rhs_vf(i)%sf)/2d0
        end do


        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => null()
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => null()
        end do
        ! ==================================================================

    end subroutine s_2nd_order_tvd_rk ! ------------------------------------

    !> 3rd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_3rd_order_tvd_rk(t_step) ! --------------------------------

        integer, intent(IN) :: t_step

        integer :: i, j !< Generic loop iterator

        ! Stage 1 of 3 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
!$acc enter data attach(q_prim_vf(i)%sf)
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
!$acc enter data attach(q_prim_vf(i)%sf)
        end do

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (any(com_wrt) .or. any(cb_wrt) .or. probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        do i = 1, sys_size
            q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) = &
                q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                + dt*rhs_vf(i)%sf
        end do


        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)

        ! ==================================================================

        ! Stage 2 of 3 =====================================================
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
!$acc enter data attach(q_prim_vf(i)%sf)
        end do

        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
!$acc enter data attach(q_prim_vf(i)%sf)
        end do

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) = &
                (3d0*q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + dt*rhs_vf(i)%sf)/4d0
        end do


        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)

        ! ==================================================================

        ! Stage 3 of 3 =====================================================
        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

        do i = 1, sys_size
            q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) = &
                (q_cons_ts(1)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + 2d0*q_cons_ts(2)%vf(i)%sf(0:m, 0:n, 0:p) &
                 + 2d0*dt*rhs_vf(i)%sf)/3d0
        end do


        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        do i = 1, cont_idx%end
!$acc exit data detach(q_prim_vf(i)%sf)
            q_prim_vf(i)%sf => null()
        end do

        do i = adv_idx%beg, adv_idx%end
!$acc exit data detach(q_prim_vf(i)%sf)
            q_prim_vf(i)%sf => null()
        end do
        ! ==================================================================

    end subroutine s_3rd_order_tvd_rk ! ------------------------------------



    !> This subroutine saves the temporary q_prim_vf vector
        !!      into the q_prim_ts vector that is then used in p_main
        !! @param t_step current time-step
    subroutine s_time_step_cycling(t_step) ! ----------------------------

        integer, intent(IN) :: t_step

        integer :: i !< Generic loop iterator

        if (t_step == t_step_start) then
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 1) then
            do i = 1, sys_size
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 2) then
            do i = 1, sys_size
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 3) then
            do i = 1, sys_size
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        else ! All other timesteps
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_ts(2)%vf(i)%sf(:, :, :)
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_ts(1)%vf(i)%sf(:, :, :)
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_ts(0)%vf(i)%sf(:, :, :)
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        end if

    end subroutine s_time_step_cycling ! -----------------------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_time_steppers_module() ! -------------------------

        integer :: i, j !< Generic loop iterators

        ! Deallocating the cell-average conservative variables
        do i = 1, num_ts

            do j = 1, sys_size
                deallocate (q_cons_ts(i)%vf(j)%sf)
            end do

            deallocate (q_cons_ts(i)%vf)

        end do

        deallocate (q_cons_ts)

        ! Deallocating the cell-average primitive ts variables
        if (any(com_wrt) .or. any(cb_wrt) .or. probe_wrt) then
            do i = 0, 3
                do j = 1, sys_size
                    deallocate (q_prim_ts(i)%vf(j)%sf)
                end do
                deallocate (q_prim_ts(i)%vf)
            end do
            deallocate (q_prim_ts)
        end if

        ! Deallocating the cell-average primitive variables
        do i = mom_idx%beg, E_idx
            deallocate (q_prim_vf(i)%sf)
        end do
        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                deallocate (q_prim_vf(i)%sf)
            end do
        end if

        deallocate (q_prim_vf)

        ! Deallocating the cell-average RHS variables
        do i = 1, sys_size
            deallocate (rhs_vf(i)%sf)
        end do

        deallocate (rhs_vf)

        ! Writing the footer of and closing the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_close_run_time_information_file()
        end if

    end subroutine s_finalize_time_steppers_module ! -----------------------

end module m_time_steppers
