!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /     
!!    / /  / / __/ / /___   
!!   /_/  /_/_/    \____/   
!!                       
!!  This file is part of MFC.
!!
!!  MFC is the legal property of its developers, whose names 
!!  are listed in the copyright file included with this source 
!!  distribution.
!!
!!  MFC is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published 
!!  by the Free Software Foundation, either version 3 of the license 
!!  or any later version.
!!
!!  MFC is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!  
!!  You should have received a copy of the GNU General Public License
!!  along with MFC (LICENSE).  
!!  If not, see <http://www.gnu.org/licenses/>.

!>
!! @file m_phase_change.f90
!! @brief Contains module m_phasechange
!! @author M. Rodriguez
!! @version 1.1
!! @date November 10, 2021

! TODO 
! 1. improve iteration method to include the Riemann method of Chiapolino
! 2. Add Pelanti's finite relaxation approach to the code
! 3. Add changes to the GPU branch using Spencer's suggestion of git diff of
! start of branch with latest commit

!> @brief This module is used to compute phase relaxation for pressure,
!         temperature and chemical interfacial relaxation
module m_phase_change

    ! Dependencies =============================================================

    use m_derived_types        !< definitions of the derived types
    use m_global_parameters    !< definitions of the global parameters
    use m_mpi_proxy            !< message passing interface (mpi) module proxy
    use m_variables_conversion !< state variables type conversion procedures

    use ieee_arithmetic
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_phasechange_module, &
                       s_relaxation_solver,             & 
                       s_relaxation_finite_solver,      & 
                       s_finite_ptg_relaxation,         &
                       s_infinite_p_relaxation,         &
                       s_infinite_p_relaxation_k,       &
                       s_infinite_pt_relaxation,        &
                       s_infinite_pt_relaxation_k,      &
                       s_infinite_ptg_relaxation,       &
                       s_infinite_ptg_relaxation_k

    !> @name Abstract interface for creating function pointers
    !> @{
    abstract interface

        !> @name abstract subroutine for the finite relaxation solver 
        !> @{
        subroutine s_abstract_relaxation_solver(q_cons_vf) ! -------
            import :: scalar_field, sys_size          
            type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        end subroutine

        !> @name abstract subroutine for the finite relaxation solver 
        !> @{
        subroutine s_abstract_relaxation_finite_solver(q_cons_vf, rhs_vf) ! -------
            import :: scalar_field, sys_size
            type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
            type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        end subroutine

    end interface

    !> @name parameters for the phase change part of the code
    !> @{
    integer,         parameter :: newton_iter       = 50        
    !< p_relaxk \alpha iter,                set to 25
    real(kind(0d0)), parameter :: pknewton_eps      = 1.d-15
    !< p_relaxk \alpha threshold,           set to 1e-15
    real(kind(0d0)), parameter :: ptsatnewton_eps   = 1.d-10
    !< saturation temperature tol,          set to 1e-10
    real(kind(0d0)), parameter :: ptgnewton_eps     = 1.d-8
    !< saturation ptg tolerance,            set to 1.d-10
    real(kind(0d0)), parameter :: pres_crith        = 22.06d6   
    !< critical water pressure              set to 22.06d6
    real(kind(0d0)), parameter :: pres_critl        = 1.d4
    !< critical water pressure              set to 1.d3
    real(kind(0d0)), parameter :: t_crit            = 648.d0    
    !< critical water temperature           set to 648
    real(kind(0d0)), parameter :: tsathv            = 900.d0   
    !< saturation temperature threshold,    set to 900
    real(kind(0d0)), parameter :: tsatlv            = 275.d0    
    !< factor for bracketing the solution,  set to 10
    real(kind(0d0)), parameter :: bracket_factor    = 5.d0
    !< maximum pressures allowed in the simulation
    real(kind(0d0)), parameter :: maxp = 1.d18, minp = -1.d9
    !> @}

    !> @name gibbs free energy phase change parameters
    !> @{
    real(kind(0d0)) :: n1, n2, pinf1, pinf2
    real(kind(0d0)) :: gibbsa, gibbsb, gibbsc, gibbsd
    !> @}

    real(kind(0d0)), allocatable, dimension(:) ::   gamma_min, pres_inf

    procedure(s_abstract_relaxation_solver), & 
    pointer :: s_relaxation_solver => null()

    procedure(s_abstract_relaxation_finite_solver), &
    pointer :: s_relaxation_finite_solver => null()


    contains

        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param q_cons_vf Cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface    
        subroutine s_initialize_phasechange_module()

            integer :: i

            n1    = 1.d0/fluid_pp(1)%gamma + 1.d0
            pinf1 = fluid_pp(1)%pi_inf/(1.d0 + fluid_pp(1)%gamma)
            n2    = 1.d0/fluid_pp(2)%gamma + 1.d0
            pinf2 = fluid_pp(2)%pi_inf/(1.d0 + fluid_pp(2)%gamma)
            gibbsa = (n1*fluid_pp(1)%cv - n2*fluid_pp(2)%cv +  & 
                     fluid_pp(2)%qvp - fluid_pp(1)%qvp) / &
                     (n2*fluid_pp(2)%cv - fluid_pp(2)%cv)
            gibbsb = (fluid_pp(1)%qv - fluid_pp(2)%qv) / &
                     (n2*fluid_pp(2)%cv - fluid_pp(2)%cv)
            gibbsc = (n2*fluid_pp(2)%cv - n1*fluid_pp(1)%cv) / &
                     (n2*fluid_pp(2)%cv - fluid_pp(2)%cv)
            gibbsd = (n1*fluid_pp(1)%cv - fluid_pp(1)%cv) / & 
                     (n2*fluid_pp(2)%cv - fluid_pp(2)%cv)

            allocate( gamma_min(1:num_fluids) )
            allocate( pres_inf(1:num_fluids) )

            do i = 1, num_fluids
                gamma_min(i) = 1.d0/fluid_pp(i)%gamma + 1.d0
                pres_inf(i)  = fluid_pp(i)%pi_inf / (1.d0+fluid_pp(i)%gamma)
            end do
            ! associating procedural pointer to the subroutine that will be
            ! utilized to calculate the solution of a given relaxation method
        
            ! opening and writing the header of the run-time information file
            if(relax_model == 0 .or. relax_model == 1) then
                s_relaxation_solver => s_infinite_p_relaxation
            elseif (relax_model == 2) then
                s_relaxation_solver => s_infinite_pt_relaxation
            elseif (relax_model == 3) then
                s_relaxation_solver => s_infinite_ptg_relaxation
            elseif (relax_model == 4) then
                s_relaxation_solver => s_infinite_p_relaxation_k
            elseif (relax_model == 5) then
                s_relaxation_solver => s_infinite_pt_relaxation_k
            elseif (relax_model == 6) then
                s_relaxation_solver => s_infinite_ptg_relaxation_k      
            else
                print '(a)', 'relaxation solver was not set!'
                call s_mpi_abort()
            end if

            !if (relax_model == 1) then
            !    s_relaxation_finite_solver => s_finite_ptg_relaxation
            !end if      

        end subroutine s_initialize_phasechange_module !-------------------------------

        !> the purpose of this procedure is to employ the inputted
        !!      cell-average conservative variables in order to compute
        !!      the cell-average rhs variables of the semidiscrete form
        !!      of the governing equations by utilizing the appropriate
        !!      riemann solver.        
        !!  @param q_cons_vf cell-average conservative variables
        !!  @param q_prim_vf cell-average primitive variables
        !!  @param rhs_vf cell-average rhs variables
        subroutine s_finite_ptg_relaxation(q_cons_vf, rhs_vf) ! -------

            type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
            type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

            real(kind(0d0))                                   ::  sum_alpha, tsat
            real(kind(0d0)), dimension(num_fluids)            ::  p_k, t_k, g_k, z_k
            real(kind(0d0))                                   ::  n_k, pinf_k
            real(kind(0d0))                                   ::  rho_k, rhoeq_k, rhoe
            real(kind(0d0))                                   ::  e_k, phi, psi
            real(kind(0d0))                                   ::  f1, f2, f3, f4
            real(kind(0d0))                                   ::  a1, b1, a2, b2
            real(kind(0d0))                                   ::  rho_i, q, kappa
            real(kind(0d0))                                   ::  deltap, p_i, mu 
            real(kind(0d0))                                   ::  e_i, mdot, nu, theta
            real(kind(0d0))                                   ::  mdotalpha, mdotrhoe
            real(kind(0d0)), dimension(2)                     ::          re
            real(kind(0d0)), dimension(num_fluids,num_fluids) ::          we
            integer :: i,j,k,l,r !< generic loop iterators
            do j = 0, m
                do k = 0, n
                    do l = 0, p
                       ! numerical correction of the volume fractions
                       if (mpp_lim) then
                       end if

                       if (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .gt. palpha_eps .or. &
                           q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .lt. 1.d0-palpha_eps) then
                       ! computing n_k, pinf_k, p_k, t_k, and g_k for finite relaxation
                       phi = 0.d0; psi = 0.d0; f1 = 0.d0; f2 = 0.d0; f3 = 0.d0; f4 = 0.d0;
                       !rhoe = 0.d0;
                       do i = 1, num_fluids
                          n_k    = 1.d0/fluid_pp(i)%gamma + 1.d0
                          pinf_k = fluid_pp(i)%pi_inf/(1.d0 + fluid_pp(i)%gamma)
                          rho_k  = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) &
                                  /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)
                          rhoeq_k = (q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l) & 
                              -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                              /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                          !rhoe = rhoe + q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l)
                          e_k = q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l) &
                               /q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)

                          p_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                          t_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf/(1.d0+fluid_pp(i)%gamma)) &
                              /(rho_k*fluid_pp(i)%cv)
                          g_k(i) = (n_k*fluid_pp(i)%cv-fluid_pp(i)%qvp)*t_k(i) &
                              -fluid_pp(i)%cv*t_k(i)*log(t_k(i)**(n_k) &
                              /((p_k(i)+pinf_k)**(n_k-1.d0)))+fluid_pp(i)%qv
                          z_k(i) = n_k*(p_k(i)+pinf_k);

                          phi = phi + 1.d0/(q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%cv)
                          psi = psi + 1.d0/(fluid_pp(i)%gamma*q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l))

                          f1 = f1 + (p_k(i)+n_k*pinf_k)/q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)
                          f2 = f2 + pinf_k/(q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%cv)
                          f3 = f3 + fluid_pp(i)%qv/(fluid_pp(i)%gamma & 
                                    *q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l))
                          f4 = f4 + (e_k-pinf_k/rho_k) &
                                    /(q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%cv)
                       end do
                       !todo improve this approach
                       tsat = f_tsat(p_k(1))
                       !print *,'prelax =',p_k(1),'trelax = ',t_k(1),', tsat = ',tsat

                       if ( ieee_is_nan(p_k(1)) .or. p_k(1) < 0.d0 ) then 
                           print *, 'code crashed' 
                           call s_mpi_abort()
                       end if

                       kappa = f1/psi
                       rho_i = (phi*f1-psi*f2)/(psi*f4-phi*f3)
                       p_i = (z_k(2)*p_k(1)+z_k(1)*p_k(2))/(z_k(1)+z_k(2));
                       e_i = f4/phi + f2/(rho_i*phi)

                       !mu = 1.d8
                       mu = 0.d0;
                       theta = 1.d8 
                       nu = 0.d0
                       if (t_k(1) .gt. tsat) nu = 1d-3
  
                       deltap = mu*(p_k(1)-p_k(2))
                       q = theta*(t_k(2)-t_k(1))
                       mdot = nu*(g_k(2)-g_k(1))
                       mdotalpha = mdot/rho_i
                       mdotrhoe = mdot*e_i

                       rhs_vf(1+adv_idx%beg-1)%sf(j,k,l) = & 
                              rhs_vf(1+adv_idx%beg-1)%sf(j,k,l) + deltap + q/kappa + mdotalpha
                       rhs_vf(2+adv_idx%beg-1)%sf(j,k,l) = & 
                              rhs_vf(2+adv_idx%beg-1)%sf(j,k,l) - deltap - q/kappa - mdotalpha
                       rhs_vf(1+cont_idx%beg-1)%sf(j,k,l) = & 
                              rhs_vf(1+cont_idx%beg-1)%sf(j,k,l) + mdot
                       rhs_vf(2+cont_idx%beg-1)%sf(j,k,l) = & 
                              rhs_vf(2+cont_idx%beg-1)%sf(j,k,l) - mdot
                       rhs_vf(1+internalenergies_idx%beg-1)%sf(j,k,l) = &
                              rhs_vf(1+internalenergies_idx%beg-1)%sf(j,k,l) - p_i*deltap + q + mdotrhoe
                       rhs_vf(2+internalenergies_idx%beg-1)%sf(j,k,l) = &
                              rhs_vf(2+internalenergies_idx%beg-1)%sf(j,k,l) + p_i*deltap - q - mdotrhoe
                       end if
                    end do
                end do
            end do
        end subroutine s_finite_ptg_relaxation ! --------------------------------------

        !>  the purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. for conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf cell-average conservative variables
        subroutine s_infinite_p_relaxation(q_cons_vf) ! ----------------
            !> @name relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume reynolds numbers and the weber numbers
            !> @{
            type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf 
            real(kind(0d0))                                   ::      rhoeq_k
            real(kind(0d0))                                   ::           a1
            real(kind(0d0)), dimension(num_fluids)            :: p_k, alpha_k
            !> @}
            integer :: i,j,k,l           !< generic loop iterators
            logical :: relax             !< relaxation procedure determination variable
            do j = 0, m
                do k = 0, n
                    do l = 0, p
                        ! p relaxation ==================================
                        relax = .false.
                        if (mpp_lim) then
                            call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        end if
                        if ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .gt. palpha_eps ) .and. &
                              q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .lt. 1.d0-palpha_eps ) relax = .true.
                        if (relax) then
                            do i = 1, num_fluids
                                 alpha_k(i) = q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                 rhoeq_k = (q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l) & 
                                          -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                          /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                 p_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                            end do
                            a1 = f_alpha1_prelax(p_k,alpha_k)
                            ! cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l) = 1.d0 - a1
                            call s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                        end if
                    end do
                end do
            end do
        end subroutine s_infinite_p_relaxation ! ----------------

        !>  the purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. for conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf cell-average conservative variables
        subroutine s_infinite_pt_relaxation(q_cons_vf) ! ----------------
            type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf 
            !> @name relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume reynolds numbers and the weber numbers
            !> @{
            real(kind(0d0)), dimension(num_fluids)            :: p_k, alpha_k 
            real(kind(0d0))                                   :: rhoeq_k, rhoe, a1
            real(kind(0d0))                                   :: rhoalpha1, rhoalpha2
            !> @}
            integer :: i, j, k, l        !< generic loop iterators
            logical :: relax             !< relaxation procedure determination variable
            do j = 0, m
                do k = 0, n
                    do l = 0, p
                        ! p relaxation ==================================
                        relax = .false.
                        if (mpp_lim) then
                            call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        end if
                        if ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .gt. palpha_eps ) .and. &
                              q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .lt. 1.d0-palpha_eps ) relax = .true.
                        if (relax) then
                            do i = 1, num_fluids
                                 alpha_k(i) = q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                 rhoeq_k = (q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l) & 
                                          -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                          /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                 p_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                            end do
                            a1 = f_alpha1_prelax(p_k,alpha_k)
                            ! cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l) = 1.d0 - a1
                        end if
                        call s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                        ! pt relaxation ==================================
                        rhoe = 0.d0
                        relax = .false.
                        if (mpp_lim) then
                            call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        end if
                        if ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .gt. palpha_eps ) .and. &
                              q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .lt. 1.d0-palpha_eps ) relax = .true.
                        if (relax) then
                            rhoalpha1 = q_cons_vf(cont_idx%beg)%sf(j,k,l)
                            rhoalpha2 = q_cons_vf(1+cont_idx%beg)%sf(j,k,l)
                            do i = 1, num_fluids
                                rhoe = rhoe + q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l)
                            end do
                            a1 = f_alpha1_ptrelax(rhoalpha1,rhoalpha2,rhoe)
                            ! cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l)  = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l)  = 1.d0 - a1
                        end if
                        call s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                    end do
                end do
            end do
        end subroutine s_infinite_pt_relaxation ! -----------------------

        !>  the purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. for conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf cell-average conservative variables
        subroutine s_infinite_ptg_relaxation(q_cons_vf) ! ----------------
            type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf 
            !> @name relaxed pressure, initial partial pressures, function f(p) and its partial
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume reynolds numbers and the weber numbers
            !> @{
            real(kind(0d0))                                   :: pres_relax, trelax
            real(kind(0d0)), dimension(num_fluids)            :: p_k, alpha_k, tk
            real(kind(0d0))                                   :: rhoalpha1, rhoalpha2
            real(kind(0d0))                                   :: rho, rhoe, rhoeq_k
            real(kind(0d0))                                   :: rho1, rho2
            real(kind(0d0))                                   :: a1, a2
            real(kind(0d0))                                   :: gamma, pi_inf, p_infk
            real(kind(0d0))                                   :: tsat
            real(kind(0d0)), dimension(2)                     :: re
            real(kind(0d0)), dimension(num_fluids,num_fluids) :: we
            !> @}
            integer :: i, j, k, l        !< generic loop iterators
            logical :: relax             !< relaxation procedure determination variable
            !< computing the constant saturation properties 
            do j = 0, m
                do k = 0, n
                    do l = 0, p
                        ! p relaxation ========================================
                        relax = .false.
                        if (mpp_lim) then
                            call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        end if
                        call s_convert_to_mixture_variables( q_cons_vf, rho, &
                                                             gamma, pi_inf,  &
                                                             re, we, j, k, l )
                        if ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .gt. palpha_eps ) .and. &
                              q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .lt. 1.d0-palpha_eps ) relax = .true.
                        if (relax) then
                            do i = 1, num_fluids
                                 alpha_k(i) = q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                 rhoeq_k = (q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l) & 
                                          -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                          /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                 p_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                            end do
                            a1 = f_alpha1_prelax(p_k,alpha_k)
                            ! cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l) = 1.d0 - a1
                        end if
                        call s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                        ! pt relaxation ==========================================
                        rhoe = 0.d0
                        relax = .false.
                        if (mpp_lim) then
                            call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        end if
                        call s_convert_to_mixture_variables( q_cons_vf, rho, &
                                                             gamma, pi_inf,  &
                                                             re, we, j, k, l )
                        if ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .gt. palpha_eps ) .and. &
                              q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .lt. 1.d0-palpha_eps ) relax = .true.
                        if (relax) then
                            rhoalpha1 = q_cons_vf(cont_idx%beg)%sf(j,k,l)
                            rhoalpha2 = q_cons_vf(1+cont_idx%beg)%sf(j,k,l)
                            do i = 1, num_fluids
                                rhoe = rhoe + q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l)
                            end do
                            a1 = f_alpha1_ptrelax(rhoalpha1,rhoalpha2,rhoe)
                            ! cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l)  = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l)  = 1.d0 - a1
                        end if
                        call s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                        ! checking if ptg relaxation is needed  =====================
                        rhoe = 0.d0
                        relax = .false.
                        if (mpp_lim) then
                            call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        end if
                        call s_convert_to_mixture_variables( q_cons_vf, rho, &
                                                             gamma, pi_inf,  &
                                                             re, we, j, k, l )
                        if ((q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .gt. ptgalpha_eps ) .and. &
                             q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .lt. 1.d0-ptgalpha_eps ) relax = .true.
                        if (relax) then
                           do i = 1, num_fluids
                               rhoe = rhoe + q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l) 
                           end do                   
                           pres_relax = (rhoe - pi_inf)/gamma
                           tsat = f_tsat(pres_relax)
                           do i = 1, num_fluids
                             tk(i) = ((q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l) & 
                                    -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                    /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) &
                                    -fluid_pp(i)%pi_inf & 
                                    /(1.d0+fluid_pp(i)%gamma)) &
                                    /(q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%cv &
                                    /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)) 
                           end do
                           if (tk(1) .lt. tsat) relax = .false.
                        end if
                        ! ptg relaxation procedure ===========================
                        if (relax) then
                            call s_compute_ptg_ptrelax(pres_relax,trelax,rho,rhoe)
                            p_infk = fluid_pp(1)%pi_inf/(1.d0+fluid_pp(1)%gamma)
                            rho1 = (pres_relax + p_infk)*fluid_pp(1)%gamma /& 
                                   (fluid_pp(1)%cv*trelax)
                            p_infk = fluid_pp(2)%pi_inf/(1.d0+fluid_pp(2)%gamma)
                            rho2 = (pres_relax + p_infk)*fluid_pp(2)%gamma /& 
                                   (fluid_pp(2)%cv*trelax)
                            ! calculate vapor and liquid volume fractions
                            a1 = (rho-rho2)/(rho1-rho2)
                            a2 = 1.d0 - a1
                            ! cell update of the volume fraction
                            q_cons_vf(cont_idx%beg)%sf(j,k,l)   = rho1*a1
                            q_cons_vf(1+cont_idx%beg)%sf(j,k,l) = rho2*a2
                            q_cons_vf(adv_idx%beg)%sf(j,k,l)    = a1
                            q_cons_vf(1+adv_idx%beg)%sf(j,k,l)  = a2
                        end if
                        call s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                    end do
                end do
            end do
        end subroutine s_infinite_ptg_relaxation ! -----------------------

        !> description: the purpose of this procedure is to infinitely relax
        !!              the pressures from the internal-energy equations to a
        !!              unique pressure, from which the corresponding volume
        !!              fraction of each phase are recomputed. for conservation
        !!              purpose, this pressure is finally corrected using the
        !!              mixture-total-energy equation.
        subroutine s_infinite_p_relaxation_k(q_cons_vf) ! ----------------        
            ! relaxed pressure, initial partial pressures, function f(p) and its partial
            ! derivative df(p), isentropic partial density, sum of volume fractions,
            ! mixture density, dynamic pressure, surface energy, specific heat ratio
            ! function, liquid stiffness function (two variations of the last two
            ! ones), shear and volume reynolds numbers and the weber numbers
            ! cell-average conservative variables
            type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf         
            real(kind(0d0)), dimension(num_fluids) ::      rho_k_s, pres_k_init
            ! generic loop iterators
            integer :: i, j, k, l
            ! relaxation procedure determination variable
            logical :: relax
            do j = 0, m
                do k = 0, n
                    do l = 0, p
                        ! numerical correction of the volume fractions
                        if (mpp_lim) then
                            call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        end if
                        ! pressures relaxation procedure ===================================
                        relax = .true.
                        do i = 1, num_fluids
                            if (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .gt. (1.d0-palpha_eps)) relax = .false.
                        end do
                        if (relax) then
                            ! calculating the initial pressure
                            do i = 1, num_fluids
                               pres_k_init(i) = &
                                ((q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l) & 
                                - q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                / q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) &
                                - fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                            !if (pres_k_init(i) .lt. 0.d0 .and. dabs(pres_k_init(i)) .lt. pres_inf(i)) &
                            !    pres_k_init(i) = sgm_eps
                            end do
                            call s_compute_p_relax_k(rho_k_s,pres_k_init,q_cons_vf,j,k,l)
                            ! cell update of the volume fraction
                            do i = 1, num_fluids
                              if ((q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .lt. 1.d0-palpha_eps) .and. & 
                                  (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .gt. palpha_eps)) &
                                    q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) = & 
                                    q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / rho_k_s(i)
                            end do
                        end if
                        call s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                    end do
                end do
            end do
        end subroutine s_infinite_p_relaxation_k ! -----------------------

        !>  the purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. for conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf cell-average conservative variables
        subroutine s_infinite_pt_relaxation_k(q_cons_vf) ! ----------------
            ! relaxed pressure, initial partial pressures, function f(p) and its partial
            ! derivative df(p), isentropic partial density, sum of volume fractions,
            ! mixture density, dynamic pressure, surface energy, specific heat ratio
            ! function, liquid stiffness function (two variations of the last two
            ! ones), shear and volume reynolds numbers and the weber numbers
            ! cell-average conservative variables
            type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf         
            real(kind(0d0)), dimension(num_fluids) ::     rho_k_s, pres_k_init
            real(kind(0d0))                        ::     rhoe, pstar, tstar
            ! generic loop iterators
            integer :: i,j,k,l
            ! relaxation procedure determination variable
            logical :: relax
            do j = 0, m
                do k = 0, n
                    do l = 0, p
                        ! numerical correction of the volume fractions
                        if (mpp_lim) then
                            call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        end if
                        ! p relaxation==============================================
                        relax = .true.
                        do i = 1, num_fluids
                            if (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .gt. (1d0-palpha_eps)) relax = .false.
                        end do
                        if (relax) then
                            do i = 1, num_fluids
                               ! calculating the initial pressure
                               pres_k_init(i) = &
                                ((q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l) & 
                                - q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                / q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) &
                                - fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                               if (pres_k_init(i) .le. 0.d0) pres_k_init(i) = 1e-2
                           end do
                           call s_compute_p_relax_k(rho_k_s,pres_k_init,q_cons_vf,j,k,l)
                           ! cell update of the volume fraction
                           do i = 1, num_fluids
                            if ((q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .lt. 1.d0-palpha_eps) .and. & 
                                (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .gt. palpha_eps)) &
                                 q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) = & 
                                 q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / rho_k_s(i)
                           end do
                        end if
                        call s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                        ! pt relaxation==============================================
                        if (mpp_lim) then
                            call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        end if
                        relax = .false.
                        rhoe = 0.d0
                        do i = 1, num_fluids
                            if ((q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .lt. 1.d0-palpha_eps) .and. & 
                                (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .gt. palpha_eps)) relax = .true.
                        end do
                        if (relax) then
                            do i = 1, num_fluids
                               rhoe = rhoe + q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l) &
                                           - q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv
                            end do                   
                            call s_compute_pt_relax_k(pstar,tstar,rhoe,q_cons_vf,j,k,l)
                            do i = 1, num_fluids
                               q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) = & 
                                    (gamma_min(i)-1.d0)*q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) & 
                                    *fluid_pp(i)%cv*tstar/(pstar+pres_inf(i))
                            end do
                        end if
                        call s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                    end do
                end do
            end do
        end subroutine s_infinite_pt_relaxation_k ! -----------------------

        !>  the purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. for conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf cell-average conservative variables
        subroutine s_infinite_ptg_relaxation_k(q_cons_vf) ! ----------------
            !> @name relaxed pressure, initial partial pressures, function f(p) and its partial
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume reynolds numbers and the weber numbers
            !> @{
            type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf 
            real(kind(0d0))                                   :: pres_relax, trelax
            real(kind(0d0)), dimension(num_fluids)            :: pres_k_init, rho_k_s
            real(kind(0d0))                                   :: rhoalpha1, rhoalpha2
            real(kind(0d0))                                   :: rho, rhoe, rhoeq
            real(kind(0d0))                                   :: rcv, bsum, tmix, tliquid
            real(kind(0d0))                                   :: rho1, rho2
            real(kind(0d0))                                   :: a1, a2
            real(kind(0d0))                                   :: gamma, pi_inf, p_infk
            real(kind(0d0))                                   :: tsat, dyn_pres
            real(kind(0d0)), dimension(2)                     :: re
            real(kind(0d0)), dimension(num_fluids,num_fluids) :: we
            !> @}
            integer :: i, j, k, l        !< generic loop iterators
            logical :: relax             !< relaxation procedure determination variable
            !< computing the constant saturation properties 
            do j = 0, m
                do k = 0, n
                    do l = 0, p
                        ! numerical correction of the volume fractions
                        !if (mpp_lim) then
                        !    call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        !end if
                        ! p relaxation==============================================
                        relax = .true.
                        do i = 1, num_fluids
                            if (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .gt. (1.d0-palpha_eps)) relax = .false.
                            !if ((q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .lt. 1.d0-palpha_eps) .and. & 
                            !    (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .gt. palpha_eps)) relax = .true.
                        end do
                        if (relax) then
                            do i = 1, num_fluids
                              ! calculating the initial pressure
                               pres_k_init(i) = &
                                ((q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l) & 
                                - q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                / q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) &
                                - fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                            end do
                            call s_compute_p_relax_k(rho_k_s,pres_k_init,q_cons_vf,j,k,l)
                            ! cell update of the volume fraction
                            do i = 1, num_fluids
                              if (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .gt. palpha_eps) &
                              !if ((q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .lt. 1.d0-palpha_eps) .and. & 
                              !    (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .gt. palpha_eps)) &
                                    q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) = & 
                                    q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / rho_k_s(i)
                            end do
                        end if
                        call s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                        ! pt relaxation==============================================
                        !if (mpp_lim) then
                        !    call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        !end if
                        rhoeq = 0.d0
                        relax = .true.
                        do i = 1, num_fluids
                            if (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .gt. (1.d0-palpha_eps)) relax = .false.
                            !if ((q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .lt. 1.d0-palpha_eps) .and. & 
                            !    (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .gt. palpha_eps)) relax = .true.
                        end do
                        if (relax) then
                            do i = 1, num_fluids
                               rhoeq = rhoeq + q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l) &
                                           - q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv
                            end do                   
                            call s_compute_pt_relax_k(pres_relax,trelax,rhoeq,q_cons_vf,j,k,l)
                            do i = 1, num_fluids
                                q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) = & 
                                (gamma_min(i)-1.d0)*q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) & 
                                *fluid_pp(i)%cv*trelax/(pres_relax+pres_inf(i))
                            end do
                        end if
                        call s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                        ! checking if ptg relaxation is needed  =====================
                        !if (mpp_lim) then
                        !    call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        !end if
                        rhoe = 0.d0; bsum = 0.d0; rcv  = 0.d0; dyn_pres = 0.d0
                        relax = .false.
                        !if (mpp_lim) then
                        !    call s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        !end if
                        if ((q_cons_vf(adv_idx%beg)%sf(j,k,l) .gt. ptgalpha_eps) .and. &
                            (q_cons_vf(adv_idx%beg)%sf(j,k,l) .lt. 1.d0-ptgalpha_eps)) relax = .true.
                        !if ((q_cons_vf(adv_idx%beg)%sf(j,k,l) .gt. ptgalpha_eps) .and. &
                        !    (q_cons_vf(adv_idx%beg)%sf(j,k,l) .lt. 1.d0-ptgalpha_eps) .and. & 
                        !    (q_cons_vf(adv_idx%beg+1)%sf(j,k,l) .gt. ptgalpha_eps) .and. & 
                        !    (q_cons_vf(adv_idx%beg+1)%sf(j,k,l) .lt. 1.d0-ptgalpha_eps) ) relax = .true.
                        if (relax) then
                           !call s_convert_to_mixture_variables( q_cons_vf, rho, &
                           !                                     gamma, pi_inf,  &
                           !                                     re, we, j, k, l )
                           do i = 1, num_fluids
                               rhoe = rhoe + q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l) 
                               bsum = bsum + q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)*pres_inf(i)
                               rcv = rcv + q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%cv
                           end do                   
                           !do i = mom_idx%beg, mom_idx%end
                           !   dyn_pres = dyn_pres + 5d-1*q_cons_vf(i)%sf(j,k,l) * & 
                           !   q_cons_vf(i)%sf(j,k,l) / max(rho,sgm_eps)
                           !end do
                           !pres_relax = (q_cons_vf(e_idx)%sf(j,k,l) - dyn_pres - pi_inf)/gamma
                           pres_relax = (rhoe - pi_inf)/gamma
                           tmix = gamma*(pres_relax+bsum)/rcv
                           if(pres_relax .lt. pres_crith .and. pres_relax .gt. pres_critl .and. tmix .lt. t_crit) then
                             tsat = f_tsat(pres_relax)
                             tliquid = ((q_cons_vf(internalenergies_idx%beg)%sf(j,k,l) & 
                                    -q_cons_vf(cont_idx%beg)%sf(j,k,l)*fluid_pp(1)%qv) &
                                    /q_cons_vf(adv_idx%beg)%sf(j,k,l) &
                                    -fluid_pp(1)%pi_inf & 
                                    /(1.d0+fluid_pp(1)%gamma)) &
                                    /(q_cons_vf(cont_idx%beg)%sf(j,k,l)*fluid_pp(1)%cv &
                                    /q_cons_vf(adv_idx%beg)%sf(j,k,l)) 
                             if (tliquid .lt. tsat) relax = .false.
                             if (tliquid .gt. t_crit) relax = .false. ! critical temperature, originally set to 700
                           else
                             relax = .false.
                           end if
                        end if
                        ! ptg relaxation procedure ===========================
                        if (relax) then
                            call s_compute_ptg_ptrelax(pres_relax,trelax,rho,rhoe)
                            p_infk = fluid_pp(1)%pi_inf/(1.d0+fluid_pp(1)%gamma)
                            rho1 = (pres_relax + p_infk)*fluid_pp(1)%gamma /& 
                                   (fluid_pp(1)%cv*trelax)
                            p_infk = fluid_pp(2)%pi_inf/(1.d0+fluid_pp(2)%gamma)
                            rho2 = (pres_relax + p_infk)*fluid_pp(2)%gamma /& 
                                   (fluid_pp(2)%cv*trelax)
                            ! calculate vapor and liquid volume fractions
                            a1 = (rho-rho2)/(rho1-rho2)
                            a2 = 1.d0 - a1 - q_cons_vf(2+adv_idx%beg)%sf(j,k,l)
                            !do i = 2, num_fluids-1
                            !   a2 = a2 - q_cons_vf(i+adv_idx%beg)%sf(j,k,l)
                            !end do
                            ! cell update of the volume fraction
                            q_cons_vf(cont_idx%beg)%sf(j,k,l)   = rho1*a1
                            q_cons_vf(1+cont_idx%beg)%sf(j,k,l) = rho2*a2
                            q_cons_vf(adv_idx%beg)%sf(j,k,l)    = a1
                            q_cons_vf(1+adv_idx%beg)%sf(j,k,l)  = a2
                        end if
                        call s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                    end do
                end do
            end do
        end subroutine s_infinite_ptg_relaxation_k ! -----------------------

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!! subroutines subroutines subroutines !!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !> @name relaxed pressure, initial partial pressures, function f(p) and its partial
        !! derivative df(p), isentropic partial density, sum of volume fractions,
        !! mixture density, dynamic pressure, surface energy, specific heat ratio
        !! function, liquid stiffness function (two variations of the last two
        !! ones), shear and volume reynolds numbers and the weber numbers
        !> @{
        subroutine s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
            type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf 
            integer, intent(in)                                    :: j, k, l
            real(kind(0d0))                                        :: sum_alpha
            !> @}
            integer :: i           !< generic loop iterators
            sum_alpha = 0d0
            do i = 1, num_fluids
               if ((q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) .lt. 0d0) .or. &
                   (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .lt. 0d0)) then
                    q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) = sgm_eps
                    q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)  = sgm_eps
                    q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l)  = 0d0
               end if
               if (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .gt. 1d0) & 
                   q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) = 1d0
               sum_alpha = sum_alpha + q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)
            end do
            do i = 1, num_fluids
                   q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) = & 
                   q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) / sum_alpha
            end do
        end subroutine s_mixture_volume_fraction_correction

        !> @name relaxed pressure, initial partial pressures, function f(p) and its partial
        !! derivative df(p), isentropic partial density, sum of volume fractions,
        !! mixture density, dynamic pressure, surface energy, specific heat ratio
        !! function, liquid stiffness function (two variations of the last two
        !! ones), shear and volume reynolds numbers and the weber numbers
        !> @{
        subroutine s_mixture_total_energy_correction(q_cons_vf, j, k, l )
            !> @name relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume reynolds numbers and the weber numbers
            !> @{
            type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf 
            integer, intent(in)                                    :: j, k, l
            real(kind(0d0))                                        :: rho, dyn_pres, e_we
            real(kind(0d0))                                        :: gamma, pi_inf, pres_relax
            real(kind(0d0)), dimension(2)                          :: re
            real(kind(0d0)), dimension(num_fluids,num_fluids)      :: we
            !> @}
            integer :: i     !< generic loop iterators
            call s_convert_to_mixture_variables( q_cons_vf, rho, &
                                                 gamma, pi_inf,  &
                                                 re, we, j, k, l )
            dyn_pres = 0.d0
            do i = mom_idx%beg, mom_idx%end
                 dyn_pres = dyn_pres + 5d-1*q_cons_vf(i)%sf(j,k,l) * & 
                   q_cons_vf(i)%sf(j,k,l) / max(rho,sgm_eps)
            end do
            e_we = 0.d0
            pres_relax = (q_cons_vf(e_idx)%sf(j,k,l) - dyn_pres - pi_inf - e_we)/gamma
            do i = 1, num_fluids
                  q_cons_vf(i+internalenergies_idx%beg-1)%sf(j,k,l) = & 
                  q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) * & 
                  (fluid_pp(i)%gamma*pres_relax + fluid_pp(i)%pi_inf) +  & 
                  q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv
            end do
        end subroutine s_mixture_total_energy_correction

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! two phase pressure function !!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>     the purpose of this subroutine is to determine the saturation
        !!         temperature by using a newton-raphson method from the provided
        !!         equilibrium pressure and eos of the binary phase system.
        !!     @param q_cons_vf cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface    
        function f_alpha1_prelax(p_k,alpha_k)
            !> @name in-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, a-d, and iteration variables, f and df
            !> @{
            real(kind(0d0))                                       ::  pstar, f_alpha1_prelax
            real(kind(0d0)), dimension(num_fluids), intent(in)    ::  p_k, alpha_k
            real(kind(0d0))                                       ::  z1, z2, pi, c1, c2
            real(kind(0d0))                                       ::  ap, bp, dp
            ! calculating coefficients, eq. c.6, pelanti 2014
            z1 = n1*(p_k(1)+pinf1)
            z2 = n2*(p_k(2)+pinf2)
            pi = (z2*p_k(1)+z1*p_k(2))/(z1+z2)
            c1 = 2.d0*n1*pinf1+(n1-1.d0)*p_k(1)
            c2 = 2.d0*n2*pinf2+(n2-1.d0)*p_k(2)
            ap = 1.d0 + n2*alpha_k(1) + n1*alpha_k(2)
            bp = c1*alpha_k(2)+c2*alpha_k(1)-(n2+1.d0)*alpha_k(1)*p_k(1)-(n1+1.d0)*alpha_k(2)*p_k(2)
            dp = -(c2*alpha_k(1)*p_k(1) + c1*alpha_k(2)*p_k(2))
            ! calculating the tstar temperature, eq. c.7, pelanti 2014
            pstar = (-bp + dsqrt(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
            f_alpha1_prelax = alpha_k(1)*((n1-1.d0)*pstar + 2.d0*p_k(1) + c1)/((n1+1.d0)*pstar+c1)
        end function f_alpha1_prelax

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! two phase pt function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>     the purpose of this subroutine is to determine the saturation
        !!         temperature by using a newton-raphson method from the provided
        !!         equilibrium pressure and eos of the binary phase system.
        !!     @param q_cons_vf cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface    
        function f_alpha1_ptrelax(rhoalpha1,rhoalpha2,e0)
            !> @name in-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, a-d, and iteration variables, f and df
            !> @{
            real(kind(0d0))                :: pstar, f_alpha1_ptrelax
            real(kind(0d0)), intent(in)    :: rhoalpha1, rhoalpha2, e0
            real(kind(0d0))                ::  cv1, cv2, q1, q2
            real(kind(0d0))                ::        ap, bp, dp
            cv1 = fluid_pp(1)%cv; q1 = fluid_pp(1)%qv;
            cv2 = fluid_pp(2)%cv; q2 = fluid_pp(2)%qv;
            ! calculating coefficients, eq. c.6, pelanti 2014
            ap = rhoalpha1*cv1 + rhoalpha2*cv2
            bp = q1*cv1*(n1-1.d0)*rhoalpha1*rhoalpha1 + q2*cv2*(n2-1.d0)*rhoalpha2*rhoalpha2 + &
                 rhoalpha1*cv1*(n1*pinf1+pinf2) + rhoalpha2*cv2*(n2*pinf2+pinf1) + &
                 rhoalpha1*rhoalpha2*(q1*cv2*(n2-1.d0)+q2*cv1*(n1-1.d0)) - &
                 e0*(cv1*(n1-1.d0)*rhoalpha1 + cv2*(n2-1.d0)*rhoalpha2)
            dp = q1*cv1*(n1-1.d0)*pinf2*rhoalpha1*rhoalpha1 + q2*cv2*(n2-1.d0)*pinf1*rhoalpha2*rhoalpha2 + &
                 pinf1*pinf2*(rhoalpha1*cv1*n1 + rhoalpha2*cv2*n2) + & 
                 rhoalpha1*rhoalpha2*(q1*cv2*(n2-1.d0)*pinf1 + q2*cv1*(n1-1.d0)*pinf2) - &
                 e0*(cv1*(n1-1.d0)*pinf2*rhoalpha1 + cv2*(n2-1.d0)*pinf1*rhoalpha2)
            ! calculating the tstar temperature, eq. c.7, pelanti 2014
            pstar = (-bp + dsqrt(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
            f_alpha1_ptrelax = (cv1*(n1-1.d0)*(pstar+pinf2)*rhoalpha1)/&
                     (cv1*(n1-1.d0)*(pstar+pinf2)*rhoalpha1 + &
                      cv2*(n2-1.d0)*(pstar+pinf1)*rhoalpha2)
        end function f_alpha1_ptrelax

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! saturation temperature functions !!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>     the purpose of this subroutine is to determine the saturation
        !!         temperature by using a newton-raphson method from the provided
        !!         equilibrium pressure and eos of the binary phase system.
        !!     @param q_cons_vf cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface    
        subroutine s_compute_fdftsat(fp,dfdp,pstar,tstar)
            real(kind(0d0)), intent(out)        :: fp, dfdp
            real(kind(0d0)), intent(in)         :: pstar, tstar
            fp = gibbsa + gibbsb/tstar + gibbsc*dlog(tstar) - & 
                 dlog((pstar+pinf2)/(pstar+pinf1)**gibbsd)
            dfdp = -gibbsb/(tstar*tstar) + gibbsc/tstar
        end subroutine s_compute_fdftsat !-------------------------------

        !>     the purpose of this subroutine is to determine the bracket of 
        !!         the pressure by finding the pressure at which b^2=4*a*c
        !!     @param p_star equilibrium pressure at the interface    
        subroutine s_compute_tsat_bracket(tstara,tstarb,pressure)
            !> @name in-subroutine variables: vapor and liquid material properties
            !> @{
            real(kind(0d0)), intent(out)   :: tstara, tstarb
            real(kind(0d0)), intent(in)    :: pressure
            real(kind(0d0))                :: fa, fb, dfdp, factor
            ! finding lower bound, getting the bracket 
            factor = 20.d0
            tstara = tsatlv
            tstarb = tsatlv+factor
            call s_compute_fdftsat(fa,dfdp,pressure,tstara)
            call s_compute_fdftsat(fb,dfdp,pressure,tstarb)
            do while ( fa*fb .gt. 0.d0 )
               if( ieee_is_nan(fb) ) then
                   fb = fa
                   tstarb = tstara
                   factor = factor-10.d0
               else 
                   factor = 20.d0
               end if
               if (tstara .gt. tsathv) then
                    print *, 'tsat bracketing failed to find lower bound'
                    print *, 'tstara :: ',tstara,', pressure :: ',pressure
                    print *, 'fa :: ',fa,', fb :: ',fb
                    call s_mpi_abort()
               end if
               fa = fb
               tstara = tstarb
               tstarb = tstara+factor
               call s_compute_fdftsat(fb,dfdp,pressure,tstarb)
            end do
        end subroutine s_compute_tsat_bracket

        !>     the purpose of this subroutine is to determine the saturation
        !!         temperature by using a newton-raphson method from the provided
        !!         equilibrium pressure and eos of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        function f_tsat(pressure)
            !> @name in-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, a-d, and iteration variables, f and df
            !> @{
            real(kind(0d0)), intent(in)    :: pressure
            real(kind(0d0))                :: f_tsat, tstar
            real(kind(0d0))                :: delta, delta_old, fp, dfdp
            real(kind(0d0))                :: fl, fh, tstarl, tstarh, tsata, tsatb
            integer :: iter                !< generic loop iterators
            call s_compute_tsat_bracket(tsata,tsatb,pressure)
            ! computing f at lower and higher end of the bracket
            call s_compute_fdftsat(fl,dfdp,pressure,tsata)
            call s_compute_fdftsat(fh,dfdp,pressure,tsatb)
            ! establishing the direction of the descent to find zero
            if(fl < 0.d0) then
                tstarl  = tsata; tstarh  = tsatb;
            else
                tstarl  = tsata; tstarh  = tsatb;
            end if
            tstar = 0.5d0*(tstarl+tstarh)
            delta_old = dabs(tstarh-tstarl)
            delta = delta_old
            call s_compute_fdftsat(fp,dfdp,pressure,tstar)
            ! combining bisection and newton-raphson methods
            do iter = 0, newton_iter
                if ((((tstar-tstarh)*dfdp-fp)*((tstar-tstarl)*dfdp-fp) > 0.d0) & ! bisect if newton out of range,
                        .or. (dabs(2.0*fp) > dabs(delta_old*dfdp))) then         ! or not decreasing fast enough.
                   delta_old = delta
                   delta = 0.5d0*(tstarh-tstarl)
                   tstar = tstarl + delta
                   if (delta .eq. 0.d0) exit
                else                    ! newton step acceptable, take it
                   delta_old = delta
                   delta = fp/dfdp
                   tstar = tstar - delta
                   if (delta .eq. 0.d0) exit
                end if
                if (dabs(delta/tstar) < ptsatnewton_eps) exit
                call s_compute_fdftsat(fp,dfdp,pressure,tstar)           
                if (fp < 0.d0) then     !maintain the bracket on the root
                   tstarl = tstar
                else
                   tstarh = tstar
                end if
                if (iter .eq. newton_iter) then
                    print *, 'tsat : ',tstar,', iter : ',iter
                    print *, 'tsat did not converge, stopping code'
                    call s_mpi_abort()
                end if                  
            end do
            f_tsat = tstar
        end function f_tsat !-------------------------------

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!! two-phase ptg relaxation functions !!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>     the purpose of this subroutine is to determine the saturation
        !!         temperature by using a newton-raphson method from the provided
        !!         equilibrium pressure and eos of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        subroutine s_compute_ptg_fdf(fp,dfdp,pstar,tstar,rho0,e0)
            !> @name in-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, a-d, and iteration variables, f and df
            !> @{
            real(kind(0d0)), intent(in)    ::   pstar, rho0, e0
            real(kind(0d0)), intent(out)   ::   fp, dfdp, tstar
            real(kind(0d0))                ::  cv1, cv2, q1, q2
            real(kind(0d0))                ::        ap, bp, dp
            real(kind(0d0))                ::  dadp, dbdp, dddp
            real(kind(0d0))                ::              dtdp
            cv1 = fluid_pp(1)%cv; q1 = fluid_pp(1)%qv;
            cv2 = fluid_pp(2)%cv; q2 = fluid_pp(2)%qv;
            ! calculating coefficients, eq. c.6, pelanti 2014
            ap = rho0*cv1*cv2*((n2-1.d0)*(pstar+n1*pinf1)-(n1-1.d0)*(pstar+n2*pinf2))
            bp = e0*((n1-1.d0)*cv1*(pstar+pinf2) - (n2-1.d0)*cv2*(pstar+pinf1)) + &
                 rho0*((n2-1.d0)*cv2*q1*(pstar+pinf1) - (n1-1.d0)*cv1*q2*(pstar+pinf2)) + &
                 cv2*(pstar+pinf1)*(pstar+n2*pinf2) - cv1*(pstar+pinf2)*(pstar+n1*pinf1)
            dp = (q2-q1)*(pstar+pinf1)*(pstar+pinf2)
            ! calculating the tstar temperature, eq. c.7, pelanti 2014
            tstar = (-bp + dsqrt(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
            ! calculating the derivatives wrt pressure of the coefficients
            dadp = rho0*cv1*cv2*((n2-1.d0)-(n1-1.d0))
            dbdp = e0*((n1-1.d0)*cv1 - (n2-1.d0)*cv2) + &
                   rho0*((n2-1.d0)*cv2*q1 - (n1-1.d0)*cv1*q2) + &
                   cv2*((pstar+pinf1)+(pstar+n2*pinf2)) - & 
                   cv1*((pstar+pinf2)+(pstar+n1*pinf1))
            dddp = (q2-q1)*((pstar+pinf1)+(pstar+pinf2))
            ! derivative of the temperature wrt to pressure, needed for dfdp
            dtdp = (-dbdp + (0.5d0/dsqrt(bp*bp-4.d0*ap*dp))*(2.d0*bp*dbdp-&
                   4.d0*(ap*dddp+dp*dadp)))/(2.d0*ap) - (dadp/ap)*tstar
            fp   = gibbsa + gibbsb/tstar + gibbsc*dlog(tstar) + & 
                   gibbsd*dlog(pstar+pinf1) - dlog(pstar+pinf2)
            dfdp = -gibbsb/(tstar*tstar)*dtdp + gibbsc/tstar*dtdp + & 
                   gibbsd/(pstar+pinf1) - 1.d0/(pstar+pinf2)
        end subroutine s_compute_ptg_fdf !------------------------

        !>     the purpose of this subroutine is to determine the bracket of 
        !!         the pressure by finding the pressure at which b^2=4*a*c
        !!     @param p_star equilibrium pressure at the interface    
        subroutine s_compute_ptg_bracket(pstara,pstarb,pstar,rho0,e0)
            !> @name in-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, a-d, and iteration variables, f and df
            !> @{
            real(kind(0d0)), intent(out)   :: pstara, pstarb
            real(kind(0d0)), intent(in)    :: pstar, rho0, e0
            real(kind(0d0))                :: fa, fb, dfdp, tstar
            real(kind(0d0))                :: ps, factor
            pstara = 0.7d0*pstar
            pstarb = pstar
            call s_compute_ptg_fdf(fa,dfdp,pstara,tstar,rho0,e0)
            call s_compute_ptg_fdf(fb,dfdp,pstarb,tstar,rho0,e0)
            factor = 1.05d0
            do while ( fa*fb .gt. 0.d0 )
                  if (pstara .gt. 1.d12) then
                         print *, 'ptg bracketing failed to find lower bound'
                         print *, 'pstara :: ',pstara
                         call s_mpi_abort()
                  end if
                  fa = fb
                  pstara = pstarb
                  pstarb = pstara*factor
                  call s_compute_ptg_fdf(fb,dfdp,pstarb,tstar,rho0,e0)
                  if( ieee_is_nan(fb) ) then
                        fb = fa
                        pstarb = pstara
                        factor = factor*0.95d0
                  else 
                        factor = 1.05d0
                  end if
            end do
        end subroutine s_compute_ptg_bracket

        !>     the purpose of this subroutine is to determine the saturation
        !!         temperature by using a newton-raphson method from the provided
        !!         equilibrium pressure and eos of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        subroutine s_compute_ptg_ptrelax(pstar,tstar,rho0,e0)
            !> @name in-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, a-d, and iteration variables, f and df
            !> @{
            real(kind(0d0)), intent(inout) :: pstar
            real(kind(0d0)), intent(out)   :: tstar
            real(kind(0d0)), intent(in)    :: rho0, e0
            real(kind(0d0))                ::  pstara, pstarb, dtstar
            real(kind(0d0))                ::  delta, delta_old, fp, dfdp
            real(kind(0d0))                ::  fl, fh, pstarl, pstarh
            integer :: iter                !< generic loop iterators
            ! computing the bracket of the root solution
            call s_compute_ptg_bracket(pstara,pstarb,pstar,rho0,e0)
            ! computing f at lower and higher end of the bracket
            call s_compute_ptg_fdf(fl,dfdp,pstara,dtstar,rho0,e0)
            call s_compute_ptg_fdf(fh,dfdp,pstarb,dtstar,rho0,e0)
            ! establishing the direction of the descent to find zero
            if(fl < 0.d0) then
                pstarl  = pstara; pstarh  = pstarb;
            else
                pstarl  = pstarb; pstarh  = pstara;
            end if
            pstar = 0.5d0*(pstara+pstarb)
            delta_old = dabs(pstarb-pstara)
            delta = delta_old
            call s_compute_ptg_fdf(fp,dfdp,pstar,dtstar,rho0,e0)
            ! combining bisection and newton-raphson methods
            do iter = 0, newton_iter
                if ((((pstar-pstarh)*dfdp-fp)*((pstar-pstarl)*dfdp-fp) > 0.d0) & ! bisect if newton out of range,
                        .or. (dabs(2.0*fp) > dabs(delta_old*dfdp))) then         ! or not decreasing fast enough.
                   delta_old = delta
                   delta = 0.5d0*(pstarh-pstarl)
                   pstar = pstarl + delta
                   if (delta .eq. 0.d0) exit                    ! change in root is negligible
                else                                              ! newton step acceptable, take it
                   delta_old = delta
                   delta = fp/dfdp
                   pstar = pstar - delta
                   if (delta .eq. 0.d0) exit
                end if
                if (dabs(delta/pstar) < ptgnewton_eps) exit           ! stopping criteria
                ! updating to next iteration
                call s_compute_ptg_fdf(fp,dfdp,pstar,tstar,rho0,e0)
                if (fp < 0.d0) then !maintain the bracket on the root
                   pstarl = pstar
                else
                   pstarh = pstar
                end if
            end do
        end subroutine s_compute_ptg_ptrelax !-------------------------------

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!  k-phase p relaxation functions !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>     the purpose of this subroutine is to determine the saturation
        !!         temperature by using a newton-raphson method from the provided
        !!         equilibrium pressure and eos of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        subroutine s_compute_pk_fdf(fp,dfdp,pstar,rho_k_s,pres_k_init,q_cons_vf,j,k,l)
            !> @name in-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, a-d, and iteration variables, f and df
            !> @{
            ! cell-average conservative variables
            type(scalar_field), dimension(sys_size), intent(in)  :: q_cons_vf
            real(kind(0d0)), intent(out)                         :: fp, dfdp
            real(kind(0d0)), intent(in)                          :: pstar
            real(kind(0d0)), dimension(num_fluids), intent(out)  :: rho_k_s
            real(kind(0d0)), dimension(num_fluids), intent(in)   :: pres_k_init
            real(kind(0d0))                                      :: numerator, denominator
            real(kind(0d0))                                      :: drhodp
            integer, intent(in)                                  :: j, k, l
            integer                                              :: i
            fp = -1.d0; dfdp = 0.d0;
            do i = 1, num_fluids
                  numerator   = gamma_min(i)*(pstar+pres_inf(i))
                  denominator = numerator + pres_k_init(i)-pstar
                  rho_k_s(i)  = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)/&
                     max(q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l),sgm_eps)*numerator/denominator
                  drhodp      = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / & 
                     max(q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l),sgm_eps) * & 
                     gamma_min(i)*(pres_k_init(i)+pres_inf(i)) / (denominator*denominator)
                  fp          = fp  + q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / rho_k_s(i)
                  dfdp        = dfdp - q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) * & 
                     drhodp / (rho_k_s(i)*rho_k_s(i))
            end do
        end subroutine s_compute_pk_fdf !------------------------

        !>     the purpose of this subroutine is to determine the bracket of 
        !!         the pressure by finding the pressure at which b^2=4*a*c
        !!     @param p_star equilibrium pressure at the interface    
        subroutine s_compute_pk_bracket(pstara,pstarb,rho_k_s,pres_k_init,q_cons_vf,j,k,l)
            !> @name in-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, a-d, and iteration variables, f and df
            !> @{
            type(scalar_field), dimension(sys_size), intent(in)  :: q_cons_vf
            real(kind(0d0)), intent(out)   :: pstara, pstarb
            real(kind(0d0))                :: fa, fb, dfdp
            real(kind(0d0)), dimension(num_fluids), intent(in)   :: pres_k_init
            real(kind(0d0)), dimension(num_fluids), intent(out)  :: rho_k_s
            integer, intent(in)            :: j, k, l
            integer                        :: i
            pstara = 1.d-2; pstarb = 1.d5;
            call s_compute_pk_fdf(fa,dfdp,pstara,rho_k_s,pres_k_init,q_cons_vf,j,k,l)
            call s_compute_pk_fdf(fb,dfdp,pstarb,rho_k_s,pres_k_init,q_cons_vf,j,k,l)
            do while ( fa*fb .gt. 0.d0 )
               !if ( ieee_is_nan(fb) .or. ieee_is_nan(fa) .or. dabs(fb) .gt. 20.0d0 ) then
               !   pstarb = pstara*bracket_factor*0.5d0
               !   call s_compute_pk_fdf(fb,dfdp,pstarb,rho_k_s,pres_k_init,q_cons_vf,j,k,l)
               if ( pstara .gt. maxp ) then
                  pstara = -1.d1; pstarb = -1.d5;
                  call s_compute_pk_fdf(fa,dfdp,pstara,rho_k_s,pres_k_init,q_cons_vf,j,k,l)
                  call s_compute_pk_fdf(fb,dfdp,pstarb,rho_k_s,pres_k_init,q_cons_vf,j,k,l)
               else if ( (pstara .lt. minp) .or. ieee_is_nan(fb) ) then
                  print *, 'p-k bracketing failed to find lower bound'
                  print *, 'location j ',j,', k ',k,', l', l,', pstara :: ',pstara
                  do i = 1, num_fluids
                     print *, 'alpha ',i,' :: ',q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)
                     print *, 'rhoal ',i,' :: ',q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)
                     print *, 'p_k ',i,' :: ',pres_k_init(i)
                  end do
                  call s_mpi_abort()                
               else 
                  fa = fb
                  pstara = pstarb
                  pstarb = pstara*bracket_factor
                  call s_compute_pk_fdf(fb,dfdp,pstarb,rho_k_s,pres_k_init,q_cons_vf,j,k,l)               
               end if
            end do
        end subroutine s_compute_pk_bracket

        !>     the purpose of this subroutine is to determine the saturation
        !!         temperature by using a newton-raphson method from the provided
        !!         equilibrium pressure and eos of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        subroutine s_compute_p_relax_k(rho_k_s,pres_k_init,q_cons_vf,j,k,l)
            !> @name in-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, a-d, and iteration variables, f and df
            !> @{
            type(scalar_field), dimension(sys_size), intent(in)  :: q_cons_vf
            real(kind(0d0))                                      :: pstar, pstara, pstarb
            real(kind(0d0))                                      :: delta, delta_old, fp, dfdp
            real(kind(0d0))                                      :: fl, fh, pstarl, pstarh
            real(kind(0d0)), dimension(num_fluids), intent(in)   :: pres_k_init
            real(kind(0d0)), dimension(num_fluids), intent(out)  :: rho_k_s 
            integer, intent(in)                                  :: j,k,l
            integer :: iter                !< generic loop iterators
            ! computing the bracket of the root solution
            call s_compute_pk_bracket(pstara,pstarb,rho_k_s,pres_k_init,q_cons_vf,j,k,l)
            ! computing f at lower and higher end of the bracket
            call s_compute_pk_fdf(fl,dfdp,pstara,rho_k_s,pres_k_init,q_cons_vf,j,k,l)
            call s_compute_pk_fdf(fh,dfdp,pstarb,rho_k_s,pres_k_init,q_cons_vf,j,k,l)
            ! establishing the direction of the descent to find zero
            if(fl < 0.d0) then
                pstarl  = pstara; pstarh  = pstarb;
            else
                pstarl  = pstarb; pstarh  = pstara;
            end if
            pstar = 0.5d0*(pstara+pstarb)
            delta_old = dabs(pstarb-pstara)
            delta = delta_old
            call s_compute_pk_fdf(fp,dfdp,pstar,rho_k_s,pres_k_init,q_cons_vf,j,k,l)
            ! combining bisection and newton-raphson methods
            do iter = 0, newton_iter
                if ((((pstar-pstarh)*dfdp-fp)*((pstar-pstarl)*dfdp-fp) > 0.d0) & ! bisect if newton out of range,
                        .or. (dabs(2.0*fp) > dabs(delta_old*dfdp))) then         ! or not decreasing fast enough.
                   delta_old = delta
                   delta = 0.5d0*(pstarh-pstarl)
                   pstar = pstarl + delta
                   if (delta .eq. 0.d0) exit                    ! change in root is negligible
                else                                              ! newton step acceptable, take it
                   delta_old = delta
                   delta = fp/dfdp
                   pstar = pstar - delta
                   if (delta .eq. 0.d0) exit
                end if
                if (dabs(delta/pstar) < pknewton_eps) exit           ! stopping criteria
                ! updating to next iteration
                call s_compute_pk_fdf(fp,dfdp,pstar,rho_k_s,pres_k_init,q_cons_vf,j,k,l)
                if (fp < 0.d0) then !maintain the bracket on the root
                   pstarl = pstar
                else
                   pstarh = pstar
                end if
            end do
        end subroutine s_compute_p_relax_k !-------------------------------

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!  k-phase pt relaxation functions !!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>     the purpose of this subroutine is to determine the saturation
        !!         temperature by using a newton-raphson method from the provided
        !!         equilibrium pressure and eos of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        subroutine s_compute_ptk_fdf(fp,dfdp,pstar,tstar,rhoe,q_cons_vf,j,k,l)
            !> @name in-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, a-d, and iteration variables, f and df
            !> @{
            ! cell-average conservative variables
            type(scalar_field), dimension(sys_size), intent(in)  :: q_cons_vf
            real(kind(0d0)), intent(out)                         :: fp, dfdp, tstar
            real(kind(0d0)), intent(in)                          :: rhoe, pstar
            real(kind(0d0))                                      :: tdem, dtdp, fa
            real(kind(0d0))                                      :: df, fk, tdemk
            integer, intent(in)                                  :: j, k, l
            integer                                              :: i
            tdem = 0.d0; dtdp = 0.d0; fa = 0.d0; df = 0.d0;
            do i = 1, num_fluids
                  tdemk = gamma_min(i)*q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) & 
                        *fluid_pp(i)%cv
                  fk = (q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%cv & 
                        *(gamma_min(i)-1.d0))/(pstar+pres_inf(i))
                  tdem = tdem + tdemk
                  fa = fa + fk
                  df = df - fk/(pstar+pres_inf(i))
            end do
            tstar = (rhoe + pstar)/tdem
            fp = 1.d0 - tstar*fa
            dfdp = -tstar*df-fa/tdem
        end subroutine s_compute_ptk_fdf !------------------------

        !>     the purpose of this subroutine is to determine the bracket of 
        !!         the pressure by finding the pressure at which b^2=4*a*c
        !!     @param p_star equilibrium pressure at the interface    
        subroutine s_compute_ptk_bracket(pstara,pstarb,rhoe,q_cons_vf,j,k,l)
            !> @name in-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, a-d, and iteration variables, f and df
            !> @{
            type(scalar_field), dimension(sys_size), intent(in)  :: q_cons_vf
            real(kind(0d0)), intent(out)   :: pstara, pstarb
            real(kind(0d0)), intent(in)    :: rhoe
            real(kind(0d0))                :: fa, fb, dfdp
            real(kind(0d0))                :: factor, tstar
            integer, intent(in)            :: j, k, l
            integer                        :: i
            pstara = 1.d1; pstarb = 1.d5;
            call s_compute_ptk_fdf(fa,dfdp,pstara,tstar,rhoe,q_cons_vf,j,k,l)
            call s_compute_ptk_fdf(fb,dfdp,pstarb,tstar,rhoe,q_cons_vf,j,k,l)
            do while ( fa*fb .gt. 0.d0 )
               !if ( ieee_is_nan(fb) .or. ieee_is_nan(fa) .or. dabs(fb) .gt. 20.d0 ) then
               !   pstarb = pstara*bracket_factor*0.5d0
               !   call s_compute_ptk_fdf(fb,dfdp,pstarb,tstar,rhoe,q_cons_vf,j,k,l)
               !if ( pstara .gt. maxp ) then
               !   pstara = -1.d1; pstarb = -1.d5;
               !   call s_compute_ptk_fdf(fa,dfdp,pstara,tstar,rhoe,q_cons_vf,j,k,l)
               !   call s_compute_ptk_fdf(fb,dfdp,pstarb,tstar,rhoe,q_cons_vf,j,k,l)
               if ( (pstara .gt. maxp) .or. ieee_is_nan(fb) ) then
                  print *, 'pt-k bracketing failed to find lower bound'
                  print *, 'location j ',j,', k ',k,', l', l
                  print *, 'rhoe :: ',rhoe
                  do i = 1, num_fluids
                     print *, 'alpha ',i,' :: ',q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)
                     print *, 'rhoal ',i,' :: ',q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)
                  end do
                  call s_mpi_abort()
               else
                  fa = fb
                  pstara = pstarb
                  pstarb = pstara*bracket_factor
                  call s_compute_ptk_fdf(fb,dfdp,pstarb,tstar,rhoe,q_cons_vf,j,k,l)
               end if
            end do
        end subroutine s_compute_ptk_bracket

        !>     the purpose of this subroutine is to determine the saturation
        !!         temperature by using a newton-raphson method from the provided
        !!         equilibrium pressure and eos of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        subroutine s_compute_pt_relax_k(pstar,tstar,rhoe,q_cons_vf,j,k,l)
            !> @name in-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, a-d, and iteration variables, f and df
            !> @{
            type(scalar_field), dimension(sys_size), intent(in)  :: q_cons_vf
            real(kind(0d0)), intent(out)                         :: pstar, tstar
            real(kind(0d0)), intent(in)                          :: rhoe
            real(kind(0d0))                                      :: pstara, pstarb
            real(kind(0d0))                                      :: delta, delta_old, fp, dfdp
            real(kind(0d0))                                      :: fl, fh, pstarl, pstarh
            integer, intent(in)                                  :: j, k, l
            integer :: iter                !< generic loop iterators
            ! computing the bracket of the root solution
            call s_compute_ptk_bracket(pstara,pstarb,rhoe,q_cons_vf,j,k,l)
            ! computing f at lower and higher end of the bracket
            call s_compute_ptk_fdf(fl,dfdp,pstara,tstar,rhoe,q_cons_vf,j,k,l)
            call s_compute_ptk_fdf(fh,dfdp,pstarb,tstar,rhoe,q_cons_vf,j,k,l)
            ! establishing the direction of the descent to find zero
            if(fl < 0.d0) then
                pstarl  = pstara; pstarh  = pstarb;
            else
                pstarl  = pstarb; pstarh  = pstara;
            end if
            pstar = 0.5d0*(pstara+pstarb)
            delta_old = dabs(pstarb-pstara)
            delta = delta_old
            call s_compute_ptk_fdf(fp,dfdp,pstar,tstar,rhoe,q_cons_vf,j,k,l)
            ! combining bisection and newton-raphson methods
            do iter = 0, newton_iter
                if ((((pstar-pstarh)*dfdp-fp)*((pstar-pstarl)*dfdp-fp) > 0.d0) & ! bisect if newton out of range,
                        .or. (dabs(2.0*fp) > dabs(delta_old*dfdp))) then         ! or not decreasing fast enough.
                   delta_old = delta
                   delta = 0.5d0*(pstarh-pstarl)
                   pstar = pstarl + delta
                   if (delta .eq. 0.d0) exit                    ! change in root is negligible
                else                                              ! newton step acceptable, take it
                   delta_old = delta
                   delta = fp/dfdp
                   pstar = pstar - delta
                   if (delta .eq. 0.d0) exit
                end if
                if (dabs(delta/pstar) < pknewton_eps) exit           ! stopping criteria
                ! updating to next iteration
                call s_compute_ptk_fdf(fp,dfdp,pstar,tstar,rhoe,q_cons_vf,j,k,l)
                if (fp < 0.d0) then !maintain the bracket on the root
                   pstarl = pstar
                else
                   pstarh = pstar
                end if
            end do
        end subroutine s_compute_pt_relax_k !-------------------------------

        subroutine s_finalize_relaxation_solver_module()
            deallocate(gamma_min,pres_inf)
            s_relaxation_solver => null()
        end subroutine

end module m_phase_change
