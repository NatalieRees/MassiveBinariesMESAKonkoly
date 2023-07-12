! ***********************************************************************
!
!   Copyright (C) 2012-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! *********************************************************************** 
      module run_binary_extras 

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      use math_lib
      
      implicit none

      integer, parameter :: ilx_reached_rlo = 1
      integer, parameter :: ilx_reached_critical = 2

      ! There is probably a better-practices way to do this, but I'm starting simple and we can make it more elegant
      ! integer, parameter :: accretion_method = 1 !0, 1, 2, 3,       

      contains
      
      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! Set these function pointers to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
         b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_start_step=> extras_binary_start_step
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          b% warn_binary_extra =.false.
         
      end subroutine extras_binary_controls


      integer function how_many_extra_binary_history_header_items(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_header_items = 0
      end function how_many_extra_binary_history_header_items


      subroutine data_for_extra_binary_history_header_items( &
           binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id, n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
      end subroutine data_for_extra_binary_history_header_items


      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 0
      end function how_many_extra_binary_history_columns


      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         
      end subroutine data_for_extra_binary_history_columns
      
      
      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart
         call binary_ptr(binary_id, b, ierr)

         ! Not necessary, but just in case, if not restart, 
         ! set xtra(1) = 0, which will be used as the proxy 
         ! for total AM accreted
         if (.not. restart) then b% xtra(1)=0d0

         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         
!          b% s1% job% warn_run_star_extras = .false.
          extras_binary_startup = keep_going

      end function  extras_binary_startup
      
      integer function extras_binary_start_step(binary_id,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr

         extras_binary_start_step = keep_going
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
      
      end function  extras_binary_start_step
      
      !Return either keep_going, retry or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_check_model = keep_going
        
      end function extras_binary_check_model
      
      
      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         

         ! Declare your variables! 
         ! real(dp) :: added_angular_momentum, total_Irot, omega_crit
         integer :: k, nz

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_finish_step = keep_going

         if (b% r(1) > b% rl(1) .and. .not. b% lxtra(ilx_reached_rlo)) then 
            ! things to do upon the first RLOF is reached
            b% lxtra(ilx_reached_rlo) = .true.
            b% s1% solver_itermin_until_reduce_min_corr_coeff = 25
            b% s1% solver_max_tries_before_reject = 40
            b% s1% max_tries_for_retry = 40
            b% s1% tiny_corr_coeff_limit = 1000    
            b% s1% corr_coeff_limit = 0.2d0
            b% s1% ignore_too_large_correction = .true.
            b% s1% ignore_min_corr_coeff_for_scale_max_correction = .true.
            b% s1% use_gold_tolerances = .true.
            b% s1% use_gold2_tolerances = .false.
            b% s1% gold_solver_iters_timestep_limit = 30
            b% s1% gold_iter_for_resid_tol3 = 10
            b% s1% gold_tol_residual_norm3 = 1d-6
            b% s1% gold_tol_max_residual3 = 1d-3
            b% s1% tol_max_correction = 1d-2
            b% s1% tol_correction_norm = 1d-3
            b% s1% max_corr_jump_limit = 1d99
            b% s1% max_resid_jump_limit = 1d99

            b% s2% solver_itermin_until_reduce_min_corr_coeff = 25
            b% s2% solver_max_tries_before_reject = 40
            b% s2% max_tries_for_retry = 40
            b% s2% tiny_corr_coeff_limit = 1000    
            b% s2% corr_coeff_limit = 0.2d0
            b% s2% ignore_too_large_correction = .true.
            b% s2% ignore_min_corr_coeff_for_scale_max_correction = .true.
            b% s2% use_gold_tolerances = .true.
            b% s2% use_gold2_tolerances = .false.
            b% s2% gold_solver_iters_timestep_limit = 30
            b% s2% gold_iter_for_resid_tol3 = 10
            b% s2% gold_tol_residual_norm3 = 1d-6
            b% s2% gold_tol_max_residual3 = 1d-3
            b% s2% tol_max_correction = 1d-2
            b% s2% tol_correction_norm = 1d-3
            b% s2% max_corr_jump_limit = 1d99
            b% s2% max_resid_jump_limit = 1d99

            b% s2% scale_max_correction = 0.1d0
            b% s2% ignore_min_corr_coeff_for_scale_max_correction = .true.
            b% s2% ignore_species_in_max_correction = .true.

            ! b% s1% drag_coefficient = 1.0d0
            ! b% s1% min_q_for_drag = 0.5
            b% s1% make_gradr_sticky_in_solver_iters = .true.

            ! b% s2% drag_coefficient = 1.0d0
            ! b% s2% min_q_for_drag = 0.5
            b% s2% make_gradr_sticky_in_solver_iters = .true.
         end if
         

         ! Accretion method. 
         ! Note that b% does not have an s% x_*ctrl equivalent, so we use the b% s2% double-pointer and put the switch in inlist2
         ! We can also define this as a global variable, but it's good to call from the inlists in case of restarts. 
         accretion_method = b% s2% x_integer_ctrl(1) 

         ! For ease of switching methods, 
         select case(accretion_method)
         case(0)
            write(*,*) 'using beta = ', b% mass_transfer_beta, 'from inlists'
         case(1)
            call switch_beta_at_omega_crit(binary_id, ierr)
         case(2)
            call set_beta_based_on_thermal_timescale(binary_id, ierr)

         ! case(3)
         !    call switch_beta_at_omega_crit_with_AM_loss(binary_id, ierr)

         end select 

      end function extras_binary_finish_step
      


      subroutine switch_beta_at_omega_crit(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr

         ! Declare your variables! 
         real(dp) :: added_angular_momentum, total_Irot, omega_crit, Eddington_factor, total_J_accreted
         integer :: k, nz

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  

         total_J_accreted = b% xtra(1) ! calling and storing total_J_accreted in a global variable accessible to mesabinary 

         ! Implementation for accretion only; if star 2 is losing more mass via winds than it gains, don't add any AM
         if (b% component_mdot(2) > 0) then
            added_angular_momentum = b% component_mdot(2) * b% s2% dt * sqrt(standard_cgrav * b% s2% mstar * b% s2% photosphere_r*Rsun)
         else
            added_angular_momentum = 0d0
         end if 

         total_J_accreted = total_J_accreted + added_angular_momentum

         b% xtra(1) = total_J_accreted ! Update global variable
         
         total_Irot = 0d0
         do k = 1, b% s2% nz
            total_Irot = total_Irot + 2d0/3d0*b% s2% dm(k) * b% s2% r(k) * b% s2% r(k)
         end do 

         ! Eddington_factor = 1
         Eddington_factor = (1-(b% s2% L(1))/(4*pi*standard_cgrav*clight* b% s2% star_mass*Msun/0.34))
         write(*,*) 'Eddington factor of second star is', Eddington_factor

         omega_crit = sqrt(Eddington_factor * standard_cgrav * b% s2% mstar/(b% s2% photosphere_r * b% s2% photosphere_r * b% s2% photosphere_r*Rsun*Rsun*Rsun))
         
         if (total_J_accreted/total_Irot >= omega_crit) then
            write(*,*) "We have reached super critical rotation ! "   
            b% mass_transfer_beta = 1d0
            b% lxtra(ilx_reached_critical) = .true. ! Throw it as a flag
         else
            b% mass_transfer_beta = 0d0
         end if 

         write(*,*) 'called accretion method 1'
         write(*,*)  'beta=', b% mass_transfer_beta

      end subroutine switch_beta_at_omega_crit


      subroutine switch_beta_at_omega_crit_with_AM_loss(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr

         ! Declare your variables! 
         real(dp) :: added_angular_momentum, total_Irot, omega_crit
         integer :: k, nz

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  

         ! We can lose AM via winds... Maybe we can accept a little more mass once we blow some stuff off. 


         if (total_J_accreted >= 0) then  ! If
            added_angular_momentum = b% component_mdot(2) * b% s2% dt * sqrt(standard_cgrav * b% s2% mstar * b% s2% photosphere_r*Rsun)  ! Note that this can be negative if b% component_mdot(2) < 0
         else
            added_angular_momentum = 0d0
         end if 

         total_J_accreted = max(0d0, total_J_accreted + added_angular_momentum) ! Never make total angular momentum negative

         write(*,*) 'added J', added_angular_momentum
         write(*,*) 'total J', total_J_accreted

         ! find total moment of inertia of the star -- I'm sure there's a better way

         total_Irot = 0d0
         do k = 1, b% s2% nz
            total_Irot = total_Irot + 2d0/3d0*b% s2% dm(k) * b% s2% r(k) * b% s2% r(k) 
         end do 

         omega_crit = sqrt(standard_cgrav * b% s2% mstar/(b% s2% photosphere_r*b% s2% photosphere_r*b% s2% photosphere_r*Rsun*Rsun*Rsun))
         ! For solid body, omega = angular momentum divided by total inertia
         if (total_J_accreted/total_Irot >= omega_crit) then

            write(*,*) "We have reached super critical rotation ! "   
            b% mass_transfer_beta = 1d0
            b% mass_transfer_alpha = 0d0
            b% mass_transfer_gamma = 0d0
            b% lxtra(ilx_reached_critical) = .true. ! Throw it as a flag
         else 
            b% mass_transfer_beta = 0d0
         end if 

         write(*,*) 'called accretion method 2'
         write(*,*)  b% mass_transfer_beta
      end subroutine switch_beta_at_omega_crit_with_AM_loss

       
      subroutine set_beta_based_on_thermal_timescale(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr

         ! Declare your variables! 
         real(dp) :: CC_thermal_time, thermal_mdot 
         real(dp) :: thermal_scale = 1d0
         integer :: k, nz, pid, k_cc

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  

         
         ! We want to find the thermal time at the edge of the core, so let's first find the edge of the convective core:
         
         ! Go from core outwards, find the convective core boundary. 
         !For now, I'm going to define this as "when there's no mixing." If overshoot, then this is the outer edge of the overshooting region... 
         ! You could instead do b% s2% mixing_type(k) != 0
         do k = b% s2% nz, 1, -1 
            if (b% s2% mixing_type(k) < 1) then 
               write(*,*) 'no mixing at', k, ' of ', b% s2% nz
               k_cc = k  ! We may want to use k as a dummy variable again in another loop so let's store this as k_convective_core, k_cc
               exit ! Exit the loop
            end if 
         end do 
         
         !!! Here's one way to do it from calling what's in the profile column
         ! pid = star_get_profile_id(b% s2, 'thermal_time_to_surface')
         ! CC_thermal_time = star_get_profile_val(b% s2, pid, k_cc)
         
         !!! Here's the direct calculation -- these two should be identical  integral cp T dm / L 
         CC_thermal_time = sum(b% s2% dm(1:k_cc)*b% s2% cp(1:k_cc)*b% s2% T(1:k_cc))/b% s2% L(1)
         thermal_mdot = b% s2% star_mass * Msun / CC_thermal_time

         ! Some write statements for debugging
         ! write(*,*) 'at zone k=', k, 'thermal time is', star_get_profile_val(b% s2, pid, k)/31560000
         ! write(*,*) 'thermal mdot ', thermal_mdot !*secyer/Msun
         ! write(*,*) 'component mdot', b% component_mdot(2)

         ! Set accretion rate for the next timestep: 
         if (b% component_mdot(2) <= thermal_scale*thermal_mdot) then
            b% mass_transfer_beta = 0 
            write(*,*) 'let beta be efficient, should be 0'
            write(*,*) 'for next step, beta is:', b% mass_transfer_beta
         else
            b% mass_transfer_beta = 1 - (thermal_scale*thermal_mdot)/b% component_mdot(2) 
            write(*,*) 'let beta be inefficient at 1 - thermal mdot/star mdot'
            write(*,*) 'for next step, beta is:', b% mass_transfer_beta
         end if  
         write(*,*) 'called accretion method 2'
         write(*,*)  b% mass_transfer_beta

      end subroutine set_beta_based_on_thermal_timescale

       





      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if      
         
 
      end subroutine extras_binary_after_evolve     
      
      end module run_binary_extras
