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
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         
!          b% s1% job% warn_run_star_extras = .false.
          extras_binary_startup = keep_going

         !!!!!! If modifying extras_binary_startup, do so below:
         if (.not. restart) b% xtra(11) = 0d0

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

         !!!!! Declare your variables here !!!!!
         integer :: accretion_method ! to toggle
         integer :: k, nz

         ! For ad hoc model 1
         real(dp) :: total_AM_accreted, AM_accreted_this_step, star_Irot, ad_hoc_omega_rot, Gamma_Edd, omega_crit

         ! For bonus ad hoc model 2
         real(dp) :: thermal_time, thermal_mdot
         integer :: k_cc, pid ! Useful ints thermal time
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_finish_step = keep_going

         if (b% r(1) > b% rl(1) .and. .not. b% lxtra(ilx_reached_rlo)) then  ! This is copied over from the basic setup
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
         

         !!!!! Implement the Ad Hoc Model(s) in Midilab below !!!!!


         ! For ease of switching between ad hoc models, these solutions use cases

         ! accretion_method = 1 ! Set manually for ad hoc model 1
         accretion_method = b% s2% x_integer_ctrl(11) ! Set to x_integer_ctrl(11) in inlist2

         select case(accretion_method) ! Pick your accretion method. Cases in fortran are handy for this.


         !!!!! '*********************************************'
         !!!!! '****** NO AD HOC MODEL ******'
         !!!!! '*********************************************'
         case(0) ! Using just the mass transfer rate specified in the inlists - No Ad hoc model
            !  write(*,*) 'using beta = ', b% mass_transfer_beta, 'from inlists'
            continue  ! Do nothing and keep going
         

         !!!!! '*********************************************'
         !!!!! '****** AD HOC MODEL 1 ******'
         !!!!! '*********************************************'

         case(1) ! Terminate at critical rotation -- Ad hoc model 1
            !  write(*,*) 'using ad hoc model 1: critical rotation rate'

            ! initialize total_AM_accreted
            total_AM_accreted = b% xtra(11) ! calling and storing total_AM_accreted in a global variable accessible to mesabinary 

            ! If star 2 is losing more mass via winds than it gains, don't add any AM
            if (b% component_mdot(2) > 0) then
               AM_accreted_this_step = b% component_mdot(2) * b% s2% dt * sqrt(standard_cgrav * b% s2% mstar * b% r(2))
            else
               AM_accreted_this_step = 0d0
            end if 

            ! Add the AM accreted this step to total_AM_accreted and update s% xtra(1). 
            total_AM_accreted = total_AM_accreted + AM_accreted_this_step
            ! Update the global variable. You could skip this by directly calling b% xtra(1) instead of total_AM_accreted
            b% xtra(11) = total_AM_accreted
            
            ! Calculate the moment of inertia by looping through the star
            star_Irot = 0d0
            do k = 1, b% s2% nz
               star_Irot = star_Irot + 2d0/3d0*b% s2% dm(k) * b% s2% r(k) * b% s2% r(k)
            end do 

            ! Define the ad hoc rotation rate as moment of AM divided by the moment of inertia
            ad_hoc_omega_rot = total_AM_accreted/star_Irot

            ! Calculate eddington ratio 
            Gamma_Edd = (b% s2% L(1))/(4*pi*standard_cgrav*clight*b% s2% mstar / 0.34)
            ! calculate critical rotation rate
            omega_crit = sqrt((1 - Gamma_Edd) * standard_cgrav * b% s2% mstar/(b% s2% r(1) * b% s2% r(1) * b% s2% r(1)))

            if (ad_hoc_omega_rot >= omega_crit) then
               write(*,*) '*********************************************'
               write(*,*) '****** Terminated at critical rotation ******'
               write(*,*) '*********************************************'
               ! b% mass_transfer_beta = 1d0 
               extras_binary_finish_step = terminate
            else
               b% mass_transfer_beta = 0d0
            end if 






         !!!!! '*********************************************'
         !!!!! '****** AD HOC MODEL 2 ******'
         !!!!! '*********************************************'
               
         case(2) ! Terminate at thermal time limited accretion
            !  write(*,*) 'using ad hoc model 2: exceeding thermal timescale mdot'
            
      !    ! Go from core outwards, find the convective core boundary. 
            do k = b% s2% nz, 1, -1 
               if (b% s2% mixing_type(k) < 1) then 
                  k_cc = k  ! We may want to use k as a dummy variable again in another loop so store this as k_convective_core, k_cc
                  exit ! Exit the loop
               end if 
            end do 
      
         !!! Get the thermal time from the convective core -- these two should be identical  integral cp T dm / L 
            thermal_time = sum(b% s2% dm(1:k_cc)*b% s2% cp(1:k_cc)*b% s2% T(1:k_cc))/b% s2% L(1)
            thermal_mdot = b% s2% mstar / thermal_time
   
      !    !!! Here's another way to get this info, by calling what's in the profile column
            ! pid = star_get_profile_id(b% s2, 'thermal_time_to_surface')
            ! thermal_time = star_get_profile_val(b% s2, pid, k_cc)
         
      !    ! Check if we've reached our stopping condition
            
            if (b% component_mdot(2) >= 3*thermal_mdot) then
               write(*,*) '*********************************************'
               write(*,*) '****** Terminated from thermal runaway ******'
               write(*,*) '*********************************************'
               extras_binary_finish_step = terminate
               ! b% mass_transfer_beta = 1 - (thermal_scale*thermal_mdot)/b% component_mdot(2) if you want to keep the run going but limit to the thermal rate?
            else
               b% mass_transfer_beta = 0 
            end if  

         end select ! End case selection
         
         !!!!! Implement the Ad Hoc Model(s) in Midilab above !!!!!

      end function extras_binary_finish_step


      
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
