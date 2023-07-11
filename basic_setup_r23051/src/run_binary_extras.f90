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
         how_many_extra_binary_history_columns = 1
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

         names(1) = 'rlof_duration'
         vals(1) = b% xtra(1)
         
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
         
!        b% s1% job% warn_run_star_extras = .false.
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
         use chem_lib, only : chem_M_div_h
         use colors_lib, only : get_abs_mag_by_name
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         real(dp) :: m_div_h, abs_mag_V_1, abs_mag_V_2
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_check_model = keep_going

         ! calculate mass transfer duration
         ! adds timestep to mass transfer duration when the primary is filling its Roche Lobe
         ! due to bug with xtra this cannot be put in extras finish step
         if (b% r(1) > b% rl(1)) then 
            b% xtra(1) = b% xtra(1) + b% time_step
            write(*,*) 'Mass transfer duration (yrs) = ', b% xtra(1)
         end if

         ! calculate post-interaction lifetime
         ! if rlof has already occured and primary not filling RL then add timestep
         if ((b% r(1) < b% rl(1)) .and. (b% lxtra(ilx_reached_rlo))) then 
            b% xtra(2) = b% xtra(2) + b% time_step
            write(*,*) 'Post-interaction duration (yrs) = ', b% xtra(2)
         end if

         ! use routine from colors/public/colors_lib.f90 to calculate V band magnitude

         ! primary
         m_div_h = chem_M_div_h(b% s1% X(1), b% s1% Z(1), b% s1% job% initial_zfracs)
         abs_mag_V_1 = get_abs_mag_by_name('V', safe_log10(b% s1% Teff), &
               b% s1% photosphere_logg, m_div_h, b% s1% photosphere_L, ierr)

         ! check for extrapolation and use previous known value if happening
         if (abs_mag_V_1 >100) then 
            abs_mag_V_1 = b% xtra(3)
         else
            b% xtra(3) = abs_mag_V_1
         end if

         ! secondary
         m_div_h = chem_M_div_h(b% s2% X(1), b% s2% Z(1), b% s2% job% initial_zfracs)
         abs_mag_V_2 = get_abs_mag_by_name('V', safe_log10(b% s2% Teff), &
               b% s2% photosphere_logg, m_div_h, b% s2% photosphere_L, ierr)

         ! check for extrapolation and use previous known value if happening
         if (abs_mag_V_2 >100) then 
            abs_mag_V_2 = b% xtra(4)
         else
            b% xtra(4) = abs_mag_V_2
         end if
         
         ! print magnitude values
         write(*,*) 'primary V band mag = ', abs_mag_V_1
         write(*,*) 'secondary V band mag = ', abs_mag_V_2

         ! if post-interaction primary contributes at least 50% of the total flux add time-step to duration
         if ((b% r(1) < b% rl(1)) .and. (b% lxtra(ilx_reached_rlo)) &
               .and. (abs_mag_V_1 < abs_mag_V_2)) then 
            b% xtra(5) = b% xtra(5) + b% time_step
         end if
        
      end function extras_binary_check_model
      
      
      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
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

         
      end function extras_binary_finish_step
      
      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if 

         ! print mass transfer duration and post-interaction lifetime
         write(*,*) 'Mass transfer duration (yrs) = ', b% xtra(1)
         write(*,*) 'Fraction of total lifetime = ', b% xtra(1)/b% binary_age
         write(*,*) 'Post-interaction lifetime (yrs) = ', b% xtra(2)
         write(*,*) 'Primary dominating V band Post-interaction duration (yrs) = ', b% xtra(5)
         write(*,*) 'Thermal timescale of primary (yrs) = ', b% s1% kh_timescale
         write(*,*) 'Nuclear timescale of primary (yrs) = ', b% s1% nuc_timescale
         
      end subroutine extras_binary_after_evolve     
      
      end module run_binary_extras
