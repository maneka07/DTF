!
!  Copyright (C) 2013, Northwestern University
!  See COPYRIGHT notice in top-level directory.
!
!  $Id: runtime_m.f90 4358 2017-06-24 16:33:27Z wkliao $

      module runtime_m
      implicit none

      integer i_time        !time step counter
      integer method        ! 0: blocking APIs, 1: nonblocking
      logical restart       ! new run or restart switch:
                            ! .FALSE. for new run, .TRUE. for restart
      integer i_time_end    !ending time step

      logical io_one_species_at_a_time ! whether read/write one species at a time

      double precision time          !current time
      double precision tstep         !timestep (non-dimensional)
      double precision time_save     !time at which to write savefiles (seconds)
      double precision time_save_inc !increment by which to write savefiles (seconds)

      double precision time_ref      !reference time (s)

      character(len=20) run_title    !unique title of run

      end module runtime_m
