!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Main program
!!
!! Outline of program workflow:
!!   -# initialise needed modules
!!   -# call diag to write initial conditions to file
!!   -# loop in time
!!      - call timestep routines of modules
!!      - call advance routines
!!      - call diagnosic routine
!!   -# deinitialise modules
!!
!! @par Include Files:
!! model.h, io.h
!!
!! @par Uses:
!! vars_module, calc_lib, diag_module, domain_module,
!! dynFromFile_module, swm_module, tracer_module, time_integration_module
!------------------------------------------------------------------
PROGRAM main_program
#include "model.h"
#include "io.h"
  use app, only: AbstractApp
  implicit none

  call main()
  STOP 0

  CONTAINS

    subroutine main()
      class(AbstractApp), pointer :: application

      application => make_app()
      call application%run(100)
      
      deallocate(application)
    end subroutine main

    function make_app() result(application)
      use app, only: new_default_app_builder, AbstractAppBuilder
      use io_module, only: Io
      use io_netcdf, only: make_netcdf_io
      use domain_module, only: make_domain_component, Domain
      use vars_module, only: make_variable_register, VariableRepository
      use calc_lib, only: make_calc_component, Calc
      use time_integration_module, only: make_time_integration_component
      use swm_module, only: make_swm_component
      use dynFromIo_module, only: make_dynFromIo_component
      use tracer_module, only: make_tracer_component
      use diag_module, only: make_diag_component

      class(AbstractApp), pointer :: application
      class(AbstractAppBuilder), pointer :: app_factory
      class(Io), pointer :: io_comp
      class(Domain), pointer :: domain_comp
      class(VariableRepository), pointer :: repo
      class(Calc), pointer :: calc_comp
      
      app_factory => new_default_app_builder()

      io_comp => make_netcdf_io()
      call app_factory%add_component(io_comp)

      domain_comp => make_domain_component(io_comp)
      call app_factory%add_component(domain_comp)

      repo => make_variable_register(domain_comp)
      call app_factory%add_component(repo)

      calc_comp => make_calc_component(domain_comp)
      call app_factory%add_component(calc_comp)

      call app_factory%add_component(make_time_integration_component(repo))

#ifdef SWM
      call app_factory%add_component( &
          make_swm_component(domain_comp, repo, io_comp, calc_comp)  &
      )
#endif

#ifdef DYNFROMFILE
      call app_factory%add_component(  &
        make_dynFromIo_component(domain_comp, repo, io_comp, calc_comp)
      )
#endif

#ifdef TRACER
      cal app_factory%add_component( &
        make_tracer_component(domain_comp, repo, io_comp)
      )
#endif

      call app_factory%add_component(  &
        make_diag_component(domain_comp, repo, io_comp, calc_comp)  &
      )

      application => app_factory%build()

      deallocate(app_factory)
    end function make_app

END PROGRAM main_program
