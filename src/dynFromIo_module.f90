!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Module for reading dynamical variables (u,v,eta) from file
!! to make offline tracer modelling possible.
!! 
!! @author Martin Claus
!!
!! This model can be used to do offline tracer calculations, etc.
!! Its purpose is to load the dynamical variables, i.e. u, v and eta,
!! from disk every time step. Alternatively, it is also possible to specify
!! the geostrophic streamfunction from which zonal and meridional velocity
!! and interface displacemnt are computed. This has the advantage of having
!! a guaranteed divergence free flow.
!!
!! @par Includes:
!! io.h, model.h
!! @par Uses:
!! types \n
!! logging, only: Logger \n
!! app, only: Component \n
!! vars_module, only: VariableRepository \n
!! domain_module, only: Domain \n
!! io_module, only: Io, Reader, HandleArgs \n
!! calc_lib, only: Calc \n
!! memchunk_module, only: MemoryChunk \n
!------------------------------------------------------------------
MODULE dynFromIo_module
#include "model.h"
#include "io.h"
  use types
  use logging, only: Logger
  use app, only: Component
  use vars_module, only: VariableRepository
  use domain_module, only: Domain
  use io_module, only: Io, Reader, HandleArgs
  use calc_lib, only: Calc
  use memchunk_module, only: MemoryChunk
  implicit none
  private

  public :: make_dynFromIo_component

  type, abstract, private :: DynInput
  contains
    procedure(call_on_self), deferred :: init
    procedure(call_on_self_with_real_arg), deferred :: get_data
  end type DynInput

  abstract interface
    subroutine call_on_self(self)
      import DynInput
      class(DynInput), intent(inout) :: self
    end subroutine

    subroutine call_on_self_with_real_arg(self, arg)
      import DynInput, KDOUBLE
      class(DynInput), intent(inout) :: self
      real(KDOUBLE), intent(in) :: arg
    end subroutine call_on_self_with_real_arg
  end interface

  type, extends(DynInput) :: PsiInput
    type(MemoryChunk) :: psi_chunk
  contains
    procedure :: init => init_psi_input
    procedure :: get_data => get_data_from_psi           !< Input stream associated with a dataset containing the geostrophic streamfunction
  end type PsiInput

  type, extends(DynInput) :: UVEInput
    type(MemoryChunk) :: eta_chunk                      !< Input stream associated with a dataset containing interace displacement
    type(MemoryChunk) :: u_chunk                        !< Input stream associated with a dataset containing zonal velocity
    type(MemoryChunk) :: v_chunk                        !< Input stream associated with a dataset containing meridional
  contains
    procedure :: init => init_uveta_input
    procedure :: get_data => get_data_from_uveta        !< Input stream associated with a dataset containing the geostrophic streamfunction
  end type UVEInput

  type, extends(Component) :: DynFromIo
    private
    class(Logger), pointer             :: log => null()
    class(Domain), pointer             :: dom => null()
    class(VariableRepository), pointer :: repo => null()
    class(Io), pointer                 :: io => null()
    class(Calc), pointer               :: calc => null()
    class(DynInput), allocatable       :: input
    real(KDOUBLE), dimension(:, :), allocatable :: eta  !< Interface displacement
    real(KDOUBLE), dimension(:, :), allocatable :: u    !< Zonal velocity
    real(KDOUBLE), dimension(:, :), allocatable :: v    !< Meridional velocity
  contains
    procedure :: initialize, &
                 finalize, &
                 step => timestep, &
                 advance
    procedure, private :: get_memchunk, reader_from_fileargs, &
                          init_psi_input, init_uveta_input, &
                          allocate_buffer, register_variables, &
                          get_data, get_data_from_uveta, get_data_from_psi
  end type DynFromIo

  contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Create a DynFromIo component
    !------------------------------------------------------------------
    function make_dynFromIo_component(log, dom, repo, io_comp, calc_comp) result(comp)
      class(Logger), target, intent(in) :: log
      class(Domain), target, intent(in) :: dom
      class(VariableRepository), target, intent(in) :: repo
      class(Io), target, intent(in) :: io_comp
      class(Calc), target, intent(in) :: calc_comp
      class(DynFromIo), pointer :: comp
      allocate(comp)
      comp%log => log
      comp%dom => dom
      comp%repo => repo
      comp%io => io_comp
      comp%calc => calc_comp
    end function make_dynFromIo_component

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Initialise component
    !!
    !! Parses namelist dynFromFile, initialise memoryChunks and allocate variable
    !! buffers.
    !------------------------------------------------------------------
    SUBROUTINE initialize(self)
      class(DynFromIo), intent(inout) :: self
      ! initial conditions will be overwritten (must be called befor any module using initial conditions is initialised)
      CHARACTER(CHARLEN) :: filename_u="", filename_v="", filename_eta="", filename_psi="",  &
                            varname_eta=OVARNAMEETA, varname_u=OVARNAMEU, varname_v=OVARNAMEV, varname_psi=OVARNAMEPSI
      integer(KINT) :: chunk_size=DEF_NT_CHUNKSIZE
      namelist / dynFromFile / &
        FileName_u, FileName_v, FileName_eta, FileName_psi, & ! Filenames of input files of dynamical variables
        varname_eta, varname_u, varname_v, varname_psi, & ! Variable names in input files
        chunk_size ! Number of timesteps to read at once
      
        ! read the namelist and close again
      open(UNIT_DYNFROMFILE_NL, file = DYNFROMFILE_NL)
      read(UNIT_DYNFROMFILE_NL, nml = dynFromFile)
      close(UNIT_DYNFROMFILE_NL)
      

      self%is_psi_input = (varname_psi .ne. "")

      ! initialise file handles
      if (self%is_psi_input) then
        call self%init_psi_input(filename_psi, varname_psi, chunk_size)
      else
        call self%init_uveta_input(  &
          filename_u, varname_u, &
          filename_v, varname_v, &
          filename_eta, varname_eta, &
          chunk_size)
      end if

      call self%allocate_buffer()
      call self%register_variables()

      ! set initial conditions
      CALL self%get_data(0._KDOUBLE)
    END SUBROUTINE initialize

    subroutine init_psi_input(self, filename, varname, chunk_size)
      class(DynFromIo) :: self
      character(*), intent(in) :: filename, varname
      integer(KINT), intent(in) :: chunk_size
      self%psi_chunk = self%get_memchunk(filename, varname, chunk_size)
    end subroutine init_psi_input

    subroutine init_uveta_input(  &
      self,  &
      filename_u, varname_u, &
      filename_v, varname_v, &
      filename_eta, varname_eta, &
      chunk_size  &
    )
      class(DynFromIo) :: self
      character(*), intent(in) :: filename_u, varname_u, &
                                  filename_v, varname_v, &
                                  filename_eta, varname_eta
      integer(KINT), intent(in) :: chunk_size
      self%u_chunk = self%get_memchunk(filename_u, varname_u, chunk_size)
      self%v_chunk = self%get_memchunk(filename_v, varname_v, chunk_size)
      self%eta_chunk = self%get_memchunk(filename_eta, varname_eta, chunk_size)
    end subroutine init_uveta_input

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Allocate buffer to store dynamical variables
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine allocate_buffer(self)
      class(DynFromIo) :: self
      integer(KINT) :: Nx, Ny
      integer :: alloc_error
      
      Nx = self%dom%Nx
      Ny = self%dom%Ny

      allocate(  &
        self%u(1:Nx,1:Ny),  &
        self%v(1:Nx,1:Ny),  &
        self%eta(1:Nx, 1:Ny),  &
        stat=alloc_error)

      if (alloc_error .ne. 0) call self%log%fatal_alloc(__FILE__, __LINE__)
    end subroutine allocate_buffer

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Register dynamical variables globally
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine register_variables(self)
      class(DynFromIo), intent(in) :: self
      call self%repo%add(self%u, "DFF_U", self%dom%u_grid)
      call self%repo%add(self%v, "DFF_V", self%dom%v_grid)
      call self%repo%add(self%eta, "DFF_ETA", self%dom%eta_grid)
    end subroutine register_variables

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Deallocate memory
    !------------------------------------------------------------------
    SUBROUTINE finalize(self)
      class(DynFromIo), intent(inout) :: self
      if (allocated(self%u)) deallocate(self%u)
      if (allocated(self%v)) deallocate(self%v)
      if (allocated(self%eta)) deallocate(self%eta)
    END SUBROUTINE finalize

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Does nothing
    !------------------------------------------------------------------
    SUBROUTINE timestep(self)
      class(DynFromIo), intent(inout) :: self
      call self%get_data(self%repo%elapsed_time())
    END SUBROUTINE timestep

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Read data from input
    !------------------------------------------------------------------
    SUBROUTINE advance(self)
      class(DynFromIo), intent(inout) :: self
      ! do nothing while advancing the application to the next time step
    END SUBROUTINE advance

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Read data from input
    !------------------------------------------------------------------
    subroutine get_data(self, time)
      class(DynFromIo), intent(inout) :: self
      real(KDOUBLE), intent(in)       :: time
      if (self%is_psi_input) then
        call self%get_data_from_psi(time)
      else
        call self%get_data_from_uveta(time)
      end if      
    end subroutine get_data

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Read data from u, v, eta input
    !------------------------------------------------------------------
    subroutine get_data_from_uveta(self, time)
      class(DynFromIo), intent(inout) :: self
      real(KDOUBLE), intent(in)       :: time
      self%u = self%u_chunk%get(time)
      self%v = self%v_chunk%get(time)
      self%eta = self%eta_chunk%get(time)
    end subroutine get_data_from_uveta

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Read data from u, v, eta input
    !!
    !! TODO: Make g configurable
    !------------------------------------------------------------------
    subroutine get_data_from_psi(self, time)
      class(DynFromIo), intent(inout) :: self
      real(KDOUBLE), intent(in)       :: time
      logical :: must_compute
      real(KDOUBLE) :: g=9.80665
      integer(KINT) :: i, j

      must_compute = (time .eq. 0 .or. .not.self%psi_chunk%is_constant())

      ! return early if nothing needs to be done
      if (.not. must_compute) return

      self%eta = self%psi_chunk%get(time)

      ! compute velocity from streamfunction
      call self%calc%evaluateStreamfunction(self%eta, self%u, self%v)

      ! compute interface displacement assuming quasi-geostrophy
      ! \eta = \psi * f / g
!$OMP parallel do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, self%dom%Ny
        do i = 1, self%dom%Nx
          self%eta(:, j) = self%dom%H_grid%f(j) * self%eta(:, j) / g
        end do
      end do
!$OMP end parallel do
    end subroutine get_data_from_psi

    type(MemoryChunk) function get_memchunk(self, filename, varname, chunk_size) result(mchunk)
      class(DynFromIo) :: self
      character(*), intent(in) :: filename, varname
      integer(KINT), intent(in) :: chunk_size
      mchunk = MemoryChunk( &
        self%reader_from_fileargs(filename, varname), &
        (/self%dom%Nx, self%dom%Ny, chunk_size/)  &
      )
    end function

    class(Reader) function reader_from_fileargs(self, filename, varname) result(var_reader)
      class(DynFromIo)         :: self
      character(*), intent(in) :: filename
      character(*), intent(in) :: varname
      allocatable :: var_reader
      type(HandleArgs) :: args
      call args%add("filename", filename)
      call args%add("varname", varname)
      var_reader = self%io%get_reader(args)
    end function

END MODULE dynFromIo_module
