program mini_adios
    use adios2
    implicit none
    type(adios2_adios)    :: adios
    type(adios2_io)       :: io
    type(adios2_engine)   :: bp_writer
    integer :: ierr

    call adios2_init (adios, adios2_debug_mode_on, ierr)
    call adios2_declare_io (io, adios, 'SimulationOutput', ierr )
    call adios2_finalize (adios, ierr)
    
end program mini_adios