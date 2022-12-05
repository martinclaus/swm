module adios2
    use types
    use logging
    use generic_list, only: list_node_t, list_data, list_insert, list_next, list_get, list_free, list_init
    use vars_module, only: getVarPtrFromRegister, variable_t

#include "io.h"

    private
    public :: initAdios2, finishAdios2, stepAdios2


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief pushes a slice of data to the server
    !------------------------------------------------------------------
    interface pushSlice
        module procedure pushSlice1D
        module procedure pushSlice2D
        module procedure pushSlice3D
    end interface pushSlice

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Handle to an adios2 IO/transport/??? object
    !------------------------------------------------------------------
    type streamHandle
        integer(KINT) :: id !< Some numeric ID to the output stream
    end type streamHandle

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Stores all that is needed to stream data from a variable to the IO application
    !------------------------------------------------------------------
    type publisher
        integer(KINT)          :: id      !< Stream task identifyer
        type(streamHandle)     :: stream  !< handle to an ADIOS2 transport (??) object
        character(len=CHARLEN) :: varname !< name of the variable to stream as declared in `addToRegister` calls
        real(KDOUBLE), DIMENSION(:,:, :), pointer  :: varData3D=>null()   !< 3D data, which should be send
        real(KDOUBLE), DIMENSION(:,:), pointer     :: varData2D=>null()   !< 2D data, which should be send
        real(KDOUBLE), dimension(:), pointer       :: varData1D=>null()   !< 1D data, which should be send
    end type publisher

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Container to store streamIO pointer for generic linked lists
    !------------------------------------------------------------------
    type :: streamIO_ptr
        type(publisher), pointer   :: publisher=>null()  !< Pointer to publisher object
    end type streamIO_ptr
 
    type(list_node_t), pointer    :: publisherList=>null()  !< Head node of diagnostic task linked list


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Subroutines applicable to all elements of a list
    !------------------------------------------------------------------
    abstract interface
        subroutine mapable(self)
            import publisher
            type(publisher), pointer, intent(inout) :: self
        end subroutine mapable
    end interface 

    contains

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief  Initialize all module variables
        !------------------------------------------------------------------
        subroutine initAdios2
            implicit none
            call initPublisherList
            call printPublishSummary
        end subroutine initAdios2

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief  Tear down the module
        !------------------------------------------------------------------
        subroutine finishAdios2
            implicit none
            call finishPublisherList
        end subroutine finishAdios2

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief  Publish at every step
        !------------------------------------------------------------------
        subroutine stepAdios2
            implicit none
            call for_each_publisher(push)
        end subroutine stepAdios2

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief  Read config to create publisher objects
        !------------------------------------------------------------------
        subroutine initPublisherList
            implicit none
            integer(KINT)  :: io_stat
            TYPE(list_node_t), POINTER :: currentNode
            type(streamIO_ptr) :: publisher_ptr
            integer(KINT)  :: n_list=0

            print *, "DEBUG1"

            !< Read namelist
            open(UNIT_ADIOS2_NL, file=ADIOS2_NL, iostat=io_stat)
            do
                print *, "DEBUG2"
                io_stat = createFromNamelist(UNIT_ADIOS2_NL, n_list, publisher_ptr%publisher)
                ! No list found or IO error encountered
                if (io_stat .NE. 0) exit
                ! create and add task
                print *, "DEBUG3"
                if (.NOT.ASSOCIATED(publisherList)) then
                    call list_init(publisherList,TRANSFER(publisher_ptr, list_data))
                    currentNode => publisherList
                else
                    call list_insert(currentNode, TRANSFER(publisher_ptr,list_data))
                    currentNode => list_next(currentNode)
                end if
            end do
            close(UNIT_ADIOS2_NL)
        end subroutine initPublisherList

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief  Create a publisher object from a namelist in an opened file
        !------------------------------------------------------------------
        function createFromNamelist(unit, n_list, self) result(io_stat)
            integer, intent(in) :: unit                     !< IO Unit of file to read from
            integer(KINT), intent(inout) :: n_list          !< Number of sucessfully read lists
            type(publisher), pointer, intent(inout) :: self  !< Pointer to created publisher object
            integer :: io_stat                              !< returned status code of file read operation
            character(len=CHARLEN) :: varname               !< Name of the variable to publish
            
            namelist / publish_nl / varname
            read (unit, nml=publish_nl, iostat=io_stat)
            ! early return in case of error
            if (io_stat .ne. 0) return

            self => createPublisher(n_list, varname)

            ! increment number of read lists
            n_list = n_list + 1
        end function createFromNamelist


        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief  Create a new publisher object
        !!
        !! Allocate a new publisher, look up and set a reference to the variable
        !! data from the global register and create a streamHandle object 
        !------------------------------------------------------------------
        function createPublisher(id, varname) result(self)
            integer(KINT), intent(in) :: id
            character(len=CHARLEN), intent(in) :: varname
            type(publisher), pointer :: self
            integer :: alloc_stat
            type(variable_t), pointer :: var_obj

            allocate(self, stat=alloc_stat)
            if (alloc_stat .ne. 0) call log_alloc_fatal(__FILE__,__LINE__)

            self%id = id
            self%varname = varname

            ! get pointer to variable data
            var_obj => getVarPtrFromRegister(varname)
            if (.not. associated(var_obj)) call log_fatal("Variable '" // trim(varname) // "' not found in register")
            
            if (associated(var_obj%data3d)) then
                self%varData3D => var_obj%data3d
            else if (associated(var_obj%data2d)) then
                self%varData2D => var_obj%data2d
            else if (associated(var_obj%data1d)) then
                self%varData1D => var_obj%data1d
            end if

            ! create transport object
            self%stream = createStream(id)

        end function createPublisher

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief  Create a stream object
        !------------------------------------------------------------------
        function createStream(id) result(stream)
            integer, intent(in) :: id
            type(streamHandle) :: stream
            stream%id = id
        end function createStream

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief  Clean up the list of publishers
        !------------------------------------------------------------------
        subroutine finishPublisherList
            implicit none
            call for_each_publisher(finishPublisher)
            call list_free(publisherList)
        end subroutine finishPublisherList

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief  Disassociate all pointers of a publisher
        !------------------------------------------------------------------
        subroutine finishPublisher(self)
            type(publisher), pointer, intent(inout) :: self
            integer :: alloc_error
            nullify(self%varData3D)
            nullify(self%varData2D)
            nullify(self%varData1D)

            deallocate(self, stat=alloc_error)
            if ( alloc_error .NE. 0 ) call log_error("Deallocation failed in " // __FILE__)
            nullify(self)
        end subroutine finishPublisher


        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief  Push data to IO Manager
        !------------------------------------------------------------------
        subroutine push(self)
            type(publisher), pointer, intent(inout) :: self

            if (.not. associated(self)) then
                call log_error("Use of null pointer when publishing data")
                return
            end if

            if (associated(self%varData3D)) then
                call pushSlice(self%varname, self%varData3D, self%stream)
            else if (associated(self%varData2D)) then
                call pushSlice(self%varname, self%varData2D, self%stream)
            else if (associated(self%varData1D)) then
                call pushSlice(self%varname, self%varData1D, self%stream)
            end if   
        end subroutine push

        subroutine pushSlice1D(name, data, stream)
            implicit none
            character(len=*), intent(in)         :: name
            real(KDOUBLE), dimension(:), pointer :: data
            type(streamHandle), intent(inout)   :: stream
            character(80) :: id_str

            write(id_str, "(I80)") stream%id
            
            call log_info("Push 1D var '" // trim(name) // "' to server via stream " // adjustl(id_str))

        end subroutine pushSlice1D

        subroutine pushSlice2D(name, data, stream)
            implicit none
            character(len=*), intent(in)            :: name
            real(KDOUBLE), dimension(:, :), pointer :: data
            type(streamHandle), intent(inout)   :: stream
            character(80) :: id_str
            
            write(id_str, "(I80)") stream%id
            
            call log_info("Push 2D var '" // trim(name) // "' to server via stream " // adjustl(id_str))
        end subroutine pushSlice2D

        subroutine pushSlice3D(name, data, stream)
            implicit none
            character(len=*), intent(in)            :: name
            real(KDOUBLE), dimension(:, :, :), pointer :: data
            type(streamHandle), intent(inout)   :: stream
            character(80) :: id_str
            
            write(id_str, "(I80)") stream%id
            
            call log_info("Push 3D var '" // trim(name) // "' to server via stream " // adjustl(id_str))
        end subroutine pushSlice3D


        subroutine printPublisher(self)
            type(publisher), pointer, intent(inout) :: self
            character(len=CHARLEN) :: msg

            if (.not. associated(self)) then
                call log_error("Try to print dangling publishing task")
                return
            end if

            write(msg, '(3(A), X, I2)') "Publish variable '", trim(self%varname), "' via stream with id", self%stream%id
            call log_info(msg)   
        end subroutine printPublisher

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Summary of published variables
        !------------------------------------------------------------------
        SUBROUTINE printPublishSummary
            call log_info("Publishing Summary:")
            call for_each_publisher(printPublisher)
          END SUBROUTINE printPublishSummary

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Apply subroutine for each publisher in publisherList
        !------------------------------------------------------------------
          subroutine for_each_publisher(apply)
            procedure(mapable) :: apply
            TYPE(list_node_t), POINTER :: currentNode
            TYPE(streamIO_ptr)         :: publisher_ptr
            currentNode => publisherList
            do while (associated(currentNode))
                if (associated(list_get(currentNode))) publisher_ptr = transfer(list_get(currentNode), publisher_ptr)
                if (ASSOCIATED(publisher_ptr%publisher)) then
                    call apply(publisher_ptr%publisher)
                end if
                currentNode => list_next(currentNode)
            end do            
        end subroutine for_each_publisher

end module adios2