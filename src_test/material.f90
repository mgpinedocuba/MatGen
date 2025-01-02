module material_mod

    implicit none

    type :: input_reader
        character(len=99) :: InputFile
        contains
        procedure :: read_input_file => read_input_file_sb
    end type input_reader
    
    type :: cell_information
        real(kind=8) :: LatPar_a
        real(kind=8) :: LatPar_b
        real(kind=8) :: LatPar_c
        real(kind=8),dimension(1:3) :: TrasVec_a
        real(kind=8),dimension(1:3) :: TrasVec_b
        real(kind=8),dimension(1:3) :: TrasVec_c
        integer(kind=8) :: NumAtoms
        character(len=2),dimension(:),allocatable :: AtomsT
        integer(kind=8),dimension(:),allocatable :: AtomsL
        real(kind=8),dimension(:,:),allocatable :: AtomsC
        contains
        procedure :: read_information => read_information_sb
        procedure :: print_information => print_information_sb
        procedure :: build_unit_cell => build_unit_cell_sb
    end type cell_information

    private :: read_input_file_sb
    private :: read_information_sb
    private :: print_information_sb
    private :: build_unit_cell_sb
    
    contains

    subroutine read_input_file_sb(this)
        implicit none
        class(input_reader) :: this
        integer(kind=8) :: ierr
        if ( command_argument_count() < 1 ) then
            write(*,*) "Error: Please provide the input file as a command-line argument."
            stop
        end if
        call get_command_argument(1, this%InputFile)
        ! checking if the file exist
        open(unit=99, file=trim(this%InputFile), status="old", action="read", iostat=ierr)
        if ( ierr/=0 ) then
            write(*,*) "Error: Unable to open file '", trim(this%InputFile),"'"
            stop
        end if
        close(99)
        write(*,*) "Input File: ",trim(this%InputFile)
    end subroutine read_input_file_sb
    
    subroutine read_information_sb(this, input)
        implicit none
        class(cell_information) :: this
        character(*), intent(in) :: input
        integer(kind=8) :: ierr, allostat, i
        open(unit=99, file=trim(input), status="old", action="read", iostat=ierr)
        if ( ierr/=0 ) then
            write(*,*) "Error: Unable to open file '", trim(input),"'"
            stop
        end if
        read(99,*)
        read(99,*) this%LatPar_a, this%LatPar_b, this%LatPar_c
        read(99,*) this%TrasVec_a(1:3)
        read(99,*) this%TrasVec_b(1:3)
        read(99,*) this%TrasVec_c(1:3)
        read(99,*) this%NumAtoms
        allocate(this%AtomsT(this%NumAtoms), stat=allostat)
        if ( allostat/=0 ) then ; stop ; end if
        allocate(this%AtomsL(this%NumAtoms), stat=allostat)
        if ( allostat/=0 ) then ; stop ; end if
        allocate(this%AtomsC(this%NumAtoms,3), stat=allostat)
        if ( allostat/=0 ) then ; stop ; end if
        do i=1,this%NumAtoms
            read(99,*) this%AtomsT(i), this%AtomsL(i), this%AtomsC(i,1:3)
        end do
        close(99)
    end subroutine read_information_sb

    subroutine print_information_sb(this)
        implicit none
        class(cell_information) :: this
        integer(kind=8) :: i
        write(*,*) "Lattice Parameter A: ", this%LatPar_a
        write(*,*) "Lattice Parameter B: ", this%LatPar_b
        write(*,*) "Lattice Parameter C: ", this%LatPar_c
        write(*,*) "Translational Vector A: ", this%TrasVec_a
        write(*,*) "Translational Vector B: ", this%TrasVec_b
        write(*,*) "Translational Vector C: ", this%TrasVec_c
        write(*,*) "Number of Atoms:", this%NumAtoms
        write(*,*) "Simbol | Label | Coordinates in Crystal Reference"
        do i=1,this%NumAtoms
            write(*,*) this%AtomsT(i), this%AtomsL(i), this%AtomsC(i,1:3)
        end do
    end subroutine print_information_sb

    subroutine build_unit_cell_sb(this)
        implicit none
        class(cell_information) :: this
        real(kind=8),dimension(1:3) :: Rxyz, c1, c2, c3
        integer(kind=8) :: i, ierr

        open(unit=99, file="UnitCell.xyz", action="write", iostat=ierr)
        if ( ierr /= 0 ) then
            write(*,*) "Error: Unit Cell file not generated"
            stop
        end if
        write(99,"(I9)") this%NumAtoms
        write(99,*)
        do i=1,this%NumAtoms
            c1(:) = this%LatPar_a * this%TrasVec_a(1:3)
            c2(:) = this%LatPar_b * this%TrasVec_b(1:3)
            c3(:) = this%LatPar_c * this%TrasVec_c(1:3)
            Rxyz(:) = this%AtomsC(i,1)*c1+this%AtomsC(i,2)*c2+this%AtomsC(i,3)*c3
            write(99,*) this%AtomsT(i), this%AtomsL(i), Rxyz(1:3)
        end do
        close(99)
    end subroutine build_unit_cell_sb

end module material_mod