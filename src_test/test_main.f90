include "material.f90"
program test_main
    use material_mod
    implicit none
    type(input_reader) :: reader
    type(cell_information) :: cell

    call reader%read_input_file()
    call cell%read_information(reader%InputFile)
    call cell%print_information()
    call cell%build_unit_cell()
end program test_main