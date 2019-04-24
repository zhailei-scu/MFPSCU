program MEANFIELDNUCLEATION_IMPLANT_OLD
    use NUCLEATION_IMPLANT

    implicit none

    call InitSimu_IMPLANT()

    call NucleationSimu_IMPLANT()

    pause
    stop

end program
