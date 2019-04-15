program MEANFIELDNUCLEATION_GRUB
    use NUCLEATION_GRUB_Test

    implicit none

    call InitSimu()

    call NucleationSimu()

    pause
    stop

end program
