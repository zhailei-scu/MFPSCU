program MEANFIELDNUCLEATION_SPACEDIST
    use NUCLEATION_GRUB_Test

    implicit none

    !call InitSimu_SpaceDist()

    !call NucleationSimu_SpaceDist()

    call InitSimu()

    call NucleationSimu()



    pause
    stop

end program
