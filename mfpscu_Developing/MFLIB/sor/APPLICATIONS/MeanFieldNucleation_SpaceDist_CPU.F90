program MEANFIELDNUCLEATION_SPACEDIST_CPU
    use MF_SimBoxArray_AppShell_CPU
    !use NUCLEATION_GRUB

    implicit none

    integer::nmpi
    integer::procid
    !---Body----
    nmpi = 1
    procid = 1

    !---Exclute the main process---
    call AppShell_Main_CPU(nmpi,procid)

    !call InitSimu()

    !call NucleationSimu()

    write(*,*) "---End MEANFIELDNUCLEATION_SPACEDIST---"

    pause
    stop

end program MEANFIELDNUCLEATION_SPACEDIST_CPU
