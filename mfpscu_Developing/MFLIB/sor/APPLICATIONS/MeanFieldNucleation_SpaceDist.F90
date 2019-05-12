program MEANFIELDNUCLEATION_SPACEDIST
    use MF_SimBoxArray_AppShell_CPU
    use NUCLEATION_GRUB

    implicit none

    integer::nmpi
    integer::procid
    !---Body----
    nmpi = 1
    procid = 1

    !---Exclute the main process---
    call AppShell_Main_CPU(nmpi,procid)

    write(*,*) "---End MEANFIELDNUCLEATION_SPACEDIST---"

    pause
    stop

end program
