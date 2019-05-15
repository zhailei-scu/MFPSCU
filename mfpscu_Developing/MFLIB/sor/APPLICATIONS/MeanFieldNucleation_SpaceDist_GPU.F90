program MEANFIELDNUCLEATION_SPACEDIST_GPU
    use MF_SimBoxArray_AppShell_GPU
    !use NUCLEATION_GRUB

    implicit none

    integer::nmpi
    integer::procid
    !---Body----
    nmpi = 1
    procid = 1

    !---Exclute the main process---
    call AppShell_Main_GPU(nmpi,procid)

    write(*,*) "---End MEANFIELDNUCLEATION_SPACEDIST GPU---"

    pause
    stop

end program MEANFIELDNUCLEATION_SPACEDIST_GPU
