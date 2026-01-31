!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
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
