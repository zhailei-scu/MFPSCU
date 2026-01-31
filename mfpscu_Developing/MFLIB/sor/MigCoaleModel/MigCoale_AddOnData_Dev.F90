!*********************************************************************************!
!--- Description:
!--- Author : Lei Zhai, Insti. of Nucle. Sci. and Tech., Sichuan University
!--- Email : zhaileiytp@163.com
!--- data: From 2017-09 to 2021-12
!--- License: MIT License (There are no limitation for anyone to use/modify/sale of this code, be happy to use this code)
!--- Please ref "Lei Zhai, Chaoqiong Ma, Jiechao Cui and Qing Hou, GPU-based acceleration of Monte Carlo simulations for migration-coalescence evolution of gas bubbles in materials
!---                        Modelling Simul. Mater. Sci. Eng. 2019. 27 055008,https://iopscience.iop.org/article/10.1088/1361-651X/ab1d14
!*********************************************************************************!
module MIGCOALE_ADDONDATA_DEV
    use MIGCOALE_ADDONDATA_HOST
    implicit none

    real(kind=KMCDF),constant::dm_RNFACTOR
    real(kind=KMCDF),constant::dm_FREESURDIFPRE
    real(kind=KMCDF),constant::dm_GBSURDIFPRE

    contains

    subroutine CopyAddOnDataToDev()
        implicit none

        dm_RNFACTOR = m_RNFACTOR
        dm_FREESURDIFPRE = m_FREESURDIFPRE
        dm_GBSURDIFPRE = m_GBSURDIFPRE

        return
    end subroutine CopyAddOnDataToDev



end module MIGCOALE_ADDONDATA_DEV
