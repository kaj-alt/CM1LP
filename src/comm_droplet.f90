!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Jian Sun - 09/28/2022:                                         !
!    This module is generated following John Dennis's suggestion !
!    to separate the MPI communication of droplets from other    !
!    MPI calls in the existing comm.F code                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module comm_droplet_module

implicit none

private

public :: comm_droplet_number,comm_droplet_value
public :: setupIndexPointers, setupDepartDroplet
public :: makeContiguous
public :: CollectDropCount, CollectDropInfo

contains



end module comm_droplet_module

