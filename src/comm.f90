  MODULE comm_module

  implicit none
  private
  public :: nabor

  !GPU enabled routines 
  public :: comm_2d_start,comm_2dns_end,comm_2dew_end
  public :: comm_1s_start,comm_1s_end
  public :: comm_1w_start,comm_1w_end
  public :: comm_1p_start,comm_1p_end
  public :: comm_1t_start,comm_1t_end
  public :: comm_2sn_start,comm_2sn_end
  public :: comm_3s_start,comm_3s_end
  public :: comm_3t_start,comm_3t_end
  public :: comm_3w_start,comm_3w_end
  public :: comm_3u_start,comm_3u_end
  public :: comm_3v_start,comm_3v_end
  public :: comm_1s2d_start,comm_1s2d_end
  public :: comm_2we_start,comm_2we_end
  public :: comm_2d_corner
  public :: comm_1s_tend_halo,comm_1u_tend_halo,comm_1v_tend_halo
  public :: getcorneru,getcornerv,getcornert,getcornerw,getcorner
  public :: getcorneru3,getcornerv3,getcornerw3
  public :: prepcorners3,prepcornert,prepcorners
  public :: comm_all_s
  public :: getcorner3_2d

  public :: sync

  CONTAINS

! Synchronization routine 
     subroutine sync()
      use mpi
      integer :: ierr
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
     end subroutine

!-----------------------------------------------------------------------
!  message passing routines
!-----------------------------------------------------------------------


      integer function nabor(i,j,nx,ny)
      implicit none
      integer i,j,nx,ny
      integer newi,newj

      newi=i
      newj=j

      if ( newi .lt.  1 ) newi = nx
      if ( newi .gt.  nx) newi = 1

      if ( newj .lt.  1 ) newj = ny
      if ( newj .gt.  ny) newj = 1

      nabor = (newi-1) + (newj-1)*nx

      end function nabor

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine prepcorners(s,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,comm)
      use input, only : ib,ie,jb,je,kb,ke,imp,jmp,kmp,cmp,rmp,kmt,ni,nj,nk, &
          tbc,bbc,cgs1,cgs2,cgs3,cgt1,cgt2,cgt3
      use bc_module, only: bcs,bcs2
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: s
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      integer, intent(inout), dimension(rmp) :: reqs_p
      integer, intent(in) :: comm

      integer :: i,j
      logical, parameter :: Debug = .FALSE.

      !$acc data present(s,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2, &
      !$acc              pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2)

!--------------------------------------------!
!  This subroutine is ONLY for parcel_interp !
!--------------------------------------------!

      IF( comm.eq.1 )THEN
        call bcs(s)
      ENDIF


      IF( bbc.eq.1 .or. bbc.eq.2 .or. bbc.eq.3 )THEN
        ! extrapolate:
        !$omp parallel do default(shared) private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present)
        do j=0,nj+1
        do i=0,ni+1
          s(i,j,0) = cgs1*s(i,j,1)+cgs2*s(i,j,2)+cgs3*s(i,j,3)
        enddo
        enddo
      ENDIF

      IF( tbc.eq.1 .or. tbc.eq.2 )THEN
        ! extrapolate:
        !$omp parallel do default(shared) private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) 
        do j=0,nj+1
        do i=0,ni+1
          s(i,j,nk+1) = cgt1*s(i,j,nk)+cgt2*s(i,j,nk-1)+cgt3*s(i,j,nk-2)
        enddo
        enddo
      ENDIF

      !$acc end data

      end subroutine prepcorners

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine prepcorners3(s,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,  &
                                n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,reqs_s,comm)
      use input, only : ib,ie,jb,je,kb,ke,imp,jmp,kmp,cmp,rmp,kmt,ni,nj,nk, &
          tbc,bbc,cgs1,cgs2,cgs3,cgt1,cgt2,cgt3
      use bc_module, only: bcs
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: s
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, intent(inout), dimension(cmp,cmp,kmt+1) :: n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2
      integer, intent(inout), dimension(rmp) :: reqs_s
      integer, intent(in) :: comm
      logical, parameter :: Debug = .FALSE.
      integer :: i,j
      !$acc data present (s,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32, &
      !$acc               n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2)


      IF( comm.eq.1 )THEN
        call bcs(s)

      ENDIF

      IF( bbc.eq.1 .or. bbc.eq.2 .or. bbc.eq.3 )THEN
        ! extrapolate:
!$omp parallel do default(shared)  &
!$omp private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
        do j=jb,je
        do i=ib,ie
          s(i,j,0) = cgs1*s(i,j,1)+cgs2*s(i,j,2)+cgs3*s(i,j,3)
        enddo
        enddo
      ENDIF

      IF( tbc.eq.1 .or. tbc.eq.2 )THEN
        ! extrapolate:
!$omp parallel do default(shared)  &
!$omp private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
        do j=jb,je
        do i=ib,ie
          s(i,j,nk+1) = cgt1*s(i,j,nk)+cgt2*s(i,j,nk-1)+cgt3*s(i,j,nk-2)
        enddo
        enddo
      ENDIF

      !$acc end data

      end subroutine prepcorners3

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine prepcornert(t,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                               tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,comm)
      use input, only : ib,ie,jb,je,kb,ke,imp,jmp,kmp,cmp,rmp,kmt,ni,nj,nk, &
          tbc,bbc,cgs1,cgs2,cgs3,cgt1,cgt2,cgt3
      use bc_module, only: bcw, bct2
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: t
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(cmp,jmp,kmt) :: tkw1,tkw2,tke1,tke2
      real, intent(inout), dimension(imp,cmp,kmt) :: tks1,tks2,tkn1,tkn2
      integer, intent(inout), dimension(rmp) :: reqs_p
      integer, intent(in) :: comm
      logical, parameter :: Debug = .FALSE.
      integer :: i,j

!--------------------------------------------!
!  This subroutine is ONLY for parcel_interp !
!--------------------------------------------!

      !$acc data present (t,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2, &
      !$acc               tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2)

      if(Debug) print *,'prepcornert'
      IF( comm.eq.1 )THEN
        call bcw(t,0)
      ENDIF


      !$acc end data

      end subroutine prepcornert




      subroutine SetMsgParams(openacc,gdirect,device,gpudirect)

         logical, intent(inout) :: openacc, gdirect
         logical, optional :: device,gpudirect

      ! Set logical flags if subroutine is called on GPU or CPU resident
      ! data. The default to for GPU resident data.  If the device flag 
      ! is not set, it is assumed that it is GPU resident.
      if(present(device)) then
         if(device) then
            openacc = .true.
            if(present(gpudirect)) then
               gdirect = gpudirect
            else
               gdirect = .true.
            endif
         else
            openacc = .false.
            gdirect = .false.
         endif
      else
         openacc = .true.
         if(present(gpudirect)) then
            gdirect = gpudirect
         else
            gdirect = .true.
         endif
      endif

      end subroutine SetMsgParams


  END MODULE comm_module

