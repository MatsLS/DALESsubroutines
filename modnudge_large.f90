!> \file modnudge.f90
!!  Nudges theta_l and q_t profiles to the initial profiles on a timescale tnudge1d
!!  Also nudges theta_l and q_t (possibility for u_i) to prescribed fields on a timescale tnudge3d
!>

!>
!!  Nudges theta_l and q_t profiles to the initial profiles on a timescale tnudget1d
!!  Also nudges theta_l and q_t (possibility for u_i) to prescribed fields on a timescale tnudge3d
!>
!!  \author Thijs Heus,MPI-M, Mats Steerneman
!!  \par Revision list
!!  \todo Documentation
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!



module modnudge


implicit none
PRIVATE
PUBLIC :: initnudge, nudge,exitnudge
SAVE
  real(kind=8), dimension(:,:,:,:), allocatable :: thlnudge3D,qtnudge3D!,unudge3D,vnudge3D,wnudge3D !MS: these arrays contain the 3D fields to be nudged towards, as given in nudge3D.inp.iexpnr.
  real, dimension(:,:)  , allocatable :: tnudge,unudge1D,vnudge1D,wnudge1D,thlnudge1D,qtnudge1D !MS: these arrays contain the 1D profiles to be nudged towards, as given in nudge.inp.iexpnr (standard DALES). 
  real, dimension(:,:)  , allocatable :: t1Dnudge,t3Dnudge !MS: these arrays are the nudging time scales for 1D and 3D nudging.
  real, dimension(:)    , allocatable :: timenudge
  real(kind = 8), dimension(:) , allocatable :: tnudgearr

  real :: t1Dnudgefac = 10.                    !MS: multiplication factor for the nudging time of 1D tendency. (default = 10s)
  real :: t3Dnudgefac = 10.                    !MS: multiplication factor for the nudging time of 3D tendency. (default = 10s)
  real :: tnudgestop = 7200.                   !MS: nudging only takes place before this time. (default = 7200s)
  real :: tnudgestart = 0.                     !MS: nudging only takes place after this time. (default = 0s)

  logical :: lnudge = .false.                  !MS: boolean given in namoptions that dictates if the run has nudging (standard DALES).
  logical :: lunudge1D,lvnudge1D,lwnudge1D,lthlnudge1D,lqtnudge1D !MS: booleans that dictate if a variable should have 1D nudging (standard DALES).
  logical :: lthlnudge3D,lqtnudge3D,end_of_file!,lunudge3D,lvnudge3D,lwnudge3D !MS: booleans that dictate if a variable should have 3D nudging.
  integer :: ntnudge = 10000
  integer :: ntnudge3D = 10                    !MS: Amount of time points at which desired fields for nudging are provided. (default = 10)
  integer :: knudgestop = 1                    !MS: knudgestop indicates until what level the 3D forcing takes place. (default = 1)
  integer :: knudgestart = 1                   !MS: knudgestart indicates at what level the 3D forcing starts. (default = 1)
contains
  subroutine initnudge
    use modmpi,   only :myid,mpierr,comm3d,mpi_logical,D_MPI_BCAST
    use modglobal,only :ifnamopt,fname_options,runtime,cexpnr,ifinput,k1,kmax,checknamelisterror,itot,jtot !MS: added itot, jtot
    implicit none
    integer :: ierr,k,t,i,j
    logical :: file_exists
    real,allocatable,dimension(:) :: height
    character(1) :: chmess1
    namelist /NAMNUDGE/ &
       lnudge,knudgestop,t1Dnudgefac,t3Dnudgefac,tnudgestop,knudgestart,tnudgestart,ntnudge3D !MS: Read all (new) inputs
    allocate(tnudge(k1,ntnudge),unudge1D(k1,ntnudge),vnudge1D(k1,ntnudge),wnudge1D(k1,ntnudge),&
           thlnudge1D(k1,ntnudge),qtnudge1D(k1,ntnudge))
    allocate(t1Dnudge(k1,ntnudge),t3Dnudge(k1,ntnudge)) !MS: allocate nudging time scales
    allocate(timenudge(0:ntnudge), height(k1))
    tnudge     =0
    unudge1D   =0
    vnudge1D   =0
    wnudge1D   =0
    thlnudge1D =0
    qtnudge1D  =0
    timenudge  =0
    tnudgearr  =0
    !unudge3D   =0
    !vnudge3D   =0
    !wnudge3D   =0
    thlnudge3D =0 !MS zero by default
    qtnudge3D  =0 !MS zero by default

    if(myid==0)then

      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMNUDGE,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMNUDGE')
      write(6 ,NAMNUDGE)
      close(ifnamopt)
    end if
    !MS: Broadcast important variables
    call D_MPI_BCAST(lnudge     , 1,0,comm3d,mpierr)
    call D_MPI_BCAST(knudgestop , 1,0,comm3d,mpierr)
    call D_MPI_BCAST(tnudgestop , 1,0,comm3d,mpierr)
    call D_MPI_BCAST(tnudgestart, 1,0,comm3d,mpierr)
    call D_MPI_BCAST(knudgestart, 1,0,comm3d,mpierr)
    call D_MPI_BCAST(ntnudge3D  , 1,0,comm3d,mpierr)

    allocate(thlnudge3D(itot,jtot,130,ntnudge3D),qtnudge3D(itot,jtot,130,ntnudge3D)) !MS: allocate 3D nudging arrays (with fluctuations) hardcoded to only have values until k=130, because of file size.
    allocate(tnudgearr(ntnudge3D)) !MS: allocate tnudgearr, which contains nudging time points

    if (.not. lnudge) return
    if(myid==0) then
      t = 0
      !read 1D nudge profiles.
      open (ifinput,file='nudge.inp.'//cexpnr)

      do while (timenudge(t) < runtime)
        t = t + 1
        if (t > ntnudge) then
           write (*,*) "Too many time points in file ", 'nudge.inp.'//cexpnr, ", the limit is ntnudge = ", ntnudge 
           stop
        end if


        chmess1 = "#"
        ierr = 1 ! not zero
        !search for the next line consisting of "# time", from there onwards the profiles will be read
        do while (.not.(chmess1 == "#" .and. ierr ==0))
          read(ifinput,*,iostat=ierr) chmess1, timenudge(t)
          if (ierr < 0) then
            stop 'STOP: No time dependent nudging data for end of run'
          end if

        end do
        write(6,*) 'time', timenudge(t)
        write(6,*) ' height    t_nudge    u_nudge    v_nudge    w_nudge    thl_nudge    qt_nudge'
        do  k=1,kmax
          read (ifinput,*) &
                height    (k), & 
                tnudge    (k,t), & 
                unudge1D  (k,t), & 
                vnudge1D  (k,t), & 
                wnudge1D  (k,t), & 
                thlnudge1D(k,t), &
                qtnudge1D (k,t)
        end do

       ! do k=kmax,1,-1
       !   write (6,'(f7.1,6e12.4)') &
       !         height    (k), &
       !         tnudge    (k,t), &
       !         unudge1D  (k,t), &
       !         vnudge1D  (k,t), &
       !         wnudge1D  (k,t), &
       !         thlnudge1D(k,t), &
       !         qtnudge1D (k,t)
       ! end do
      end do
      close(ifinput)
      !MS: compute nudging time scales.
      t1Dnudge = t1Dnudgefac*tnudge
      t3Dnudge = t3Dnudgefac*tnudge

      !MS: read 3D nudge fields.
      inquire(file='nudge3D.inp.'//cexpnr, exist=file_exists) !MS: check if 3D nudge input file exists

      if (.not. file_exists) then !MS: if file does not exist, print this
        write(6,*) 'The input file for the fields to be nudged ', 'nudge3D.inp.'//cexpnr, ' was not found.'
        STOP
      end if
      
      if (file_exists) then !MS: if the file does exist, read its contents
        open(ifinput,file='nudge3D.inp.'//cexpnr,form='unformatted')
        write(*,*) kmax, jtot, itot

        read(ifinput) tnudgearr !MS: read nudging time points
        !read(ifinput, iostat=ierr) unudge3D
        !read(ifinput) vnudge3D
        !read(ifinput) wnudge3D
        read(ifinput) thlnudge3D !MS: read fluctuations in thl
        read(ifinput) qtnudge3D !MS: read fluctuations in qt

        close(ifinput)
        write(*,*) thlnudge3D(1,1,1,1) !MS: print value for checks
      end if

    end if
    !MS: broadcast nudging arrays
    call D_MPI_BCAST(timenudge ,ntnudge+1   ,0,comm3d,mpierr)
    call D_MPI_BCAST(tnudge    ,k1*ntnudge  ,0,comm3d,mpierr)
    call D_MPI_BCAST(t1Dnudge  ,k1*ntnudge  ,0,comm3d,mpierr)
    call D_MPI_BCAST(t3Dnudge  ,k1*ntnudge  ,0,comm3d,mpierr)
    call D_MPI_BCAST(unudge1D  ,k1*ntnudge  ,0,comm3d,mpierr)
    call D_MPI_BCAST(vnudge1D  ,k1*ntnudge  ,0,comm3d,mpierr)
    call D_MPI_BCAST(wnudge1D  ,k1*ntnudge  ,0,comm3d,mpierr)
    call D_MPI_BCAST(thlnudge1D,k1*ntnudge  ,0,comm3d,mpierr)
    call D_MPI_BCAST(qtnudge1D ,k1*ntnudge  ,0,comm3d,mpierr)
    !call D_MPI_BCAST(unudge3D  ,itot*jtot*kmax,0,comm3d,mpierr)
    !call D_MPI_BCAST(vnudge3D  ,itot*jtot*kmax,0,comm3d,mpierr)
    !call D_MPI_BCAST(wnudge3D  ,itot*jtot*kmax,0,comm3d,mpierr)
    call D_MPI_BCAST(tnudgearr ,ntnudge3D     ,0,comm3d,mpierr)
    call D_MPI_BCAST(thlnudge3D,itot*jtot*130*ntnudge3D,0,comm3d,mpierr)
    call D_MPI_BCAST(qtnudge3D ,itot*jtot*130*ntnudge3D,0,comm3d,mpierr)
    !MS: check nudging booleans: if arrays have any value above 1e-8, nudging is activated, otherwise, it is deactivated for this variable
    lunudge1D   = any(abs(unudge1D)>1e-8)
    lvnudge1D   = any(abs(vnudge1D)>1e-8)
    lwnudge1D   = any(abs(wnudge1D)>1e-8)
    lthlnudge1D = any(abs(thlnudge1D)>1e-8)
    lqtnudge1D  = any(abs(qtnudge1D)>1e-8)
    !lunudge3D   = any(abs(unudge3D)>1e-8)
    !lvnudge3D   = any(abs(vnudge3D)>1e-8)
    !lwnudge3D   = any(abs(wnudge3D)>1e-8)
    lthlnudge3D = any(abs(thlnudge3D)>1e-8)
    lqtnudge3D  = any(abs(qtnudge3D)>1e-8)
    !MS: broadcast booleans
    !call D_MPI_BCAST(lunudge3D   ,1,0,comm3d,mpierr)
    !call D_MPI_BCAST(lvnudge3D   ,1,0,comm3d,mpierr)
    !call D_MPI_BCAST(lwnudge3D   ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lthlnudge3D ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lqtnudge3D  ,1,0,comm3d,mpierr)

  end subroutine initnudge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!MS: completely rewritten from standard DALES
  subroutine nudge
    use modmpi,     only : myidx,myidy,myid
    use modglobal,  only : timee,rtimee,i1,j1,kmax,rdt,itot,jtot,pi,imax,jmax
    use modfields,  only : up,vp,wp,thlp,qtp,u0av,v0av,qt0av,thl0av,thl0,qt0,u0,v0,w0
    use modstartup, only : randthl
    implicit none

    integer k,t,i,j,ti,tn
    integer is,ie,js,je
    real :: dtm,dtp,currtnudge,dtmn,dtpn

    if (.not.(lnudge)) return !MS: check if lnudge = true
!     if (rk3step/=3) return
    if (timee==0) return !MS: nudging not applied if the simulation time = 0

    !MS: determine the nudging time point index
    do ti = 1,ntnudge3D-1
      if (rtimee<tnudgearr(1)) then
        tn = 1
      end if
      if (rtimee>=tnudgearr(ti) .and. rtimee<tnudgearr(ti+1)) then
        tn = ti + 1
      end if
    end do

    !MS: determine fractions for use in extrapolation of fields
    if (tn>1) then
      dtmn = ( rtimee-tnudgearr(tn-1) ) / ( tnudgearr(tn)-tnudgearr(tn-1) )
      dtpn = ( tnudgearr(tn)-rtimee)/ ( tnudgearr(tn)-tnudgearr(tn-1) )
    end if

    !MS: standard DALES nudge time point indexing, not used
    t=1
    do while(rtimee>timenudge(t))
      t=t+1
    end do
    if (rtimee>timenudge(1)) then
      t=t-1
    end if

    dtm = ( rtimee-timenudge(t) ) / ( timenudge(t+1)-timenudge(t) )
    dtp = ( timenudge(t+1)-rtimee)/ ( timenudge(t+1)-timenudge(t) )

    !MS: compute the first and last i and j index for each processor.
    is = myidx * imax + 1
    ie = is + imax - 1

    js = myidy * jmax + 1
    je = js + jmax - 1

    !MS: loop over k values
    do k=1,kmax
      currtnudge   = max(rdt,t1Dnudge(k,t)*dtp+t1Dnudge(k,t+1)*dtm)
      if (rtimee >= tnudgestart .and. rtimee < tnudgestop) then !MS: only apply nudging between the specified times.
        if (k >= knudgestart .and. k < knudgestop) then !MS: only apply 3D nudging between specified heights.
          do j=1,jtot !MS: loop over j,i values
            do i=1,itot
              if (i >= is .and. i <= ie .and. &
                j >= js .and. j <= je) then !MS: check if i,j in the range of used processor
                !MS: 3D field nudge.
                !if(lunudge3D) up (i-is+2,j-js+2,k)=up (i-is+2,j-js+2,k)-(u0 (i-is+2,j-js+2,k)-&
                !    (u0av (k)+unudge3D(i,j,k)))/t3Dnudge (k,t)
                !if(lvnudge3D) vp (i-is+2,j-js+2,k)=vp (i-is+2,j-js+2,k)-(v0 (i-is+2,j-js+2,k)-&
                !    (v0av (k)+vnudge3D(i,j,k)))/t3Dnudge (k,t)
                !if(lwnudge3D) wp (i-is+2,j-js+2,k)=wp (i-is+2,j-js+2,k)-(w0 (i-is+2,j-js+2,k)-&
                !    (         wnudge3D(i,j,k)))/ti3Dnudge (k,t)
                if (ntnudge3D==1) then !MS: no extrapolation for single nudging time point
                  if(lthlnudge3D) thlp (i-is+2,j-js+2,k)=thlp (i-is+2,j-js+2,k)-(thl0 (i-is+2,j-js+2,k)-&
                      (thl0av (k)+thlnudge3D(i,j,k,tn)))/t3Dnudge (k,tn)
                  if(lqtnudge3D) qtp (i-is+2,j-js+2,k)=qtp (i-is+2,j-js+2,k)-(qt0 (i-is+2,j-js+2,k)-&
                      (qt0av (k)+qtnudge3D(i,j,k,tn)))/t3Dnudge (k,tn)
                end if
                if (ntnudge3D>1) then !MS: multiple nudging time points
                  if (tn==1) then !MS: no extrapolation for first nudging time point
                    if(lthlnudge3D) thlp (i-is+2,j-js+2,k)=thlp (i-is+2,j-js+2,k)-(thl0 (i-is+2,j-js+2,k)-&
                        (thl0av (k)+thlnudge3D(i,j,k,tn)))/t3Dnudge (k,tn)
                    if(lqtnudge3D) qtp (i-is+2,j-js+2,k)=qtp (i-is+2,j-js+2,k)-(qt0 (i-is+2,j-js+2,k)-&
                        (qt0av (k)+qtnudge3D(i,j,k,tn)))/t3Dnudge (k,tn)
                  end if
                  if (tn>1) then !MS: Extrapolation 3D nudging
                    if(lthlnudge3D) thlp (i-is+2,j-js+2,k)=thlp (i-is+2,j-js+2,k)-(thl0 (i-is+2,j-js+2,k)-&
                        (thl0av (k)+(thlnudge3D(i,j,k,tn-1)*dtpn + thlnudge3D(i,j,k,tn)*dtmn)))/currtnudge
                    if(lqtnudge3D) qtp (i-is+2,j-js+2,k)=qtp (i-is+2,j-js+2,k)-(qt0 (i-is+2,j-js+2,k)-&
                        (qt0av (k)+(qtnudge3D(i,j,k,tn-1)*dtpn + qtnudge3D(i,j,k,tn)*dtmn)))/currtnudge
                  end if
                end if
              end if
            end do
          end do
        end if
        !MS: 1D profile nudge.
        if (ntnudge3D==1) then !MS: no extrapolation for single nudging time point
          if(lthlnudge1D) thlp (2:i1,2:j1,k)=thlp (2:i1,2:j1,k)-&
                (thl0av (k)-(thlnudge1D (k,tn)))/t1Dnudge(k,tn)
          if(lqtnudge1D) qtp (2:i1,2:j1,k)=qtp (2:i1,2:j1,k)-&
                (qt0av (k)-(qtnudge1D (k,tn)))/t1Dnudge(k,tn)
        end if
        if (ntnudge3D>1) then !MS: multiple nudging time points
          if (tn == 1) then !MS: no extrapolation for first nudging time point
            if(lthlnudge1D) thlp (2:i1,2:j1,k)=thlp (2:i1,2:j1,k)-&
                  (thl0av (k)-(thlnudge1D (k,tn)))/t1Dnudge(k,tn)
            if(lqtnudge1D) qtp (2:i1,2:j1,k)=qtp (2:i1,2:j1,k)-&
                  (qt0av (k)-(qtnudge1D (k,tn)))/t1Dnudge(k,tn)
          end if
          if (tn>1) then !MS: Extrapolation 1D nudging
            if(lthlnudge1D) thlp (2:i1,2:j1,k)=thlp (2:i1,2:j1,k)-&
                  (thl0av (k)-(thlnudge1D (k,tn-1)*dtpn+thlnudge1D (k,tn)*dtmn))/currtnudge
            if(lqtnudge1D) qtp (2:i1,2:j1,k)=qtp (2:i1,2:j1,k)-&
                  (qt0av (k)-(qtnudge1D (k,tn-1)*dtpn+qtnudge1D (k,tn)*dtmn))/currtnudge
          end if
        end if
        !MS standard DALES 1D nudging for u,v,w
        if(lunudge1D  ) up  (2:i1,2:j1,k)=up  (2:i1,2:j1,k)-&
              (u0av  (k)-(unudge1D  (k,t)*dtp+unudge1D  (k,t+1)*dtm))/currtnudge
        if(lvnudge1D  ) vp  (2:i1,2:j1,k)=vp  (2:i1,2:j1,k)-&
              (v0av  (k)-(vnudge1D  (k,t)*dtp+vnudge1D  (k,t+1)*dtm))/currtnudge
        if(lwnudge1D  ) wp  (2:i1,2:j1,k)=wp  (2:i1,2:j1,k)-&
              (         -(wnudge1D  (k,t)*dtp+wnudge1D  (k,t+1)*dtm))/currtnudge
      end if
    end do
  end subroutine nudge

  subroutine exitnudge
  deallocate(timenudge)
  end subroutine exitnudge

end module
