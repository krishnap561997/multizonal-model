!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 rpi,rp,vmag,rep,rmass,dmass,fbx,fby,fbz,fqs,fqsx,fqsy,fqsz
      integer*4 i

      integer*4 icalld
      save      icalld
      data      icalld /0/

!
! External:
!
      real*8 rmu,rhof,rg
      common /parameters/ rmu,rhof,rg

      real*8 rlx,rrx,rly,rry,rlz,rrz
      common /domainsize/ rlx,rrx,rly,rry,rlz,rrz

      real*8 dpmin
      common /dpm/ dpmin

      real*8 y_out(PPICLF_LRS         , PPICLF_LPART) ! Normal ordering
      real*8 rprop_out(PPICLF_LRP + 1 , PPICLF_LPART) ! Normal ordering
      integer*4 iprop_out(5  , PPICLF_LPART) ! Normal ordering
      integer*4 npart_out
      common /outlet/ y_out,rprop_out,iprop_out,npart_out

! Langevin
      real*8 lan_c0,lan_cy,lan_cg,lan_gij,lan_cw
      real*8 lan_cw1,lan_cw2,lan_cw3  
      common /langevin/ lan_c0,lan_cy,lan_cg,lan_gij,lan_cw
      real*8 ppiclf_ran,fun_erfinv
!
      real*8 dum_ufl(3)

      rpi  = 4.0*atan(1.0)

!     REMOVE PARTICLES TOUCHING WALLS OR OUTLET
      do i=1,ppiclf_npart
         rp = ppiclf_rprop(PPICLF_R_JDP,i)/2.0
         if(ppiclf_y(PPICLF_JX,i).gt. rrx-rp .or.
     >      ppiclf_y(PPICLF_JX,i).lt. rlx+rp .or.
     >      ppiclf_y(PPICLF_JY,i).gt. rry-rp .or.
     >      ppiclf_y(PPICLF_JY,i).lt. rly+rp .or.
     >      ppiclf_y(PPICLF_JZ,i).gt. rrz-rp .or.
     >      ppiclf_y(PPICLF_JZ,i).lt. rlz+rp     )then
            call ppiclf_solve_MarkForRemoval(i)
            call ppiclf_solve_SaveRemoved(i)
         endif
      enddo

      if (icalld.eq.3) icalld = 0
      icalld = icalld + 1

! evaluate ydot
      do i=1,ppiclf_npart

c        ##########################################################
c        Langevin drifting velocity
         if (icalld.eq.1) then       
          if (ppiclf_cycle.eq.0) then  
            ppiclf_rprop(PPICLF_R_JDUX,i)=0.d0
            ppiclf_rprop(PPICLF_R_JDUY,i)=0.d0
            ppiclf_rprop(PPICLF_R_JDUZ,i)=0.d0
          else
            lan_gij=1.d0-lan_cg*ppiclf_dt*ppiclf_rprop(PPICLF_R_JCS2S,i)
            lan_gij=min(lan_gij,1.d0)
            lan_gij=max(lan_gij,0.d0)
c            lan_cw=sqrt(lan_c0*ppiclf_dt
c     >              *max(ppiclf_rprop(PPICLF_R_JDISS,i),0.d0))
            lan_cw1=sqrt(lan_c0*ppiclf_dt)
     >              *max(ppiclf_rprop(PPICLF_R_JFX,i),0.d0)
            lan_cw2=sqrt(lan_c0*ppiclf_dt)
     >              *max(ppiclf_rprop(PPICLF_R_JFY,i),0.d0)
            lan_cw3=sqrt(lan_c0*ppiclf_dt)
     >              *max(ppiclf_rprop(PPICLF_R_JFZ,i),0.d0)


            ppiclf_rprop(PPICLF_R_JDUX,i)=
     >        lan_cw1*sqrt(2.0/3.0)*fun_erfinv(2.d0*ppiclf_ran(2)-1.d0)
     >       +lan_gij*ppiclf_rprop(PPICLF_R_JDUX,i)
            ppiclf_rprop(PPICLF_R_JDUY,i)=
     >        lan_cw2*sqrt(2.0/3.0)*fun_erfinv(2.d0*ppiclf_ran(2)-1.d0)
     >       +lan_gij*ppiclf_rprop(PPICLF_R_JDUY,i)
            ppiclf_rprop(PPICLF_R_JDUZ,i)=
     >        lan_cw3*sqrt(2.0/3.0)*fun_erfinv(2.d0*ppiclf_ran(2)-1.d0)
     >       +lan_gij*ppiclf_rprop(PPICLF_R_JDUZ,i)
          endif

c         Add perturbations from Langevin model to fluid velocity
          ppiclf_rprop(PPICLF_R_JUTX,i) = ppiclf_rprop(PPICLF_R_JUX,i)
     >                                  + ppiclf_rprop(PPICLF_R_JDUX,i)
          ppiclf_rprop(PPICLF_R_JUTY,i) = ppiclf_rprop(PPICLF_R_JUY,i)
     >                                  + ppiclf_rprop(PPICLF_R_JDUY,i)
          ppiclf_rprop(PPICLF_R_JUTZ,i) = ppiclf_rprop(PPICLF_R_JUZ,i)
     >                                  + ppiclf_rprop(PPICLF_R_JDUZ,i)
         endif

         if (ppiclf_rprop(PPICLF_R_JDP,i).le.dpmin) then ! Small particles
           ppiclf_ydot(PPICLF_JX ,i) = ppiclf_rprop(PPICLF_R_JUTX,i)
           ppiclf_ydot(PPICLF_JY ,i) = ppiclf_rprop(PPICLF_R_JUTY,i)
     >                               + ppiclf_rprop(PPICLF_R_JVS,i)
           ppiclf_ydot(PPICLF_JZ ,i) = ppiclf_rprop(PPICLF_R_JUTZ,i)
  
           ppiclf_ydot(PPICLF_JVX,i) = 0.0
           ppiclf_ydot(PPICLF_JVY,i) = 0.0
           ppiclf_ydot(PPICLF_JVZ,i) = 0.0
         else
c        Particle motion with force models
         vmag  = sqrt((ppiclf_rprop(PPICLF_R_JUTX,i)
     >                -ppiclf_y(PPICLF_JVX,i))**2
     >               +(ppiclf_rprop(PPICLF_R_JUTY,i)
     >                -ppiclf_y(PPICLF_JVY,i))**2
     >               +(ppiclf_rprop(PPICLF_R_JUTZ,i)
     >                -ppiclf_y(PPICLF_JVZ,i))**2)
         rep   = rhof*vmag*ppiclf_rprop(PPICLF_R_JDP,i)/rmu

         ! set ydot for all PPICLF_SLN number of equations
         ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
         ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
         ppiclf_ydot(PPICLF_JZ ,i) = ppiclf_y(PPICLF_JVZ,i)

         ! gravity
         rmass = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *ppiclf_rprop(PPICLF_R_JRHOP,i)
         dmass = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *(ppiclf_rprop(PPICLF_R_JRHOP,i)-rhof)
         fbx  = 0.d0
         fby  = dmass*rg
         fbz  = 0.d0

         ! quasi-steady
         fqs   = 3.d0*rpi*rmu*ppiclf_rprop(PPICLF_R_JDP,i)
     >          *(1.d0+0.15d0*rep**0.687d0)
         fqsx  = fqs
     >   *(ppiclf_rprop(PPICLF_R_JUTX,i)-ppiclf_y(PPICLF_JVX,i))
         fqsy  = fqs
     >   *(ppiclf_rprop(PPICLF_R_JUTY,i)-ppiclf_y(PPICLF_JVY,i))
         fqsz  = fqs
     >   *(ppiclf_rprop(PPICLF_R_JUTZ,i)-ppiclf_y(PPICLF_JVZ,i))

         ppiclf_ydot(PPICLF_JVX,i) = (fbx+fqsx)/rmass
         ppiclf_ydot(PPICLF_JVY,i) = (fby+fqsy)/rmass
         ppiclf_ydot(PPICLF_JVZ,i) = (fbz+fqsz)/rmass

         endif

      enddo 
! evaluate ydot

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop)
!
      implicit none
!
! Input:
!
      real*8 y    (PPICLF_LRS)
      real*8 ydot (PPICLF_LRS)
      real*8 ydotc(PPICLF_LRS)
      real*8 rprop(PPICLF_LRP)
!
! Output:
!
      real*8 map  (PPICLF_LRP_PRO)
!
! Internal:
!
      real*8 dp_norm
!

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_EvalNearestNeighbor
     >                                        (i,j,yi,rpropi,yj,rpropj)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 i
      integer*4 j
      real*8 yi    (PPICLF_LRS)
      real*8 rpropi(PPICLF_LRP)
      real*8 yj    (PPICLF_LRS)
      real*8 rpropj(PPICLF_LRP)
!
! Internal:
!
      real*8 ksp,erest
      common /ucollision/ ksp,erest
#ifdef PPICLC
      BIND(C, name="ucollision") :: /ucollision/ ! c binding
#endif


      return
      end

!-----------------------------------------------------------------------
      subroutine ppiclf_solve_SaveRemoved(i)
!
      implicit none
!
      include "PPICLF"
!
      real*8 rlx,rrx,rly,rry,rlz,rrz
      common /domainsize/ rlx,rrx,rly,rry,rlz,rrz

      real*8 rlx_out,rrx_out,rly_out,rry_out,rlz_out,rrz_out
      common /outletsize/rlx_out,rrx_out,rly_out,rry_out,rlz_out,rrz_out

      real*8 y_out(PPICLF_LRS         , PPICLF_LPART) ! Normal ordering
      real*8 rprop_out(PPICLF_LRP + 1 , PPICLF_LPART) ! Normal ordering
      integer*4 iprop_out(5  , PPICLF_LPART) ! Normal ordering
      integer*4 npart_out
      common /outlet/ y_out,rprop_out,iprop_out,npart_out



! Input:
!
      integer*4 i,j
      real*8 rp
!
      npart_out = npart_out + 1

!     SAVE PARTICLE POSITION WHEN OUT OF DOMIAN
      do j=1,PPICLF_LRS
        y_out(j,npart_out) = ppiclf_y(j,i)
      enddo

      do j=1,PPICLF_LRP
        rprop_out(j,npart_out) = ppiclf_rprop(j,i)
      enddo

      rprop_out(PPICLF_LRP+1,npart_out) = ppiclf_time
      
      do j=1,3
        iprop_out(j,npart_out) = ppiclf_iprop(j+4,i)
      enddo
      
      iprop_out(4,npart_out) = ppiclf_cycle
!     
      rp = ppiclf_rprop(PPICLF_R_JDP,i)/2.0

      if (ppiclf_y(PPICLF_JX,i).gt. rrx-rp) iprop_out(5,npart_out) = 2
      if (ppiclf_y(PPICLF_JX,i).lt. rlx+rp) iprop_out(5,npart_out) = 4

      if (ppiclf_y(PPICLF_JY,i).gt. rry-rp) then
        if (ppiclf_y(PPICLF_JX,i).ge.rlx_out .and. 
     >      ppiclf_y(PPICLF_JX,i).le.rrx_out .and.
     >      ppiclf_y(PPICLF_JZ,i).ge.rlz_out .and.
     >      ppiclf_y(PPICLF_JZ,i).le.rrz_out      ) then
          iprop_out(5,npart_out) = 7
        else
          iprop_out(5,npart_out) = 3
        endif
      endif

      if (ppiclf_y(PPICLF_JY,i).lt. rly+rp) iprop_out(5,npart_out) = 1
      if (ppiclf_y(PPICLF_JZ,i).gt. rrz-rp) iprop_out(5,npart_out) = 6
      if (ppiclf_y(PPICLF_JZ,i).lt. rlz+rp) iprop_out(5,npart_out) = 5

      return
      end

!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteParticleOutletVTU(filein1)
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      character (len = *) filein1
!
! Internal:
!
      real*4  rout_pos(3      *PPICLF_LPART) 
     >       ,rout_sln(PPICLF_LRS*PPICLF_LPART)
     >       ,rout_lrp((PPICLF_LRP+1)*PPICLF_LPART)
     >       ,rout_lip(5      *PPICLF_LPART)
      character*6 filein
      character*15 vtufile
      character*6  prostr
      integer*4 icalld1
      save      icalld1
      data      icalld1 /0/
      integer*4 vtu,pth,prevs(2,ppiclf_np)
      integer*8 idisp_pos,idisp_sln,idisp_lrp,idisp_lip,stride_len
      integer*4 iint, nnp, nxx, npt_total, jx, jy, jz, if_sz, isize,
     >          iadd, if_pos, if_sln, if_lrp, if_lip, ic_pos, ic_sln,
     >          ic_lrp, ic_lip, i, j, ie, nps, nglob, nkey, ndum,
     >          icount_pos, icount_sln, icount_lrp, icount_lip, iorank,
     >          ierr, ivtu_size
      integer*4 ppiclf_iglsum
      external ppiclf_iglsum
!
! External:
!
      real*8 y_out(PPICLF_LRS           , PPICLF_LPART) ! Normal ordering
      real*8 rprop_out(PPICLF_LRP+1     , PPICLF_LPART) ! Normal ordering
      integer*4 iprop_out(5  , PPICLF_LPART) ! Normal ordering
      integer*4 npart_out
      common /outlet/ y_out,rprop_out,iprop_out,npart_out
!
!
!
      call ppiclf_printsi(' *Begin WriteParticleOutletVTU$'
     >                   ,ppiclf_cycle)

      icalld1 = icalld1+1

      nnp   = ppiclf_np
      nxx   = npart_out

      npt_total = ppiclf_iglsum(nxx,1)

      jx    = 1
      jy    = 2
      jz    = 1
      if (ppiclf_ndim .eq. 3)
     >jz    = 3

c      if_sz = len(filein1)
c      if (if_sz .lt. 3) then
c         filein = 'par'
c      else 
         write(filein,'(A6)') filein1
c      endif

! --------------------------------------------------
! COPY PARTICLES TO OUTPUT ARRAY
! --------------------------------------------------

      isize = 4

      iadd = 0
      if_pos = 3*isize*npt_total
      if_sln = 1*isize*npt_total
      if_lrp = 1*isize*npt_total
      if_lip = 1*isize*npt_total

      ic_pos = iadd
      ic_sln = iadd
      ic_lrp = iadd
      ic_lip = iadd
 
      if (nxx.ne.0) then
      do i=1,nxx
         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = sngl(y_out(jx,i))
         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = sngl(y_out(jy,i))
         ic_pos = ic_pos + 1
         if (ppiclf_ndim .eq. 3) then
            rout_pos(ic_pos) = sngl(y_out(jz,i))
         else
            rout_pos(ic_pos) = 0.0
         endif
      enddo
      do j=1,PPICLF_LRS
      do i=1,nxx
         ic_sln = ic_sln + 1
         rout_sln(ic_sln) = sngl(y_out(j,i))
      enddo
      enddo
      do j=1,PPICLF_LRP+1
      do i=1,nxx
         ic_lrp = ic_lrp + 1
         rout_lrp(ic_lrp) = sngl(rprop_out(j,i))
      enddo
      enddo
      do j=1,5
      do i=1,nxx
         ic_lip = ic_lip + 1
         rout_lip(ic_lip) = real(iprop_out(j,i))
      enddo
      enddo
      endif ! nxx.ne.0

! --------------------------------------------------
! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
! --------------------------------------------------
      do i=1,nnp
         prevs(1,i) = i-1
         prevs(2,i) = nxx
      enddo

      nps   = 1 ! index of new proc for doing stuff
      nglob = 1 ! unique key to sort by
      nkey  = 1 ! number of keys (just 1 here)
      ndum = 2
      call pfgslib_crystal_ituple_transfer(ppiclf_cr_hndl,prevs,
     >                 ndum,nnp,nnp,nps)
      call pfgslib_crystal_ituple_sort(ppiclf_cr_hndl,prevs,
     >                 ndum,nnp,nglob,nkey)

      stride_len = 0
      if (ppiclf_nid .ne. 0) then
      do i=1,ppiclf_nid
         stride_len = stride_len + prevs(2,i)
      enddo
      endif

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A6,I5.5,A4)') filein,icalld1,'.vtu'

      if (ppiclf_nid .eq. 0) then !---------------------------

      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (ppiclf_iendian .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (ppiclf_iendian .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') ppiclf_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') ppiclf_cycle
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='yes') '  </FieldData>'
      write(vtu,'(A)',advance='no') '  <Piece '
      write(vtu,'(A)',advance='no') 'NumberOfPoints="'
      write(vtu,'(I0)',advance='no') npt_total
      write(vtu,'(A)',advance='yes') '" NumberOfCells="0"> '

! -----------
! COORDINATES 
! -----------
      iint = 0
      write(vtu,'(A)',advance='yes') '   <Points>'
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3,iint)
      iint = iint + 3*isize*npt_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'


      do ie=1,PPICLF_LRS
         write(prostr,'(A1,I2.2)') "y",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      do ie=1,PPICLF_LRP+1
         write(prostr,'(A4,I2.2)') "rprop",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      do ie=1,5
         write(prostr,'(A3,I2.2)') "tag",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      write(vtu,'(A)',advance='yes') '   </PointData> '

! ----------
! END MATTER
! ----------
      write(vtu,'(A)',advance='yes') '   <Cells> '
      write(vtu,'(A)',advance='no')  '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="connectivity" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="types" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='yes') '   </Cells> '
      write(vtu,'(A)',advance='yes') '  </Piece> '
      write(vtu,'(A)',advance='yes') ' </UnstructuredGrid> '

! -----------
! APPEND DATA  
! -----------
      write(vtu,'(A)',advance='no') ' <AppendedData encoding="raw">'
      close(vtu)

      open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >    ,position='append')
      write(vtu) '_'
      close(vtu)

      inquire(file=vtufile,size=ivtu_size)
      endif ! ------------ nid .eq. 0 -------------------------

      call ppiclf_bcast(ivtu_size, isize)

      ! byte-displacements
      idisp_pos = ivtu_size + isize*(3*stride_len + 1)

      ! how much to write
      icount_pos = 3*nxx
      icount_sln = 1*nxx
      icount_lrp = 1*nxx
      icount_lip = 1*nxx

      iorank = -1

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      ! write
        call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
        call ppiclf_byte_set_view(idisp_pos,pth)
        call ppiclf_byte_write_mpi(rout_pos,icount_pos,iorank,pth,ierr)
        call ppiclf_byte_close_mpi(pth,ierr)

      call mpi_barrier(ppiclf_comm,ierr)

      do i=1,PPICLF_LRS
         idisp_sln = ivtu_size + isize*(3*npt_total 
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + i)

         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_sln
           close(vtu)
         endif
   
         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*npart_out + 1
   
         ! write
c         if (nxx.ne.0) then
           call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
           call ppiclf_byte_set_view(idisp_sln,pth)
           call ppiclf_byte_write_mpi(rout_sln(j),icount_sln,iorank,pth
     >                             ,ierr)
           call ppiclf_byte_close_mpi(pth,ierr)
c        endif ! nxx.ne.0
      enddo

      do i=1,PPICLF_LRP+1
         idisp_lrp = ivtu_size + isize*(3*npt_total  
     >                         + PPICLF_LRS*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + i)

         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_lrp
           close(vtu)
         endif
   
         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*npart_out + 1
   
         ! write
c        if (nxx.ne.0) then
           call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
           call ppiclf_byte_set_view(idisp_lrp,pth)
           call ppiclf_byte_write_mpi(rout_lrp(j),icount_lrp,iorank,pth
     >                             ,ierr)
           call ppiclf_byte_close_mpi(pth,ierr)
c         endif
      enddo

      do i=1,5
         idisp_lip = ivtu_size + isize*(3*npt_total
     >                         + PPICLF_LRS*npt_total
     >                         + (PPICLF_LRP+1)*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + PPICLF_LRP + 1 + i)
         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_lip
           close(vtu)
         endif

         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*npart_out + 1
   
         ! write
c         if (nxx.ne.0) then
           call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
           call ppiclf_byte_set_view(idisp_lip,pth)
           call ppiclf_byte_write_mpi(rout_lip(j),icount_lip,iorank,pth
     >                             ,ierr)
           call ppiclf_byte_close_mpi(pth,ierr)
c         endif
      enddo

      if (ppiclf_nid .eq. 0) then
      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      call ppiclf_printsi(' *End WriteParticleOutletVTU$',ppiclf_cycle)

      return
      end

!-----------------------------------------------------------------------
      REAL*8 FUNCTION ppiclf_ran(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     $        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     $        IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
c Long period (> 2 ! 1018 ) random number generator of L’Ecuyer with
c Bays-Durham shuffle and added safeguards. Returns a uniform random deviate
c between 0.0 and 1.0 (exclusive of the endpoint values).
c Call with idum a negative integer to initialize; thereafter, do not alter
c idum between successive deviates in a sequence. RNMX should approximate the
c largest floating value that is less than 1.
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum1=max(-idum,1)
         idum2=idum1
         do j=NTAB+8,1,-1
            k=idum1/IQ1
            idum1=IA1*(idum1-k*IQ1)-k*IR1
            if (idum1.lt.0) idum1=idum1+IM1
            if (j.le.NTAB) iv(j)=idum1
         enddo
         iy=iv(1)
      endif
      k=idum1/IQ1
      idum1=IA1*(idum1-k*IQ1)-k*IR1
      if (idum1.lt.0) idum1=idum1+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum1
      if(iy.lt.1)iy=iy+IMM1
      ppiclf_ran=min(AM*iy,RNMX)
      return
      END

!-----------------------------------------------------------------------
      real*8 function fun_erfinv(x)

      implicit none
      real*8 x,a
      real*8 :: pi=4.*atan(1.d0)

      if (x.eq.0.d0)then
         fun_erfinv=0.d0
      else
         a=8.0*(pi-3.0)/(3.0*pi*(4.0-pi))
         fun_erfinv=x/abs(x)*sqrt(sqrt((2.0/pi/a+0.5*log(1.0-x**2))**2
     $             -log(1.0-x**2)/a)-(2.0/pi/a+0.5*log(1.0-x**2)))
      endif

      return
      end
