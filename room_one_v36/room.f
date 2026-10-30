#include "PPICLF_USER.h"
#include "PPICLF_STD.h"

#define IOPART uparam(1)                         /* IOSTEP FOR PARTICLES */
#define srcin_len   uparam( 2)
#define npart_inj  uparam( 3)
#define usr_endinj  uparam( 4)
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
c
c     ############################################################
c     NEK5000
      real*8 buffthick,buffmag

      common /cdsmag/ ediff(lx1,ly1,lz1,lelv)
      integer ie,ix,iy,iz,eg
      
c     ############################################################
c     PPICLF
      real*8 rlx,rrx,rly,rry,rlz,rrz,dum1,dum2
      common /domainsize/ rlx,rrx,rly,rry,rlz,rrz

      real*8 rlx_out,rrx_out,rly_out,rry_out,rlz_out,rrz_out
      common /outletsize/rlx_out,rrx_out,rly_out,rry_out,rlz_out,rrz_out

c     ##################################################################
c     SMAGORINSKY

      
      ie  = gllel(eg)

      buffmag    = 20.0

      buffthick  = rry_out - rly_out  !size (in y) of top element in domain

      udiff = ediff(ix,iy,iz,ie)
      utrans= param(1)

      if (abs(rry_out-y).le.buffthick) then      
      if (x.ge.rlx_out.and.x.le.rrx_out) then
      if (z.ge.rlz_out.and.z.le.rrz_out) then
        udiff =udiff 
     >      +((ediff(ix,iy,iz,ie)*buffmag-ediff(ix,iy,iz,ie))/buffthick)
     >      *(y-rry)
      endif
      endif
      endif

      if(ifield.eq.2) then
         utrans = param(7)
         udiff = param(8)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'ZPER'  ! for nelx,nely,nelz

c     ############################################################
c     NEK5000
      real*8 buffthick,buffmag,bufftarget,buff,yy,beta
      integer ix,iy,iz,ie,eg

c     ############################################################
c     PPICLF
      real*8 rlx,rrx,rly,rry,rlz,rrz,dum1,dum2
      common /domainsize/ rlx,rrx,rly,rry,rlz,rrz

      real*8 rlx_out,rrx_out,rly_out,rry_out,rlz_out,rrz_out
      common /outletsize/rlx_out,rrx_out,rly_out,rry_out,rlz_out,rrz_out
c
c
c
c     ############################################################
      buffmag    = 1.d2
      bufftarget = 1.d-3
      beta = 2.0
c     ############################################################

      buffthick  = rry_out - rly_out  !size (in y) of top element in domain

      ie = gllel(eg)

      ffx = 0.0
      ffy = 0.0
      ffz = 0.0

      if (abs(rry_out-y).le.buffthick) then      
      if (x.ge.rlx_out.and.x.le.rrx_out) then
      if (z.ge.rlz_out.and.z.le.rrz_out) then
        if (vy(ix,iy,iz,ie).lt.0.d0) then
          yy = (y-rry)/(rry_out-rry)
          buff = buffmag * yy**beta
          ffy = ffy - buff*(vy(ix,iy,iz,ie)-bufftarget)
     >        / vtrans(ix,iy,iz,ie,1)         
        endif 
      endif
      endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg
c     e = gllel(eg)
 
      qvol   = 0.0
      source = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk
      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'  ! for nelx,nely,nelz

      include "PPICLF"

c     ############################################################
c     PPICLF
      real*8 y_out(PPICLF_LRS, PPICLF_LPART) ! Normal ordering
      real*8 rprop_out(PPICLF_LRP + 1 , PPICLF_LPART) ! Normal ordering
      integer*4 iprop_out(5  , PPICLF_LPART) ! Normal ordering
      integer*4 npart_out
      common /outlet/ y_out,rprop_out,iprop_out,npart_out

      ! Langevin
      real*8  lan_diss(lx1,ly1,lz1,lelt)
      real*8  lan_cs2s(lx1,ly1,lz1,lelt)
      common /lan/ lan_diss,lan_cs2s

      !Particle Injection
      real*8 y(PPICLF_LRS    , PPICLF_LPART) ! Normal ordering
      real*8 rprop(PPICLF_LRP, PPICLF_LPART) ! Normal ordering
      real*8 npart 

      integer i,j
c     ############################################################
c     NEK5000

      integer n,e

c     SMAGORINSKY
      common /cdsmag/ ediff(lx1,ly1,lz1,lelv)

      common /dynsmg/ sij (lx1*ly1*lz1,ldim,ldim)
     $              , mij (lx1*ly1*lz1,3*ldim-3)
     $              , lij (lx1*ly1*lz1,3*ldim-3)
     $              , dg2 (lx1*ly1*lz1,lelv)
     $              , num (lx1*ly1*lz1,lelv)
     $              , den (lx1*ly1*lz1,lelv)
     $              , snrm(lx1*ly1*lz1,lelv)
      real sij,mij,lij,dg2,num,den,snrm

c     TIME AVERAGE
      common /avgtime/ 
     $        uavg(lx1,ly1,lz1,lelv)
     $    ,   vavg(lx1,ly1,lz1,lelv)
     $    ,   wavg(lx1,ly1,lz1,lelv)
     $    ,   tavg(lx1,ly1,lz1,lelt)
     $    ,   ediffavg(lx1,ly1,lz1,lelv)
      integer icalled
      save    icalled
      data    icalled /0/

      real atime,timel,dtime
      save atime,timel
      real alpha,beta
      logical ifverbose

c     LANGEVIN MODEL - KRISHNA 02/26/2022
      real*8 rlx,rrx,rly,rry,rlz,rrz,dum1,dum2
      common /domainsize/ rlx,rrx,rly,rry,rlz,rrz

      integer exy, eyz, ezx, eg, k, lxyz, ex, ey, ez
      integer lelxy, lelyz,lelzx, mxy, myz, mzx
      real*8 rnu,xp,yp,zp,f_par,f_perp,ylim
     >       wrx, wlx, wry, wly, wrz, wlz,tmp1,tmp2 
      real*8 xl(lx1,ly1,lz1,lelv), xp1(lx1,ly1,lz1,lelv),
     >       xl_f(lx1,ly1,lz1,lelv),yl_f(lx1,ly1,lz1,lelv), 
     >       zl_f(lx1,ly1,lz1,lelv), 
     >       yl(lx1,ly1,lz1,lelv), yp1(lx1,ly1,lz1,lelv),
     >       zl(lx1,ly1,lz1,lelv), zp1(lx1,ly1,lz1,lelv),
     >       f1(lx1,ly1,lz1,lelv), f2(lx1,ly1,lz1,lelv),
     >       f3(lx1,ly1,lz1,lelv),
     >       vxx(lx1,ly1,lz1,lelv), vxy(lx1,ly1,lz1,lelv),
     >       vxz(lx1,ly1,lz1,lelv), vyx(lx1,ly1,lz1,lelv),
     >       vyy(lx1,ly1,lz1,lelv), vyz(lx1,ly1,lz1,lelv),
     >       vzx(lx1,ly1,lz1,lelv), vzy(lx1,ly1,lz1,lelv),
     >       vzz(lx1,ly1,lz1,lelv),
     >       wrx_grad(ly1,lz1,lelz,lely), wlx_grad(ly1,lz1,lelz,lely),
c      real*8 wx(ly1,lz1,lelz*lely)
     >       wry_grad(lz1,lx1,lelz,lelx),wly_grad(lz1,lx1,lelz,lelx),
c      real*8 wy(lz1,lx1,lelz*lelx)
     >       wrz_grad(lx1,ly1,lelx,lely),wlz_grad(lx1,ly1,lelx,lely),
c      real*8 wz(lx1,ly1,lelx*lely)
     >       f1_diss(lx1,ly1,lz1,lelv),f2_diss(lx1,ly1,lz1,lelv),
     >       f3_diss(lx1,ly1,lz1,lelv)

c     RANDOM NUMBER PERTURBATIONS - KRISHNA 02/26/2022
      real*8 rnum_arr(lx1,lz1,lelx,lelz)
      real*8 rdum
      common /randomarray/ rnum_arr
      external ran2

c     ############################################################
c     NEK5000

c     SMAGORINSKY
c     Compute eddy viscosity using dynamic smagorinsky model
      if (istep.eq.0) then
         ifuservp  = .true.
         ifexplvis = .true.
         param(30) = 1
      endif

      if(ifuservp) then
        if(nid.eq.0) print*,'Calculating eddy viscosity'
        call set_grid_spacing
     
        do e=1,lelv
           call eddy_visc(ediff,e)
        enddo

        call eddy_visc_01(ediff)
C       call filter_s0mine(ediff,0.95,2,'ediff')


        n=nx1*ny1*nz1*nelv
        call copy(t(1,1,1,1,2),ediff,n)
c        call copy(t(1,1,1,1,3),dg2,n)
      endif

c     OUTPUT INITIAL CONDITION
      if (istep.eq.0) then
c        call gfldr("FRONT.ini") 
        ifxyo = .true.
        call prepost(.true.,'   ' )
        ifxyo = .false.
      endif

c     -----------------------------------
c     LANGEVIN MODEL - KRISHNA 02/26/2022
      !Part 1 - Storing Xl, Yl and Zl divided by nu for the mesh

      lxyz = lx1*ly1*lz1
      rnu = param(2)/param(1)
      if (istep.eq.0) then
        do i=1,n
       
        !IF LEFT WALL IS CLOSER, FLAG SET TO 1
          wrx = abs(rrx - xm1(i,1,1,1))
          wry = abs(rry - ym1(i,1,1,1))
          wrz = abs(rrz - zm1(i,1,1,1))
          wlx = abs(xm1(i,1,1,1)-rlx)
          wly = abs(ym1(i,1,1,1)-rly)
          wlz = abs(zm1(i,1,1,1)-rlz)
  
          if (wrx.ge.wlx) then
            xl_f(i,1,1,1) = 1
            xl(i,1,1,1) = wlx/rnu
          else
            xl_f(i,1,1,1) = 0
            xl(i,1,1,1) = wrx/rnu
          endif
  
          if (wry.ge.wly) then
            yl_f(i,1,1,1) = 1
            yl(i,1,1,1) = wly/rnu
          else
            yl_f(i,1,1,1) = 0
            yl(i,1,1,1) = wry/rnu
          endif
    
          if (wrz.ge.wlz) then
            zl_f(i,1,1,1) = 1
            zl(i,1,1,1) = wlz/rnu
          else
            zl_f(i,1,1,1) = 0
            zl(i,1,1,1) = wrz/rnu
          endif
        enddo
      endif

      !Part 2 - Calculating U* values at different walls
      call gradm1(vxx,vxy,vxz,vx)
      call gradm1(vyx,vyy,vyz,vy)
      call gradm1(vzx,vzy,vzz,vz)

      lelxy = lelx*lely
      lelyz = lely*lelz
      lelzx = lelz*lelx
      mxy = lx1*ly1*lelxy
      myz = ly1*lz1*lelyz
      mzx = lz1*lx1*lelzx

      call rzero(wrx_grad, myz)
      call rzero(wlx_grad, myz)  
      call rzero(wry_grad, mzx)
      call rzero(wly_grad, mzx)
      call rzero(wrz_grad, mxy)
      call rzero(wlz_grad, mxy)

      do e=1,lelv
        do k=1,lz1
          do j=1,ly1
            do i=1,lx1
          
              eg = lglel(e)
              call get_exyz(ex,ey,ez,eg,lelx,lely,lelz)
              if (ex.eq.1.and.i.eq.1) then
                wlx_grad(j,k,ey,ez) = sqrt(rnu*sqrt(vyx(i,j,k,e)**2 
     >                                + vzx(i,j,k,e)**2))
              endif
              if (ex.eq.lelx.and.i.eq.lx1) then
                wrx_grad(j,k,ey,ez) = sqrt(rnu*sqrt(vyx(i,j,k,e)**2
     >                                + vzx(i,j,k,e)**2))
              endif
              if (ey.eq.1.and.j.eq.1) then
                wly_grad(k,i,ez,ex) = sqrt(rnu*sqrt(vzy(i,j,k,e)**2 
     >                                + vxy(i,j,k,e)**2))
              endif
              if (ey.eq.lely.and.j.eq.ly1) then
                wry_grad(k,i,ez,ex) = sqrt(rnu*sqrt(vzy(i,j,k,e)**2 
     >                                + vxy(i,j,k,e)**2))
              endif
              if (ez.eq.1.and.k.eq.1) then
                wlz_grad(i,j,ex,ey) = sqrt(rnu*sqrt(vxz(i,j,k,e)**2 
     >                                + vyz(i,j,k,e)**2))
              endif
              if (ez.eq.lelz.and.k.eq.lz1) then
                wrz_grad(i,j,ex,ey) = sqrt(rnu*sqrt(vxz(i,j,k,e)**2
     >                                + vyz(i,j,k,e)**2))
              endif
            enddo
          enddo
        enddo
      enddo
      call rzero(wx, myz)
      call gop(wrx_grad, wx,'+  ', myz)
      call rzero(wx, myz)
      call gop(wlx_grad, wx,'+  ', myz)
      call rzero(wy, mzx)
      call gop(wly_grad, wy,'+  ', mzx)
      call rzero(wy, mzx)
      call gop(wry_grad, wy,'+  ', mzx)
      call rzero(wz, mxy)
      call gop(wlz_grad, wz,'+  ', mxy)
      call rzero(wz,mxy)
      call gop(wrz_grad, wz,'+  ', mxy)

      !Part 3 - Assigning f1,f2 & f3 based on x+,y+ and z+ values
      ylim   = 100
      f_perp = 0.0116*ylim**2/(1 + 0.203*ylim + 0.0014*(ylim**2.421))
      f_par  = 0.3*ylim/(1 + 0.03*(ylim**1.4))
      do e = 1,lelv
        do k = 1,lz1
          do j = 1,ly1
            do i = 1,lx1
              eg = lglel(e)
              call get_exyz(ex,ey,ez,eg,lelx,lely,lelz)
              if (xl_f(i,j,k,e).eq.1) then 
                 xp = xl(i,j,k,e)*wlx_grad(j,k,ey,ez)
              else
                 xp = xl(i,j,k,e)*wrx_grad(j,k,ey,ez)
              endif
            
              if (yl_f(i,j,k,e).eq.1) then
                 yp = yl(i,j,k,e)*wly_grad(k,i,ez,ex)
              else
                 yp = yl(i,j,k,e)*wry_grad(k,i,ez,ex)
              endif
              
              if (zl_f(i,j,k,e).eq.1) then
                 zp = zl(i,j,k,e)*wlz_grad(i,j,ex,ey)
              else
                 zp = zl(i,j,k,e)*wrz_grad(i,j,ex,ey)
              endif
              xp1(i,j,k,e) = xp;
              yp1(i,j,k,e) = yp;
              zp1(i,j,k,e) = zp;
              
             !Assigning f1,f2,f3 based on proximity to the walls.
              if (xp.gt.100 .and. yp.gt.100 .and. zp.gt.100) then
                f1(i,j,k,e) = 1
                f2(i,j,k,e) = 1
                f3(i,j,k,e) = 1
              elseif (yp.gt.100) then
                if (xp.gt.100 .and. zp.le.100) then
                  f1(i,j,k,e) = 0.3*zp/(1 + 0.03*(zp**1.4))/f_par
                  f2(i,j,k,e) = 0.3*zp/(1 + 0.03*(zp**1.4))/f_par
                  f3(i,j,k,e) = 0.0116*zp**2/(1 + 0.203*zp + 0.0014
     >                          *(zp**2.421))/f_perp
                elseif (xp.le.100 .and. zp.gt.100) then
                  f1(i,j,k,e) = 0.0116*xp**2/(1 + 0.203*xp + 0.0014
     >                          *(xp**2.421))/f_perp
                  f2(i,j,k,e) = 0.3*xp/(1 + 0.03*(xp**1.4))/f_par
                  f3(i,j,k,e) = 0.3*xp/(1 + 0.03*(xp**1.4))/f_par
                else 
                  tmp1 = 0.3*xp/(1 + 0.03*(xp**1.4))/f_par
                  tmp2 = 0.3*zp/(1 + 0.03*(zp**1.4))/f_par
                  f1(i,j,k,e) = 0.0116*xp**2/(1 + 0.203*xp + 0.0014
     >                          *(xp**2.421))/f_perp
                  f2(i,j,k,e) = AMIN1(tmp1,tmp2)
                  f3(i,j,k,e) = 0.0116*zp**2/(1 + 0.203*zp + 0.0014
     >                          *(zp**2.421))/f_perp
                endif
              else
                if (xp.gt.100) then
                  if (zp.gt.100) then
                    f1(i,j,k,e) = 0.3*yp/(1 + 0.03*(yp**1.4))/f_par
                    f2(i,j,k,e) = 0.0116*yp**2/(1 + 0.203*yp + 0.0014 
     >                            *(yp**2.421))/f_perp
                    f3(i,j,k,e) = 0.3*yp/(1 + 0.03*(yp**1.4))/f_par
                  else
                    tmp1 = 0.3*yp/(1 + 0.03*(zp**1.4))/f_par
                    tmp2 = 0.3*zp/(1 + 0.03*(zp**1.4))/f_par
                    f1(i,j,k,e) = AMIN1(tmp1,tmp2)
                    f2(i,j,k,e) = 0.0116*yp**2/(1 + 0.203*yp + 0.0014
     >                            *(yp**2.421))/f_perp
                    f3(i,j,k,e) = 0.0116*zp**2/(1 + 0.203*zp + 0.0014
     >                            *(zp**2.421))/f_perp
                  endif
                else
                  if (zp.gt.100) then
                    tmp1 = 0.3*xp/(1 + 0.03*(xp**1.4))/f_par
                    tmp2 = 0.3*yp/(1 + 0.03*(yp**1.4))/f_par
                    f1(i,j,k,e) = 0.0116*xp**2/(1 + 0.203*xp + 0.0014
     >                            *(xp**2.421))/f_perp
                    f2(i,j,k,e) = 0.0116*yp**2/(1 + 0.203*yp + 0.0014
     >                            *(yp**2.421))/f_perp
                    f3(i,j,k,e) = AMIN1(tmp1,tmp2)
                  else
                    f1(i,j,k,e) = 0.0116*xp**2/(1 + 0.203*xp + 0.0014
     >                            *(xp**2.421))/f_perp
                    f2(i,j,k,e) = 0.0116*yp**2/(1 + 0.203*yp + 0.0014
     >                            *(yp**2.421))/f_perp 
                    f3(i,j,k,e) = 0.0116*zp**2/(1 + 0.203*zp + 0.0014
     >                            *(zp**2.421))/f_perp
                  endif
                endif
              endif
            enddo
          enddo
        enddo
      enddo

c     Multiplying f1,f2,f3 with dissipation for interpolation in ppiclf
      do i=1,n
        f1_diss(i,1,1,1) = sqrt(lan_diss(i,1,1,1))*f1(i,1,1,1)
        f2_diss(i,1,1,1) = sqrt(lan_diss(i,1,1,1))*f2(i,1,1,1)
        f3_diss(i,1,1,1) = sqrt(lan_diss(i,1,1,1))*f3(i,1,1,1)
      enddo
c
c
c     ---------------------------------------------------------
c     RANDOM NUMBERS ARRAY FOR FLUID INLET - KRISHNA 02/26/2022
      call rzero(wy,mzx)
      call rzero(rnum_arr,mzx)
      if (nid.eq.0) then
        if (istep.eq.0) rdum = ran2(-1) ! init random numbers
          do i=1,mzx
            rnum_arr(i,1,1,1) =2.0*ran2(2)-1.0 
          enddo
      endif
      call rzero(wy, mzx)
      call gop(rnum_arr, wy,'+  ', mzx)
c
c
c     ############################################################
c     PPICLF

c     --------------------------------------------
c     PARTICLE INJECTION AT INLET - KRISHNA 02/26/2022
      !Conditions: 1) t = time of sampling  2)t>=0 
      ! 3) t < t_{end of injection} 
      
c      if ((mod(istep,IDNINT(IOPART)).eq.0).and.
c      if ((mod(istep,IDNINT(IOPART)).eq.(IDNINT(IOPART)-1)).and.
      if   ((istep.gt.0).and.
     >     (time.le.usr_endinj)) then
         call my_place_particle(npart,y,rprop) ! All outputs
         ! Add them into solver
         call ppiclf_solve_AddParticles(npart,y,rprop) ! All inputs
      endif

c     -------------------------------------------
c     INTERPOLATE FIELDS TO PARTICLE POSITIONS
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JUX,vx(1,1,1,1))
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JUY,vy(1,1,1,1))
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JUZ,vz(1,1,1,1))
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JFX,f1_diss)
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JFY,f2_diss)
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JFZ,f3_diss)
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JCS2S,lan_cs2s )

c     ------------------------------------------
c     INTEGRATE PARTICLES
      call ppiclf_solve_IntegrateParticle(istep ,
     >                               INT(IOPART),
     >                                    dt    ,
     >                                    time  )

c     -----------------------------------------
c     OUTPUT PAROUT FILES
      if (ppiclf_iostep .gt.0)then
       if (mod(ppiclf_cycle,ppiclf_iostep) .eq. 0) then
        call ppiclf_io_WriteParticleOutletVTU('parout')

        if (.true.) then
          npart_out = 0 !outputs particles out only in between [iostep-1,iostep]
          do i=1,PPICLF_LPART
            do j=1,PPICLF_LRS
              y_out(j,i) = 0.0
            enddo
            do j=1,PPICLF_LRP+1
              rprop_out(j,i) = 0.0
            enddo
            do j=1,5
              iprop_out(j,i) = 0.0
            enddo
          enddo
        endif

       endif
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)

      include 'SIZE'
      include 'TSTEP'
      include 'PARALLEL'
      include 'NEKUSE'


c     ############################################################
c     NEK5000
      real*8 amp,uin,rnum,amp_pert,rdum
      real*8 width_ac,width_in,length_in
      real*8 rlx_in1,rrx_in1,rlz_in1,rrz_in1
     $       rlx_in2,rrx_in2,rlz_in2,rrz_in2
     $       rlx_in3,rrx_in3,rlz_in3,rrz_in3
     $       rlx_in4,rrx_in4,rlz_in4,rrz_in4
      real*8 shp,deg_angle,rad_angle,rpi
      common /inletparam/rlx_in1,rlx_in2, rlx_in3,
     >       rlx_in4,rlz_in1,rlz_in2,rlz_in3,rlz_in4
      common /inletparam2/ uin, rad_angle, width_in, length_in

c     ############################################################
c     PPICLF
      real*8 rlx,rrx,rly,rry,rlz,rrz,dum1,dum2
      common /domainsize/ rlx,rrx,rly,rry,rlz,rrz

c     RANDOM PERTURBATIONS
      real*8 rnum_arr(lx1,lz1,lelx,lelz)
      common /randomarray/ rnum_arr
      integer ex,ey,ez


c     ############################################################
c     USER DEFINES INLET PARAMETERS
c      amp = 1.0d-2
      amp = 3.33d-3
      uin = 7.1428 ! Inlet velocity (m/s)
      amp_pert = 2.5d-1*uin

      ! Inlet size (mts)
      width_ac  = 0.907
      width_in  = 0.055
      length_in = 0.44

      ! Inlet angle (degrees)
      deg_angle = 40.0
c     ############################################################
c
c
c
      uy   = 0.
      uz   = 0.
      temp = 0.

      rpi=4.d0*atan(1.d0)
      rad_angle = rpi*deg_angle/180.0

      ! Vein 1 (closer to x=0)
      rlx_in1 = rrx/2 - width_ac/2 - width_in/2
      rrx_in1 = rrx/2 - width_ac/2 + width_in/2
      rlz_in1 = rrz/2 - length_in/2
      rrz_in1 = rrz/2 + length_in/2

      ! Vein 2 (closer to z=rrz)
      rlx_in2 = rrx/2 - length_in/2
      rrx_in2 = rrx/2 + length_in/2
      rlz_in2 = rrz/2 + width_ac/2 - width_in/2
      rrz_in2 = rrz/2 + width_ac/2 + width_in/2

      ! Vein 3 (closer to x=rrx)
      rlx_in3 = rrx/2 + width_ac/2 - width_in/2
      rrx_in3 = rrx/2 + width_ac/2 + width_in/2
      rlz_in3 = rrz/2 - length_in/2
      rrz_in3 = rrz/2 + length_in/2

      ! Vein 4 (closer to z=0)
      rlx_in4 = rrx/2 - length_in/2
      rrx_in4 = rrx/2 + length_in/2
      rlz_in4 = rrz/2 - width_ac/2 - width_in/2
      rrz_in4 = rrz/2 - width_ac/2 + width_in/2

      if (ifield .eq. 1) then               ! velocity
        if(iside.eq.3) then                 ! Top wall
          shp  =
     >          (0.5*(tanh((x-rlx_in1)/amp) - tanh((x-rrx_in1)/amp))) 
     >        * (0.5*(tanh((z-rlz_in1)/amp) - tanh((z-rrz_in1)/amp)))
     >        + (0.5*(tanh((x-rlx_in2)/amp) - tanh((x-rrx_in2)/amp))) 
     >        * (0.5*(tanh((z-rlz_in2)/amp) - tanh((z-rrz_in2)/amp)))
     >        + (0.5*(tanh((x-rlx_in3)/amp) - tanh((x-rrx_in3)/amp))) 
     >        * (0.5*(tanh((z-rlz_in3)/amp) - tanh((z-rrz_in3)/amp)))
     >        + (0.5*(tanh((x-rlx_in4)/amp) - tanh((x-rrx_in4)/amp))) 
     >        * (0.5*(tanh((z-rlz_in4)/amp) - tanh((z-rrz_in4)/amp)))
          call get_exyz(ex,ey,ez,ieg,lelx,lely,lelz)
          rnum = rnum_arr(ix,iz,ex,ez)
         ! Vein 1
          if (x.ge.rlx_in1.and.x.le.rrx_in1) then
            if (z.ge.rlz_in1.and.z.le.rrz_in1) then
              call get_exyz(ex,ey,ez,ieg,lelx,lely,lelz)
              rnum = rnum_arr(ix,iz,ex,ez)
              uy  = -sin(rad_angle)*uin*shp+(amp_pert*rnum)
              rnum = rnum_arr(ix,iz,ex-1,ez)
              ux  = -cos(rad_angle)*uin*shp+(amp_pert*rnum)
            endif
          endif       

          ! Vein 2
          if (x.ge.rlx_in2.and.x.le.rrx_in2) then
            if (z.ge.rlz_in2.and.z.le.rrz_in2) then
              call get_exyz(ex,ey,ez,ieg,lelx,lely,lelz)
              rnum = rnum_arr(ix,iz,ex,ez)
              uy  = -sin(rad_angle)*uin*shp+(amp_pert*rnum)
              rnum = rnum_arr(ix,iz,ex,ez+1)
              uz  = cos(rad_angle)*uin*shp+(amp_pert*rnum)
            endif
          endif       

          ! Vein 3
          if (x.ge.rlx_in3.and.x.le.rrx_in3) then
            if (z.ge.rlz_in3.and.z.le.rrz_in3) then
              call get_exyz(ex,ey,ez,ieg,lelx,lely,lelz)
              rnum = rnum_arr(ix,iz,ex,ez)
              uy  = -sin(rad_angle)*uin*shp+(amp_pert*rnum)
              rnum = rnum_arr(ix,iz,ex+1,ez)
              ux  = cos(rad_angle)*uin*shp+(amp_pert*rnum)
            endif
          endif     

          ! Vein 4
          if (x.ge.rlx_in4.and.x.le.rrx_in4) then
            if (z.ge.rlz_in4.and.z.le.rrz_in4) then
              call get_exyz(ex,ey,ez,ieg,lelx,lely,lelz)
              rnum = rnum_arr(ix,iz,ex,ez)
              uy  = -sin(rad_angle)*uin*shp+(amp_pert*rnum)
              rnum = rnum_arr(ix,iz,ex,ez-1)
              uz  = -cos(rad_angle)*uin*shp+(amp_pert*rnum)
            endif
          endif     

        endif
      elseif (ifield .eq. 2) then     !temperature
c        temp = 0.0
c        if (iside.eq.5) then  
c          temp  = 0.0
c        endif
      endif
c     ############################################################

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ux = 0.0
      uy = 0.0
      uz = 0.0
      temp = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat
      include 'SIZE'
      include 'TOTAL'
      

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      include 'SIZE'
      include 'TOTAL'
c
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2  !  Modify geometry 

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

c     ############################################################
c     NEK5000
      integer ifc,iel,eg,ex,ey,ez,ntot,nelx_out,nely_out,nelz_out

c     ############################################################
c     PPICLF
      integer i,j,k

      integer*4 ndiam
      parameter(ndiam = 10)                            ! USER DEFINED PARAMETER
      real*8 rhop, dp(ndiam)
      common /dpinfo/ dp, rhop
      data rhop /1000.0/                               ! USER DEFINED PARAMETER
      data dp /0.2d-6,1.0d-6,5.0d-6,10.0e-6,15.0e-6
     >        ,20.0e-6,25.0d-6,30.0e-6,40.0e-6,50.0d-6/! USER DEFINED PARAMETER

      real*8 dpmin
      common /dpm/ dpmin

      real*8 rmu,rhof,rg
      common /parameters/ rmu,rhof,rg

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      real*8 y(PPICLF_LRS    , PPICLF_LPART) ! Normal ordering
      real*8 rprop(PPICLF_LRP, PPICLF_LPART) ! Normal ordering

      real*8 rlx,rrx,rly,rry,rlz,rrz,dum1,dum2
      common /domainsize/ rlx,rrx,rly,rry,rlz,rrz

      real*8 rlx_out,rrx_out,rly_out,rry_out,rlz_out,rrz_out
      common /outletsize/rlx_out,rrx_out,rly_out,rry_out,rlz_out,rrz_out

      integer*4 imethod,iendian,npart

      real*8 rdum, ran2
      external ran2

      real*8 y_out(PPICLF_LRS         , PPICLF_LPART) ! Normal ordering
      real*8 rprop_out(PPICLF_LRP + 1 , PPICLF_LPART) ! Normal ordering
      integer*4 iprop_out(5  , PPICLF_LPART) ! Normal ordering
      integer*4 npart_out
      common /outlet/ y_out,rprop_out,iprop_out,npart_out

      ! Langevin
      real*8 lan_c0,lan_cy,lan_cg,lan_gij,lan_cw
      common /langevin/ lan_c0,lan_cy,lan_cg,lan_gij,lan_cw

c     ############################################################
c     NEK5000
      ntot = lelx*lely*lelz ! Number of elements in BOX 1

c     USER DEFINED PARAMETERS
c     SPECIFIC TO CURRENT CASE
c     MODIFIY IF GEOMETRY CHANGES
      nelx_out = 4          ! Number of elements in BOX 2
      nely_out = 3
      nelz_out = 4
      
c     SET BC BETWEEN TWO BOXES
      do iel=1,nelt

        eg = lglel(iel)

        if (eg.le.ntot) then
c         CHANGE BC IN BOTTOM BOX (BOX 1)
          ifc=3 !top face, if it was generated by genbox
          call get_exyz(ex,ey,ez,eg,lelx,lely,lelz)        
          if (ex.ge.29.and.ex.le.32) then
            if (ez.ge.29.and.ez.le.32) then
              if(cbc(ifc,iel,1) .eq. 'v  ') cbc(ifc,iel,1) = 'E  '
              if(cbc(ifc,iel,2) .eq. 't  ') cbc(ifc,iel,2) = 'E  '
            endif
          endif
          if(cbc(ifc,iel,2) .eq. 't  ') cbc(ifc,iel,2) = 'I  '
        else
c         CHANGE BC IN TOP BOX (BOX 2)
          ifc=1 !bottom face, if it was generated by genbox
          call get_exyz(ex,ey,ez,eg-ntot,nelx_out,nely_out,nelz_out) 
          if (ey.eq.1) then
            if(cbc(ifc,iel,1) .eq. 'v  ') cbc(ifc,iel,1) = 'E  '
            if(cbc(ifc,iel,2) .eq. 't  ') cbc(ifc,iel,2) = 'E  '
          endif

        endif

      enddo

c     ############################################################
c     PPICLF

!     Zero all arrays to track outlet particles
      npart_out = 0
      do i=1,PPICLF_LPART

        do j=1,PPICLF_LRS
          y_out(j,i) = 0.0
        enddo

        do j=1,PPICLF_LRS+1
          rprop_out(j,i) = 0.0
        enddo

        do j=1,5
          iprop_out(j,i) = 0.0
        enddo

      enddo

      ! Pass to library to Init MPI
      call ppiclf_comm_InitMPI(nekcomm,
     >                         nid    , ! nid already defined in Nek5000
     >                         np     ) ! np already defined in Nek5000

c     USER DEFINED PARAMETERS
      ! Set initial conditions and parameters for particles
      imethod = 1
      ndim    = 3
      iendian = 0
      npart   = 4000 ! NUMBER OF PARTICLES PER RANK PER DIAMETER ( NTOT=npart*size(dp)*np )
      dpmin = 11.0d-6

      ! Domain size
      rlx     =  0.0d0
      rrx     =  10.0d0
      rly     =  0.0d0
      rry     =  3.2d0
      rlz     =  0.0d0
      rrz     =  10.0d0

      ! Outlet size
      rlx_out = 4.7
      rrx_out = 5.3
      rly_out = 3.2
      rry_out = 3.7
      rlz_out = 4.7
      rrz_out = 5.3



      ! Fluid and physics parameters
      rhof = param(1)
      rmu  = param(2)
      rg   = -9.8

      ! Langevin parameters
      lan_c0  = 2.1d0
      lan_cy  = 0.039d0
      lan_cg=(5.d-1+7.5d-1*lan_c0)/(2.d0*lan_cy)

      k = 0
      rdum    = ran2(-1-nid) ! init random numbers
      do i=1,npart
        do j=1,ndiam
          k = k + 1

          y(PPICLF_JX,k)  = rlx + dp(j)*2.0 + (rrx - rlx - dp(j)*4.0) 
     $                     * ran2(2)                                 ! NEED TO CHECK SEED FOR RAN2
c          y(PPICLF_JY,k)  = 1.0 
          y(PPICLF_JY,k)  = rly + dp(j)*2.0 + (rry - rly - dp(j)*4.0) 
     $                     * ran2(2)
          y(PPICLF_JZ,k)  = rlz + dp(j)*2.0 + (rrz - rlz - dp(j)*4.0) 
     $                     * ran2(2)

          rprop(PPICLF_R_JRHOP,k) = rhop
          rprop(PPICLF_R_JDP  ,k) = dp(j)
          rprop(PPICLF_R_JVOLP,k) = pi/6.0D0*rprop(PPICLF_R_JDP,k)**3
          rprop(PPICLF_R_JVS  ,k) = rhop*dp(j)**2*rg/18./rmu
c          rprop(PPICLF_R_JVS  ,k) = vsettling(j) 

          y(PPICLF_JVX,k) = 0.0
          y(PPICLF_JVY,k) = 0.0 !rprop(PPICLF_R_JVS  ,k)
          y(PPICLF_JVZ,k) = 0.0
        enddo
      enddo

c      call ppiclf_io_ReadParticleVTU('parIC.vtu')

      call ppiclf_solve_InitParticle(imethod   ,
     >                               ndim      ,
     >                               iendian   ,
     >                               k         ,
     >                               y(1,1)    ,
     >                               rprop(1,1))

c      call ppiclf_solve_InitTargetBins('x',2,1)
c      call ppiclf_solve_InitTargetBins('y',2,1)


      ! Specify Overlap Mesh
      call ppiclf_comm_InitOverlapMesh(nelt,lx1,ly1,lz1,xm1,ym1,zm1)

      return
      end


!-----------------------------------------------------------------------
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
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
      ran2=min(AM*iy,RNMX)
      return
      END

c ######################################################################
c     SMAGORINSKY ROUTINES

c-----------------------------------------------------------------------
      subroutine eddy_visc(ediff,e)
c
c     Compute eddy viscosity using dynamic smagorinsky model
c
      include 'SIZE'
      include 'TOTAL'

      real ediff(nx1*ny1*nz1,nelv)
      integer e

      common /dynsmg/ sij (lx1*ly1*lz1,ldim,ldim)
     $              , mij (lx1*ly1*lz1,3*ldim-3)
     $              , lij (lx1*ly1*lz1,3*ldim-3)
     $              , dg2 (lx1*ly1*lz1,lelv)
     $              , num (lx1*ly1*lz1,lelv)
     $              , den (lx1*ly1*lz1,lelv)
     $              , snrm(lx1*ly1*lz1,lelv)
      real sij,mij,lij,dg2,num,den,snrm

      parameter(lxyz=lx1*ly1*lz1)
      common /xzmp0/ ur (lxyz) , us (lxyz) , ut (lxyz)
      real           vr (lxyz) , vs (lxyz) , vt (lxyz)
     $     ,         wr (lxyz) , ws (lxyz) , wt (lxyz)
      common /xzmp1/ w1(lx1*lelv),w2(lx1*lelv)

      !! NOTE CAREFUL USE OF EQUIVALENCE HERE !!
      equivalence (vr,lij(1,1)),(vs,lij(1,2)),(vt,lij(1,3))
     $          , (wr,lij(1,4)),(ws,lij(1,5)),(wt,lij(1,6))

      common /sgsflt/ fh(lx1*lx1),fht(lx1*lx1),diag(lx1)

      integer nt
      save    nt
      data    nt / -9 /

      ntot = nx1*ny1*nz1

      if (nt.lt.0) call
     $   set_ds_filt(fh,fht,nt,diag,nx1)! dyn. Smagorinsky filter

      call comp_gije(sij,vx(1,1,1,e),vy(1,1,1,e),vz(1,1,1,e),e)
      call comp_sije(sij)

      call mag_tensor_e(snrm(1,e),sij)
      call cmult(snrm(1,e),2.0,ntot)

      call comp_mij   (mij,sij,dg2,ur,us,fh,fht,nt,e)
      call comp_lij   (lij,vx,vy,vz,ur,us,ut,fh,fht,e)

c     Compute numerator (ur) & denominator (us) for Lilly contraction

      n = nx1*ny1*nz1
      do i=1,n
         ur(i) = mij(i,1)*lij(i,1)+mij(i,2)*lij(i,2)+mij(i,3)*lij(i,3)
     $      + 2*(mij(i,4)*lij(i,4)+mij(i,5)*lij(i,5)+mij(i,6)*lij(i,6))
         us(i) = mij(i,1)*mij(i,1)+mij(i,2)*mij(i,2)+mij(i,3)*mij(i,3)
     $      + 2*(mij(i,4)*mij(i,4)+mij(i,5)*mij(i,5)+mij(i,6)*mij(i,6))
      enddo

c     smoothing numerator and denominator in time
      call copy (vr,ur,nx1*nx1*nx1)
      call copy (vs,us,nx1*nx1*nx1)

      beta1 = 0.0                   ! Temporal averaging coefficients
      if (istep.gt.1) beta1 = 0.9   ! Retain 90 percent of past
      beta2 = 1. - beta1

      do i=1,n
         num (i,e) = beta1*num(i,e) + beta2*vr(i)
         den (i,e) = beta1*den(i,e) + beta2*vs(i)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine eddy_visc_01(ediff)
c
c     Compute eddy viscosity using dynamic smagorinsky model
c
      include 'SIZE'
      include 'TOTAL'

      real ediff(nx1*ny1*nz1,nelv)!,uni(lx1*ly1*lz1,nelv)

      common /dynsmg/ sij (lx1*ly1*lz1,ldim,ldim)
     $              , mij (lx1*ly1*lz1,3*ldim-3)
     $              , lij (lx1*ly1*lz1,3*ldim-3)
     $              , dg2 (lx1*ly1*lz1,lelv)
     $              , num (lx1*ly1*lz1,lelv)
     $              , den (lx1*ly1*lz1,lelv)
     $              , snrm(lx1*ly1*lz1,lelv)
      real sij,mij,lij,dg2,num,den,snrm

      real dumr
      integer n

      ! Langevin
      real*8  lan_diss(lx1,ly1,lz1,lelt)
      real*8  lan_cs2s(lx1,ly1,lz1,lelt)
      common /lan/ lan_diss,lan_cs2s

      do i=1,lx1*ly1*lz1*lelv
         cdyn = 0
         if (den(i,1).gt.0) cdyn = 0.5*num(i,1)/den(i,1)
c         cdyn = cdyn*uni(i,1)
c         cdyn = max(cdyn,0.)
         ediff(i,1) = param(2)+max(cdyn*dg2(i,1)*snrm(i,1),0.d0)

         lan_diss(i,1,1,1)=max(cdyn*dg2(i,1)*snrm(i,1),0.d0)
     $                    *snrm(i,1)**2
c         lan_cs2s(i,1,1,1)= cdyn*snrm(i,1)
         lan_cs2s(i,1,1,1)= snrm(i,1)*0.18**2.0
      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine set_ds_filt(fh,fht,nt,diag,nx) ! setup test filter

      INCLUDE 'SIZE'

      real fh(nx*nx),fht(nx*nx),diag(nx)

c Construct transfer function
      call rone(diag,nx)

      diag(nx-0) = 0.01
      diag(nx-1) = 0.10
      diag(nx-2) = 0.50
      diag(nx-3) = 0.90
      diag(nx-4) = 0.99
      nt = nx - 2

c      diag(nx-0) = 0.05
c      diag(nx-1) = 0.50
c      diag(nx-2) = 0.95
c      nt = nx - 1

      call build_1d_filt(fh,fht,diag,nx,nid)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_grid_spacing
c
c     Compute D^2, the grid spacing used in the DS sgs model.
c
      include 'SIZE'
      include 'TOTAL'


      common /dynsmg/ sij (lx1*ly1*lz1,ldim,ldim)
     $              , mij (lx1*ly1*lz1,3*ldim-3)
     $              , lij (lx1*ly1*lz1,3*ldim-3)
     $              , dg2 (lx1*ly1*lz1,lelv)
     $              , num (lx1*ly1*lz1,lelv)
     $              , den (lx1*ly1*lz1,lelv)
     $              , snrm(lx1*ly1*lz1,lelv)
      real sij,mij,lij,dg2,num,den,snrm

      real dxmax, dymax, dzmax, rdum
      real damping,eps

      integer im,ip

      real*8 rlx,rrx,rly,rry,rlz,rrz
      common /domainsize/ rlx,rrx,rly,rry,rlz,rrz

      integer e,eg,ex,ey,ez,n,m

      gamma = 1.
      gamma = gamma/ndim

      n = nx1*ny1*nz1*nelv
      call rone(dg2,n)
c      call cmult(dg2,1.0d-3,n)
      return          ! Comment this line for a non-trivial Delta defn
c
c      eps = 0.12 ! From Özgökmen, Fischer 2009
c      dxmax = 1.0
c      dymax = 1.0
c      dzmax = 1.0
c      do e=1,nelv
cC          in the gll mesh the largest gridspace ins at the middle of the element
c        im = int(nx1/2)
c        ip = int(nx1/2) + 1
c        dxmax = xm1(ip,0,0,e) - xm1(im,0,0,e) 
c        dymax = ym1(0,ip,0,e) - ym1(0,im,0,e) 
c        dzmax = zm1(0,0,ip,e) - zm1(0,0,im,e) 
c        do m=1,nx1*ny1*nz1
C              damping = (1.0 - exp(-(xm1(m,1,1,e) - rlx)/eps)**2)
C     >                * (1.0 - exp(-(xm1(m,1,1,e) - rrx)/eps)**2)
C     >                * (1.0 - exp(-(ym1(m,1,1,e) - rly)/eps)**2)
C     >                * (1.0 - exp(-(ym1(m,1,1,e) - rry)/eps)**2)
C     >                * (1.0 - exp(-(zm1(m,1,1,e) - rlz)/eps)**2)
C     >                * (1.0 - exp(-(zm1(m,1,1,e) - rrz)/eps)**2)
c              damping = 1.0         
c              dg2(m,e) = damping * (dxmax*dymax*dzmax)**gamma
c              dg2(m,e) =  dg2(m,e)**2
c        enddo
c      enddo

      call dsavg(dg2)  ! average neighboring elements

      return
      end

c-----------------------------------------------------------------------
      subroutine comp_lij(lij,u,v,w,fu,fv,fw,fh,fht,e)
c
c     Compute Lij for dynamic Smagorinsky model:
c                    _   _      _______
c          L_ij  :=  u_i u_j  - u_i u_j
c
      include 'SIZE'
c
      integer e
c
      real lij(lx1*ly1*lz1,3*ldim-3)
      real u  (lx1*ly1*lz1,lelv)
      real v  (lx1*ly1*lz1,lelv)
      real w  (lx1*ly1*lz1,lelv)
      real fu (1) , fv (1) , fw (1)
     $   , fh (1) , fht(1)

      call tens3d1(fu,u(1,e),fh,fht,nx1,nx1)  ! fh x fh x fh x u
      call tens3d1(fv,v(1,e),fh,fht,nx1,nx1)
      call tens3d1(fw,w(1,e),fh,fht,nx1,nx1)

      n = nx1*ny1*nz1
      do i=1,n
         lij(i,1) = fu(i)*fu(i)
         lij(i,2) = fv(i)*fv(i)
         lij(i,3) = fw(i)*fw(i)
         lij(i,4) = fu(i)*fv(i)
         lij(i,5) = fv(i)*fw(i)
         lij(i,6) = fw(i)*fu(i)
      enddo

      call col3   (fu,u(1,e),u(1,e),n)    !  _______
      call tens3d1(fv,fu,fh,fht,nx1,nx1)  !  u_1 u_1
      call sub2   (lij(1,1),fv,n)

      call col3   (fu,v(1,e),v(1,e),n)    !  _______
      call tens3d1(fv,fu,fh,fht,nx1,nx1)  !  u_2 u_2
      call sub2   (lij(1,2),fv,n)

      call col3   (fu,w(1,e),w(1,e),n)    !  _______
      call tens3d1(fv,fu,fh,fht,nx1,nx1)  !  u_3 u_3
      call sub2   (lij(1,3),fv,n)

      call col3   (fu,u(1,e),v(1,e),n)    !  _______
      call tens3d1(fv,fu,fh,fht,nx1,nx1)  !  u_1 u_2
      call sub2   (lij(1,4),fv,n)

      call col3   (fu,v(1,e),w(1,e),n)    !  _______
      call tens3d1(fv,fu,fh,fht,nx1,nx1)  !  u_2 u_3
      call sub2   (lij(1,5),fv,n)

      call col3   (fu,w(1,e),u(1,e),n)    !  _______
      call tens3d1(fv,fu,fh,fht,nx1,nx1)  !  u_3 u_1
      call sub2   (lij(1,6),fv,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_mij(mij,sij,dg2,fs,fi,fh,fht,nt,e)
c
c     Compute Mij for dynamic Smagorinsky model:
c
c                     2 _  ____     _______
c          M_ij  :=  a  S  S_ij  -  S  S_ij
c
      include 'SIZE'
c
      integer e
c
      real mij(lx1*ly1*lz1,3*ldim-3)
      real dg2(lx1*ly1*lz1,lelv)
      real fs (1) , fi (1) , fh (1) , fht(1)

      real magS(lx1*ly1*lz1)
      real sij (lx1*ly1*lz1*ldim*ldim)

      integer imap(6)
      data imap / 0,4,8,1,5,2 /

      integer n
      real dumr

      n = nx1*ny1*nz1

      call mag_tensor_e(magS,sij)
      call cmult(magS,2.0,n)

c     Filter S
      call tens3d1(fs,magS,fh,fht,nx1,nx1)  ! fh x fh x fh x |S|

c     a2 is the test- to grid-filter ratio, squared

      a2 = nx1-1       ! nx1-1 is number of spaces in grid
      a2 = a2 /(nt-1)  ! nt-1 is number of spaces in filtered grid

      do k=1,6
         jj = n*imap(k) + 1
         call col3   (fi,magS,sij(jj),n)
         call tens3d1(mij(1,k),fi,fh,fht,nx1,nx1) ! fh x fh x fh x (|S|S_ij)
         call tens3d1(fi,sij(jj),fh,fht,nx1,nx1)  ! fh x fh x fh x S_ij
         do i=1,n
            mij(i,k) = (a2**2 * fs(i)*fi(i) - mij(i,k))*dg2(i,e)
         enddo
      enddo


c      n = nx1*ny1*nz1*lelv
c      dumr = glmin(dg2,n)
c      if (nid.eq.0.and.istep.eq.0) write(999,*) 'Min', dumr
c
c      dumr = glmax(dg2,n)
c      if (nid.eq.0.and.istep.eq.0) write(999,*) 'Max', dumr

      return
      end

c-----------------------------------------------------------------------
      subroutine smooth_fld(u,ddir)

      include 'SIZE'
      include 'TOTAL'

      real*8 u(lx1,ly1,lz1,nelv),udum(lx1,ly1,lz1,nelv)
      integer i,j,k,im,jm,km,ie,ddir

      nxyz = nx1*ny1*nz1
      n    = nxyz*nelv
      call copy(udum,u,n)

      if (ddir.eq.1) then
c       x derivative
        do ie=1,nelv
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            if (i.eq.1) then
              im = i+1
              udum(i,j,k,ie) = u(im,j,k,ie)
            elseif (i.eq.nx1) then
              im = i-1
              udum(i,j,k,ie) = u(im,j,k,ie)
            endif
         enddo
         enddo
         enddo
        enddo
      elseif (ddir.eq.2) then
c       y derivative
        do ie=1,nelv
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            if (j.eq.1) then
              jm = j+1
              udum(i,j,k,ie) = u(i,jm,k,ie)
            elseif (j.eq.ny1) then
              jm = j-1
              udum(i,j,k,ie) = u(i,jm,k,ie)
            endif
         enddo
         enddo
         enddo
        enddo
      elseif (ddir.eq.3) then
c       z derivative
        do ie=1,nelv
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            if (k.eq.1) then
              km = k+1
              udum(i,j,k,ie) = u(i,j,km,ie)
            elseif (k.eq.nz1) then
              km = k-1
              udum(i,j,k,ie) = u(i,j,km,ie)
            endif
         enddo
         enddo
         enddo
        enddo

      endif

      call copy(u,udum,n)
      call dsavg(u)

      return
      end

c-----------------------------------------------------------------------
      subroutine filter_s0mine(scalar,wght,ncut,name5) ! filter scalar field

      include 'SIZE'
      include 'TOTAL'

      real scalar(1)
      character*5 name5

      parameter (l1=lx1*lx1)
      real intdv(l1),intuv(l1),intdp(l1),intup(l1),intv(l1),intp(l1)
      save intdv    ,intuv    ,intdp    ,intup    ,intv    ,intp

      common /ctmp0/ intt
      common /screv/ wk1,wk2
      common /scrvh/ zgmv,wgtv,zgmp,wgtp,tmax(100)

      real intt (lx1,lx1)
      real wk1  (lx1,lx1,lx1,lelt)
      real wk2  (lx1,lx1,lx1)
      real zgmv (lx1),wgtv(lx1),zgmp(lx1),wgtp(lx1)


      integer icall
      save    icall
      data    icall /0/

      logical ifdmpflt

      imax = nid
      imax = iglmax(imax,1)
      jmax = iglmax(imax,1)

      if (icall.eq.0) call build_new_filter(intv,zgm1,lx1,ncut,wght,nio)
C       call build_new_filter(intv,zgm1,lx1,ncut,wght,nio)

      icall = 1

      call filterq(scalar,intv,lx1,lz1,wk1,wk2,intt,if3d,fmax)
      fmax = glmax(fmax,1)

      if (nio.eq.0) write(6,1) istep,fmax,name5
    1 format(i8,' sfilt:',1pe12.4,a10)

      return
      end

c     ###############################################################
c     Particle Injection at inlet
      subroutine my_place_particle(k,y,rprop)
      include 'SIZE'
      include 'TOTAL'

c     Need to add rad_angle, rry, rlx_in1,..., rlz_in1,..., length_in,
c     width_in,particle diameters and uin to list of variables/commons.


      integer i,j,k,l,k2
      real*8  y(*)
      real*8  rprop(*)
      real*8  y1,v1,tmp,rdum,wid1,len1
      real*8 uin, rad_angle,rlx_in1,rlx_in2, rlx_in3,
     >       rlx_in4,rlz_in1,rlz_in2,rlz_in3,rlz_in4,
     >       length_in,width_in
      common /inletparam/ rlx_in1,rlx_in2, rlx_in3,
     >       rlx_in4,rlz_in1,rlz_in2,rlz_in3,rlz_in4
      common /inletparam2/ uin, rad_angle, width_in, length_in
      real*8 rlx,rrx,rly,rry,rlz,rrz
      common /domainsize/ rlx,rrx,rly,rry,rlz,rrz
      integer*4 ndiam
      parameter(ndiam=10)
      real*8 rhop, dp(ndiam)
      common /dpinfo/ dp, rhop
 
      real*8 rmu,rhof,rg
      common /parameters/ rmu,rhof,rg
           
      if (istep.le.IDNINT(IOPART)) rdum = ran2(-1-nid)
     
      k=0
      v1 = 0.
      tmp = 0.
      k2 = 0

      if (nid.eq.0) then
        k = 0
        v1 = -uin*sin(rad_angle)
        tmp = uin*cos(rad_angle)

        !Inlet 1
        do j=1,ndiam
          do i=1,npart_inj
            k = k+1
            !Randomize y1 - distance from top wall
            y1 = (rry - dp(j)*2.0) - sin(rad_angle)*srcin_len*ran2(2)
            len1 = length_in*ran2(2)         !Varies from 0 to length_in
            wid1 = width_in*ran2(2)          !Varies from 0 to width_in
            x1 = (rlx_in1 - (rry-y1)/tan(rad_angle)) + wid1
            z1 = rlz_in1 + len1
            k2 = PPICLF_LRP*(k-1)
            rprop(k2+PPICLF_R_JRHOP) = rhop
            rprop(k2+PPICLF_R_JDP) = dp(j)
            rprop(k2+PPICLF_R_JVOLP) = pi/6.0D0*dp(j)**3
            rprop(k2+PPICLF_R_JVS) = rhop*dp(j)**2*rg
     >                               /18./rmu
              
            k2 = PPICLF_LRS*(k-1)
            y(k2+PPICLF_JX)  = (rlx_in1 - (rry-y1)/tan(rad_angle))
     >                         + wid1
            y(k2+PPICLF_JY)  = y1
            y(k2+PPICLF_JZ)  = rlz_in1 + len1
            y(k2+PPICLF_JVX) = -1.*tmp
            y(k2+PPICLF_JVY) = v1
            y(k2+PPICLF_JVZ) = 0.
          enddo
        enddo

        !Inlet 2
        do j=1,ndiam
          do i=1,npart_inj
            k = k+1

            y1 = (rry - dp(j)*2.0) - sin(rad_angle)*srcin_len*ran2(2)
            len1 = length_in*ran2(2)         !Varies from 0 to length_in
            wid1 = width_in*ran2(2)          !Varies from 0 to width_in
  
            k2 = PPICLF_LRP*(k-1)
            rprop(k2+PPICLF_R_JRHOP) = rhop
            rprop(k2+PPICLF_R_JDP) = dp(j)
            rprop(k2+PPICLF_R_JVOLP) = pi/6.0D0*dp(j)**3
            rprop(k2+PPICLF_R_JVS) = rhop*dp(j)**2*rg
     >                               /18./rmu            
            k2 = PPICLF_LRS*(k-1)
            y(k2+PPICLF_JX)  = rlx_in2 + len1
            y(k2+PPICLF_JY)  = y1
            y(k2+PPICLF_JZ)  = (rlz_in2 + (rry-y1)/tan(rad_angle))
     >                         + wid1
            y(k2+PPICLF_JVX) = 0.
            y(k2+PPICLF_JVY) = v1
            y(k2+PPICLF_JVZ) = tmp
          enddo
        enddo
        
        !Inlet 3
        do j=1,ndiam
          do i=1,npart_inj
            k = k+1
            y1 = (rry - dp(j)*2.0) - sin(rad_angle)*srcin_len*ran2(2)
            len1 = length_in*ran2(2)         !Varies from 0 to length_in
            wid1 = width_in*ran2(2)          !Varies from 0 to width_in
  
            k2 = PPICLF_LRP*(k-1)
            rprop(k2+PPICLF_R_JRHOP) = rhop
            rprop(k2+PPICLF_R_JDP) = dp(j)
            rprop(k2+PPICLF_R_JVOLP) = pi/6.0D0*dp(j)**3
            rprop(k2+PPICLF_R_JVS) = rhop*dp(j)**2*rg
     >                               /18./rmu
              
            k2 = PPICLF_LRS*(k-1)
            y(k2+PPICLF_JX)  = (rlx_in3 + (rry-y1)/tan(rad_angle))
     >                         + wid1
            y(k2+PPICLF_JY)  = y1
            y(k2+PPICLF_JZ)  = rlz_in3 + len1
            y(k2+PPICLF_JVX) = tmp
            y(k2+PPICLF_JVY) = v1
            y(k2+PPICLF_JVZ) = 0.
          enddo
        enddo

        !Inlet 4
        do j=1,ndiam
          do i=1,npart_inj
            k = k+1
            y1 = (rry - dp(j)*2.0) - sin(rad_angle)*srcin_len*ran2(2)
           !Randomize the y location        
            len1 = length_in*ran2(2)         !Varies from 0 to length_in
            wid1 = width_in*ran2(2)          !Varies from 0 to width_in
  
            k2 = PPICLF_LRP*(k-1)
            rprop(k2+PPICLF_R_JRHOP) = rhop
            rprop(k2+PPICLF_R_JDP) = dp(j)
            rprop(k2+PPICLF_R_JVOLP) = pi/6.0D0*dp(j)**3
            rprop(k2+PPICLF_R_JVS) = rhop*dp(j)**2*rg
     >                               /18./rmu            
            k2 = PPICLF_LRS*(k-1)
            y(k2+PPICLF_JX)  = rlx_in4 + len1
            y(k2+PPICLF_JY)  = y1
            y(k2+PPICLF_JZ)  = (rlz_in4 - (rry-y1)/tan(rad_angle))
     >                         + wid1
            y(k2+PPICLF_JVX) = 0.
            y(k2+PPICLF_JVY) = v1
            y(k2+PPICLF_JVZ) = -1.*tmp
          enddo
        enddo
      endif

      return
      end


c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
