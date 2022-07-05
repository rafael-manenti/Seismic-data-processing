!==================================================================================
! Rafael Rodrigues Manenti
! Subroutine used to obtain the radial transform panel
! nt -> number of traces
! ns -> number of samples
! x ->  original data panel (usually in X-T domain)
! ntheta -> number of radial traces
! v ->  data in the radial domain (radial panel)
! x0 -> x coordinate for the origin of the radial traces
! t0 -> t coordinate for the origin of the radial traces
! thetai -> initial angle
! dtheta -> angle aperture
!==================================================================================

      Subroutine get_rd_inverse(ns,nt,ntheta,x0,t0,thetai,dtheta,v,x)
      complex x,v
      dimension x(ns,nt),v(ns,ntheta)
      allocatable rtheta(:)
      allocate (rtheta(ntheta))

      i0=t0
      a_nt=nt
      j0=x0
      pi=acos(-1.0)
      ti=thetai*pi/180
      dt=(dtheta*pi/180)/(ntheta-1)

      do i=1,ntheta
         rtheta(i)=ti+dt*(i-1)
      end do

         do i=1,ns
            do j=1,ntheta-1
               do k= 1,nt
                  a=tan(pi*0.5-rtheta(j))
                  b=tan(pi*0.5-rtheta(j+1))
                  x1=x0+a*(i-t0) !  radial trace length
                  x2=x0+b*(i-t0)
                  m=int(x1)
                  n=int(x2)
                  if (m.ne.k) go to 10 
                  alpha=x1-m
                  beta=m-x2
                  c=alpha**2+beta**2
                  x(i,k)=(beta**2)/c*v(i,j)+(alpha**2)/c*v(i,j+1)
 10               continue
               end do
            end do
         end do
      x(i0,j0)=0.
      deallocate (rtheta)

      return
      end  
