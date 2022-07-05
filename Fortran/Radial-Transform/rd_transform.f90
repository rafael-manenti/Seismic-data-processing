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
      Subroutine get_rd_sec(ns,nt,x,x0,t0,thetai,dtheta,ntheta,v)
      complex x,v
      dimension x(ns,nt),v(ns,ntheta)

      pi=acos(-1.0)
      ti=thetai*pi/180 ! initial angle in rad
      dt=(dtheta*pi/180)/(ntheta-1) ! delta angle between each trace
      v=0.0
      i0=t0 
      a_nt=nt
 
      do j=1,ntheta
      
         theta=ti+dt*(j-1) ; a=tan(pi*0.5-theta)
         do i=i0,ns
         
            x1=x0 + a*(i-i0)
            
            if (x1.gt.a_nt.or.x1.lt.0.0) exit
            
            j1=int(x1)
            alpha=x1-j1
            beta=1.0-alpha
            c=(alpha**2)+(beta**2) ! Shepard interpolation with p = 2 using only two samples
            v(i,j)=((beta**2)/c)*x(i,j1)+((alpha**2)/c)*x(i,j1+1)

         end do
      end do

      return
      end

