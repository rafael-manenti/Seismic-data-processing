
!============================================================
!	Subroutine for power method SVD
!	It was created to obtain only the first eigen image, but it can be easily extended for other eigen images
!	09/February/2013
!	Rafael Rodrigues Manenti
!	m and n are the matrix indices
!	A is input matrix
!	Lambda is dominant eigenvalue, or first eigenvalue
!	x is the dominant eigenvector or first eigenvector
!	n_iter is number of iterations for calculating lambda and x.
!	n_iter is set to 50
!-----------------------------------------------------------------------
! 	modified by Milton J. Porsani 12/02/2013 (terça de carnaval)
!============================================================
     subroutine svd_power(m,n,X,sigma,U,V)
      dimension  U(m), V(n), X(m,n)
      allocatable XTX(:,:),x1(:),xaux(:)
      niter_max=50 ; tol=1.e-15     
      allocate (XTX(n,n),x1(n),xaux(n))
!------------------------------------------------------
      call xtx_xty(m,n,X,xtx)
 
      x1=0.; xaux=1
      call matvetmult(n,n,xtx,xaux,x1)
      Q0=dot_product(xaux,x1)
      xaux=x1/sqrt(Q0)
!=========================================================================================
     ikey=0
     k=0
     do while(ikey.eq.0)
        k=k+1                          
         call matvetmult(n,n,XTX,xaux,x1)
         Q1=dot_product(xaux,x1)      
         xaux=x1/sqrt(Q1)
!                                   write(*,*)k,Q1,Q0,Q1-q0,sqrt(Q1)                        
         if(abs(Q1-Q0).lt.tol.or.k.eq.niter_max)ikey=1  
         Q0=Q1                              
      end do    
                             
!=========================================================================================        
      xx=dot_product(xaux,xaux)                                                                 !
      slambda=Q1/xx              ! sp positive                                                   
      sigma=sqrt(slambda) 

!=========================================================================================
      V(1:n)=xaux(1:n)/sqrt(xx)      !normalizing for unitary energy p/ ter energia unitária                        
      call matvetmult (m,n,X,v,u)    ! generating eigenvector
      u=u/sigma                                         !                !
!=========================================================================================
      deallocate (XTX,x1,xaux)
      return
      end
