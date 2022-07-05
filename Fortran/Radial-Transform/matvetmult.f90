!=========================================================================================
!
!	            SUBROUTINE TO MULTIPLY MATRIX BY VECTOR
!	
!=========================================================================================
      Subroutine matvetmult (m,n,a,b,c)
      real a(m,n),b(n),c(m)
      c=0.
      do i=1,m
         do j=1,n
            c(i)=c(i)+a(i,j)*b(j)
         end do
      end do
      return
      end
