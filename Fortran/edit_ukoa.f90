 !==========================================================================================
 !
 !	Program used to edit the ukoa information of the land seismic data
 !	Used to upload info to header in ProMAX SeisSpace
 !
 !==========================================================================================

 character*264 arq1,arq2,arq3,arq4 
allocatable x(:),y(:),z(:),xedit(:),yedit(:),zedit(:)
m=966
n=m*2
n1=1
arq1='ukoa.edit'
arq2='ukoa-interp'
allocate (x(m),y(m),z(m),xedit(n),yedit(n),zedit(n))

    open(10,file=arq1)
    open(20,file=arq2)

    do i=1,m
        read(10,*)x(i),y(i),z(i)
    end do
   
    do i=1,m
       xedit(n1)=x(i)
       yedit(n1)=y(i)
       zedit(n1)=z(i)
       n1=n1+1
       xedit(n1)=0.5*(x(i)+x(i+1))
       yedit(n1)=0.5*(y(i)+y(i+1))
       zedit(n1)=0.5*(z(i)+z(i+1))
       n1=n1+1
    end do

    xedit(n)=x(m)
    yedit(n)=y(m)
    zedit(n)=z(m)

   do i=1,n
      write(20,*)xedit(i),yedit(i),zedit(i)
   end do

   close(10)
   close(20)

   deallocate(x,y,z,xedit,yedit,zedit)
   end  
