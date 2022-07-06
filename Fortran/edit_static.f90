 !==========================================================================================
 !
 !	Program used to edit the statics information of the land seismic data
 !	Used to upload info to header in ProMAX SeisSpace
 !
 !==========================================================================================
     character*264 arq1,arq2,arq3,arq4 
     integer, allocatable:: sx(:),gx(:),sstat(:),gstat(:),sstat2(:),gstat2(:),sx2(:),gx2(:)

arq1='sstat_full_edit_part1'
arq2='gstat_full_edit_part1'
arq3='hdrfile-rl008'
arq4='hdrfile-update-rl008'

m=47434 !vetor header
m1=16718 !vetor fonte
m2=16876 !vetor receptor

     open(10,file=arq1)
     open(20,file=arq2)
     open(30,file=arq3)
     open(40,file=arq4)     

    allocate (sx(m1),gx(m2),sstat(m1),gstat(m2),sstat2(m),gstat2(m),sx2(m),gx2(m))

       do i=1,m
          read(30,*)sx2(i),gx2(i),sstat2(i),gstat2(i)
       end do

    do i=1,m1
       read(10,*)sx(i),sstat(i)
    end do

    do i=1,m2
       read(20,*)gx(i),gstat(i)
    end do    



    do i=1,m
   
       do j=1,m2
          if(gx(j).eq.gx2(i)) gstat2(i)=gstat(j)
       end do
      
       do j=1,m1
          if(sx(j).eq.sx2(i)) sstat2(i)=sstat(j)
       end do

          write(40,*)sx2(i),gx2(i),sstat2(i),gstat2(i)

    end do

    deallocate (sx,gx,sstat,gstat,sstat2,gstat2,sx2,gx2)
    end
    
  
