!==================================================================
!
!
!	Program to apply the forward and inverse radial-wavelet transform
!	Wavelet transform was all written by Lucas Andrade de Almeida
!
!==================================================================	
	
	allocatable x(:,:), XL(:,:,:), XH(:,:,:), HL(:,:,:), HH(:,:,:)
      allocatable g(:), h(:), HX(:,:), a(:,:), g1(:), h1(:)
      allocatable b(:,:), xdec(:), HLaux(:,:), HHaux(:,:), HLfilt(:,:)
      allocatable XHaux(:,:), dif(:,:), T(:,:), HXaux(:,:),XHfilt(:,:)
      allocatable HHfilt(:,:), vx(:,:),vi(:,:), HLaux2(:),xaux2(:,:) 
      allocatable xaux(:,:),vaux(:,:),vf(:,:),xf(:,:),vdec(:),d(:,:)
      ndec=3
      n=1024
      nt=1024
      ip=nt/2
      igmin=0
      igmax=11
      ihmin=0
      ihmax=11
      nshift=0
      !ntiros=362
      ntiros=1
      it=96
      ia=1001
      iat=ia+nshift

      ncf=2
      !ncfLH=2
      lw=20
      t0=1 ! amostra onde ocorre a primeira chegada da onda direta
      x0=48 ! traço onde ocorre a primeira chegada da onda direta
      ntheta=1000 ! n. de ângulos que se deseja mapear
      thetai=76. ! ângulo inicial
      thetaf=102. ! ângulo final
      dtheta=thetaf-thetai ! abertura do ângulo que se deseja varrer
      len=4*nft
      allocate(xaux(ia,it),vaux(iat,ntheta),vf(iat,ntheta),xf(ia,it),vdec(ia),d(ia,it))
      vaux=0.
      allocate(XL(0:n-1,0:ip-1,ndec), h(ihmin:ihmax), x(0:2*n-1,0:nt-1))
      allocate (XH(0:n-1,0:ip-1,ndec), HL(0:n-1,0:ip-1,ndec),HLaux2(n))
      allocate(g(igmin:igmax), HH(0:n-1,0:ip-1,ndec), h1(ihmin:ihmax))
      allocate (xdec(n), HLaux(1:n,1:ip), HHaux(1:n,1:ip), g1(igmin:igmax)) 
      allocate (HX(0:2*n-1,0:nt-1), XHaux(1:n,1:ip), dif(0:2*n-1,0:nt-1))
      allocate (T(0:n-1,0:ip-1), HXaux(1:1001,1:96), HLfilt(1:n,1:ip))
      allocate (XHfilt(1:n,1:ip),HHfilt(1:n,1:ip),vx(1:n,1:ip),vi(1:n,1:ip))
      allocate (xaux2(iat,it))
 
      h(0)=0.015404 ; h(1)=0.003491 ; h(2)=-0.117990 ; h(3)=-0.048312
      h(4)=0.491056 ; h(5)=0.787641 ; h(6)=0.337929; h(7)=-0.072638
      h(8)=-0.021060 ; h(9)=0.044725 ; h(10)=0.001768 ; h(11)=-0.007801
 
      g(0)=-h(11) ; g(1)=h(10) ; g(2)=-h(9) ; g(3)=h(8) ; g(4)=-h(7) 
      g(5)=h(6) ; g(6)=-h(5) ; g(7)=h(4) ; g(8)=-h(3) ; g(9)=h(2) 
      g(10)=-h(1) ; g(11)=h(0) 
       
   call system("rm l5090-tiro1-recon.ad rd-500.ad rd-rec-500.ad orig.ad")
       

       open(10,file='l5090-mute.ad',form='unformatted',status='unknown',recl=4*ia, access='direct')
          open(20,file='l5090-tiro1-recon.ad',form='unformatted',status='unknown',recl=4*ia, access='direct')
      open(30, file='orig.ad',form='unformatted',status='unknown',recl=4*iat, access='direct') 
      !open(40,file='HH1prenmo.ad',form='unformatted',status='unknown',recl=4*501, access='direct')
      open(50,file='HL1l5090-500.ad',form='unformatted',status='unknown',recl=4*n, access='direct')
      open(60,file='HH1l5090-500.ad',form='unformatted',status='unknown',recl=4*n, access='direct')
      open(70,file='LH1l5090-500.ad',form='unformatted',status='unknown',recl=4*n, access='direct')
      open(80,file='LL1l5090-500.ad',form='unformatted',status='unknown',recl=4*n, access='direct')
       open(90,file='rd-500.ad',form='unformatted',status='unknown',recl=4*(ia+nshift), access='direct')
       open(100,file='rd-rec-500.ad',form='unformatted',status='unknown',recl=4*ia, access='direct')
      
      x=0.0
      ig=1
      k1=(ig-1)*it 

       xaux2=0.
      
       do j=1,it
            read(10,rec=k1+j)(xaux(i,j),i=1,ia)
       end do

       write(*,*)'leitura'
    
       xaux2(nshift+1:iat,1:it)=xaux(1:ia,1:it)

       write(*,*)'shift'

      call get_rd_sec(iat,it,xaux2,x0,t0,thetai,dtheta,ntheta,vaux) ! chama a subrotina para obter a transformada direta

       write(*,*)'radial direta'

      xaux2=0.0

       do m=1, ntheta
          write(90,rec=k1+m) (vaux(i,m), i=1, ia+nshift)
       end do

       write(*,*)'gravacao radial'

      x(0:ia+nshift-1,0:ntheta-1)=vaux(1:ia+nshift,1:ntheta)
 
      
      XL=0.0 ; HL=0.0 ; XH=0.0 ; HH=0.0 ; HX=0.0
  
      j=ndec
      call F(j,n,nt,ip,ihmax,ihmin,igmax,igmin,x,g,h,XL,HL,XH,HH)
       write(*,*)'wavelet direta'

      HX=0.0
      do io=0, ip-1
      do l=0, n-1
      HX(l,io)=XL(l,io,ndec)
      end do
      end do
     
      do m=1, ip
          write(50,rec=m) (HL(i,m-1,1), i=0, 511)
      end do
      do m=1, ip
          write(60,rec=m) (HH(i,m-1,1), i=0, 511)
      end do
      do m=1, ip
          write(70,rec=m) (XH(i,m-1,1), i=0, 511)
      end do
      do m=1, ip
          write(80,rec=m) (XL(i,m-1,1), i=0, 511)
      end do
       HH=0.
        HL=0.
!        XH=0.
!        XL=0.
       write(*,*)'gravacao wavelet'


      j=ndec
      call Inv(j,n,nt,ip,g,h,igmax,igmin,ihmax,ihmin,XH,HL,HH,HX)
       
       write(*,*)'wavelet inversa'

       do m=1, ntheta
          write(100,rec=k1+m) (HX(i,m-1), i=0, ia-1)
       end do

        xaux=0.
        vf(1:ia+nshift,1:ntheta)=HX(0:ia+nshift-1,0:ntheta-1)
       
        call get_rd_inverse(iat,it,ntheta,x0,t0,thetai,dtheta,vf,xaux2)

       write(*,*)'radial inversa'

       xaux(1:ia,1:it)=xaux2(1+nshift:iat,1:it)

       write(*,*)'troca de vetor'
      
       do m=1, it
          write(20,rec=m) (xaux(i,m), i=1, ia)
       end do
       write(*,*)'gravacao do dado'
     
      close (10)
      close (20)
     close (30)
     close (40)
      close (50)
     close (60)
      close (70)
      close (80)
      close (90)
      close (100) 
      deallocate (x,XL,XH,HL,HH,g,h,HX,xdec,HHaux,XHaux,HLaux,dif,T)
      deallocate (g1,h1, HXaux,HHfilt,HLfilt,XHfilt)
      deallocate (xaux,vaux,xf,vf,vdec,d,vx,xaux2)
      end
      


      subroutine F(j,n,nt,ip,ihmax,ihmin,igmax,igmin,x,g,h,XL,HL,XH,HH)
      allocatable a(:,:), b(:,:)
      dimension XL(0:n-1,0:ip-1,j), h(ihmin:ihmax), x(0:2*n-1,0:nt-1)
      dimension  XH(0:n-1,0:ip-1,j), HL(0:n-1,0:ip-1,j)
      dimension  g(igmin:igmax), HH(0:n-1,0:ip-1,j)
      allocate (a(0:2*n-1,0:ip-1), b(0:2*n-1,0:ip-1))
      
      a=0.0 ; b=0.0
      do k=1,j
      !*****************************************************************
      !TransformaÆo Horizontal: Fixa uma linha, e varia as colunas
      do l=0, 2*n-1
      do io= 0, ip-1
         do i=0, ihmax
         ic=mod(2*io+i,2*ip)
         a(l,io)=a(l,io)+h(i)*x(l,ic)
         end do
      end do
      end do
      do l=0, 2*n-1
      do io=0, ip-1
         do i=0, igmax
         ic=mod(2*io+i,2*ip)
         b(l,io)=b(l,io)+g(i)*x(l,ic)
         end do
      end do
      end do
      !*****************************************************************
      !TransformaÆo Vertical: Fixa uma coluna, e varia as linhas
      do io=0,ip-1
      do l=0,n-1
          do i=0, ihmax
          ic=mod(2*l+i,2*n)
         XL(l,io,k)=XL(l,io,k)+h(i)*a(ic,io)
         HL(l,io,k)=HL(l,io,k)+h(i)*b(ic,io)
         end do
      end do
      end do
      do io=0, ip-1
      do l=0, n-1
         do i=0, igmax
         it=mod(2*l+i,2*n)
         XH(l,io,k)=XH(l,io,k)+g(i)*a(it,io)
         HH(l,io,k)=HH(l,io,k)+g(i)*b(it,io)
         end do
      end do
      end do
      !*****************************************************************
      !Zera o vetor x, e passa o conte£do de XL do n¡vel pra ele
      x=0.0
      do io=0,ip-1
      do l=0,n-1
      x(l,io)=XL(l,io,k)
      end do
      end do
      !*****************************************************************
      !Zera a e b para usarmos de novo
      a=0.0
      b=0.0
      end do
      !*****************************************************************
      deallocate (a, b)
      return
      end
      
      subroutine Inv(j,n,nt,ip,g,h,igmax,igmin,ihmax,ihmin,XH,HL,HH,HX)
      allocatable s(:,:), t(:,:), u(:,:)

      dimension HX(0:2*n-1,0:nt-1), HH(0:n-1,0:ip-1,j)
      dimension  XH(0:n-1,0:ip-1,j), HL(0:n-1,0:ip-1,j)
      dimension g(igmin:igmax), h(ihmin:ihmax)
      allocate (s(0:2*n-1,0:nt-1), t(0:2*n-1,0:nt-1))
      allocate (u(0:2*n-1,0:nt-1))


      !*****************************************************************
      !ReconstruÆo vertical: Fixa uma coluna, e varia as linhas
      do k= j,1,-1
          do io=0, nt-1
           do l=0, 2*n-1
            do i=0, ihmax
            iu=mod(l,2)
            iy=mod(i,2)
            if ((iy.eq.0.and.iu.eq.0).or.(iy.ne.0.and.iu.ne.0)) then
            ie=(l-i)/2
            iv= mod(ie,n)
            !if (iv.lt.0) then
            !iv=0
            !end if 
            fa=fa+h(i)*HX(iv,io)
            else
            fa=fa
            end if
           end do

            do i=0, igmax
            iu=mod(l,2)
            iy=mod(i,2)
            if ((iy.eq.0.and.iu.eq.0).or.(iy.ne.0.and.iu.ne.0)) then
            ie=(l-i)/2
            iv= mod(ie,n)
            !if (iv.lt.0) then
            !iv=0
            !end if 
            fb=fb+g(i)*XH(iv,io,k)
            else
            fb=fb
            end if
           end do
           s(l,io)=fa+fb
           fa=0
           fb=0
          end do
         end do

         do io=0, nt-1
           do l=0, 2*n-1
          do i=0, ihmax
          iu=mod(l,2)
          iy=mod(i,2)
          if ((iy.eq.0.and.iu.eq.0).or.(iy.ne.0.and.iu.ne.0)) then
          ie=(l-i)/2
          iv= mod(ie,n)
          !if (iv.lt.0) then
          !iv=0
          !end if 
          fa=fa+h(i)*HL(iv,io,k)
          else
          fa=fa
          end if
          end do

          do i=0, igmax
          iu=mod(l,2)
          iy=mod(i,2)
          if ((iy.eq.0.and.iu.eq.0).or.(iy.ne.0.and.iu.ne.0)) then
          ie=(l-i)/2
          iv= mod(ie,n)
          !if (iv.lt.0) then
          !iv=0
          !HX(iv,io)=0
          !end if 
          fb=fb+g(i)*HH(iv,io,k)
          else
          fb=fb
          end if
          end do
          t(l,io)=fa+fb
          fa=0
          fb=0
          end do
         end do
      !*****************************************************************
      !ReconstruÆo Horizontal: Fixa uma linha, e varia as  colunas
          do l=0, 2*n-1
           do io=0, nt-1
          do i=0, ihmax
          iu=mod(io,2)
          iy=mod(i,2)
          if ((iy.eq.0.and.iu.eq.0).or.(iy.ne.0.and.iu.ne.0)) then
          ie=(io-i)/2
          iv= mod(ie,ip)
          if (iv.lt.0) then
          iv=ip+iv      
          end if 
          fa=fa+h(i)*s(l,iv)
          else
          fa=fa
          end if
          end do

          do i=0, igmax
          iu=mod(io,2)
          iy=mod(i,2)
          if ((iy.eq.0.and.iu.eq.0).or.(iy.ne.0.and.iu.ne.0)) then
          ie=(io-i)/2
          iv= mod(ie,ip)
          if (iv.lt.0) then
          iv=ip+iv       
          end if 
          fb=fb+g(i)*t(l,iv)
          else
          fb=fb
          end if
          end do
          u(l,io)=fa+fb
          fa=0
          fb=0
          end do
         end do
      !*****************************************************************
      !Passa o conte£do de u para HX para usarmos novamente, e zera u,
      ! s,t para usarmos novamente.
         HX=0.0
         HX=u
         s=0.0
         t=0.0
         u=0.0
      end do
      !*****************************************************************
      deallocate(s,t,u)
      return
      end

 !=========================================================================================
!                    SUBROTINA MULTIPLICAÇÃO DE MATRIZES POR VETORES
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
!=========================================================================================
!                           SUBROTINA CÁLCULO DO PRODUTO INTERNO
!=========================================================================================
      Subroutine dot(m,a,b,c)
      real a(m),b(m),c
      c=0.
      do i=1,m
         c=c+a(i)*b(i)
      end do
      return
      end


!==================================================================================
! Rafael Rodrigues Manenti
! Subrotina utilizada para obter a transformada radial do conjunto de traços
! nt -> número de traços
! ns -> número de amostragem por traço
! x -> conjunto do dado original
! ntheta -> número de traços da transformada radial
! v -> conjunto de traços no domínio RT
! x0 -> coordenada da distância de localização do foco
! t0 -> coordenada do tempo de localização do foco
! thetai -> ângulo inicial
! dtheta -> ângulo de mapeamento do dado
!==================================================================================
      Subroutine get_rd_sec(ns,nt,x,x0,t0,thetai,dtheta,ntheta,v)
      dimension x(ns,nt),v(ns,ntheta)

      pi=acos(-1.0)
      ti=thetai*pi/180 ! ângulo inicial em radiano
      dt=(dtheta*pi/180)/(ntheta-1) ! ângulo entre os traços no domínio radial
      v=0.0
      i0=t0 
      a_nt=nt
      write(*,*)'1'
      do j=1,ntheta
         theta=ti+dt*(j-1) ; a=tan(pi*0.5-theta)
         do i=i0,ns
            x1=x0 + a*(i-i0)
            if (x1.gt.a_nt.or.x1.lt.0.0) exit
            j1=int(x1)
            alpha=x1-j1 ; beta=1.0-alpha
            c=(alpha**2)+(beta**2)
          
            v(i,j)=((beta**2)/c)*x(i,j1)+((alpha**2)/c)*x(i,j1+1)
!      write(*,*)'3'
      write(110,*)i,j,v(i,j)

         end do
!        write(*,*)j
      end do

      return
      end

!=======================================================================================
! Rafael Rodrigues Manenti
! Subrotina utilizada para obter a transformada radial inversa do conjunto de traços
!=======================================================================================

      Subroutine get_rd_inverse(ns,nt,ntheta,x0,t0,thetai,dtheta,v,x)
      dimension x(ns,nt),v(ns,ntheta)
      allocatable rtheta(:)
      allocate (rtheta(ntheta))
      x=0.
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
                  x1=x0+a*(i-t0) ! obtém o alcance do traço radial
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

