c$Id:$
      subroutine elmt07(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c....  Sara Grbcic - 3D Micropolar continuum with IM  23/11/2017
c....				 - trilinear hexaedron (8 nodes)
c....              - linear interpolation for displacements and microrotations

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    23/11/2017
c-----[--.----+----.----+----.-----------------------------------------]

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
c         ix(*)     - Global nodal connections
c         tl(*)     - Nodal temp vector
c         ndf       - Degree of freedoms/node
c         ndm       - Mesh coordinate dimension
c         nst       - Element array dimension
c         isw       - Solution option switch

c      Outputs:
c         s(nst,*)  - Element array
c         p(nst)      - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none
      
      include  'cdata.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'eqsym.h'
      include  'evdata.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'part1.h'
      include  'pmod2d.h'
      include  'strnum.h'
      include  'comblk.h'

      character text2(1)*3
      logical pcomp
	integer l,ndf,ndm,nst,isw,nhv,tdof,pdof,lint,i,j
      integer ix(*)
      real*8  d(*),ul(ndf,nen,1), xl(3,8), tl(*), s(nst,nst), p(*)
      real*8  ElastModulus,PoissonRatio
      real*8  mi,ni,lambda,alfa,beta,gamma
      real*8  dd1(9,9),dd2(9,9),sg(4,8)
      real*8  xsj,shp(4,8),dv
      real*8  Bu(9,nst),Bfi(9,nst),BuT(nst,9),BfiT(nst,9)
	real*8  IfiNfi(9,nst),IfiNfiT(nst,9)
      real*8  gpxL,gpyL,gpzL,gpx,gpy,gpz 
      real*8  shpin3d(4,3),G(9,9),GT(9,9),Greg(9,9)
      real*8  kk(nst,nst),ff(9,nst),hh(9,9),ffT(nst,9),hhin(9,9)
	real*8  p1(nst,9),p2(nst,nst),p3(nst,9),p4(nst,nst),p5(nst,nst)
	real*8	pp6(nst,nst),p7(nst,9),p8(nst,nst),p9(9,9),p10(9,nst),p11(9,nst)
      real*8  p12(9,9),p13(9,9),p14(nst,9),p15(nst,nst),p16(9,nst)
      real*8  p17(9,1),p18(9,1),p19(9,1),p21(9,1),eldof(nst,1)
      real*8  alphaIM3D(9,1),epsi(9,1),sig(9,1),mistress(9,1),p20(9,nst)
      real*8  load(6),Nu(3,nst)
      real*8  Nfi(3,nst),NuT(nst,3),NfiT(nst,3)
      real*8  geomx,geomy,geomz
      real*8  modif(3,3),omega,romega
      
      
c     Output element type

      if(isw.eq.0 .and. ior.lt.0) then
        write(*,*) 'Elmt  07: 3-D Micropolar Linear Elastic Material'
        write(*,*) 'Displacement based approach'
        write(*,*) 'Lagrange interpolation + INCOMPATIBLE MODES Hex8IM'
        
c     Input material properties

      elseif(isw.eq.1) then
		write(iow,2001)
          write(iow,'(A)') 'trilinear hexaedron (8 nodes per element)'
		if(ior.lt.0) write(*,2001)
          call dinput(d,6)
          call dinput(load,6) !BODY LOAD
          
          !FOR PROVIDAS PATCH TEST 3
          text2(1)='NO'
          
          WRITE(iow,'(A)') text2(1)   
          if(pcomp(text2(1),'YES',3)) then
          write(iow,'(A)') '      PROVIDAS TEST 3'
          else
          write(iow,'(A)') '  NOT PROVIDAS TEST 3'   
          endif
          
          write(iow,'(A)')    '   MICROPOLAR MATERIAL PARAMETERS:'
          write(iow,'(A,E10.3)') '    Mi (LAMES PARAM)        = ',d(1)
          write(iow,'(A,E10.3)') '    Ni                      = ',d(2)
          write(iow,'(A,E10.3)') '    Lambda (LAMES PARAM)    = ',d(3)
          write(iow,'(A,E10.3)') '    Alpha                   = ',d(4)
          write(iow,'(A,E10.3)') '    Beta                    = ',d(5)
          write(iow,'(A,E10.3)') '    Gamma                   = ',d(6)
          write(iow,'(A)')    '   -------------------------------------'
          ElastModulus=d(1)*(3*d(3)+2.0d0*d(1))/(d(3)+d(1))
          PoissonRatio=d(3)/(2.0d0*(d(3)+d(1)))
          
          write(iow,'(A)')    '          CORRESPONDS TO:'
          write(iow,'(A,E10.3)') '          Elasticity Modulus   =   ',
     &ElastModulus
          write(iow,'(A,E10.3)') '          Poisson ratio        =   ',
     &PoissonRatio
          
      elseif(isw.eq.3 .or. isw.eq.6) then  !Form stiffness and residual

		mi		=	d(1)
		ni		=	d(2)
		lambda	=	d(3)
		alfa	=	d(4)
		beta	=	d(5)
          gamma   =   d(6)
		
          
c		Form constitutive matrix dd1(9,9) and dd2(9,9)

          call zeromatrix(dd1,9,9)
          call zeromatrix(dd2,9,9)
          call zeromatrix(kk,nst,nst)
          call zeromatrix(ff,9,nst)
          call zeromatrix(ffT,nst,9)
          call zeromatrix(hh,9,9)
          call zeromatrix(hhin,9,9)
          call zeromatrix(alphaIM3D,9,1)
          
		dd1(1,1)	=	lambda + 2.0d0*mi
		dd1(1,5)	=	lambda
          dd1(1,9)	=	lambda
		dd1(2,2)	=	mi + ni
		dd1(2,4)	=	mi - ni
          dd1(3,3)	=	mi + ni
		dd1(3,7)	=	mi - ni
		dd1(4,2)	=	mi - ni
		dd1(4,4)	=	mi + ni
          dd1(5,1)	=	lambda
		dd1(5,5)	=	lambda + 2.0d0*mi
          dd1(5,9)	=	lambda
          dd1(6,6)	=	mi + ni
		dd1(6,8)	=	mi - ni          
		dd1(7,3)	=	mi - ni
          dd1(7,7)	=	mi + ni 
          dd1(8,6)	=	mi - ni
          dd1(8,8)	=	mi + ni 
          dd1(9,1)	=	lambda 
		dd1(9,5)	=	lambda
          dd1(9,9)	=	lambda+ 2.0d0*mi
          
		dd2(1,1)	=	gamma + 2.0d0*alfa
		dd2(1,5)	=	gamma
          dd2(1,9)	=	gamma
		dd2(2,2)	=	alfa + beta
		dd2(2,4)	=	alfa - beta
          dd2(3,3)	=	alfa + beta
		dd2(3,7)	=	alfa - beta
		dd2(4,2)	=	alfa - beta
		dd2(4,4)	=	alfa + beta
          dd2(5,1)	=	gamma
		dd2(5,5)	=	gamma + 2.0d0*alfa
          dd2(5,9)	=	gamma
          dd2(6,6)	=	alfa + beta
		dd2(6,8)	=	alfa - beta      
		dd2(7,3)	=	alfa - beta
          dd2(7,7)	=	alfa + beta
          dd2(8,6)	=	alfa - beta
          dd2(8,8)	=	alfa + beta
          dd2(9,1)	=	gamma 
		dd2(9,5)	=	gamma
          dd2(9,9)	=	gamma + 2.0d0*alfa
          
c			Compute Gauss quadrature points and weights

			l=2  !number of int points per direction

			call int3d(l,lint,sg) !sg(4,8) - array of points and weights
                                    !sg(1,i) - x coordinate of point GP_i
                                    !sg(2,i) - y coordinate of point GP_i
                                    !sg(3,i) - z coordinate of point GP_i
                                    !sg(4,i) = w_i (weight factor)
              
          do i=1,nst
           p(i) = 0.0d0
          enddo
          
c			Compute shape functions and derivatives
				
          do l=1,lint !loop Gauss points
              
              gpxL=0.0d0
              gpyL=0.0d0
              gpzL=0.0d0
                  
              gpxL=sg(1,l)
              gpyL=sg(2,l)
              gpzL=sg(3,l)
              
c		Fill matrices with zero values
                  
c...      compatible part
              call zeromatrix(Bu,9,nst)
              call zeromatrix(BuT,nst,9)
              call zeromatrix(Bfi,9,nst)
              call zeromatrix(BfiT,nst,9)
              call zeromatrix(IfiNfi,9,nst)
              call zeromatrix(IfiNfiT,nst,9)
              call zeromatrix(Nu,3,nst)
              call zeromatrix(NuT,nst,3)
              call zeromatrix(Nfi,3,nst)
              call zeromatrix(NfiT,nst,3)
              
c...      incompatible part              
                call zeromatrix(G,9,9)
                call zeromatrix(Greg,9,9)
                call zeromatrix(GT,9,9)

              call shp3d(sg(1,l),xsj,shp,xl,ndm,nel)
                      
                  dv = 0.0d0
                  dv = xsj*sg(4,l)
                  
				
			!shp(1,i) = Ni,x in Gauss point l, shp(2,i) = Ni,y in Gauss point l
			!shp(3,i) = Ni,z in Gauss point l, shp(4,i) = Ni in Gauss point l
                  !i=1..8
                  
c          ---   GAUSS POINTS GLOBAL CS   ---    
                  gpx=0.0d0
                  gpy=0.0d0
                  gpz=0.0d0
                  
               do i=1,8
                  gpx=gpx+shp(4,i)*xl(1,i)
                  gpy=gpy+shp(4,i)*xl(2,i)
                  gpz=gpz+shp(4,i)*xl(3,i)
               end do !i
               
         !     write(iow,*)			!write blank line
         !     write(iow,100) 'GAUSS POINT=',l	!A-string format, I-integer format
         !     write(iow,'(A,E15.6,A,E15.6,A,E15.6,A)') '(',gpxL,',',gpyL,','
         !& ,gpzL,')'
         !     write(iow,*)			!write blank line                   
         !     write(iow,100) 'GAUSS POINT Global =',l	!A-string format, I-integer format
         !     write(iow,'(A,E15.6,A,E15.6,A,E15.6,A)') '(',gpx,',',gpy,','
         !&,gpz,')' !F, real 3dec pl
         !     write(iow,*)			!write blank line
                
                
              !do i=1,8      
              !    write(iow,*) i
              !    write(iow,'(A,E10.3)') '	dN_i/dx=',shp(1,i)
              !    write(iow,'(A,E10.3)') '	dN_i/dy=',shp(2,i)
              !    write(iow,'(A,E10.3)') '	dN_i/dz=',shp(3,i)
              !    write(iow,'(A,E10.3)') '	N_i=',shp(4,i)
              !end do !i
      
c...      incompatible part
          
              call shapefinm3d(gpxL,gpyL,gpzL,xl,shpin3d,ndm)
              
                  !write(iow,'(A,E10.3)') '	dM_1/dx=',shpin3d(1,1)
                  !write(iow,'(A,E10.3)') '	dM_1/dy=',shpin3d(2,1)
                  !write(iow,'(A,E10.3)') '	dM_1/dz=',shpin3d(3,1)
                  !write(iow,'(A,E10.3)') '	dM_2/dx=',shpin3d(1,2)
                  !write(iow,'(A,E10.3)') '	dM_2/dy=',shpin3d(2,2)
                  !write(iow,'(A,E10.3)') '	dM_2/dz=',shpin3d(3,2)
                  !write(iow,'(A,E10.3)') '	dM_3/dx=',shpin3d(1,3)
                  !write(iow,'(A,E10.3)') '	dM_3/dy=',shpin3d(2,3)
                  !write(iow,'(A,E10.3)') '	dM_3/dz=',shpin3d(3,3)
                  !write(iow,'(A,E10.3)') '	M_1=',shpin3d(4,1)
                  !write(iow,'(A,E10.3)') '	M_2=',shpin3d(4,2)
                  !write(iow,'(A,E10.3)') '	M_3=',shpin3d(4,3)
              
                          


      
c    !!!     FORM MATRICES
      
c              Matrix Nu[3,nst] 
                  
              do i=1,8
                Nu(1,3*i-2)=shp(4,i)
                Nu(2,3*i-1)=shp(4,i)
                Nu(3,3*i)=shp(4,i)
              enddo !i
          
c              Matrix Nfi[1,nst]          
              do i=1,8
                Nfi(1,3*i-2)=shp(4,i)
                Nfi(2,3*i-1)=shp(4,i)
                Nfi(3,3*i)=shp(4,i)
              enddo !i
               
c             Matrix Bu[9,nst]              
              do i=1,8
                  Bu(1,i*6-5)=shp(1,i)
                  Bu(2,i*6-5)=shp(2,i)
                  Bu(3,i*6-5)=shp(3,i)
                  Bu(4,i*6-4)=shp(1,i)
                  Bu(5,i*6-4)=shp(2,i)
                  Bu(6,i*6-4)=shp(3,i)
                  Bu(7,i*6-3)=shp(1,i)
                  Bu(8,i*6-3)=shp(2,i)
                  Bu(9,i*6-3)=shp(3,i)
              end do !i
              
              !write (iow,*)
              !write(iow,'(A)') '-------	Matrix Bu	-------'
              !do i=1,9
              !    write (iow,*) (Bu(i,j),j=1,nst)
              !    write(iow,'(A)') '	new row'
              !end do

c             Matrix Bfi[9,nst]              
              do i=1,8
                  Bfi(1,i*6-2)=shp(1,i)
                  Bfi(2,i*6-2)=shp(2,i)
                  Bfi(3,i*6-2)=shp(3,i)
                  Bfi(4,i*6-1)=shp(1,i)
                  Bfi(5,i*6-1)=shp(2,i)
                  Bfi(6,i*6-1)=shp(3,i)
                  Bfi(7,i*6)=shp(1,i)
                  Bfi(8,i*6)=shp(2,i)
                  Bfi(9,i*6)=shp(3,i)
              end do !i
              
              !write (iow,*)
              !write(iow,'(A)') '-------	Matrix Bfi	-------'
              !do i=1,9
              !    write (iow,*) (Bfi(i,j),j=1,nst)
              !    write(iow,'(A)') '	new row'
              !end do
             
c             Matrix IfiNfi[9,nst]              
              do i=1,8
                  IfiNfi(2,i*6)=shp(4,i)
                  IfiNfi(3,i*6-1)=-shp(4,i)
                  IfiNfi(4,i*6)=-shp(4,i)
                  IfiNfi(6,i*6-2)=shp(4,i)
                  IfiNfi(7,i*6-1)=shp(4,i)
                  IfiNfi(8,i*6-2)=-shp(4,i)
              end do !i
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  Matrix IfiNfi	-------'
              !do i=1,9
              !    write (iow,*) (IfiNfi(i,j),j=1,nst)
              !    write(iow,'(A)') '	new row'
              !end do
              
              
c             Matrix G[9,9]  !!NOTE: this matrix has to be modified
              
              do i=1,3
                  Greg(1,i*3-2) = shpin3d(1,i)
                  Greg(2,i*3-2) = shpin3d(2,i)
                  Greg(3,i*3-2) = shpin3d(3,i)
                  Greg(4,i*3-1) = Greg(1,i*3-2)
                  Greg(5,i*3-1) = Greg(2,i*3-2)
                  Greg(6,i*3-1) = Greg(3,i*3-2)
                  Greg(7,i*3) = Greg(1,i*3-2)
                  Greg(8,i*3) = Greg(2,i*3-2)
                  Greg(9,i*3) = Greg(3,i*3-2)
              end do !i
              
              !write (iow,*)
              !write(iow,'(A)') '-------	Matrix G NOT modified	-------'
              !do i=1,9
              !    write (iow,*) (Greg(i,j),j=1,9)
              !    write(iow,'(A)') '	new row'
              !end do
              
           call matrixmodificationg3d(sg,xl,ndm,omega,modif)
          
          romega=1.d0/omega
          
          
          !do i=1,3
          !   write(iow,*) i
          !   write (iow,*)  'integral(1,i) po x=', modif(1,i)
          !   write (iow,*)'integral(2,i) po y= ', modif(2,i)
          !   write (iow,*) 'integral(3,i) po z = ',modif(3,i)
          !end do
              
              do i=1,3
                  G(1,i*3-2) = shpin3d(1,i)-romega*modif(1,i)
                  G(2,i*3-2) = shpin3d(2,i)-romega*modif(2,i)
                  G(3,i*3-2) = shpin3d(3,i)-romega*modif(3,i)
                  G(4,i*3-1) = G(1,i*3-2)
                  G(5,i*3-1) = G(2,i*3-2)
                  G(6,i*3-1) = G(3,i*3-2)
                  G(7,i*3) = G(1,i*3-2)
                  G(8,i*3) = G(2,i*3-2)
                  G(9,i*3) = G(3,i*3-2)
              end do !i
              
              
              !write (iow,*)
              !write(iow,'(A)') '-------	Matrix G modified	-------'
              !do i=1,9
              !    write (iow,*) (G(i,j),j=1,9)
              !    write(iow,'(A)') '	new row'
              !end do
              
          do i=1,3 !loop Nu -> NuT
			do j=1,nst
			NuT(j,i)=Nu(i,j)
			end do !j
          end do !i
          
          do i=1,3 !loop Nu -> NuT
			do j=1,nst
			NfiT(j,i)=Nfi(i,j)
			end do !j
          end do !i 

		do i=1,9 !loop Bu -> BuT
			do j=1,nst
			BuT(j,i)=Bu(i,j)
			end do !j
          end do !i

          do i=1,9 !loop Bfi -> BfiT
			do j=1,nst
			BfiT(j,i)=Bfi(i,j)
			end do !j
          end do !i
                  
		do i=1,9 !loop IfiNfi -> IfiNfiT
			do j=1,nst
			IfiNfiT(j,i)=IfiNfi(i,j)
			end do !j
          end do !i
          
          do i=1,9 !loop G -> GT
			do j=1,9
			GT(j,i)=G(i,j)
			end do !j
          end do !i
          
c  -----     MAIN MATRICES TO FORM STIFFNESS MATRIX     -----     
              
c  ....      MATRIX [K] = 48x48      ....       
              
		call matmultip(BuT,dd1,p1,nst,9,9) !BuT[nst,9].C1[9,9]=p1[nst,9]
		call matmultip(p1,Bu,p2,nst,9,nst) !p1[nst,9].Bu[9,nst]=p2[9,9]
		call matmultip(IfiNfiT,dd1,p3,nst,9,9) !IfiNfiT[nst,9].C1[9,9]=p3[nst,9]
		call matmultip(p3,Bu,p4,nst,9,nst) !p3[nst,9].Bu[9,nst]=p4[nst,nst]
		call matmultip(p1,IfiNfi,p5,nst,9,nst) !p1[nst,9].IfiNfi[9,nst]=p5[nst,nst]
		call matmultip(p3,IfiNfi,pp6,nst,9,nst) !p3[nst,9].IfiNfi[9,nst]=pp6[nst,nst]
		call matmultip(BfiT,dd2,p7,nst,9,9) !BfiT[nst,9].C2[9,9]=p7[nst,9]
		call matmultip(p7,Bfi,p8,nst,9,nst) !p7[nst,9].Bfi[9,nst]=p8[nst,nst]

		do i=1,nst
			do j=1,nst
		kk(i,j)=kk(i,j)+(p2(i,j)+p4(i,j)+p5(i,j)+pp6(i,j)+p8(i,j))
     &		      *dv
			end do !j
          end do !i
          
          !write (iow,*)
          !    write(iow,'(A)') '	Volume load'
          !    write(iow,'(A,i5)') 'ELEMENT = ',n
          !    do i=1,nst
          !        write (iow,*) (p(i))
          !        write (iow,*)
          !        write(iow,'(A)') '	new row'
          !    end do
          
          
          
c  ....      MATRIX [F] = 9x48      ....
              
              call matmultip(GT,dd1,p9,9,9,9)        !GT[9,9].C1[9,9]=p9[9,9]
              call matmultip(p9,Bu,p10,9,9,nst)      !p9[9,9].Bu[9,nst]=p10[9,nst]
              call matmultip(p9,IfiNfi,p11,9,9,nst)  !p9[9,9].IfiNfi[9,nst]=p11[9,nst]
              
              
              do i=1,9
			    do j=1,nst
                  ff(i,j)=ff(i,j)+(p10(i,j)+p11(i,j))*dv
			    end do !j
              end do !i
              
c  ....      MATRIX [H] = 9 x 9      ....         
            
              call matmultip(GT,dd1,p12,9,9,9)        !GT[9,9].C1[9,9]=p12[9,9]
              call matmultip(p12,G,p13,9,9,9)         !p12[9,9].G[9,9]=p13[9,9]
              
              do i=1,9
			    do j=1,9
                  hh(i,j)=hh(i,j)+p13(i,j)*dv
			    end do !j
              end do !i
              
c ....        RESIDUAL
              
          geomx = 0.0d0
          geomy = 0.0d0
          geomz = 0.0d0
          
      do i=1,nen
          geomx = geomx + shp(4,i)*xl(1,i)
          geomy = geomy + shp(4,i)*xl(2,i)
          geomz = geomz + shp(4,i)*xl(3,i)
      enddo
          
          if(pcomp(text2(1),'YES',3)) then
              do i=1,nen
          p(6*i-5)=p(6*i-5)+NuT(3*i-2,1)*load(1)*dv
          p(6*i-4)=p(6*i-4)+NuT(3*i-1,2)*load(2)*dv
          p(6*i-3)=p(6*i-3)+NuT(3*i,3)*load(3)*dv
          p(6*i-2)=p(6*i-2)+NfiT(3*i-2,1)*2.0d0*(geomx-geomy-geomz)*dv
          p(6*i-1)=p(6*i-1)+NfiT(3*i-1,2)*2.0d0*(geomx-geomy-geomz)*dv
          p(6*i)=p(6*i)+NfiT(3*i,3)*2.0d0*(geomx-geomy-geomz)*dv
              enddo !i
          else
              do i=1,nen
              p(6*i-5)=p(6*i-5)+NuT(3*i-2,1)*load(1)*dv
              p(6*i-4)=p(6*i-4)+NuT(3*i-1,2)*load(2)*dv
              p(6*i-3)=p(6*i-3)+NuT(3*i,3)*load(3)*dv
              p(6*i-2)=p(6*i-2)+NfiT(3*i-2,1)*load(4)*dv
              p(6*i-1)=p(6*i-1)+NfiT(3*i-1,2)*load(5)*dv
              p(6*i)=p(6*i)+NfiT(3*i,3)*load(6)*dv
              enddo !i
          endif
              
          end do !l   END LOOP GP
          
c ....        FORM THE CONDENSED STIFNESS MATRIX
          
              !write (iow,*)
              !write(iow,'(A)') '	Stiffness Matrix'
              !write(iow,'(A,i5)') 'ELEMENT = ',n
              !do i=1,nst
              !    write (iow,*) (kk(i,j),j=1,nst)
              !    write (iow,*)
              !    write(iow,'(A)') '	new row'
              !end do
          
              do i=1,9 !loop F[9,nst] -> FT[nst,9]
			    do j=1,nst
			    ffT(j,i)=ff(i,j)
			    end do !j
              end do !i
              
              call invert(hh,9,9)
              
              do i=1,9
			    do j=1,9
                  hhin(i,j)=hh(i,j)
			    end do !j
              end do !i
              
c  ....      MATRIX [S] = [K]-[FT].[Hin].[F]  48x48      ....         
              
              call matmultip(ffT,hhin,p14,nst,9,9)     !FT[48,9].Hin[9,9]=p14[48,9]
              call matmultip(p14,ff,p15,nst,9,nst)     !p14[48,9].F[9,48]=p15[48,48]
              
              do i=1,nst
			    do j=1,nst
                  s(i,j)=kk(i,j)-p15(i,j)
			    end do !j
              end do !i
              
      elseif(isw.eq.4) then  !Output stresses and strains at GP
          
          
          mi		=	d(1)
		ni		=	d(2)
		lambda	=	d(3)
		alfa	=	d(4)
		beta	=	d(5)
          gamma   =   d(6)
		
          
c		Form constitutive matrix dd1(9,9) and dd2(9,9)

          call zeromatrix(dd1,9,9)
          call zeromatrix(dd2,9,9)
          call zeromatrix(kk,nst,nst)
          call zeromatrix(ff,9,nst)
          call zeromatrix(hh,9,9)
          call zeromatrix(ffT,nst,9)
          call zeromatrix(hhin,9,9)
          call zeromatrix(eldof,nst,1)
          call zeromatrix(alphaIM3D,9,1)
          
		dd1(1,1)	=	lambda + 2.0d0*mi
		dd1(1,5)	=	lambda
          dd1(1,9)	=	lambda
		dd1(2,2)	=	mi + ni
		dd1(2,4)	=	mi - ni
          dd1(3,3)	=	mi + ni
		dd1(3,7)	=	mi - ni
		dd1(4,2)	=	mi - ni
		dd1(4,4)	=	mi + ni
          dd1(5,1)	=	lambda
		dd1(5,5)	=	lambda + 2.0d0*mi
          dd1(5,9)	=	lambda
          dd1(6,6)	=	mi + ni
		dd1(6,8)	=	mi - ni          
		dd1(7,3)	=	mi - ni
          dd1(7,7)	=	mi + ni 
          dd1(8,6)	=	mi - ni
          dd1(8,8)	=	mi + ni 
          dd1(9,1)	=	lambda 
		dd1(9,5)	=	lambda
          dd1(9,9)	=	lambda+ 2.0d0*mi
          
		dd2(1,1)	=	gamma + 2.0d0*alfa
		dd2(1,5)	=	gamma
          dd2(1,9)	=	gamma
		dd2(2,2)	=	alfa + beta
		dd2(2,4)	=	alfa - beta
          dd2(3,3)	=	alfa + beta
		dd2(3,7)	=	alfa - beta
		dd2(4,2)	=	alfa - beta
		dd2(4,4)	=	alfa + beta
          dd2(5,1)	=	gamma
		dd2(5,5)	=	gamma + 2.0d0*alfa
          dd2(5,9)	=	gamma
          dd2(6,6)	=	alfa + beta
		dd2(6,8)	=	alfa - beta      
		dd2(7,3)	=	alfa - beta
          dd2(7,7)	=	alfa + beta
          dd2(8,6)	=	alfa - beta
          dd2(8,8)	=	alfa + beta
          dd2(9,1)	=	gamma 
		dd2(9,5)	=	gamma
          dd2(9,9)	=	gamma + 2.0d0*alfa
          
c			Compute Gauss quadrature points and weights

			l=2  !number of int points per direction

			call int3d(l,lint,sg) !sg(4,8) - array of points and weights
                                    !sg(1,i) - x coordinate of point GP_i
                                    !sg(2,i) - y coordinate of point GP_i
                                    !sg(3,i) - z coordinate of point GP_i
                                    !sg(4,i) = w_i (weight factor)
              
c			Compute shape functions and derivatives 
				
          do l=1,lint !loop Gauss points
              
              gpxL=0.0d0
              gpyL=0.0d0
              gpzL=0.0d0
                  
              gpxL=sg(1,l)
              gpyL=sg(2,l)
              gpzL=sg(3,l)
              
c		Fill matrices with zero values
                  
c...      compatible part
              call zeromatrix(Bu,9,nst)
              call zeromatrix(BuT,nst,9)
              call zeromatrix(Bfi,9,nst)
              call zeromatrix(BfiT,nst,9)
              call zeromatrix(IfiNfi,9,nst)
              call zeromatrix(IfiNfiT,nst,9)
              call zeromatrix(alphaIM3D,9,1) !PAZI

              call shp3d(sg(1,l),xsj,shp,xl,ndm,nel)
                      
                  dv = 0.0d0
                  dv = xsj*sg(4,l)
                  
				
			!shp(1,i) = Ni,x in Gauss point l, shp(2,i) = Ni,y in Gauss point l
			!shp(3,i) = Ni,z in Gauss point l, shp(4,i) = Ni in Gauss point l
                  !i=1..8
                  
c          ---   GAUSS POINTS GLOBAL CS   ---    
                  gpx=0.0d0
                  gpy=0.0d0
                  gpz=0.0d0
                  
               do i=1,8
                  gpx=gpx+shp(4,i)*xl(1,i)
                  gpy=gpy+shp(4,i)*xl(2,i)
                  gpz=gpz+shp(4,i)*xl(3,i)
               end do !i
               
         ! write(iow,*)			!write blank line
         ! write(iow,100) 'GAUSS POINT=',l	!A-string format, I-integer format
         ! write(iow,'(A,E15.6,A,E15.6,A,E15.6,A)') '(',gpxL,',',gpyL,','
         !& ,gpzL,')'
         ! write(iow,*)			!write blank line                   
         ! write(iow,100) 'GAUSS POINT Global =',l	!A-string format, I-integer format
         ! write(iow,'(A,E15.6,A,E15.6,A,E15.6,A)') '(',gpx,',',gpy,',',gpz,
         !&     ')' !F, real 3dec pl
         ! write(iow,*)			!write blank line
                
                
              !do i=1,8      
              !    write(iow,*) i
              !    write(iow,'(A,E10.3)') '	dN_i/dx=',shp(1,i)
              !    write(iow,'(A,E10.3)') '	dN_i/dy=',shp(2,i)
              !    write(iow,'(A,E10.3)') '	dN_i/dz=',shp(3,i)
              !    write(iow,'(A,E10.3)') '	N_i=',shp(4,i)
              !end do !i
      
c...      incompatible part
          
              call shapefinm3d(gpxL,gpyL,gpzL,xl,shpin3d,ndm)
              
                  !write(iow,'(A,E10.3)') '	dM_1/dx=',shpin3d(1,1)
                  !write(iow,'(A,E10.3)') '	dM_1/dy=',shpin3d(2,1)
                  !write(iow,'(A,E10.3)') '	dM_1/dz=',shpin3d(3,1)
                  !write(iow,'(A,E10.3)') '	dM_2/dx=',shpin3d(1,2)
                  !write(iow,'(A,E10.3)') '	dM_2/dy=',shpin3d(2,2)
                  !write(iow,'(A,E10.3)') '	dM_2/dz=',shpin3d(3,2)
                  !write(iow,'(A,E10.3)') '	dM_3/dx=',shpin3d(1,3)
                  !write(iow,'(A,E10.3)') '	dM_3/dy=',shpin3d(2,3)
                  !write(iow,'(A,E10.3)') '	dM_3/dz=',shpin3d(3,3)
                  !write(iow,'(A,E10.3)') '	M_1=',shpin3d(4,1)
                  !write(iow,'(A,E10.3)') '	M_2=',shpin3d(4,2)
                  !write(iow,'(A,E10.3)') '	M_3=',shpin3d(4,3)

              
              
c     incompatible part
              
                call zeromatrix(G,9,9)
                call zeromatrix(GT,9,9)
      
c    !!!     FORM MATRICES
               
c             Matrix Bu[9,nst]              
              do i=1,8
                  Bu(1,i*6-5)=shp(1,i)
                  Bu(2,i*6-5)=shp(2,i)
                  Bu(3,i*6-5)=shp(3,i)
                  Bu(4,i*6-4)=shp(1,i)
                  Bu(5,i*6-4)=shp(2,i)
                  Bu(6,i*6-4)=shp(3,i)
                  Bu(7,i*6-3)=shp(1,i)
                  Bu(8,i*6-3)=shp(2,i)
                  Bu(9,i*6-3)=shp(3,i)
              end do !i
              
              !write (iow,*)
              !write(iow,'(A)') '-------	Matrix Bu	-------'
              !do i=1,9
              !    write (iow,*) (Bu(i,j),j=1,nst)
              !    write(iow,'(A)') '	new row'
              !end do

c             Matrix Bfi[9,nst]              
              do i=1,8
                  Bfi(1,i*6-2)=shp(1,i)
                  Bfi(2,i*6-2)=shp(2,i)
                  Bfi(3,i*6-2)=shp(3,i)
                  Bfi(4,i*6-1)=shp(1,i)
                  Bfi(5,i*6-1)=shp(2,i)
                  Bfi(6,i*6-1)=shp(3,i)
                  Bfi(7,i*6)=shp(1,i)
                  Bfi(8,i*6)=shp(2,i)
                  Bfi(9,i*6)=shp(3,i)
              end do !i
              
              !write (iow,*)
              !write(iow,'(A)') '-------	Matrix Bfi	-------'
              !do i=1,9
              !    write (iow,*) (Bfi(i,j),j=1,nst)
              !    write(iow,'(A)') '	new row'
              !end do
             
c             Matrix IfiNfi[9,nst]              
              do i=1,8
                  IfiNfi(2,i*6)=shp(4,i)
                  IfiNfi(3,i*6-1)=-shp(4,i)
                  IfiNfi(4,i*6)=-shp(4,i)
                  IfiNfi(6,i*6-2)=shp(4,i)
                  IfiNfi(7,i*6-1)=shp(4,i)
                  IfiNfi(8,i*6-2)=-shp(4,i)
              end do !i
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  Matrix IfiNfi	-------'
              !do i=1,9
              !    write (iow,*) (IfiNfi(i,j),j=1,nst)
              !    write(iow,'(A)') '	new row'
              !end do
              
              
c             Matrix G[9,9]  !!NOTE: this matrix has to be modified 
              
           call matrixmodificationg3d(sg,xl,ndm,omega,modif)
          
          romega=1.d0/omega
              
              do i=1,3
                  G(1,i*3-2) = shpin3d(1,i)-romega*modif(1,i)
                  G(2,i*3-2) = shpin3d(2,i)-romega*modif(2,i)
                  G(3,i*3-2) = shpin3d(3,i)-romega*modif(3,i)
                  G(4,i*3-1) = G(1,i*3-2)
                  G(5,i*3-1) = G(2,i*3-2)
                  G(6,i*3-1) = G(3,i*3-2)
                  G(7,i*3) = G(1,i*3-2)
                  G(8,i*3) = G(2,i*3-2)
                  G(9,i*3) = G(3,i*3-2)
              end do !i
              
		do i=1,9 !loop Bu -> BuT
			do j=1,nst
			BuT(j,i)=Bu(i,j)
			end do !j
          end do !i

          do i=1,9 !loop Bfi -> BfiT
			do j=1,nst
			BfiT(j,i)=Bfi(i,j)
			end do !j
          end do !i
                  
		do i=1,9 !loop IfiNfi -> IfiNfiT
			do j=1,nst
			IfiNfiT(j,i)=IfiNfi(i,j)
			end do !j
          end do !i
          
          do i=1,9 !loop G -> GT
			do j=1,9
			GT(j,i)=G(i,j)
			end do !j
          end do !i
          
c  -----     MAIN MATRICES TO FORM STIFNESS MATRIX     -----     
              
c  ....      MATRIX [K] = 48x48      ....       
              
		call matmultip(BuT,dd1,p1,nst,9,9) !BuT[nst,9].C1[9,9]=p1[nst,9]
		call matmultip(p1,Bu,p2,nst,9,nst) !p1[nst,9].Bu[9,nst]=p2[9,9]
		call matmultip(IfiNfiT,dd1,p3,nst,9,9) !IfiNfiT[nst,9].C1[9,9]=p3[nst,9]
		call matmultip(p3,Bu,p4,nst,9,nst) !p3[nst,9].Bu[9,nst]=p4[nst,nst]
		call matmultip(p1,IfiNfi,p5,nst,9,nst) !p1[nst,9].IfiNfi[9,nst]=p5[nst,nst]
		call matmultip(p3,IfiNfi,pp6,nst,9,nst) !p3[nst,9].IfiNfi[9,nst]=pp6[nst,nst]
		call matmultip(BfiT,dd2,p7,nst,9,9) !BfiT[nst,9].C2[9,9]=p7[nst,9]
		call matmultip(p7,Bfi,p8,nst,9,nst) !p7[nst,9].Bfi[9,nst]=p8[nst,nst]

		do i=1,nst
			do j=1,nst
		kk(i,j)=kk(i,j)+(p2(i,j)+p4(i,j)+p5(i,j)+pp6(i,j)+p8(i,j))
     &		      *dv
			end do !j
          end do !i
          
          
          
c  ....      MATRIX [F] = 9x48      ....
              
              call matmultip(GT,dd1,p9,9,9,9)        !GT[9,9].C1[9,9]=p9[9,9]
              call matmultip(p9,Bu,p10,9,9,nst)      !p9[9,9].Bu[9,nst]=p10[9,nst]
              call matmultip(p9,IfiNfi,p11,9,9,nst)  !p9[9,9].IfiNfi[9,nst]=p11[9,nst]
              
              
              do i=1,9
			    do j=1,nst
                  ff(i,j)=ff(i,j)+(p10(i,j)+p11(i,j))*dv
			    end do !j
              end do !i
              
c  ....      MATRIX [H] = 9 x 9      ....         
            
              call matmultip(GT,dd1,p12,9,9,9)        !GT[9,9].C1[9,9]=p12[9,9]
              call matmultip(p12,G,p13,9,9,9)         !p12[9,9].G[9,9]=p13[9,9]
              
              do i=1,9
			    do j=1,9
                  hh(i,j)=hh(i,j)+p13(i,j)*dv
			    end do !j
              end do !i              
              
          end do !l   END LOOP GP
          
              do i=1,9 !loop F[9,nst] -> FT[nst,9]
			    do j=1,nst
			    ffT(j,i)=ff(i,j)
			    end do !j
              end do !i
              
              call invert(hh,9,9)
              
              do i=1,9
			    do j=1,9
                  hhin(i,j)=hh(i,j)
			    end do !j
              end do !i
              
          write(iow,*)
		write(iow,'(A,i4,A)') '- - - - - - - - - E L E M E N T = ',n,'- - - ' 
		write(iow,*)				!write blank line
          
c     put element DOFs in a column vector
          
      do i=1,8
          eldof(6*i-5,1)=ul(1,i,1)
          eldof(6*i-4,1)=ul(2,i,1)
          eldof(6*i-3,1)=ul(3,i,1)
          eldof(6*i-2,1)=ul(4,i,1)
          eldof(6*i-1,1)=ul(5,i,1)
          eldof(6*i,1) = ul(6,i,1)
      end do !i
      
c  ....      alphaIM3D = -[Hin].[F].eldis       ....     
                  
              call matmultip(hhin,ff,p16,9,9,nst)      !Hin[9,9].F[9,nst]=p16[9,nst]
              call matmultip(p16,eldof,p17,9,nst,1)    !p16[9,nst].eldof[nst,1]=p17[9,1]
              
              do i=1,9
                  alphaIM3D(i,1)=-p17(i,1)
              end do !i
              
c  ....      strains = [Bu].eldof+IfiNfi.eldof+[G].alphaIM3D       ....               
          
           write(iow,*) ' < < < < < < <    STRESSES AND STRAINS    > > >
     & > > > >'
      
          do l=1,lint !loop Gauss points
              
              gpxL=0.0d0
              gpyL=0.0d0
              gpzL=0.0d0
                  
              gpxL=sg(1,l)
              gpyL=sg(2,l)
              gpzL=sg(3,l)
              
c		Fill matrices with zero values
                  
c...      compatible part
              call zeromatrix(Bu,9,nst)
              call zeromatrix(Bfi,9,nst)
              call zeromatrix(IfiNfi,9,nst)

              call shp3d(sg(1,l),xsj,shp,xl,ndm,nel)
                      
                  dv = 0.0d0
                  dv = xsj*sg(4,l)
                  
			!shp(1,i) = Ni,x in Gauss point l, shp(2,i) = Ni,y in Gauss point l
			!shp(3,i) = Ni,z in Gauss point l, shp(4,i) = Ni in Gauss point l
                  !i=1..8
                  
c          ---   GAUSS POINTS GLOBAL CS   ---    
                  gpx=0.0d0
                  gpy=0.0d0
                  gpz=0.0d0
                  
               do i=1,8
                  gpx=gpx+shp(4,i)*xl(1,i)
                  gpy=gpy+shp(4,i)*xl(2,i)
                  gpz=gpz+shp(4,i)*xl(3,i)
               end do !i
               
          write(iow,*)			!write blank line
          write(iow,100) 'GAUSS POINT=',l	!A-string format, I-integer format
          write(iow,'(A,E15.6,A,E15.6,A,E15.6,A)') '(',gpxL,',',gpyL,','
     & ,gpzL,')'
          write(iow,*)			!write blank line                   
          write(iow,100) 'GAUSS POINT Global =',l	!A-string format, I-integer format
          write(iow,'(A,E15.6,A,E15.6,A,E15.6,A)') '(',gpx,',',gpy,','
     & ,gpz,')' 
          write(iow,*)			!write blank line
                
                
              !do i=1,8      
              !    write(iow,*) i
              !    write(iow,'(A,E10.3)') '	dN_i/dx=',shp(1,i)
              !    write(iow,'(A,E10.3)') '	dN_i/dy=',shp(2,i)
              !    write(iow,'(A,E10.3)') '	dN_i/dz=',shp(3,i)
              !    write(iow,'(A,E10.3)') '	N_i=',shp(4,i)
              !end do !i
      
c...      incompatible part
          
              call shapefinm3d(gpxL,gpyL,gpzL,xl,shpin3d,ndm)
              
                  !write(iow,'(A,E10.3)') '	dM_1/dx=',shpin3d(1,1)
                  !write(iow,'(A,E10.3)') '	dM_1/dy=',shpin3d(2,1)
                  !write(iow,'(A,E10.3)') '	dM_1/dz=',shpin3d(3,1)
                  !write(iow,'(A,E10.3)') '	dM_2/dx=',shpin3d(1,2)
                  !write(iow,'(A,E10.3)') '	dM_2/dy=',shpin3d(2,2)
                  !write(iow,'(A,E10.3)') '	dM_2/dz=',shpin3d(3,2)
                  !write(iow,'(A,E10.3)') '	dM_3/dx=',shpin3d(1,3)
                  !write(iow,'(A,E10.3)') '	dM_3/dy=',shpin3d(2,3)
                  !write(iow,'(A,E10.3)') '	dM_3/dz=',shpin3d(3,3)
                  !write(iow,'(A,E10.3)') '	M_1=',shpin3d(4,1)
                  !write(iow,'(A,E10.3)') '	M_2=',shpin3d(4,2)
                  !write(iow,'(A,E10.3)') '	M_3=',shpin3d(4,3)
              

                call zeromatrix(G,9,9)
      
c    !!!     FORM MATRICES
               
c             Matrix Bu[9,nst]              
              do i=1,8
                  Bu(1,i*6-5)=shp(1,i)
                  Bu(2,i*6-5)=shp(2,i)
                  Bu(3,i*6-5)=shp(3,i)
                  Bu(4,i*6-4)=shp(1,i)
                  Bu(5,i*6-4)=shp(2,i)
                  Bu(6,i*6-4)=shp(3,i)
                  Bu(7,i*6-3)=shp(1,i)
                  Bu(8,i*6-3)=shp(2,i)
                  Bu(9,i*6-3)=shp(3,i)
              end do !i
              
              !write (iow,*)
              !write(iow,'(A)') '-------	Matrix Bu	-------'
              !do i=1,9
              !    write (iow,*) (Bu(i,j),j=1,nst)
              !    write(iow,'(A)') '	new row'
              !end do

c             Matrix Bfi[9,nst]              
              do i=1,8
                  Bfi(1,i*6-2)=shp(1,i)
                  Bfi(2,i*6-2)=shp(2,i)
                  Bfi(3,i*6-2)=shp(3,i)
                  Bfi(4,i*6-1)=shp(1,i)
                  Bfi(5,i*6-1)=shp(2,i)
                  Bfi(6,i*6-1)=shp(3,i)
                  Bfi(7,i*6)=shp(1,i)
                  Bfi(8,i*6)=shp(2,i)
                  Bfi(9,i*6)=shp(3,i)
              end do !i
              
              !write (iow,*)
              !write(iow,'(A)') '-------	Matrix Bfi	-------'
              !do i=1,9
              !    write (iow,*) (Bfi(i,j),j=1,nst)
              !    write(iow,'(A)') '	new row'
              !end do
             
c             Matrix IfiNfi[9,nst]              
              do i=1,8
                  IfiNfi(2,i*6)=shp(4,i)
                  IfiNfi(3,i*6-1)=-shp(4,i)
                  IfiNfi(4,i*6)=-shp(4,i)
                  IfiNfi(6,i*6-2)=shp(4,i)
                  IfiNfi(7,i*6-1)=shp(4,i)
                  IfiNfi(8,i*6-2)=-shp(4,i)
              end do !i
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  Matrix IfiNfi	-------'
              !do i=1,9
              !    write (iow,*) (IfiNfi(i,j),j=1,nst)
              !    write(iow,'(A)') '	new row'
              !end do
              
              
c             Matrix G[9,9]              
              
              do i=1,3
                  G(1,i*3-2) = shpin3d(1,i)
                  G(2,i*3-2) = shpin3d(2,i)
                  G(3,i*3-2) = shpin3d(3,i)
                  G(4,i*3-1) = shpin3d(1,i)
                  G(5,i*3-1) = shpin3d(2,i)
                  G(6,i*3-1) = shpin3d(3,i)
                  G(7,i*3) = shpin3d(1,i)
                  G(8,i*3) = shpin3d(2,i)
                  G(9,i*3) = shpin3d(3,i)
              end do !i
              
          call matmultip(Bu,eldof,p18,9,nst,1)      !Bu[9,nst].eldis[nst,1]=p18[9,1]
          call matmultip(G,alphaIM3D,p19,9,9,1)     !G[9,9].alphaIM3D[9,1]=p19[9,1]
          call matmultip(IfiNfi,eldof,p21,9,nst,1)  !IfiNfi[9,nst].eldis[nst,1]=p21[9,1]
          
          call zeromatrix(epsi,9,1)
          
          do i=1,9
               epsi(i,1)=p18(i,1)+p19(i,1)+p21(i,1)
          end do
          
          write(iow,*)
          write(iow,100) '    STRAINS AT GAUSS POINT',l
		write(iow,*) 
		write(iow,'(A,e15.6)') '				EPSILON (1,1)=',epsi(1,1)
		write(iow,'(A,e15.6)') '				EPSILON (1,2)=',epsi(2,1)
          write(iow,'(A,e15.6)') '				EPSILON (1,3)=',epsi(3,1)
          write(iow,'(A,e15.6)') '				EPSILON (2,1)=',epsi(4,1)
		write(iow,'(A,e15.6)') '				EPSILON (2,2)=',epsi(5,1)
          write(iow,'(A,e15.6)') '				EPSILON (2,3)=',epsi(6,1)
          write(iow,'(A,e15.6)') '				EPSILON (3,1)=',epsi(7,1)
		write(iow,'(A,e15.6)') '				EPSILON (3,2)=',epsi(8,1)
          write(iow,'(A,e15.6)') '				EPSILON (3,3)=',epsi(9,1)
		write(iow,*)
          
          call zeromatrix(sig,9,1)
          
          call matmultip(dd1,epsi,sig,9,9,1)  !C1[9,9].epsi[9,1]=sig[9,1]
          
          call zeromatrix(mistress,9,1)
          
          call matmultip(dd2,Bfi,p20,9,9,nst) !C2[9,9].Bfi[9,nst]=p20[9,nst]
          call matmultip(p20,eldof,mistress,9,nst,1) !p20[9,nst].eldof[nst,1]=mistress[9,1]
						
          write(iow,100) '    STRESSES AT GAUSS POINT',l
          write(iow,*)			!write blank line
          write(iow,*) '			NORMAL STRESSES SIGMA'
		write(iow,*) 
		write(iow,'(A,e15.6)') '				SIGMA (1,1)=',sig(1,1)
		write(iow,'(A,e15.6)') '				SIGMA (1,2)=',sig(2,1)
          write(iow,'(A,e15.6)') '				SIGMA (1,3)=',sig(3,1)
          write(iow,'(A,e15.6)') '				SIGMA (2,1)=',sig(4,1)
		write(iow,'(A,e15.6)') '				SIGMA (2,2)=',sig(5,1)
          write(iow,'(A,e15.6)') '				SIGMA (2,3)=',sig(6,1)
          write(iow,'(A,e15.6)') '				SIGMA (3,1)=',sig(7,1)
		write(iow,'(A,e15.6)') '				SIGMA (3,2)=',sig(8,1)
          write(iow,'(A,e15.6)') '				SIGMA (3,3)=',sig(9,1)
		write(iow,*)
		write(iow,'(A,e15.6)') '			COUPLE STRESSES MI'
		write(iow,*)
		write(iow,'(A,e15.6)') '				MI (1,1)=', mistress(1,1)
		write(iow,'(A,e15.6)') '				MI (1,2)=', mistress(2,1)
          write(iow,'(A,e15.6)') '				MI (1,3)=', mistress(3,1)
		write(iow,'(A,e15.6)') '				MI (2,1)=', mistress(4,1)
		write(iow,'(A,e15.6)') '				MI (2,2)=', mistress(5,1)
          write(iow,'(A,e15.6)') '				MI (2,3)=', mistress(6,1)
          write(iow,'(A,e15.6)') '				MI (3,1)=', mistress(7,1)
		write(iow,'(A,e15.6)') '				MI (3,2)=', mistress(8,1)
          write(iow,'(A,e15.6)') '				MI (3,3)=', mistress(9,1)
          
          end do !l  
      
	endif

100	format(A,i5)
2001  format(
     & /5x,'T h r e e   D i m e n s i o n a l   M i c r o p o l a r   
     &  S o l i d   E l e m e n t'/)
      
      end
