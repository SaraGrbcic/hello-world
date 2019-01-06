c$Id:$
      subroutine elmt08(d,ul,xl,ix,tl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c....  Sara Grbcic - 3D GEOMETRICALLY NONLINEAR MICROPOLAR CONTINUUM
c....				 - trilinear hexaedron (8 nodes)
c....              - linear interpolation for displacements and microrotations

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    17/05/2018
c-----[--.----+----.----+----.-----------------------------------------]

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
c         ix(*)     - Global nodal coj1ections
c         tl(*)     - Nodal temp vector
c         ndf       - Degree of freedoms/node
c         ndm       - Mesh coordinate dimension
c         nst       - Element array dimension
c         isw       - Solution option switch

c      Outputs:
c         s(nst,*)  - Element array
c         r(ndf,nen)  - Residual Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none
      
      include  'cdata.h'  !numnp,numel,nummat,nen,neq,ipr
      include  'eldata.h' !dm,n,ma,mct,iel,nel,pstyp,eltyp,elty2,elty3
      include  'eltran.h' !bpr,ctan,psil
      include  'eqsym.h' !neqs
      include  'evdata.h'!imtyp,mf,mq,evtyp 
      include  'hdata.h'  !history variables
      include  'iofile.h' !ior,iow,ilg
      include  'part0.h' !npart,ndfp,ndfo,ndfg,ndfl
      include  'part1.h'  !nqp,nqs,nqr,nittyp
      include  'pmod2d.h'
      include  'strnum.h'
      include  'comblk.h'
      include  'counts.h' !count variables niter ect.
      include  'eldatp.h'
      include  'fdata.h'

      character text2(1)*3,text3(1)*3
      logical pcomp
	integer l,ndf,ndm,nst,isw,nhv,tdof,pdof,lint,i,j,m,j1
      integer ix(*),nn1,ngp,fn
      
      real*8  d(*),load(6),xl(3,nel),tl(*), s(nst,nst)
      real*8  r(ndf,nen),mi,ni,lambda,alfa,beta,gamma
      real*8  ElastModulus,PoissonRatio
      real*8  Bstress(3,3),Gstress(3,3),Estrain(3,3),Knew(3,3)
      real*8  Iden(3,3),F(3,3),FT(3,3),Q(3,3),QT(3,3),prod1(3,3)
      real*8  gpxL,gpyL,gpzL,gpx,gpy,gpz
      real*8  xsj,shp(4,nen),dv,sg(4,8),ul(ndf,nen,*),DeltaFi(3) !sg(4,ngp)
      real*8  quaDeltaFi(4),quaNEW(4),quaOLD(4)
      real*8  QB(3,3),Ninabla(3),vect1(3),skewvect1(3,3),term1(3,3)
      real*8  FQT(3,3),skewFQT(3,3),axskewFQT(3),vect2(3),matterm2(3,3)
      real*8  term2(3,3),Njnabla(3),term3(3,3),vect3(3),skewvect3(3,3)
      real*8  term4(3,3),vect4(3),matterm4(3,3),QG(3,3),vect5(3)
      real*8  skewvect5(3,3),term5(3,3),matterm6(3,3),term6(3,3)
      real*8  BstressT(3,3),matterm7(3,3),matterm7T(3,3),trmatterm7
      real*8  term7(3,3),nulmatrix(3,3),blokKg(6,6)
      real*8  Kg(nst,nst),matterm8a(3,3),matterm8b(3,3),term8(3,3)
      real*8  scnabla,term9(3,3),blokKgpom(6,6,nen,nen)
      real*8  term10(3,3),matterm10a(3,3),matterm10b(3,3)
      real*8  term11(3,3),vect11(3),skewvect11(3,3)
      real*8  term12(3,3),vect12(3),skewvect12(3,3),matterm12(3,3)
      real*8  term13(3,3),vect13(3),skewvect13(3,3)
      real*8  term14(3,3),vect14(3),skewvect14(3,3),matterm14(3,3)
      real*8  term15(3,3)
      real*8  term16(3,3)
      real*8  term17(3,3)
      real*8  term18(3,3),matterm18a(3,3),trmatterm18a
      real*8  term19(3,3),Ricci1(3,3),Ricci2(3,3),Ricci3(3,3)
      real*8  matterm19a(3,3),matterm19b1(3,3),matterm19b2(3,3)
      real*8  matterm19b3(3,3),dskewmatterm19b1(3,3)
      real*8  dskewmatterm19b2(3,3),dskewmatterm19b3(3,3)
      real*8  axmatterm19b1(3),axmatterm19b2(3),axmatterm19b3(3)
      real*8  matterm19(3,3),blokKm(6,6),Km(nst,nst),res1(3)
      real*8  res2(3),res2a(3,3),skewres2a(3,3),axres2a(3)
      real*8  res3(3),blr1(6),blr2(6)
      real*8  newvect(3),newvect2(3),newmatrix(3,3)
      real*8  Kold(3,3),skewDeltaFi(3,3),normDeltaFi,const1,const2
      real*8  skewDeltaFiSQ(3,3),Hmat(3,3),K1old(3),K2old(3),K3old(3)
      real*8  pDeltaFipx(3),pDeltaFipy(3),pDeltaFipz(3),QTHmat(3,3)
      real*8  K1new(3),K2new(3),K3new(3)
      real*8  K1vect(3),K2vect(3),K3vect(3)
      real*8  extRPT2(6),extRPT3(6)
      real*8  pomocna(3,3)
      real*8  Q3T(3,3),Q4T(3,3),Q7T(3,3),Q8T(3,3),normFi3,normFi4
      real*8  normFi7,normFi8,surfaceload(3)
      real*8  folload3(3),folload4(3),folload7(3),folload8(3)
      real*8  forceF1,numsteps,skewfolload3(3,3),skewfolload4(3,3)
      real*8  skewfolload7(3,3),skewfolload8(3,3),surloadF1(3)

      
c     Output element type

      if(isw.eq.0 .and. ior.lt.0) then
        write(*,*) 'Elmt  08: 3-D Micropolar Geom NonLinear, Elastic'
        
c     Input material properties

      elseif(isw.eq.1) then
          
c.... set history storage for curvatures + quaternions
          
      ngp = 8 !number of gauss integration points
          
      nh1  = nh1 + ngp*4 +ngp*9 !8-numGPs*4-quaternion parameters + 8-numGPs*9-curvature tensor
      hplmax = nh1 !for paraview, nh1 + nh2 + nh3
      
		write(iow,2001)
          write(iow,'(A)') 'trilinear hexaedron (8 nodes per element)'
		if(ior.lt.0) write(*,2001)
          call dinput(d,6)
          call dinput(load,6) !BODY LOAD
          call dinput(surfaceload,3) !SURFACE LOAD
          
          write(iow,'(A)')    '   LOAD:'
          write(iow,'(A,E10.3)') '    px        = ',load(1)
          write(iow,'(A,E10.3)') '    py        = ',load(2)
          write(iow,'(A,E10.3)') '    pz        = ',load(3)
          write(iow,'(A,E10.3)') '    mx        = ',load(4)
          write(iow,'(A,E10.3)') '    my        = ',load(5)
          write(iow,'(A,E10.3)') '    mz        = ',load(6)
          
          write(iow,'(A)')    ' SURFACE LOAD:'
          write(iow,'(A,E10.3)') '    psx        = ',surfaceload(1)
          write(iow,'(A,E10.3)') '    psy        = ',surfaceload(2)
          write(iow,'(A,E10.3)') '    psz        = ',surfaceload(3)
          
          !FOR FOLLOWER LOAD
          text3(1)='YES'
          
          WRITE(iow,'(A)') text3(1)   
          if(pcomp(text3(1),'YES',3)) then
          write(iow,'(A)') 'WE HAVE FOLLOWER LOAD'
          else
          write(iow,'(A)') '  NO FOLLOWER LOAD'   
          endif
          
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
          
      elseif(isw.eq.14) then
          
      ngp = 8 !number of gauss integration points
      
      nn1 = 0
      
      do j1 = 1,ngp
          
c   Initialize quaternion parameters - specific for the 0th iteration of the 1st step
             
             do i=1,3
             quaOLD(i) = 0.0d0
             enddo
             quaOLD(4) = 1.0d0
             
c    Store data to history variables (quaternions)        
             
             do i=1,4
                  !hr (nh2+nn1-1+i) = quaOLD(i)
                  hr (nh1+nn1-1+i) = quaOLD(i)
             enddo !
             
c   Initialize curvature parameters
             
             call pzero(K1old,3)
             call pzero(K2old,3)
             call pzero(K3old,3)
             
             call zeromatrix(Knew,3,3)

             
c    Store data to history variables (curvatures)        

              do i=1,3
                  !hr(nh2+nn1+3+i) = K1old(i)
                  hr(nh1+nn1+3+i) = K1old(i)
              enddo !i
              
              do i=1,3
                  !hr(nh2+nn1+6+i) = K2old(i)
                  hr(nh1+nn1+6+i) = K2old(i)
              enddo !i
              
              do i=1,3
                  !hr(nh2+nn1+9+i) = K3old(i)
                  hr(nh1+nn1+9+i) = K3old(i)
              enddo !i
              
c     Advance history variable pointers
              
              nn1 = nn1 + 13
             
              !write (iow,*)
              !write(iow,'(A)') '-------	  Matrix Q	-------'
              !do i=1,3
              !    write (iow,*) (Q(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
              !!
              !write (iow,*)
              !write(iow,'(A)') '-------	  quaOLD nit=0	-------'
              !do i=1,4
              !    write (iow,*) quaOLD(i)
              !    write(iow,'(A)') '	new row'
              !end do
              !
              !write (iow,*)
              !write(iow,'(A)') '-------	  K1old nit=0	-------' 
              !do i=1,3
              !    write (iow,*) K1old(i)
              !    write(iow,'(A)') '	new row'
              !end do
              
      enddo !j1
          
          
      elseif(isw.eq.3 .or. isw.eq.6) then  !Form stiffness and residual
          
         write(iow,*)'TOTAL ELEMENT NO =',numel 
         write(iow,*)'ELEMENT NO =',n
         write(iow,*)'NSTEP =',nstep
         write(iow,*)'NITER =',niter
         
         !if(niter.eq.0) then
         !    write (iow,*) 'matrix Q has to be identity' 
         ! endif

		mi		=	d(1)
		ni		=	d(2)
		lambda	=	d(3)
		alfa	=	d(4)
		beta	=	d(5)
          gamma   =   d(6)
		
          call zeromatrix(s,nst,nst)
          call zeromatrix(r,ndf,nen)
          call pzero(blr1,6)
          call pzero(blr2,6)
          call zeromatrix(Iden,3,3)
          call zeromatrix(Bstress,3,3)
          call zeromatrix(BstressT,3,3)
          call zeromatrix(Gstress,3,3)
          call zeromatrix(Estrain,3,3)
          call zeromatrix(Knew,3,3)
          call zeromatrix(F,3,3) !deformation gradient
          call zeromatrix(FT,3,3)
          call zeromatrix(Q,3,3) !rotation matrix
          call zeromatrix(nulmatrix,3,3)
          call zeromatrix(term1,3,3)
          call zeromatrix(term2,3,3)
          call zeromatrix(term3,3,3)
          call zeromatrix(term4,3,3)
          call zeromatrix(term5,3,3)
          call zeromatrix(term6,3,3)
          call zeromatrix(term7,3,3)
          call zeromatrix(term8,3,3)
          call zeromatrix(term9,3,3)
          call zeromatrix(term10,3,3)
          call zeromatrix(term11,3,3)
          call zeromatrix(term12,3,3)
          call zeromatrix(term13,3,3)
          call zeromatrix(term14,3,3)
          call zeromatrix(term15,3,3)
          call zeromatrix(term16,3,3)
          call zeromatrix(term17,3,3)
          call zeromatrix(term18,3,3)
          call zeromatrix(term19,3,3)
          call zeromatrix(Ricci1,3,3)
          call zeromatrix(Ricci2,3,3)
          call zeromatrix(Ricci3,3,3)
          call zeromatrix(matterm19,3,3)
          call pzero(Ninabla,3)  !vector of interpolation function derivatives
          call pzero(Njnabla,3)  !vector of interpolation function derivatives, j-th node summation term
          call zeromatrix(matterm7T,3,3)
          call zeromatrix(blokKg,6,6)
          call zeromatrix(Kg,nst,nst)
          call pzero(blokKgpom,nst*nst)
          call zeromatrix(blokKm,6,6) 
          call zeromatrix(Km,nst,nst)
          call zeromatrix(Hmat,3,3)
          
c         Form Ricci submatrices
          
          Ricci1(2,3) =  1.0d0
          Ricci1(3,2) = -1.0d0
          Ricci2(1,3) = -1.0d0
          Ricci2(3,1) =  1.0d0
          Ricci3(1,2) =  1.0d0
          Ricci3(2,1) = -1.0d0
              
          
c			Compute Gauss quadrature points and weights

			l=2  !number of int points per direction 

			call int3d(l,lint,sg) !sg(4,8) - array of points and weights
                                    !sg(1,i) - x coordinate of point GP_i
                                    !sg(2,i) - y coordinate of point GP_i
                                    !sg(3,i) - z coordinate of point GP_i
                                    !sg(4,i) = w_i (weight factor)
              
c         Set history variables pointer
              
              nn1 = 0
              
              !write(iow,100) 'nn1=',nn1
              
c			Compute shape functions and derivatives 
				
          do l=1,lint !loop Gauss points
              
          call zeromatrix(F,3,3) !deformation gradient
          call zeromatrix(Q,3,3) !rotation matrix
          call pzero(pDeltaFipx,3)
          call pzero(pDeltaFipy,3)
          call pzero(pDeltaFipz,3)
          call pzero(DeltaFi,3)
              
              gpxL=0.0d0
              gpyL=0.0d0
              gpzL=0.0d0
                  
              gpxL=sg(1,l)
              gpyL=sg(2,l)
              gpzL=sg(3,l)
              
          call shp3d(sg(1,l),xsj,shp,xl,ndm,nel)
                      
                  dv = 0.0d0
                  dv = xsj*sg(4,l)
                  
				
			!shp(1,i) = Ni,x in Gauss point l, shp(2,i) = Ni,y in Gauss point l
			!shp(3,i) = Ni,z in Gauss point l, shp(4,i) = Ni in Gauss point l

                  
c          ---   GAUSS POINTS GLOBAL CS   ---    
                  gpx=0.0d0
                  gpy=0.0d0
                  gpz=0.0d0
                  
               do i=1,nel
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
         !&  ,gpz,     ')' !F, real 3dec pl
         !     write(iow,*)			!write blank line 
                
                
              !do i=1,8      
              !    write(iow,*) i
              !    write(iow,'(A,E10.3)') '	dN_i/dx=',shp(1,i)
              !    write(iow,'(A,E10.3)') '	dN_i/dy=',shp(2,i)
              !    write(iow,'(A,E10.3)') '	dN_i/dz=',shp(3,i)
              !    write(iow,'(A,E10.3)') '	N_i=',shp(4,i)
              !end do !i
               
c     Form the increment of microrotation DeltaFi = Summ N_i*DeltaFi_i (summation over nodes)
               
              
          call pzero(DeltaFi,3)
              
          !write(iow,'(A)') '-------	  niter.ne.0	----'
          !
          !    do i=1,3      
          !        write(iow,*) i
          !        write(iow,'(A,E10.3)') ' prije	DeltaFi(i)=',DeltaFi(i)
          !    end do !i
               
          do j=1,3
               do i=1,nel
                  DeltaFi(j) = DeltaFi(j) + shp(4,i)*ul(3+j,i,3) !ul(4,i,3) - inkrement of the j-th component Fi_j in node i
               enddo !i
          enddo !j
          
              !do i=1,3      
              !    write(iow,*) i
              !    write(iow,'(A,E10.3)') ' new DeltaFi(i)=',DeltaFi(i)
              !end do !i
          
          call rotqua (DeltaFi, quaDeltaFi) !convert increment of rotation vector DeltaFi to 4x1 UNIT quaternion quaDeltaFi
          
              !do i=1,4      
              !    write(iow,*) i
              !    write(iow,'(A,E10.3)') '	quaDeltaFi(i)=',quaDeltaFi(i)
              !end do !i
              
c   Recover data from history variables
              
              if (niter.eq.0) then
              
              do i=1,4
                  !quaOLD(i) = hr (nh2+nn1-1+i)
                  quaOLD(i) = hr (nh1+nn1-1+i)
              enddo !
              
              else
                  
              do i=1,4
                  quaOLD(i) = hr (nh2+nn1-1+i)
                  !quaOLD(i) = hr (nh1+nn1-1+i)
              enddo !
              
              endif
              
              write (iow,*)
              write(iow,'(A)') '-------	  quaOLD read from history	----'
              do i=1,4
                  write(iow,*) i
                  write(iow,'(A,E16.5)') '	quaOLD(i)=',quaOLD(i)
              end do
              
          call quamul (quaDeltaFi, quaOLD, quaNEW) !multiply 2 quaternions quaDeltaFi*quaOLD into quaNEW
          call quanrm ( quaNEW ) !normalize the new quaternion
          
              !do i=1,4      
              !    write(iow,*) i
              !    write(iow,'(A,E16.5)') '	quaOLD(i)=',quaOLD(i)
              !end do !i
              !
              !do i=1,4      
              !    write(iow,*) i
              !    write(iow,'(A,E16.5)') '	quaNEW(i)=',quaNEW(i)
              !end do !i
          
              do i=1,4
                  quaOLD(i) = quaNEW(i) !I HAVE TO STORE quaOLD(4)
              enddo !i
              
c    Store data to history variables (quaternions)
              
             do i=1,4
                  hr (nh2+nn1-1+i) = quaOLD(i)
              enddo !
              
          call quamatnew (quaNEW, Q) !convert quaternion to matrix
          
              !write (iow,*)
              !write(iow,'(A)') '----Matrix Q converted from quaNEW----'
              !do i=1,3
              !    write (iow,'(E13.6)') (Q(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
                    
               !endif !niter.ne.0
               
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
          
c     Form the Identity Matrix
          
          do i=1,3
              Iden(i,i) = 1.0d0
          enddo !i
          
c     Form the F matrix (deformation gradient)
          
c .........First the part without the identity matrix (only GRAD(u) - summ over nodes)
          
          !write(iow,'(A,E10.3)') '	ul(1,5,1)=',ul(1,5,1)
          
          do i=1,nel
              F(1,1) = F(1,1) + ul(1,i,1)*shp(1,i)
              F(1,2) = F(1,2) + ul(1,i,1)*shp(2,i)
              F(1,3) = F(1,3) + ul(1,i,1)*shp(3,i)
              F(2,1) = F(2,1) + ul(2,i,1)*shp(1,i)
              F(2,2) = F(2,2) + ul(2,i,1)*shp(2,i)
              F(2,3) = F(2,3) + ul(2,i,1)*shp(3,i)
              F(3,1) = F(3,1) + ul(3,i,1)*shp(1,i)
              F(3,2) = F(3,2) + ul(3,i,1)*shp(2,i)
              F(3,3) = F(3,3) + ul(3,i,1)*shp(3,i)
          enddo !i
          
c .........Add the identity matrix (F = I + GRAD(u))
          
          do i=1,3
              F(i,i) = F(i,i) + 1.0d0
          enddo !i
          
              write (iow,*)
              write(iow,'(A)') '-------	  Matrix F	-------'
              do i=1,3
                  write (iow,*) (F(i,j),j=1,3)
                  write(iow,'(A)') '	new row'
              end do
          
c .........Transpose of the deformation gradient matrix
          
		do i=1,3 !loop F -> FT
			do j=1,3
			FT(j,i)=F(i,j)
			end do !j
          end do !i
              
          
c     Form the E matrix (Biot-like strain tensor)
          
c .........Transpose of the rotation matrix
          
		do i=1,3 !loop Q -> QT
			do j=1,3
			QT(j,i)=Q(i,j)
			end do !j
          end do !i
          
c .........Form the matrix product QT.F
          
          call matmultip(QT,F,prod1,3,3,3)
          
c .........Form the strain as E = QT.F - I
          
          do i=1,3
              do j=1,3
                  Estrain(i,j) = prod1(i,j) - Iden(i,j)
              enddo !j
          enddo !i
          
c     Form the B matrix (Biot-like stress tensor)
              
          Bstress(1,1) = (lambda + 2.0d0*mi)*Estrain(1,1) + 
     &                    lambda*Estrain(2,2) + lambda*Estrain(3,3)
          Bstress(1,2) = (mi + ni)*Estrain(1,2) + (mi - ni)*Estrain(2,1)
          Bstress(1,3) = (mi + ni)*Estrain(1,3) + (mi - ni)*Estrain(3,1)
          Bstress(2,1) = (mi - ni)*Estrain(1,2) + (mi + ni)*Estrain(2,1)
          Bstress(2,2) = (lambda + 2.0d0*mi)*Estrain(2,2) + 
     &                    lambda*Estrain(1,1) + lambda*Estrain(3,3)
          Bstress(2,3) = (mi + ni)*Estrain(2,3) + (mi - ni)*Estrain(3,2)
          Bstress(3,1) = (mi - ni)*Estrain(1,3) + (mi + ni)*Estrain(3,1)
          Bstress(3,2) = (mi - ni)*Estrain(2,3) + (mi + ni)*Estrain(3,2)
          Bstress(3,3) = (lambda + 2.0d0*mi)*Estrain(3,3) + 
     &                    lambda*Estrain(1,1) + lambda*Estrain(2,2)
          
              write (iow,*)
              write(iow,'(A)') '-------	  Matrix Bstress	-------'
              do i=1,3
                  write (iow,*) (Bstress(i,j),j=1,3)
                  write(iow,'(A)') '	new row'
              end do
          
c .........Transpose of the Biot-like stress tensor
          
		do i=1,3 !loop Bstress -> BstressT
			do j=1,3
			BstressT(j,i)=Bstress(i,j)
			end do !j
          end do !i
          
c   Form the K matrix (curvature tensor)
          
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
          !if(niter.ne.0) then
          
c   Recover data from history variables (curvatures)
          
          if (niter.eq.0) then

              do i=1,3
                K1old(i) = hr(nh1+nn1+3+i)
              enddo !i
              
              do i=1,3
                K2old(i) = hr(nh1+nn1+6+i)
              enddo !i
              
              do i=1,3
                K3old(i) = hr(nh1+nn1+9+i)
              enddo !i
              
          else
              
              do i=1,3
                K1old(i) = hr(nh2+nn1+3+i)
              enddo !i
              
              do i=1,3
                K2old(i) = hr(nh2+nn1+6+i)
              enddo !i
              
              do i=1,3
                K3old(i) = hr(nh2+nn1+9+i)
              enddo !i
              
          endif
             
              !write (iow,*)
              !write(iow,'(A)') '-------	  K1old read from history	----'
              !do j1=1,3
              !    write (iow,'(E16.5)') K1old(j1)
              !    write(iow,'(A)') '	new row'
              !end do
              !
              !write (iow,*)
              !write(iow,'(A)') '-------	  K2old read from history	----'
              !do j1=1,3
              !    write (iow,'(E16.5)') K2old(j1)
              !    write(iow,'(A)') '	new row'
              !end do
              !
              !write (iow,*)
              !write(iow,'(A)') '-------	  K3old read from history	----'
              !do j1=1,3
              !    write (iow,'(E16.5)') K3old(j1)
              !    write(iow,'(A)') '	new row'
              !end do
             
c Compute Knew = Kold + .... 
             
c     Form H matrix
c         Norm of DeltaFi
             
              !do i=1,3      
              !    write(iow,*) i
              !    write(iow,'(A,E10.3)') ' again DeltaFi(i)=',DeltaFi(i)
              !end do !i
             
          normDeltaFi = sqrt(DeltaFi(1)**2+DeltaFi(2)**2+DeltaFi(3)**2) !ovo je ok

          
          !write(iow,'(A,E10.3)') '	normDeltaFi=',normDeltaFi
          
c         Form skew matrix from vector
          call hatmatrixfromvector(DeltaFi,skewDeltaFi)
c         Form matrix skewDeltaFi**2
          call matmultip(skewDeltaFi,skewDeltaFi,skewDeltaFiSQ,3,3,3)
c         Form constants
          
          if (normDeltaFi.eq.0.0d0) then
              const1=0.0d0
              const2=0.0d0
          else
          const1 = (1.0d0-cos(normDeltaFi))/((normDeltaFi)**2)
          const2 = (normDeltaFi-sin(normDeltaFi))/((normDeltaFi)**3)
          endif
          
          !write(iow,'(A,E10.3)') 'const1=',const1 !ok
          !write(iow,'(A,E10.3)') 'const2=',const2 !ok
          
             do j1=1,3
              do i=1,3
              Hmat(j1,i) = Iden(j1,i) + const1*skewDeltaFi(j1,i) +
     &                     const2*skewDeltaFiSQ(j1,i)
              enddo !i
             enddo !j1        
             
c         Form matrix QT.Hmat
          call matmultip(QT,Hmat,QTHmat,3,3,3)
          
c             Form partialDeltaFipartialX
          
          do j1=1,3
               do i=1,8
               pDeltaFipx(j1) = pDeltaFipx(j1) + shp(1,i)*ul(3+j1,i,3) !ul(4,i,3) - inkrement of the j-th component Fi_j in node i
               enddo !i
          enddo !j
          
c             Form partialDeltaFipartialY
          
          do j1=1,3
               do i=1,8
               pDeltaFipy(j1) = pDeltaFipy(j1) + shp(2,i)*ul(3+j1,i,3) !ul(4,i,3) - inkrement of the j-th component Fi_j in node i
               enddo !i
          enddo !j
          
c             Form partialDeltaFipartialZ
          
          do j1=1,3
               do i=1,8
               pDeltaFipz(j1) = pDeltaFipz(j1) + shp(3,i)*ul(3+j1,i,3) !ul(4,i,3) - inkrement of the j-th component Fi_j in node i
               enddo !i
          enddo !j
          
c         Form K1new
          
          call matvectmultip(QTHmat,pDeltaFipx,K1vect,3,3)
          
          do i=1,3
              K1new(i) = K1old(i) + K1vect(i)
          enddo !i
          
c         Form K2new
          
          call matvectmultip(QTHmat,pDeltaFipy,K2vect,3,3)
          
          do i=1,3
              K2new(i) = K2old(i) + K2vect(i)
          enddo !i
          
c         Form K3new
          
          call matvectmultip(QTHmat,pDeltaFipz,K3vect,3,3)
          
          do i=1,3
              K3new(i) = K3old(i) + K3vect(i)
          enddo !i
          
              !do i=1,3      
              !    write(iow,*) i
              !    write(iow,'(A,E16.5)') ' K1new(i)=',K1new(i)
              !end do !i
              !do i=1,3      
              !    write(iow,*) i
              !    write(iow,'(A,E16.5)') ' K2new(i)=',K2new(i)
              !end do !i
              !do i=1,3      
              !    write(iow,*) i
              !    write(iow,'(A,E16.5)') ' K3new(i)=',K3new(i)
              !end do !i
          
          
c     FORM CURVATURE TENSOR
          
          do i=1,3
              Knew(i,1) = K1new(i)
              Knew(i,2) = K2new(i)
              Knew(i,3) = K3new(i)
          enddo !i
             
             
c  K1old := K1new, K2old := K2new, K3old := K3new
             
          do i=1,3
              K1old(i) = K1new(i)
          enddo!i
          
          do i=1,3
              K2old(i) = K2new(i) 
          enddo!i
          
          do i=1,3
              K3old(i) = K3new(i)
          enddo!i
             
c    Store data to history variables (curvatures) UPDATE
            

              do i=1,3
                  hr(nh2+nn1+3+i) = K1old(i) 
              enddo !i
              
              do i=1,3
                  hr(nh2+nn1+6+i) = K2old(i) 
              enddo !i
              
              do i=1,3
                  hr(nh2+nn1+9+i) = K3old(i) !l-Gauss point, nh3+10,nh3+11,nh3+12
              enddo !i
              
              
c     Advance history variable pointers
              
              nn1 = nn1+13
                

          !endif !niter.ne.0 2nd time
          
          !write(iow,100) 'nn1=',nn1
          
          
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
c     Form the G matrix (Biot-like couple stress tensor)
          
          Gstress(1,1) = (alfa + 2.0d0*beta)*Knew(1,1) + 
     &                    alfa*Knew(2,2) + alfa*Knew(3,3)
          Gstress(1,2) = (beta + gamma)*Knew(1,2) + 
     &                   (beta - gamma)*Knew(2,1)
          Gstress(1,3) = (beta + gamma)*Knew(1,3) + 
     &                   (beta - gamma)*Knew(3,1)
          Gstress(2,1) = (beta - gamma)*Knew(1,2) + 
     &                   (beta + gamma)*Knew(2,1)
          Gstress(2,2) = (alfa + 2.0d0*beta)*Knew(2,2) + 
     &                    alfa*Knew(1,1) + alfa*Knew(3,3)
          Gstress(2,3) = (beta + gamma)*Knew(2,3) + 
     &                   (beta - gamma)*Knew(3,2)
          Gstress(3,1) = (beta - gamma)*Knew(1,3) + 
     &                   (beta + gamma)*Knew(3,1)
          Gstress(3,2) = (beta - gamma)*Knew(2,3) + 
     &                   (beta + gamma)*Knew(3,2)
          Gstress(3,3) = (alfa + 2.0d0*beta)*Knew(3,3) + 
     &                    alfa*Knew(1,1) + alfa*Knew(2,2)
          
              !write (iow,*)
              !write(iow,'(A)') '-------	  Matrix Gstress	-------'
              !do i=1,3
              !    write (iow,*) (Gstress(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
          

c  -----     TERMS NEEDED IN THE RESIDUAL r[6,8]   -----

c ->->->->->-> LOOP i- for i-th nodal residual

          do i=1,8
              
              !write(iow,*) i
              
c  >>>>>>>>>>  TERM R1: QB(Ni \nabla)  -------- 
              
c .............. Form the matrix product Q.Bstress = QB
          
                  call matmultip(Q,Bstress,QB,3,3,3)
              
                  Ninabla(1) = shp(1,i)    !i = i-th node residual
                  Ninabla(2) = shp(2,i)
                  Ninabla(3) = shp(3,i)
                  
                  call matvectmultip(QB,Ninabla,res1,3,3)
                  
              !write (iow,*)
              !write(iow,'(A)') '-------	  res1	-------'
              !do m=1,3
              !    write (iow,*) res1(m)
              !    write(iow,'(A)') '	new row'
              !end do
                  
c  >>>>>>>>>>  TERM R2: 2*Ni*ax ( skew (F.BT.QT) )  -------- 
                  
                  call triplematmultip(F,BstressT,QT,res2a,3,3,3,3)
                  call skewmatrixpart(res2a,skewres2a,3,3)
                  call axialofskew(skewres2a,axres2a)
                  
                  do m=1,3
                  res2(m) = 2.0d0*shp(4,i)*axres2a(m)
                  enddo !m
                  
              !write (iow,*)
              !write(iow,'(A)') '-------	  res2	-------'
              !do m=1,3
              !    write (iow,*) res2(m)
              !    write(iow,'(A)') '	new row'
              !end do
                  
c  >>>>>>>>>>  TERM R3: QG(Ni \nabla)  -------- 
              
c .............. Form the matrix product Q.Gstress = QG
          
                  call matmultip(Q,Gstress,QG,3,3,3)
                  call matvectmultip(QG,Ninabla,res3,3,3)
                  
              !write (iow,*)
              !write(iow,'(A)') '-------	  res3	-------'
              !do m=1,3
              !    write (iow,*) res3(m)
              !    write(iow,'(A)') '	new row'
              !end do
                  
c  >>>>>>>>>>  TERMs extRPT2 and extRPT3:  -------- 
              
              
              if(pcomp(text2(1),'YES',3)) then
                  !write(iow,'(A)') '	PT 3'
                  do m=1,3
                  extRPT3(m) = shp(4,i)*load(m)
                  enddo !m
                  do m=1,3
                  extRPT3(m+3) = shp(4,i)*2.0d0*(gpx-gpy-gpz)
                  enddo !m
              else
                  !write(iow,'(A)') '	PT 2'
                  do m=1,6
                  extRPT2(m) = shp(4,i)*load(m)
                  enddo !m
              endif
                  
c  xxxxxxxxxxxxxxxxxx     FORM THE ELEMENT RESIDUAL VECTOR     xxxxxxxxxxxxxxxxxx
                  
c ............. first form the blocks blr1(6) and blr2(6) ..............
                  
          call pzero(blr1,6)
          call pzero(blr2,6)
                  
          do m=1,3
              blr1(m)   = res1(m)
              blr2(m+3) = res2(m) + res3(m)
          enddo !m
          
          do m=1,6
		r(m,i)=r(m,i)-(blr1(m)+blr2(m)-extRPT2(m)-extRPT3(m))*dv
          enddo !m

          enddo !i = i-th node residual
          
              !write (iow,*)
              !write(iow,'(A)') '-------	  residual	-------'
              !do m=1,6
              !    write (iow,*) (r(m,i),i=1,8)
              !    write(iow,'(A)') '	new row'
              !end do
                  
          
          
c  -----     TERMS NEEDED IN THE GEOMETRIC STIFFNESS MATRIX Kg[6,6]   -----
          
c ->->->->->-> LOOP i- for i-th nodal residual and j-for j-th element node
          
          do i=1,8
              do j=1,8

c  >>>>>>>>>>  1st TERM: -hat(QB(Ni \nabla))*Nj  -------- 
          
c .............. Form the matrix product Q.Bstress = QB
          
                  call matmultip(Q,Bstress,QB,3,3,3)
                  
c .............. Form the vector of shape function derivatives
                  
                  Ninabla(1) = shp(1,i)    !i = i-th node residual
                  Ninabla(2) = shp(2,i)
                  Ninabla(3) = shp(3,i)
                  
c .............. Form the matrix.vector product QB(3,3).Ninabla(3) = vect1(3)
                  
                  call matvectmultip(QB,Ninabla,vect1,3,3)
                  
c .............. Form the skew matrix from vect1(3)
                  
                  call hatmatrixfromvector(vect1,skewvect1)
                  
c .............. Form the term -hat(vect1)(3,3)*Nj = term1(3,3)
                  
              do m=1,3
                  do j1=1,3
                  term1(m,j1) = -skewvect1(m,j1)*shp(4,j) !j = j-th node summation term
                  enddo !n
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  term1	-------'
              !do i=1,3
              !    write (iow,*) (term1(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
              
c  >>>>>>>>>>  2nd TERM: +2*lambda*Nj*Q(Ni \nabla).(vecL)T --------
          
c .............. Form the matrix product F.QT = FQT
          
                  call matmultip(F,QT,FQT,3,3,3)
                  
c .............. Form the skew matrix part skew(FQT)
                  
                  call skewmatrixpart(FQT,skewFQT,3,3)
                  
c .............. Form the axial vector of a skew matrix vecL = axskewFQT (written as L in the notes)
                  
                  call axialofskew(skewFQT,axskewFQT)
                  
c .............. Form the matrix.vector product Q(3,3).Ninabla(3) = vect2(3)
                  
                  call matvectmultip(Q,Ninabla,vect2,3,3)
                  
c .............. Form the tensor product between vectors vect2(3) \otimes axskewFQT(3)               
                  
                  call tensorprodvectors(vect2,axskewFQT,matterm2)
                  
c .............. Form the term 2*lambda*Nj*Q(3,3).Ninabla(3) \otimes axskew(FQT)  =  term2(3,3) MNOZENJE MATTRICE SKALAROM 2*lambda*Nj
                  
              do m=1,3
                  do j1=1,3
                  term2(m,j1) = 2.0d0*lambda*shp(4,j)*matterm2(m,j1) !j = j-th node summation term
                  enddo !n
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  term2	-------'
              !do i=1,3
              !    write (iow,*) (term2(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
                                   
c  >>>>>>>>>>  3rd TERM: hat(QB(Nj \nabla))*Ni --------
c ...............QB already formed                 
c .............. Form the vector of shape function derivatives
                  
                  Njnabla(1) = shp(1,j)    !j = j-th node summation term
                  Njnabla(2) = shp(2,j)
                  Njnabla(3) = shp(3,j)
                  
c .............. Form the matrix.vector product QB(3,3).Njnabla(3) = vect3(3)
                  
                  call matvectmultip(QB,Njnabla,vect3,3,3)
                  
c .............. Form the skew matrix from vect3(3)
                  
                  call hatmatrixfromvector(vect3,skewvect3)
                  
c .............. Form the term hat(vect3)(3,3)*Ni = term3(3,3)
                  
              do m=1,3
                  do j1=1,3
                  term3(m,j1) = skewvect3(m,j1)*shp(4,i) !i = i-th node residual
                  enddo !n
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  term3	-------'
              !do i=1,3
              !    write (iow,*) (term3(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
                  
              
c  >>>>>>>>>>  4th TERM: +2*lambda*Ni*(vecL).(Q(Nj \nabla))T = 2*lambda*Ni*(vecL) \otimes (Q(Nj \nabla)) --------
              
c................(vecL)= axskewFQT already formed
c................ Njnabla(3) already formed
                  
c .............. Form the matrix.vector product Q(3,3).Njnabla(3) = vect4(3)
                  
                  call matvectmultip(Q,Njnabla,vect4,3,3)
                  
c .............. Form the tensor product between vectors axskewFQT \otimes vect4             
                  
                  call tensorprodvectors(axskewFQT,vect4,matterm4)
                  
c .............. Form the term 2*lambda*Ni*axskewFQT \otimes Q(3,3).Njnabla(3)  =  term4(3,3) MNOZENJE MATTRICE SKALAROM 2*lambda*Ni
                  
              do m=1,3
                  do j1=1,3
                  term4(m,j1) = 2.0d0*lambda*shp(4,i)*matterm4(m,j1) !i = i-th node residual
                  enddo !j1
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  term4	-------'
              !do i=1,3
              !    write (iow,*) (term1(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
              
c  >>>>>>>>>>  5th TERM: -hat(QG(Ni \nabla))*Nj  --------
          
c .............. Form the matrix product Q.Gstress = QG
          
                  call matmultip(Q,Gstress,QG,3,3,3)
                  
c................ Ninabla(3) already formed
                  
c .............. Form the matrix.vector product QG(3,3).Ninabla(3) = vect5(3)
                  
                  call matvectmultip(QG,Ninabla,vect5,3,3)
                  
c .............. Form the skew matrix from vect5(3)
                  
                  call hatmatrixfromvector(vect5,skewvect5)
                  
c .............. Form the term -hat(vect5)(3,3)*Nj = term1(3,3)
                  
              do m=1,3
                  do j1=1,3
                  term5(m,j1) = -skewvect5(m,j1)*shp(4,j) !j = j-th node summation term
                  enddo !j1
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  term5	-------'
              !do i=1,3
              !    write (iow,*) (term5(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
          
c  >>>>>>>>>>  6th TERM: 4*lambda*Ni*Nj* axskewFQT \otimes axskewFQT   --------
              
c .............. Form the tensor product between vectors axskewFQT \otimes axskewFQT             
                  
                  call tensorprodvectors(axskewFQT,axskewFQT,matterm6)
                  
c .............. Form the term 4*lambda*Ni*Nj* axskewFQT \otimes axskewFQT  =  term6(3,3) MNOZENJE MATTRICE SKALAROM 4*lambda*Ni*Nj
                  
              do m=1,3
                  do j1=1,3
                  term6(m,j1) = 4.0d0*lambda*shp(4,i)*shp(4,j)*
     &                          matterm6(m,j1) !i = i-th node residual, j = j-th node summation term
                  enddo !n
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  term6	-------'
              !do i=1,3
              !    write (iow,*) (term6(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
              
c  >>>>>>>>>>  7th TERM: Ni*Nj*[(FBTQT)T-tr(FBTQT)*Iden]   --------
              
              call triplematmultip(F,BstressT,QT,matterm7,3,3,3,3)
              
c ............Transpose of the matrix F.BstressT.QT = matterm7
          
		    do m=1,3 !loop matterm7 -> matterm7T
			    do j1=1,3
			    matterm7T(j1,m)=matterm7(m,j1)
			    end do !j
              end do !i
              
c ............Find the trace of matrix F.BstressT.QT = matterm7
              
              call traceofmatrix(matterm7,trmatterm7)
              
c .............. Form the term Ni*Nj*[(FBTQT)T-tr(FBTQT)*Iden]  =  term7(3,3)
                  
              do m=1,3
                  do j1=1,3
                  term7(m,j1) = shp(4,i)*shp(4,j)*(matterm7T(m,j1)-
     &                         trmatterm7*Iden(m,j1)) !i = i-th node residual
                  enddo !n
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  term7	-------'
              !do i=1,3
              !    write (iow,*) (term7(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
              
              
c  -----     FORMING THE BLOCK MATRIX blokKg(6,6)
              
              do m=1,3
                  do j1=1,3
                      blokKg(m,j1) = nulmatrix(m,j1)
                      blokKgpom(m,j1,i,j) = nulmatrix(m,j1)
                  enddo !j1
                  do j1=4,6
                      blokKg(m,j1) = term1(m,j1-3)+term2(m,j1-3)
                      blokKgpom(m,j1,i,j) = term1(m,j1-3)+term2(m,j1-3)
                  enddo !j1
              enddo !m - this fills the first 3rows and all the 6 columns
              
              do m=4,6
                  do j1=1,3
                      blokKg(m,j1) = term3(m-3,j1)+term4(m-3,j1)
                      blokKgpom(m,j1,i,j) = term3(m-3,j1)+term4(m-3,j1)
                  enddo !j1
                  do j1=4,6
                      blokKg(m,j1) =term5(m-3,j1-3)+term6(m-3,j1-3)+
     &                              term7(m-3,j1-3)
                      blokKgpom(m,j1,i,j) =term5(m-3,j1-3)+
     &                              term6(m-3,j1-3)+term7(m-3,j1-3)
                  enddo !j1
              enddo !m - this fills the last 3rows and all the six columns
              
              
c  ....      MATRIX [Kg] = 48x48      ....       

              do m=1,6
                  do j1=1,6
                      Kg(6*i-6+m,6*j-6+j1) = blokKg(m,j1)
                  enddo !j1
              enddo !m
              
              enddo !i -----> FROM ABOVE
          enddo !j
              
              !write (iow,*)
              !write(iow,'(A)') '-------	Matrix Kg	-------'
              !do m=1,nst
              !    write (iow,*) (Kg(m,j1),j1=1,nst)
              !    write(iow,'(A)') '	new row'
              !end do
              
c  -----     TERMS NEEDED IN THE MATERIAL STIFFNESS MATRIX Km[i,j]   -----
              
              
          
c ->->->->->-> LOOP i- for i-th nodal residual and j-for j-th element node
          
          do i=1,8
              do j=1,8
                  
                  Ninabla(1) = shp(1,i)    !i = i-th node residual
                  Ninabla(2) = shp(2,i)
                  Ninabla(3) = shp(3,i)
                  
                  Njnabla(1) = shp(1,j)    !j = j-th node summation term
                  Njnabla(2) = shp(2,j)
                  Njnabla(3) = shp(3,j)

c  >>>>>>>>>>  8th TERM: lambda*Q.(Ni \nabla) \otimes (Nj \nabla) . QT  -------- 
              
              
       call tensorprodvectors(Ninabla,Njnabla,matterm8a)
       call triplematmultip(Q,matterm8a,QT,matterm8b,3,3,3,3)
       
              do m=1,3
                  do j1=1,3
                  term8(m,j1) = matterm8b(m,j1)*lambda
                  enddo !j1
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	term8	-------'
              !do m=1,3
              !    write (iow,*) (term8(m,j1),j1=1,3)
              !    write(iow,'(A)') '	new row'
              !end do

c  >>>>>>>>>>  9th TERM: (mi+ni)*scalarnabla.Iden  --------

              scnabla=0.0d0
          
              do m=1,3
                  scnabla = scnabla + Ninabla(m)*Njnabla(m)
              enddo !m
          
          
              do m=1,3
                  do j1=1,3
                  term9(m,j1) = (mi+ni)*scnabla*Iden(m,j1)
                  enddo !j1
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	term9	-------'
              !do m=1,3
              !    write (iow,*) (term9(m,j1),j1=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
              
c  >>>>>>>>>>  10th TERM: (mi-ni)*Q.(Nj \nabla) \otimes (Ni \nabla) . QT  -------- 
              
              
       call tensorprodvectors(Njnabla,Ninabla,matterm10a)
       call triplematmultip(Q,matterm10a,QT,matterm10b,3,3,3,3)
       
              do m=1,3
                  do j1=1,3
                  term10(m,j1) = matterm10b(m,j1)*(mi-ni)
                  enddo !j1
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	term10	-------'
              !do m=1,3
              !    write (iow,*) (term10(m,j1),j1=1,3)
              !    write(iow,'(A)') '	new row'
              !end do

c  >>>>>>>>>>  11st TERM: (mi+ni)*hat(F(Ni \nabla))*Nj  -------- 
                  
c .............. Form the matrix.vector product F(3,3).Ninabla(3) = vect11(3)
                  
                  call matvectmultip(F,Ninabla,vect11,3,3)
                  
c .............. Form the skew matrix from vect11(3)
                  
                  call hatmatrixfromvector(vect11,skewvect11)
                  
c .............. Form the term (mi+ni)*hat(F(Ni \nabla))*Nj = term11(3,3)
                  
              do m=1,3
                  do j1=1,3
                  term11(m,j1) = (mi+ni)*shp(4,j)*skewvect11(m,j1) !j = j-th node summation term
                  enddo !n
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	term11	-------'
              !do m=1,3
              !    write (iow,*) (term11(m,j1),j1=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
              
          
c  >>>>>>>>>>  12th TERM: -(mi-ni)*Nj*Q.FT.hat(Q(Ni \nabla))  -------- 
                  
c .............. Form the matrix.vector product Q(3,3).Ninabla(3) = vect12(3)
                  
                  call matvectmultip(Q,Ninabla,vect12,3,3)
                  
c .............. Form the skew matrix from vect12(3)
                  
                  call hatmatrixfromvector(vect12,skewvect12)
                  
c .............. Form the triple matrix product Q.FT.skewvect12 =matterm12
                  
                 call triplematmultip(Q,FT,skewvect12,matterm12,3,3,3,3)
                  
c .............. Form the term -(mi-ni)*Nj*Q.FT.hat(Q(Ni \nabla)) = term12(3,3)
                  
              do m=1,3
                  do j1=1,3
                  term12(m,j1) = -(mi-ni)*shp(4,j)*matterm12(m,j1) !j = j-th node summation term
                  enddo !n
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	term12	-------'
              !do m=1,3
              !    write (iow,*) (term12(m,j1),j1=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
              
c  >>>>>>>>>>  13th TERM: -(mi+ni)*Ni*hat(F(Nj \nabla))  -------- 
                  
c .............. Form the matrix.vector product F(3,3).Njnabla(3) = vect13(3)
                  
                  call matvectmultip(F,Njnabla,vect13,3,3)
                  
c .............. Form the skew matrix from vect13(3)
                  
                  call hatmatrixfromvector(vect13,skewvect13)
                  
c .............. Form the term -(mi+ni)*Ni*hat(F(Nj \nabla)) = term13(3,3)
                  
              do m=1,3
                  do j1=1,3
                  term13(m,j1) = -(mi+ni)*shp(4,i)*skewvect13(m,j1) !i = i-th node residual
                  enddo !n
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  term13	-------'
              !do i=1,3
              !    write (iow,*) (term13(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
          
c  >>>>>>>>>>  14th TERM: (mi-ni)*Ni*hat(Q(Nj \nabla)).F.QT  --------  
                  
c .............. Form the matrix.vector product Q(3,3).Njnabla(3) = vect14(3)
                  
                  call matvectmultip(Q,Njnabla,vect14,3,3) 
                  
c .............. Form the skew matrix from vect14(3)
                  
                  call hatmatrixfromvector(vect14,skewvect14)
                  
c .............. Form the triple matrix product skewvect14.F.QT =matterm14
                  
                 call triplematmultip(skewvect14,F,QT,matterm14,3,3,3,3)
                  
c .............. Form the term (mi-ni)*Ni*hat(Q(Nj \nabla)).F.QT = term14(3,3)
                  
              do m=1,3
                  do j1=1,3
                  term14(m,j1) = (mi-ni)*shp(4,i)*matterm14(m,j1) !i = i-th node residual
                  enddo !n
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  term14	-------'
              !do i=1,3
              !    write (iow,*) (term1(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
              
c  >>>>>>>>>>  15th TERM: alpha*Q.(Ni \nabla) \otimes (Nj \nabla) . QT  -------- 
       
              do m=1,3
                  do j1=1,3
                  term15(m,j1) = matterm8b(m,j1)*alfa
                  enddo !j1
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  term15	-------'
              !do i=1,3
              !    write (iow,*) (term15(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
              
c  >>>>>>>>>>  16th TERM: (beta+gamma)*scalarnabla.Iden  --------

              do m=1,3
                  do j1=1,3
                  term16(m,j1) = (beta+gamma)*scnabla*Iden(m,j1)
                  enddo !j1
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  term16	-------'
              !do i=1,3
              !    write (iow,*) (term16(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
              
c  >>>>>>>>>>  17th TERM: (beta-gamma)*Q.(Nj \nabla) \otimes (Ni \nabla) . QT  -------- 
       
              do m=1,3
                  do j1=1,3
                  term17(m,j1) = matterm10b(m,j1)*(beta-gamma)
                  enddo !j1
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  term17	-------'
              !do i=1,3
              !    write (iow,*) (term17(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
          
c  >>>>>>>>>>  18th TERM: -(mi+ni)*Ni*Nj*[(F.FT)-tr(F.FT).Iden]  --------
          
c ............Find the product F.FT = matterm18a
          
              call matmultip(F,FT,matterm18a,3,3,3)
              
c ............Find the trace of matrix F.FT
              
              call traceofmatrix(matterm18a,trmatterm18a)
              
c .............. Form the term -(mi+ni)*Ni*Nj*[(F.FT)-tr(F.FT).Iden]  =  term17(3,3)
                  
              do m=1,3
                  do j1=1,3
                  term18(m,j1) = -(mi+ni)*shp(4,i)*shp(4,j)*
     &             (matterm18a(m,j1)-trmatterm18a*Iden(m,j1))
                  enddo !n
              enddo !m
              
              !write (iow,*)
              !write(iow,'(A)') '-------	  term18	-------'
              !do i=1,3
              !    write (iow,*) (term18(i,j),j=1,3)
              !    write(iow,'(A)') '	new row'
              !end do
              
c  >>>>>>>>>>  19th TERM: (mi-ni)*Ni*Nj*[m1 m2 m3]  -------- where m1 = ax ( 2 skew (F.QT.Ricci1.F.QT) )
c                                                                  m2 = ax ( 2 skew (F.QT.Ricci2.F.QT) )
c                                                                  m3 = ax ( 2 skew (F.QT.Ricci3.F.QT) )
              
c ............Find the product F.QT = matterm19a
              
              call matmultip(F,QT,matterm19a,3,3,3)
              
c ............ Form the triple matrix products matterm19a.Ricci1.matterm19a = matterm19b1
c ............ Form the triple matrix products matterm19a.Ricci2.matterm19a = matterm19b2                  
c ............ Form the triple matrix products matterm19a.Ricci3.matterm19a = matterm19b3
              
              call triplematmultip(matterm19a,Ricci1,matterm19a,
     &                                matterm19b1,3,3,3,3)
              call triplematmultip(matterm19a,Ricci2,matterm19a,
     &                                matterm19b2,3,3,3,3)
              call triplematmultip(matterm19a,Ricci3,matterm19a,
     &                                matterm19b3,3,3,3,3)
              
c ............ Form the double skew symm matrix of matterm19b1 = dskewmatterm19b1
c ............ Form the double skew symm matrix of matterm19b2 = dskewmatterm19b2
c ............ Form the double skew symm matrix of matterm19b3 = dskewmatterm19b3              
              
             call doubleskewmatrixpart(matterm19b1,dskewmatterm19b1,3,3)
             call doubleskewmatrixpart(matterm19b2,dskewmatterm19b2,3,3)
             call doubleskewmatrixpart(matterm19b3,dskewmatterm19b3,3,3)
             
c ............ Form the axial vector of dskewmatterm19b1 = axmatterm19b1
c ............ Form the axial vector of dskewmatterm19b2 = axmatterm19b2
c ............ Form the axial vector of dskewmatterm19b3 = axmatterm19b3
             
             call axialofskew(dskewmatterm19b1,axmatterm19b1)
             call axialofskew(dskewmatterm19b2,axmatterm19b2)
             call axialofskew(dskewmatterm19b3,axmatterm19b3)
             
c ............ Store the axial vectors in a matrix = matterm19(3,3)
             
             do m=1,3
                 matterm19(m,1) = axmatterm19b1(m)
                 matterm19(m,2) = axmatterm19b2(m)
                 matterm19(m,3) = axmatterm19b3(m)
             enddo !m
             
              do m=1,3
                do j1=1,3
                term19(m,j1) = (mi-ni)*shp(4,i)*shp(4,j)*matterm19(m,j1)
                enddo !n
              enddo !m
          
          
c  -----     FORMING THE BLOCK MATRIX blokKm(6,6)
              
              do m=1,3
                  do j1=1,3
                     blokKm(m,j1) = term8(m,j1)+term9(m,j1)+term10(m,j1)
                  enddo !j1 this fills the first 3rows and the first 3 columns
                  do j1=4,6
                     blokKm(m,j1) = term11(m,j1-3)+term12(m,j1-3)
                  enddo !j1 j1 this fills the first 3rows and the second 3 columns
              enddo !m - this fills the first 3rows and all the 6 columns
              
              do m=4,6
                  do j1=1,3
                      blokKm(m,j1) = term13(m-3,j1)+term14(m-3,j1)
                  enddo !j1 this fills the second 3rows and the first 3 columns
                  do j1=4,6
                      blokKm(m,j1) =term15(m-3,j1-3)+term16(m-3,j1-3)+
     &                              term17(m-3,j1-3)+term18(m-3,j1-3)+
     &                              term19(m-3,j1-3)
                  enddo !j1
              enddo !m - this fills the last 3rows and all the 6 columns
              
              !write (iow,*)
              !write(iow,'(A)') '-------	block Matrix blKm	-------'
              !write(iow,'(A,i5)') 'i = ',i
              !write(iow,'(A,i5)') 'j = ',j
              !do m=1,6
              !    write (iow,*) (blokKm(m,j1),j1=1,6)
              !    write(iow,'(A)') '	new row'
              !end do
              
              
c  ....      MATRIX [Km] = 48x48      ....       

              do m=1,6
                  do j1=1,6
                      Km(6*i-6+m,6*j-6+j1) = blokKm(m,j1)
                  enddo !j1
              enddo !m
              
              enddo !i -----> FROM ABOVE
          enddo !j
              
              !write (iow,*)
              !write(iow,'(A)') '-------	Matrix Km	-------'
              !do m=1,nst
              !    write (iow,*) (Km(m,j1),j1=1,nst)
              !    write(iow,'(A)') '	new row'
              !end do
          
c  xxxxxxxxxxxxxxxxxx     FORM THE ELEMENT STIFFNESS MATRIX     xxxxxxxxxxxxxxxxxx
          
		do i=1,nst
			do j=1,nst
		s(i,j)=s(i,j)+(Kg(i,j) + Km(i,j))*dv
			end do !j
          end do !i
          
              !write (iow,*)
              !write(iow,'(A)') '	Stiffness Matrix'
              !write(iow,'(A,i5)') 'GP = ',l
              !do i=1,nst
              !    write (iow,*) (s(i,j),j=1,nst)
              !    write (iow,*)
              !    write(iow,'(A)') '	new row'
              !end do
          
          end do !l   END LOOP GP
          
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
c CONTRIBUTION OF THE SURFACE LOADING INTO THE RESIDUAL 
          
          
      if(n.eq.numel) then
          
      
      forceF1 = 0.0d0
       
      numsteps = 1.0d0  
      !
      forceF1 = 0.0833333333d0*surfaceload(1)/numsteps*nstep
      !
      !write(iow,'(A,E10.3)') '    forceF1        = ',forceF1
          
	call zeromatrix(Q3T,3,3)
	call zeromatrix(Q4T,3,3)
	call zeromatrix(Q7T,3,3)
	call zeromatrix(Q8T,3,3)
      
      call pzero(surloadF1,3)
	
c Form orientation matrices

c 		Find the norm of nodal microrotations

	normFi3 = sqrt(ul(4,3,1)**2+ul(5,3,1)**2+ul(6,3,1)**2)
	normFi4 = sqrt(ul(4,4,1)**2+ul(5,4,1)**2+ul(6,4,1)**2)
	normFi7 = sqrt(ul(4,7,1)**2+ul(5,7,1)**2+ul(6,7,1)**2)
	normFi8 = sqrt(ul(4,8,1)**2+ul(5,8,1)**2+ul(6,8,1)**2)
	
c       Form orientation matrices NOT TRANSPOSED!!!

	Q3T(1,1)=cos(normFi3)
	Q3T(1,2)=-sin(normFi3)
	Q3T(2,1)=sin(normFi3)
	Q3T(2,2)=cos(normFi3)
      Q3T(3,3)=1.0d0
	
	Q4T(1,1)=cos(normFi4)
	Q4T(1,2)=-sin(normFi4)
	Q4T(2,1)=sin(normFi4)
	Q4T(2,2)=cos(normFi4)
      Q4T(3,3)=1.0d0
	
	Q7T(1,1)=cos(normFi7)
	Q7T(1,2)=-sin(normFi7)
	Q7T(2,1)=sin(normFi7)
	Q7T(2,2)=cos(normFi7)
      Q7T(3,3)=1.0d0
	
	Q8T(1,1)=cos(normFi8)
	Q8T(1,2)=-sin(normFi8)
	Q8T(2,1)=sin(normFi8)
	Q8T(2,2)=cos(normFi8)
      Q8T(3,3)=1.0d0
	
	        !      write (iow,*)
         !         write(iow,'(A)') '----Matrix Q3T ----'
         !         do i=1,3
         !             write (iow,'(E13.6)') (Q3T(i,j),j=1,3)
         !             write(iow,'(A)') '	new row'
         !         end do
			      !
			      !write (iow,*)
         !         write(iow,'(A)') '----Matrix Q4T ----'
         !         do i=1,3
         !             write (iow,'(E13.6)') (Q4T(i,j),j=1,3)
         !             write(iow,'(A)') '	new row'
         !         end do
			      !
			      !write (iow,*)
         !         write(iow,'(A)') '----Matrix Q7T ----'
         !         do i=1,3
         !             write (iow,'(E13.6)') (Q7T(i,j),j=1,3)
         !             write(iow,'(A)') '	new row'
         !         end do
			      !
			      !write (iow,*)
         !         write(iow,'(A)') '----Matrix Q8T ----'
         !         do i=1,3
         !             write (iow,'(E13.6)') (Q8T(i,j),j=1,3)
         !             write(iow,'(A)') '	new row'
         !         end do
			  
	
c  Contribution of external surface loading to the residual vector
      
      call pzero(surloadF1,3)
      
      surloadF1(1) = forceF1

	call matvectmultip(Q3T,surloadF1,folload3,3,3)
	call matvectmultip(Q4T,surloadF1,folload4,3,3)
	call matvectmultip(Q7T,-surloadF1,folload7,3,3)
	call matvectmultip(Q8T,-surloadF1,folload8,3,3)
	
	
	     !      write (iow,*)
      !         write(iow,'(A)') '-------	folload3	-------'
      !         do m=1,3
      !            write (iow,'(E16.5)') folload3(m)
      !            write(iow,'(A)') '	new row'
      !         end do
			   !
	     !      write (iow,*)
      !         write(iow,'(A)') '-------	folload4	-------'
      !         do m=1,3
      !            write (iow,'(E16.5)') folload4(m)
      !            write(iow,'(A)') '	new row'
      !         end do
			   !
	     !      write (iow,*)
      !         write(iow,'(A)') '-------	folload7	-------'
      !         do m=1,3
      !            write (iow,'(E16.5)') folload7(m)
      !            write(iow,'(A)') '	new row'
      !         end do
			   !
			   !write (iow,*)
      !         write(iow,'(A)') '-------	folload8	-------'
      !         do m=1,3
      !            write (iow,'(E16.5)') folload8(m)
      !            write(iow,'(A)') '	new row'
      !         end do
               
c    Contribute to the element residual
               
          do m=1,3
		r(m,3)=r(m,3)+folload3(m)
          r(m,4)=r(m,4)+folload4(m)
          r(m,7)=r(m,7)+folload7(m)
          r(m,8)=r(m,8)+folload8(m)
          enddo !m
          
          
      if(niter.ne.0) then
          
c    Contribute to the element stiffness matrix
          
          call hatmatrixfromvector(folload3,skewfolload3)
          call hatmatrixfromvector(folload4,skewfolload4)
          call hatmatrixfromvector(folload7,skewfolload7)
          call hatmatrixfromvector(folload8,skewfolload8)
          

c    ....For the 3rd node
          
      fn=3    
          
      do i=1,3
          do j=1,3
          s(6*fn-6+i,6*fn-3+j)=s(6*fn-6+i,6*fn-3+j)+skewfolload3(i,j)
          enddo !i
      enddo !j
      
c    ....For the 4th node
      
      fn=4    
          
      do i=1,3
          do j=1,3
          s(6*fn-6+i,6*fn-3+j)=s(6*fn-6+i,6*fn-3+j)+skewfolload4(i,j)
          enddo !i
      enddo !j
      
c    ....For the 7th node
      
      fn=7    
          
      do i=1,3
          do j=1,3
          s(6*fn-6+i,6*fn-3+j)=s(6*fn-6+i,6*fn-3+j)+skewfolload7(i,j)
          enddo !i
      enddo !j
      
c    ....For the 8th node
      
      fn=8    
          
      do i=1,3
          do j=1,3
          s(6*fn-6+i,6*fn-3+j)=s(6*fn-6+i,6*fn-3+j)+skewfolload8(i,j)
          enddo !i
      enddo !j
          
      endif !(niter.ne.0)
               
      endif !(n.eq.numel)
      
			   
			   
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
          
          
          !write(iow,100) 'nn1=',nn1
          !    !
              !write (iow,*)
              !write(iow,'(A)') '	Residual'
              !do m=1,6
              !    write (iow,*) (r(m,j),j=1,8)
              !    write (iow,*)
              !    write(iow,'(A)') '	new row'    
              !end do
          
!          write(iow,'(A)') 'Element displacements'
!		write(iow,*)				!write blank line
!          do m=1,8
!          write(iow,100) 'node =',m
!          write(iow,'(A,E12.4,E12.4,E12.4)') '(x,y,z) = ',xl(1,m), 
!     &      xl(2,m),xl(3,m)
!          write(iow,'(A,e18.5)') '	displ u1 ',ul(1,m,1)
!		write(iow,'(A,e18.5)') '	displ u2 ',ul(2,m,1)
!          write(iow,'(A,e18.5)') '	displ u3 ',ul(3,m,1)
!          write(iow,'(A,e18.5)') '	displ fi1 ',ul(4,m,1)
!		write(iow,'(A,e18.5)') '	displ fi2 ',ul(5,m,1)
!          write(iow,'(A,e18.5)') '	displ fi3 ',ul(6,m,1)
!          enddo!m
!          
!              
!          write(iow,'(A)') 'Element displacement INCREMENTS'
!		write(iow,*)				!write blank line
!          do m=1,nen
!          write(iow,100) 'node =',m
!          write(iow,'(A,E12.4,E12.4,E12.4)') '(x,y,z) = ',xl(1,m), 
!     &      xl(2,m),xl(3,m)
!          write(iow,'(A,e18.5)') '	DELTA u1 ',ul(1,m,3)
!		write(iow,'(A,e18.5)') '	DELTA u2 ',ul(2,m,3)
!          write(iow,'(A,e18.5)') '	DELTA u3 ',ul(3,m,3)
!          write(iow,'(A,e18.5)') '	DELTA fi1 ',ul(4,m,3)
!		write(iow,'(A,e18.5)') '	DELTA fi2 ',ul(5,m,3)
!          write(iow,'(A,e18.5)') '	DELTA fi3 ',ul(6,m,3) 
!          enddo!m
!
!!		write(iow,*)				!write blank line
!!		write(iow,*)				!write blank line
          
          
	endif

100	format(A,i5)
2001  format(
     & /5x,'T h r e e   D i m e n s i o n a l   M i c r o p o l a r   
     &  S o l i d   E l e m e n t'/)
      
      end
