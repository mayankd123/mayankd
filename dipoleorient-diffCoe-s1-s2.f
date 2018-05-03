c       program to calculate local dipole-orientational of solvent molecule around ions
c       This program is made Mayank Dixit at 1/12/2014
        implicit real*8(a-h,o-z)
c
        parameter (ns1=297,ns2=891,nmol=1190,itm=10001,jtmm=1)
        parameter (dionp=2.6,delt=0.001,nx=201,ny=101)
        parameter (box=5.01848,hbox=box/2.0)
        dimension pos(6,3,nmol,itm)
        parameter (ncor=(itm-1)/2)
        parameter (ncora=ncor-1000,mnt=500)  
        parameter (gpnaow=0.600,gpclow=0.600)        
        dimension dnawtt(ns1,ncora),dclwtt(ns1,ncora)
        dimension dnamht(ns2,ncora),dclmht(ns2,ncora)
        dimension rsdwt(ns1,ncora),rsdmh(ns2,ncora)
        dimension idxnaw(ns1,itm),idxclw(ns1,itm)
        dimension idxnam(ns2,itm),idxclm(ns2,itm)
        dimension pna(nx),pcl(nx),ppna(nx),ppcl(nx)
        dimension phina(nx),phicl(nx),dnas1t(ncora,jtmm)
        dimension dcls1t(ncora,jtmm),rnas1t(ncora,jtmm)
        dimension rcls1t(ncora,jtmm),conas1(ncora),cocls1(ncora)
        dimension dfnas1(ncora),dfcls1(ncora)
c
        dimension dnas2t(ncora,jtmm),rs1t(ncora,jtmm)
        dimension dcls2t(ncora,jtmm),rnas2t(ncora,jtmm)
        dimension rcls2t(ncora,jtmm),conas2(ncora),cocls2(ncora)
        dimension dfnas2(ncora),dfcls2(ncora),dfs1(ncora)
        dimension rs2t(ncora,jtmm),dfs2(ncora)
c
c input files
        open (21,file='1-10ps.dat')
c output files
        open (22,file='output1')
        open (51,file='dipole-na-s1.dat')
        open (52,file='dipole-cl-s1.dat')
        open (53,file='dipole-na-s2.dat')
        open (54,file='dipole-cl-s2.dat')
        open (55,file='difsion-na-s1.dat')
        open (56,file='difsion-cl-s1.dat')
        open (57,file='difsion-na-s2.dat')
        open (58,file='difsion-all.dat')
c end of files
c box= box length of simulated box and hbox=half of the box length of simulated box
c end of parameters
c
 80     format (I5,A5,A5,I5,3f8.3)
        do 9 j=1,itm
c
        do 10  m=1,ns1
        do 11 i=1,3
c
        read(21,80) kk,y,z,kk1,pos(i,1,m,j),pos(i,2,m,j),pos(i,3,m,j)
        write(22,80)kk,y,z,kk1,pos(i,1,m,j),pos(i,2,m,j),pos(i,3,m,j)
 11      continue
 10      continue
        do 26  m=ns1+1,nmol-2
        do 25 i=1,4
c   
        read(21,80) kk,y,z,kk1,pos(i,1,m,j),pos(i,2,m,j),pos(i,3,m,j)
        write(22,80)kk,y,z,kk1,pos(i,1,m,j),pos(i,2,m,j),pos(i,3,m,j)
 25      continue
 26      continue
c   
        do 33  m=nmol-1,nmol
c
        read(21,80) kk,y,z,kk1,pos(1,1,m,j),pos(1,2,m,j),pos(1,3,m,j)
        write(22,80)kk,y,z,kk1,pos(1,1,m,j),pos(1,2,m,j),pos(1,3,m,j)
 33      continue
c
        read(21,*)
        read(21,*)
        read(21,*)
c
        write(22,*)
        write(22,*)
        write(22,*)
c
 9       continue
c
c        call aorient_s1s2(pos,phina,phicl)
c***************************************************
c=======================================================================
        call dip_diff(ncora,dfnas1,dfnas2,dfcls1,dfcls2,dfs1,dfs2)
        do 78 i=1,ncora
        time=i*delt
        write(58,*)time,dfnas1(i),dfnas2(i),dfcls1(i),dfcls2(i),dfs1(i),
     6  dfs2(i)
c
 78    continue
c
      call coe_err(avgns1,stans1,avgns2,stans2,avgcs1,stacs1,
     6 avgcs2,stacs2)
c      call diff_coe_err 
      write(59,*) "avgns1,stans1,avgns2,stans2,avgcs1,stacs1,avgcs2,
     6 stacs2" 
      write(59,*) avgns1,stans1,avgns2,stans2,avgcs1,stacs1,avgcs2,
     6 stacs2

c
         stop
         end
c*******************************************
      subroutine s1_di_ornt(pos,idxnaw,idxclw)
        implicit real*8(a-h,o-z)
        parameter (ns1=297,ns2=891,nmol=1190,itm=10001)
        parameter (delt=0.001,box=5.01848,hbox=box/2.0)
        parameter (c1=-0.70,c2=0.265,c3=0.435)
        parameter (gpnaow=0.600,gpclow=0.600,nx=201,ny=101)
        dimension pos(6,3,nmol,itm)
        dimension idxnaw(ns1,itm),idxclw(ns1,itm)
        dimension rnacos(ns1,itm),rclcos(ns1,itm)
        parameter (ncor=(itm-1)/2)
        parameter (ncora=ncor-1000)
c
      dimension dna0(ns1),dnaag0(ns1),dnat(ns1,ncora),dnawtt(ns1,ncora)
      dimension dcl0(ns1),dclag0(ns1),dclt(ns1,ncora),dclwtt(ns1,ncora)
      dimension pna(nx),pcl(nx)
c
        aechrg=1.609*(10**-19.0)
        Deby=3.33564*10**(-30.0) 
c
        do 12 j=1,itm
c
c
        dxce=(pos(1,1,nmol,j)+pos(1,1,nmol-1,j))/2
        dyce=(pos(1,2,nmol,j)+pos(1,2,nmol-1,j))/2
        dzce=(pos(1,3,nmol,j)+pos(1,3,nmol-1,j))/2
c
         if(abs(dxce).gt.hbox)then
         dxce=dxce-sign(box,dxce)
         endif
         if(abs(dyce).gt.hbox)then
         dyce=dyce-sign(box,dyce)
         endif
         if(abs(dzce).gt.hbox)then
         dzce=dzce-sign(box,dzce)
         endif
c*********end of PBC********************************************
c
c distance and unit vector of centerofmass-Clion
c
        dxcecl=pos(1,1,nmol,j)-dxce
        dycecl=pos(1,2,nmol,j)-dyce
        dzcecl=pos(1,3,nmol,j)-dzce
c
         if(abs(dxcecl).gt.hbox)then
         dxcecl=dxcecl-sign(box,dxcecl)
         endif
         if(abs(dycecl).gt.hbox)then
         dycecl=dycecl-sign(box,dycecl)
         endif
         if(abs(dzcecl).gt.hbox)then
         dzcecl=dzcecl-sign(box,dzcecl)
         endif
c
        rcecl=sqrt((dxcecl)**2+(dycecl)**2+(dzcecl)**2)
        rcecl=abs(rcecl)
c
        uxcecl=dxcecl/rcecl
        uycecl=dycecl/rcecl
        uzcecl=dzcecl/rcecl
c end of  distance and unit vector of centerofmass-Clion
c
c distance between ions and water molecules
c
        do 13 mm=1,ns1
c***********************************************************
        dxnaow=pos(1,1,mm,j)-pos(1,1,nmol-1,j)
        dynaow=pos(1,2,mm,j)-pos(1,2,nmol-1,j)
        dznaow=pos(1,3,mm,j)-pos(1,3,nmol-1,j)
c
         if(abs(dxnaow).gt.hbox)then
         dxnaow=dxnaow-sign(box,dxnaow)
         endif
         if(abs(dynaow).gt.hbox)then
         dynaow=dynaow-sign(box,dynaow)
         endif
         if(abs(dznaow).gt.hbox)then
         dznaow=dznaow-sign(box,dznaow)
         endif
c*********end of PBC********************************************
        rnaow=sqrt((dxnaow)**2+(dynaow)**2+(dznaow)**2)
        rnaow=abs(rnaow)
        uxnaow=dxnaow/rnaow
        uynaow=dynaow/rnaow
        uznaow=dznaow/rnaow
c**********************************************************
c**********************************************************
        dxclow=pos(1,1,mm,j)-pos(1,1,nmol,j)
        dyclow=pos(1,2,mm,j)-pos(1,2,nmol,j)
        dzclow=pos(1,3,mm,j)-pos(1,3,nmol,j)
c
         if(abs(dxclow).gt.hbox)then
         dxclow=dxclow-sign(box,dxclow)
         endif
         if(abs(dyclow).gt.hbox)then
         dyclow=dyclow-sign(box,dyclow)
         endif
         if(abs(dzclow).gt.hbox)then
         dzclow=dzclow-sign(box,dzclow)
         endif
c*********end of PBC********************************************
        rclow=sqrt((dxclow)**2+(dyclow)**2+(dzclow)**2)
        rclow=abs(rclow)
c**********************************************************
         if (rnaow.lt.gpnaow)then
         idxnaw(mm,j)=mm
         endif
         if (rclow.lt.gpclow)then
         idxclw(mm,j)=mm
         endif
c

 13     continue
c
 12     continue
c
c
        return
        end
c****************************************************
c
c*****************************************************
      subroutine s2_di_ornt(pos,idxnam,idxclm)
        implicit real*8(a-h,o-z)
        parameter (ns1=297,ns2=891,nmol=1190,itm=10001)
        parameter (nx=201,ny=101)
        parameter (box=5.01848,hbox=box/2.0)
        parameter (c1=0.1390,c2=-0.459,c3=0.160,c4=0.160)
        parameter (gpnaom=0.600,gpclom=0.600)
        dimension pos(6,3,nmol,itm)
        dimension idxnam(ns2,itm),idxclm(ns2,itm)
        dimension rnacos(ns2,itm),rclcos(ns2,itm)
        parameter (ncor=(itm-1)/2)
        parameter (ncora=ncor-1000)
c
      dimension dna0(ns2),dnaag0(ns2),dnat(ns2,ncora)
      dimension dnamht(ns2,ncora)
      dimension dcl0(ns2),dclag0(ns2),dclt(ns2,ncora)
      dimension dclmht(ns2,ncora),ppna(nx),ppcl(nx)
c
        aechrg=1.609*(10**-19.0)
        Deby=3.33564*10**(-30.0)
c
        do 12 j=1,itm
c
c
        dxce=(pos(1,1,nmol,j)+pos(1,1,nmol-1,j))/2
        dyce=(pos(1,2,nmol,j)+pos(1,2,nmol-1,j))/2
        dzce=(pos(1,3,nmol,j)+pos(1,3,nmol-1,j))/2
c        write(2,*) dxce,dyce,dzce,j
c
         if(abs(dxce).gt.hbox)then
         dxce=dxce-sign(box,dxce)
         endif
         if(abs(dyce).gt.hbox)then
         dyce=dyce-sign(box,dyce)
         endif
         if(abs(dzce).gt.hbox)then
         dzce=dzce-sign(box,dzce)
         endif
c        write(2,*) dxce,dyce,dzce,j
c*********end of PBC********************************************
c
c
        dxcecl=pos(1,1,nmol,j)-dxce
        dycecl=pos(1,2,nmol,j)-dyce
        dzcecl=pos(1,3,nmol,j)-dzce
c
         if(abs(dxcecl).gt.hbox)then
         dxcecl=dxcecl-sign(box,dxcecl)
         endif
         if(abs(dycecl).gt.hbox)then
         dycecl=dycecl-sign(box,dycecl)
         endif
         if(abs(dzcecl).gt.hbox)then
         dzcecl=dzcecl-sign(box,dzcecl)
         endif
c
        rcecl=sqrt((dxcecl)**2+(dycecl)**2+(dzcecl)**2)
        rcecl=abs(rcecl)
c
        uxcecl=dxcecl/rcecl
        uycecl=dycecl/rcecl
        uzcecl=dzcecl/rcecl
c
        do 13 mm=ns1+1,nmol-2
c*****************************************************************
        dxnaow=pos(1,1,mm,j)-pos(1,1,nmol-1,j)
        dynaow=pos(1,2,mm,j)-pos(1,2,nmol-1,j)
        dznaow=pos(1,3,mm,j)-pos(1,3,nmol-1,j)
c
         if(abs(dxnaow).gt.hbox)then
         dxnaow=dxnaow-sign(box,dxnaow)
         endif
         if(abs(dynaow).gt.hbox)then
         dynaow=dynaow-sign(box,dynaow)
         endif
         if(abs(dznaow).gt.hbox)then
         dznaow=dznaow-sign(box,dznaow)
         endif
c*********end of PBC********************************************
        rnaow=sqrt((dxnaow)**2+(dynaow)**2+(dznaow)**2)
        rnaow=abs(rnaow)
        uxnaow=dxnaow/rnaow
        uynaow=dynaow/rnaow
        uznaow=dznaow/rnaow
c        write(2,*) dxnaow,dynaow,dznaow,rnaow,mm
c******************************************************************
c******************************************************************
        dxclow=pos(1,1,mm,j)-pos(1,1,nmol,j)
        dyclow=pos(1,2,mm,j)-pos(1,2,nmol,j)
        dzclow=pos(1,3,mm,j)-pos(1,3,nmol,j)
c
         if(abs(dxclow).gt.hbox)then
         dxclow=dxclow-sign(box,dxclow)
         endif
         if(abs(dyclow).gt.hbox)then
         dyclow=dyclow-sign(box,dyclow)
         endif
         if(abs(dzclow).gt.hbox)then
         dzclow=dzclow-sign(box,dzclow)
         endif
c*********end of PBC********************************************
        rclow=sqrt((dxclow)**2+(dyclow)**2+(dzclow)**2)
        rclow=abs(rclow)
        uxclow=dxclow/rclow
        uyclow=dyclow/rclow
        uzclow=dzclow/rclow
c        write(2,*) dxclow,dyclow,dzclow,rclow,mm
c*******************************************************************
c
c

         if (rnaow.lt.gpnaom)then
c
         idxnam(mm-ns1,j)=mm-ns1
         endif
         if (rclow.lt.gpclom)then
         idxclm(mm-ns1,j)=mm-ns1
         endif
c
 13     continue
c
 12     continue
c
        return
        end
c===========================================
c subroutine to calculate diffusion constant
c===========================================
c********************************************
        subroutine diffusion(pos,rsdwt,rsdmh)
        implicit real*8(a-h,o-z)
        parameter (ns1=297,ns2=891,nmol=1190,itm=10001)
        parameter (box=5.01848,hbox=box/2.0)
        parameter (delt=0.001)
        dimension pos(6,3,nmol,itm)
        dimension idxnaw(ns1,itm),idxclw(ns1,itm)
c
        parameter (ncor=(itm-1)/2)
        parameter (ncora=ncor-1000)
        dimension rsdwt(ns1,ncora),rsdmh(ns2,ncora)
c============================================================
c calculation of diffusion coefficient of water
        do mm=1,ns1
c
        do i=1,ncora
        sum=0.0

c
        do j=1001,ncor
        disp1=(pos(1,1,mm,i+j)-pos(1,1,mm,j))**2 
        disp2=(pos(1,2,mm,i+j)-pos(1,2,mm,j))**2
        disp3=(pos(1,3,mm,i+j)-pos(1,3,mm,j))**2 
        rsq=disp1+disp2+disp3 
        sum=sum+rsq
        enddo
        rsdwt(mm,i)=sum/(ncora)
        enddo
        enddo
c
c       do mm=1,ns1
c       do i=10,ncor/2,10
c       time=float(i)*delt
c       write(71,*)time,rsdwt(mm,i)
c       enddo
c       enddo
c end of calculation of diffusion coefficeint of water 
c=============================================================
c calculation of diffusion coefficient of water
        do mm=ns1+1,nmol-2
c
        do i=1,ncora
        sum=0.0
c
        do j=1001,ncor
        disp1=(pos(1,1,mm,i+j)-pos(1,1,mm,j))**2
        disp2=(pos(1,2,mm,i+j)-pos(1,2,mm,j))**2
        disp3=(pos(1,3,mm,i+j)-pos(1,3,mm,j))**2
        rsq=disp1+disp2+disp3
        sum=sum+rsq
        enddo
        rsdmh(mm-ns1,i)=sum/(ncora)
        enddo
        enddo
c
c       do mm=1,nmeoh
c       do i=10,ncor/2,10
c       time=float(i)*delt
c       write(72,*)time,rsdmh(mm,i)
c       enddo
c       enddo
c end of calculation of diffusion coefficeint of water
c==============================================================
c==============================================================
       return
       end 
c
c        subroutine aorient_s1s2(pos,phina,phicl)
c        implicit real*8(a-h,o-z)
c        parameter (ns1=297,ns2=891,nmol=1190,itm=10001)
c        parameter (box=5.01848,hbox=box/2.0)
c        parameter (gpnaow=0.600,gpclow=0.600,nx=201,ny=101)
c        dimension pos(6,3,nmol,itm)
c
c      dimension uxnaow(nmol-2,itm),uynaow(nmol-2,itm),uznaow(nmol-2,itm)
c      dimension uxclow(nmol-2,itm),uyclow(nmol-2,itm),uzclow(nmol-2,itm)
c      dimension phina(nx),phicl(nx)
c      dimension rnaow(nmol-2,itm),rclow(nmol-2,itm)
c
c
c        do i=0,nx
c        phina(i)=0.0
c        phicl(i)=0.0
c        enddo
c
c        do 12 j=1,itm
c
c
c        dxce=(pos(1,1,nmol,j)+pos(1,1,nmol-1,j))/2
c        dyce=(pos(1,2,nmol,j)+pos(1,2,nmol-1,j))/2
c        dzce=(pos(1,3,nmol,j)+pos(1,3,nmol-1,j))/2
c
c         if(abs(dxce).gt.hbox)then
c         dxce=dxce-sign(box,dxce)
c         endif
c         if(abs(dyce).gt.hbox)then
c         dyce=dyce-sign(box,dyce)
c         endif
c         if(abs(dzce).gt.hbox)then
c         dzce=dzce-sign(box,dzce)
c         endif
c*********end of PBC********************************************
c
c distance and unit vector of centerofmass-Clion
c
c        dxcecl=pos(1,1,nmol,j)-dxce
c        dycecl=pos(1,2,nmol,j)-dyce
c        dzcecl=pos(1,3,nmol,j)-dzce
c
c         if(abs(dxcecl).gt.hbox)then
c         dxcecl=dxcecl-sign(box,dxcecl)
c         endif
c         if(abs(dycecl).gt.hbox)then
c         dycecl=dycecl-sign(box,dycecl)
c         endif
c         if(abs(dzcecl).gt.hbox)then
c         dzcecl=dzcecl-sign(box,dzcecl)
c         endif
c
c        rcecl=sqrt((dxcecl)**2+(dycecl)**2+(dzcecl)**2)
c        rcecl=abs(rcecl)
c
c        uxcecl=dxcecl/rcecl
c        uycecl=dycecl/rcecl
c        uzcecl=dzcecl/rcecl
c end of  distance and unit vector of centerofmass-Clion
c
c distance between ions and water molecules
c
c        do 13 mm=1,nmol-2
c***********************************************************
c        dxnaow=pos(1,1,mm,j)-pos(1,1,nmol-1,j)
c        dynaow=pos(1,2,mm,j)-pos(1,2,nmol-1,j)
c        dznaow=pos(1,3,mm,j)-pos(1,3,nmol-1,j)
c
c         if(abs(dxnaow).gt.hbox)then
c         dxnaow=dxnaow-sign(box,dxnaow)
c         endif
c         if(abs(dynaow).gt.hbox)then
c         dynaow=dynaow-sign(box,dynaow)
c         endif
c         if(abs(dznaow).gt.hbox)then
c         dznaow=dznaow-sign(box,dznaow)
c         endif
c*********end of PBC********************************************
c        rnaow(mm,j)=sqrt((dxnaow)**2+(dynaow)**2+(dznaow)**2)
c        rnaow(mm,j)=abs(rnaow(mm,j))
c        uxnaow(mm,j)=dxnaow/rnaow(mm,j)
c        uynaow(mm,j)=dynaow/rnaow(mm,j)
c        uznaow(mm,j)=dznaow/rnaow(mm,j)
c**********************************************************
c**********************************************************
c        dxclow=pos(1,1,mm,j)-pos(1,1,nmol,j)
c        dyclow=pos(1,2,mm,j)-pos(1,2,nmol,j)
c        dzclow=pos(1,3,mm,j)-pos(1,3,nmol,j)
c
c         if(abs(dxclow).gt.hbox)then
c         dxclow=dxclow-sign(box,dxclow)
c         endif
c         if(abs(dyclow).gt.hbox)then
c         dyclow=dyclow-sign(box,dyclow)
c         endif
c         if(abs(dzclow).gt.hbox)then
c         dzclow=dzclow-sign(box,dzclow)
c         endif
c*********end of PBC********************************************
c        rclow(mm,j)=sqrt((dxclow)**2+(dyclow)**2+(dzclow)**2)
c        rclow(mm,j)=abs(rclow(mm,j))
c        uxclow(mm,j)=dxclow/rclow(mm,j)
c        uyclow(mm,j)=dyclow/rclow(mm,j)
c        uzclow(mm,j)=dzclow/rclow(mm,j)
c**********************************************************
c
c
c 13     continue
c
c 12     continue
c
c          do j=1,itm   
c          do mm=1,nmol-3
c          do nn=mm+1,nmol-2
c         if ((rnaow(mm,j).lt.gpnaow).and.(rnaow(nn,j).lt.gpnaow))then
c         theta2=acos(uxnaow(mm,j)*uxnaow(nn,j)+uynaow(mm,j)*uynaow(nn,j)
c     1    +uznaow(mm,j)*uznaow(nn,j))
c         ina=int(cos(theta2)/0.01)+1
c         ina=ina+ny
c         phina(ina)=phina(ina)+1
c         write(41,*)"ina,cos(theta2),phina(ina)"
c         write(41,*)ina,cos(theta2),phina(ina)
c         endif
c         if ((rclow(mm,j).lt.gpclow).and.(rclow(nn,j).lt.gpclow))then
c         theta3=acos(uxclow(mm,j)*uxclow(nn,j)+uyclow(mm,j)*uyclow(nn,j)
c     1    +uzclow(mm,j)*uzclow(nn,j))
c         icl=int(cos(theta3)/0.01)+1
c         icl=icl+ny
c         phicl(icl)=phicl(icl)+1
c         write(42,*)"icl,cos(theta3),phicl(icl)"
c         write(42,*)icl,cos(theta3),phicl(icl)
c
c         endif
c         enddo
c         enddo
c         enddo 
c
cangular distribution function
c        do i=0,nx
c        phina(i)=phina(i)/itm
c        phicl(i)=phicl(i)/itm
c        enddo
c        write(43,*)"i,phina(ina),phicl(ina)"
c        do i=0,nx
c         ii=i-ny
c         ai=real(ii/100.0)
c        write(43,*)ai,phina(i),phicl(i)
c        enddo
c
c
c        return
c        end

c=======================
c subroutine to calculate the dipole-corelation and diffusional behaviour 
c=============================
      subroutine dip_diff(ncora,dfnas1,dfnas2,dfcls1,dfcls2,dfs1,dfs2)
        implicit real*8(a-h,o-z)
c
        parameter (ns1=297,ns2=891,nmol=1190,itm=10001,jtmm=1)
        parameter (dionp=2.6,delt=0.001,nx=201,ny=101)
        parameter (box=5.01848,hbox=box/2.0)
        dimension pos(6,3,nmol,itm)
        parameter (ncor=(itm-1)/2)
        parameter (gpnaow=0.600,gpclow=0.600)
        dimension dnawtt(ns1,ncora),dclwtt(ns1,ncora)
        dimension dnamht(ns2,ncora),dclmht(ns2,ncora)
        dimension rsdwt(ns1,ncora),rsdmh(ns2,ncora)
        dimension idxnaw(ns1,itm),idxclw(ns1,itm)
        dimension idxnam(ns2,itm),idxclm(ns2,itm)
        dimension pna(nx),pcl(nx),ppna(nx),ppcl(nx)
        dimension phina(nx),phicl(nx),dnas1t(ncora,jtmm)
        dimension dcls1t(ncora,jtmm),rnas1t(ncora,jtmm)
        dimension rcls1t(ncora,jtmm)
        dimension dfnas1(ncora),dfcls1(ncora)
c
        dimension dnas2t(ncora,jtmm),rs1t(ncora,jtmm)
        dimension dcls2t(ncora,jtmm),rnas2t(ncora,jtmm)
        dimension rcls2t(ncora,jtmm)
        dimension dfnas2(ncora),dfcls2(ncora),dfs1(ncora)
        dimension rs2t(ncora,jtmm),dfs2(ncora)
c
       call s1_di_ornt(pos,idxnaw,idxclw)
       call s2_di_ornt(pos,idxnam,idxclm)
        call diffusion(pos,rsdwt,rsdmh)
c***************************************************
c=======================================================================
c
        do k=1,ncora
        dfnas1(i)=0.0
        dfcls1(i)=0.0
        dfs1(i)=0.0
        do jj=1,jtmm
        rnas1t(i,jj)=0.0
        rcls1t(i,jj)=0.0
        rs1t(i,jj)=0.0
        enddo
        enddo
c
        do 77 i=1,ncora
        do jj=1,jtmm
         ina=0
         icl=0
         do mm=1,ns1
        rs1t(i,jj)=rs1t(i,jj)+rsdwt(mm,i)
        if(mm.eq.idxnaw(mm,jj))then
        ina=ina+1
        rnas1t(i,jj)=rnas1t(i,jj)+rsdwt(mm,i)
        endif
        if(mm.eq.idxclw(mm,jj))then
        icl=icl+1
        rcls1t(i,jj)=rcls1t(i,jj)+rsdwt(mm,i)
        endif
c
        enddo
c
        rnas1t(i,jj)=rnas1t(i,jj)/ina
        rcls1t(i,jj)=rcls1t(i,jj)/icl
c
        rs1t(i,jj)=rs1t(i,jj)/ns1
c
        dfnas1(i)=dfnas1(i)+rnas1t(i,jj)
        dfcls1(i)=dfcls1(i)+rcls1t(i,jj)
        dfs1(i)=dfs1(i)+rs1t(i,jj)
        enddo
c
        dfnas1(i)=dfnas1(i)/jtmm
        dfcls1(i)=dfcls1(i)/jtmm
        dfs1(i)=dfs1(i)/jtmm
c

 77    continue
c
        do k=1,ncora
        dfnas2(i)=0.0
        dfcls2(i)=0.0
        dfs2(i)=0.0
        do jj=1,jtmm
        rnas2t(i,jj)=0.0
        rcls2t(i,jj)=0.0
        rs2t(i,jj)=0.0
        enddo
        enddo
        do 78 i=1,ncora
        do jj=1,jtmm
         ina=0
         icl=0
         do mm=1,ns2
        rs2t(i,jj)=rs2t(i,jj)+rsdmh(mm,i)
        if(mm.eq.idxnam(mm,jj))then
        ina=ina+1
        rnas2t(i,jj)=rnas2t(i,jj)+rsdmh(mm,i)
        endif
        if(mm.eq.idxclm(mm,jj))then
        icl=icl+1
        rcls2t(i,jj)=rcls2t(i,jj)+rsdmh(mm,i)
        endif
c
        enddo
        rnas2t(i,jj)=rnas2t(i,jj)/ina
        rcls2t(i,jj)=rcls2t(i,jj)/icl
        rs2t(i,jj)=rs2t(i,jj)/ns2
        dfnas2(i)=dfnas2(i)+rnas2t(i,jj)
        dfcls2(i)=dfcls2(i)+rcls2t(i,jj)
        dfs2(i)=dfs2(i)+rs2t(i,jj)
        enddo
        dfnas2(i)=dfnas2(i)/jtmm
        dfcls2(i)=dfcls2(i)/jtmm
        dfs2(i)=dfs2(i)/jtmm
      write(61,*) i,dfnas1(i),dfcls1(i),dfnas2(i),dfcls2(i),dfs1(i)
     1   ,dfs2(i)
c
   78    continue

c
       return
       end 
c==========================================================
c subroutine to calculate the diffusion coefficients and error bars
c==================================================================
         subroutine coe_err(avgns1,stans1,avgns2,stans2,avgcs1,stacs1,
     1 avgcs2,stacs2)
c       diff_coe_err
        implicit real*8(a-h,o-z)
c
        parameter (ns1=297,ns2=891,nmol=1190,itm=10001,jtmm=1)
        parameter (dionp=2.6,delt=0.001,nx=201,ny=101,mnt=500)
        parameter (box=5.01848,hbox=box/2.0)
        dimension pos(6,3,nmol,itm)
        parameter (ncor=(itm-1)/2)
        parameter (ncora=ncor-1000)
c        parameter (ncora=ncor/2)
        parameter (gpnaow=0.600,gpclow=0.600)
        dimension dnawtt(ns1,ncora),dclwtt(ns1,ncora)
        dimension dnamht(ns2,ncora),dclmht(ns2,ncora)
        dimension rsdwt(ns1,ncora),rsdmh(ns2,ncora)
        dimension idxnaw(ns1,itm),idxclw(ns1,itm)
        dimension idxnam(ns2,itm),idxclm(ns2,itm)
        dimension pna(nx),pcl(nx),ppna(nx),ppcl(nx)
        dimension phina(nx),phicl(nx),dnas1t(ncora,jtmm)
        dimension dcls1t(ncora,jtmm),rnas1t(ncora,jtmm)
        dimension rcls1t(ncora,jtmm),conas1(ncora),cocls1(ncora)
        dimension dfnas1(ncora),dfcls1(ncora)
c
        dimension dnas2t(ncora,jtmm),rs1t(ncora,jtmm)
        dimension dcls2t(ncora,jtmm),rnas2t(ncora,jtmm)
        dimension rcls2t(ncora,jtmm),conas2(ncora),cocls2(ncora)
        dimension dfnas2(ncora),dfcls2(ncora),dfs1(ncora)
        dimension rs2t(ncora,jtmm),dfs2(ncora)
c
         call dip_diff(ncora,dfnas1,dfnas2,dfcls1,dfcls2,dfs1,dfs2)
c
         ad1ns1=(dfnas1(ncora)-dfnas1(ncora-mnt))/(6.0*mnt*delt)
         ad2ns1=(dfnas1(ncora-mnt)-dfnas1(ncora-2*mnt))/(6.0*mnt*delt)
         ad3ns1=(dfnas1(ncora-2*mnt)-dfnas1(ncora-3*mnt))/(6.0*mnt*delt)
         ad4ns1=(dfnas1(ncora-3*mnt)-dfnas1(ncora-4*mnt))/(6.0*mnt*delt)
         ad5ns1=(dfnas1(ncora-4*mnt)-dfnas1(ncora-5*mnt))/(6.0*mnt*delt)
         ad6ns1=(dfnas1(ncora-5*mnt)-dfnas1(ncora-6*mnt))/(6.0*mnt*delt)
c
         avgns1= (ad1ns1+ad2ns1+ad3ns1+ad4ns1+ad5ns1+ad6ns1)/6.0
         varns1=((ad1ns1-avgns1)**2+(ad2ns1-avgns1)**2+(ad3ns1-avgns1)**2+
     6  (ad4ns1-avgns1)**2+(ad5ns1-avgns1)**2+(ad6ns1-avgns1)**2)/6.0
         stans1=sqrt(varns1)
c         write(59,*) "avgns1,tand"
c         write(59,*) avgns1,stand
c
         ad1ns2=(dfnas2(ncora)-dfnas2(ncora-mnt))/(6.0*mnt*delt)
         ad2ns2=(dfnas2(ncora-mnt)-dfnas2(ncora-2*mnt))/(6.0*mnt*delt)
         ad3ns2=(dfnas2(ncora-2*mnt)-dfnas2(ncora-3*mnt))/(6.0*mnt*delt)
         ad4ns2=(dfnas2(ncora-3*mnt)-dfnas2(ncora-4*mnt))/(6.0*mnt*delt)
         ad5ns2=(dfnas2(ncora-4*mnt)-dfnas2(ncora-5*mnt))/(6.0*mnt*delt)
         ad6ns2=(dfnas2(ncora-5*mnt)-dfnas2(ncora-6*mnt))/(6.0*mnt*delt)
c
         avgns2= (ad1ns2+ad2ns2+ad3ns2+ad4ns2+ad5ns2+ad6ns2)/6.0
         varns2=((ad1ns2-avgns2)**2+(ad2ns2-avgns2)**2+(ad3ns2-avgns2)**2+
     6  (ad4ns2-avgns2)**2+(ad5ns2-avgns2)**2+(ad6ns2-avgns2)**2)/6.0
         stans2=sqrt(varns2)
c
c
         ad1cs1=(dfcls1(ncora)-dfcls1(ncora-mnt))/(6.0*mnt*delt)
         ad2cs1=(dfcls1(ncora-mnt)-dfcls1(ncora-2*mnt))/(6.0*mnt*delt)
         ad3cs1=(dfcls1(ncora-2*mnt)-dfcls1(ncora-3*mnt))/(6.0*mnt*delt)
         ad4cs1=(dfcls1(ncora-3*mnt)-dfcls1(ncora-4*mnt))/(6.0*mnt*delt)
         ad5cs1=(dfcls1(ncora-4*mnt)-dfcls1(ncora-5*mnt))/(6.0*mnt*delt)
         ad6cs1=(dfcls1(ncora-5*mnt)-dfcls1(ncora-6*mnt))/(6.0*mnt*delt)
c
         avgcs1= (ad1cs1+ad2cs1+ad3cs1+ad4cs1+ad5cs1+ad6cs1)/6.0
         varcs1=((ad1cs1-avgcs1)**2+(ad2cs1-avgcs1)**2+(ad3cs1-avgcs1)**2+
     6  (ad4cs1-avgcs1)**2+(ad5cs1-avgcs1)**2+(ad6cs1-avgcs1)**2)/6.0
         stacs1=sqrt(varcs1)
         write(60,*) "avgns1,tand"
         write(60,*) avgns1,stand
     
         ad1cs2=(dfcls2(ncora)-dfcls2(ncora-mnt))/(6.0*mnt*delt)
         ad2cs2=(dfcls2(ncora-mnt)-dfcls2(ncora-2*mnt))/(6.0*mnt*delt)
         ad3cs2=(dfcls2(ncora-2*mnt)-dfcls2(ncora-3*mnt))/(6.0*mnt*delt)
         ad4cs2=(dfcls2(ncora-3*mnt)-dfcls2(ncora-4*mnt))/(6.0*mnt*delt)
         ad5cs2=(dfcls2(ncora-4*mnt)-dfcls2(ncora-5*mnt))/(6.0*mnt*delt)
         ad6cs2=(dfcls2(ncora-5*mnt)-dfcls2(ncora-6*mnt))/(6.0*mnt*delt)
c
         avgcs2= (ad1cs2+ad2cs2+ad3cs2+ad4cs2+ad5cs2+ad6cs2)/6.0
         varcs2=((ad1cs2-avgcs2)**2+(ad2cs2-avgcs2)**2+(ad3cs2-avgcs2)**2+
     6  (ad4cs2-avgcs2)**2+(ad5cs2-avgcs2)**2+(ad6cs2-avgcs2)**2)/6.0
         stacs2=sqrt(varcs2)
         write(60,*) "avgns1,tand"
         write(60,*) avgns1,stand
c

      return 
      end
