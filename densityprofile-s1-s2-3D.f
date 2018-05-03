c       program to calculate density profile of solvent molecules around ion pair
c       This program is made by Mayank Dixit 9/9/2014
        implicit real*8(a-h,o-z)
c
        parameter (ns1=297,ns2=891,nmol=1190,itm=1001)
        parameter (nx=201,ny=101,dionp=2.6)
        parameter (box=5.03006,hbox=box/2.0)
        parameter (aclpop=2.0,anapop=2.0)
        dimension pos(4,3,nmol,itm)
        dimension bins1(nx,ny),bins2(nx,ny)
c******************************************************
c diop is the distance between ions
c ns1 is the number of solvent-1 molecule
c ns2 is the number of solvent-2 molecule
c ns1+ns2+2=nmol, 2=1=cation+1=anion
c itm is the number of time frames
c box is the box length of the box, hbox is the half of the box 
c aionpop population of cation and anion 
c pos is position array of all the sites
c bins1 is density array of solvent-1
c bins2 is density array of solvent-2
c******************************************************
c input files
        open (21,file='input1000.dat')
c output files
        open (22,file='output')
        open (51,file='s1-density-3D.dat')
        open (52,file='s2-density-3D.dat')
c end of files
c
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
        hdionp=dionp/2.0
        hdionp=hdionp/0.1
c ice is the center of the box
        ice=int(hdionp)+1
        jcent=(nx+1)/2
c icepls is the position of the cation and icemin is the position of the anion
c
        icepls=jcent+ice
        icemin= jcent-ice
        write(*,*)"dionp,hdionp,ice,icepls,icemin"
        write(*,*)dionp,hdionp,ice,jcent,icepls,icemin  
c
        call dense_s1(pos,bins1)
        call dense_s2(pos,bins2)
c************************************************
c assigning the position of cation and anion for solvent-1
        jj=1
          do 60 i=icepls-2, icepls+2
        bins1(i,jj)=aclpop
 60    continue
          ii=icepls
          do 61 j=1,3 
        bins1(ii,j)=aclpop
 61    continue
c        
          do 62 i=icemin-2,icemin+2
        bins1(i,jj)=anapop
 62    continue
c
c end of assigning the position of cation and anion
c*****************************************************
        do i=0,nx
        do j=-ny,ny
c the density profile of solvent-1 are wrriten in the file 51
        if (j.lt.0.0)then
        k=-j
c
        write(51,*) i,j,bins1(i,k)
        endif
        if (j.gt.0.0)then
        write(51,*) i,j,bins1(i,j)
        endif
        enddo
        enddo
c
c*********************************************************
c assigning the position of cation and anion for solvent-2
          do 63 i=icepls-2, icepls+2
        bins2(i,jj)=aclpop
 63    continue
          ii=icepls
          do 64 j=1,3
        bins2(ii,j)=aclpop
 64    continue
c
          do 65 i=icemin-2,icemin+2
        bins2(i,jj)=anapop
 65    continue
c end of assigning the position of cation and anion for solvent-2
c****************************************
        do i=0,nx    
        do j=-ny,ny
c the density profile of solvent-2 are wrriten in the file 52
        if (j.lt.0.0)then
        k=-j
        write(52,*) i,j,bins2(i,k)
        endif
c
        if (j.gt.0.0)then
        write(52,*) i,j,bins2(i,j)
        endif
        enddo
        enddo
c
         stop
         end
c*******************************************
c dense_s1 is the subroutine to calculate the density of solvent-1 around the ion pair
c******************************************** 
        subroutine dense_s1(pos,bins1)
        implicit real*8(a-h,o-z)
        parameter (ns1=297,ns2=891,nmol=1190,itm=1001)
        parameter (nx=201,ny=101)
        parameter (box=5.03006,hbox=box/2.0)
        dimension pos(4,3,nmol,itm)
        dimension rcoscl(ns1),rsincl(ns1)
        dimension bins1(nx,ny)
c
        do ii=0,nx
        do jj=0,ny
c ii refers to rcos(theta) and jj refers to rsin(theta) of the distance of solvent-1 from center of mass of ion pair (cylindrical coordinate)
        bins1(ii,jj)=0.0
        enddo
        enddo
c
         write(2,*) "mm,rceow,theta1,rcoscl(mm),rsincl(mm),ik1,jj"
        do 12 j=1,itm
c
c      coordinates of center of mass
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
c        write(2,*) dxce,dyce,dzce,j
c*********end of PBC********************************************
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
c        write(2,*) dxcecl,dycecl,dzcecl,rcecl,j
c
        uxcecl=dxcecl/rcecl
        uycecl=dycecl/rcecl
        uzcecl=dzcecl/rcecl
c
        do 13 mm=1,ns1
        dxceow=pos(1,1,mm,j)-dxce
        dyceow=pos(1,2,mm,j)-dyce
        dzceow=pos(1,3,mm,j)-dzce
c
         if(abs(dxceow).gt.hbox)then
         dxceow=dxceow-sign(box,dxceow)
         endif
         if(abs(dyceow).gt.hbox)then
         dyceow=dyceow-sign(box,dyceow)
         endif
         if(abs(dzceow).gt.hbox)then
         dzceow=dzceow-sign(box,dzceow)
         endif
c*********end of PBC********************************************
        rceow=sqrt((dxceow)**2+(dyceow)**2+(dzceow)**2)
        rceow=abs(rceow)
        uxceow=dxceow/rceow
        uyceow=dyceow/rceow
        uzceow=dzceow/rceow
c
         theta1=acos(uxceow*uxcecl+uyceow*uycecl+uzceow*uzcecl)
c
         rcoscl(mm)=rceow*cos(theta1)
         rsincl(mm)=rceow*sin(theta1)
          if ((rcoscl(mm).lt.1.0).and.(rcoscl(mm).gt.-1.0)) then
          if (rsincl(mm).lt.1.0)then
c         write(2,*) mm,rceow,theta1,rcoscl(mm),rsincl(mm)
c
         arcscl=rcoscl(mm)/0.01
         arsicl=rsincl(mm)/0.01
c
         ii=int(arcscl)+1
         jj=int(arsicl)+1
         ik1=ny-ii
c
         bins1(ik1,jj)=bins1(ik1,jj)+1.0
         write(2,*) mm,rceow,theta1,rcoscl(mm),rsincl(mm),ik1,jj
         endif
         endif
 13     continue
c
 12     continue
        do i=0,nx
        do j=0,ny
        bins1(i,j)=bins1(i,j)/itm
c bins1=bins1/itm is the time avrage of bins1
c
c bins1 is multiplied by 1000 to increase the intensity the user can change this number
c
        bins1(i,j)=bins1(i,j)*1000
        enddo
        enddo
c binwt is the number density of water around ion pair
        return
        end
c****************************************************
c end of dense_s1 
c*****************************************************
c*******************************************
c dense_s2 is the subroutine to calculate the density of solvent-2 around the ion pair
c********************************************
         subroutine dense_s2(pos,bins2)
        implicit real*8(a-h,o-z)
        parameter (ns1=297,ns2=891,nmol=1190,itm=1001)
        parameter (nx=201,ny=101)
        parameter (box=5.03006,hbox=box/2.0)
        dimension pos(4,3,nmol,itm)
        dimension rcoscl(ns2),rsincl(ns2)
        dimension bins2(nx,ny)
c
        do ii=0,nx
        do jj=0,ny
        bins2(ii,jj)=0.0
        enddo
        enddo
        do 12 j=1,itm
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
c        write(2,*) dxcecl,dycecl,dzcecl,rcecl,j
c
        uxcecl=dxcecl/rcecl
        uycecl=dycecl/rcecl
        uzcecl=dzcecl/rcecl
c
c        write(2,*) "mm,dxceow,dyceow,dzceow,rceow,thelta1,
c     1  rcoscl(mm),rsincl(mm)"
c
        do 13 mm=ns1+1,nmol-2
        dxceow=pos(1,1,mm,j)-dxce
        dyceow=pos(1,2,mm,j)-dyce
        dzceow=pos(1,3,mm,j)-dzce
c
         if(abs(dxceow).gt.hbox)then
         dxceow=dxceow-sign(box,dxceow)
         endif
         if(abs(dyceow).gt.hbox)then
         dyceow=dyceow-sign(box,dyceow)
         endif
         if(abs(dzceow).gt.hbox)then
         dzceow=dzceow-sign(box,dzceow)
         endif
c*********end of PBC********************************************
        rceow=sqrt((dxceow)**2+(dyceow)**2+(dzceow)**2)
        rceow=abs(rceow)
        uxceow=dxceow/rceow
        uyceow=dyceow/rceow
        uzceow=dzceow/rceow
c        write(2,*) dxceow,dyceow,dzceow,rceow,mm
c
         theta1=acos(uxceow*uxcecl+uyceow*uycecl+uzceow*uzcecl)
c
         rcoscl(mm-ns1)=rceow*cos(theta1)
         rsincl(mm-ns1)=rceow*sin(theta1)
       if((rcoscl(mm-ns1).lt.1.0).and.(rcoscl(mm-ns1).gt.-1.0))then
       if (rsincl(mm-ns1).lt.1.0)then
c         write(2,*) mm-ns1,rceow,theta1,rcoscl(mm-nwt),rsincl(mm-nwt)
c
         arcscl=rcoscl(mm-ns1)/0.01
         arsicl=rsincl(mm-ns1)/0.01
c
         ii=int(arcscl)+1
         jj=int(arsicl)+1
         k1=ny-ii
c
c        write(2,*) mm-ns1,dxceow,dyceow,dzceow,rceow,theta1,
c     1  rcoscl(mm-ns1),rsincl(mm-ns1)
         bins2(k1,jj)=bins2(k1,jj)+1.0
         endif
         endif
 13     continue
c
 12     continue
        do i=0,nx
        do j=0,ny
        bins2(i,j)=bins2(i,j)/itm
        bins2(i,j)=bins2(i,j)*1000
c bins2 is the number density of solvent-2 around ion pair
        enddo
        enddo
        return
        end
