c-----------------------------------------------------------------------
c ********************* JET-FINDING ROUTINES ***************************
c-----------------------------------------------------------------------
c
      subroutine jetalg(Pin)
      implicit DOUBLE PRECISION (a-h,o-z)
      common/rcone/rcone
      common/imtt/imtt
      common/imtt2/imtt2
      common/imerge/imerge
      COMMON/JETS/POUT,zkt,zeta,zphi,n1,n2,n3
      COMMON/JETSCALES/SCALEREN,SCALEFACT,Q2
      COMMON/JETPARAMS/R,IALG
      DOUBLE PRECISION P(4,7),Pin(4,7),Pout(4,7),xmin,xmax,SCALEJET
     $                 ,QQ,Y,Q2
      integer ialg,i,j,imerge,n1,n2,n3,FRAMEj
      integer imtt2(3:5,3:5),imtt(3:4,3:5)
      dimension ykt(1:3),yeta(1:3),yphi(1:3)
      dimension zkt(1:3),zeta(1:3),zphi(1:3)     
c
c
c The matrix IMTT(I,J), with 3<=I<=4 and 3<=J<=5, is defined such
c that {IMTT(I,J_0)| I=3,4, J_0 fixed}={3,4,5}-{J_0}
      imtt(3,3)=4
      imtt(4,3)=5
      imtt(3,4)=3
      imtt(4,4)=5
      imtt(3,5)=3
      imtt(4,5)=4
c The matrix IMTT2(I,J), with 3<=I,J<=5, is defined such
c that IMTT2(I_0,J_0)={3,4,5}-{I_0,J_0}
      imtt2(3,4)=5
      imtt2(3,5)=4
      imtt2(4,3)=5
      imtt2(4,5)=3
      imtt2(5,3)=4
      imtt2(5,4)=3

      DO J=1,7
          DO I=1,4 
          POUT(I,J)=PIN(I,J)
          ENDDO
      ENDDO    
c define variables for partons
      DO J=2,4
         YKT(J-1)=(PIN(1,J)**2+PIN(2,J)**2)**0.5
         IF (YKT(J-1).EQ.0) THEN 
c         IF (abs(YKT(J-1)).LT.1d-5) THEN !tocar aca
c         YKT(J-1)=0d0
         YETA(J-1)=1.d8
         YPHI(J-1)=0d0       
         ELSE       
         YETA(J-1)=0.5D0*DLOG((PIN(4,J)+PIN(3,J))/(PIN(4,J)-PIN(3,J)))
         YPHI(J-1)=ATAN2(PIN(2,J),PIN(1,J))
         ENDIF
      ENDDO

      if (ialg.eq.1) then
          call antikt(ykt,yeta,yphi,zkt,zeta,zphi,R,imerge,0)
      elseif (ialg.eq.2) then
          call antikt(ykt,yeta,yphi,zkt,zeta,zphi,R,imerge,1)     
      elseif (ialg.eq.3) then
          call kt(ykt,yeta,yphi,zkt,zeta,zphi,R,imerge)
      else 
          write(6,*)'Wrong jet algorithm parameter ialg= ',ialg
          stop
      endif
C imerge=0 means no merging  --> partons are jets
C imerge=1 means merging  --> one parton disappears
C imerge=2 means double-merging  --> two parton disappears

c Order jets according to pT      
      xmax=max(zkt(1),zkt(2),zkt(3))
      xmin=min(zkt(1),zkt(2),zkt(3))
      if (xmax.eq.0d0) then  ! no jets
             n1=1
             n2=2
             n3=3
      else       
         do i=1,3
           if(zkt(i).eq.xmax)n1=i
           if(zkt(i).eq.xmin)n3=i
         enddo
         n2=imtt2(n1+2,n3+2)-2
      endif
c c now et(n3)<=et(n2)<=et(n1)  HARDEST JET IS N1, SECOND HARDEST JET IS N2
C THEREFORE  ZKT(N1) IS THE TRANSVERSE MOMENTUM OF THE HARDEST JET
C            ZETA(N1) IS THE RAPIDITY OF THE HARDEST JET 
C            ZPHI(N1) IS THE AZIMUTHAL ANGLE OF THE HARDEST JET
C AND SO ON..

       IF(IMERGE.EQ.0) THEN  ! PARTONS ARE JETS, JUST ORDER MOMENTA (2,3,4)
          DO I=1,4
          POUT(I,N1+1) = PIN(I,N1+1) !debug
          POUT(I,N2+1) = PIN(I,N2+1)
          POUT(I,N3+1) = PIN(I,N3+1)
          ENDDO
       ELSEIF(IMERGE.EQ.1) THEN  ! TWO JETS, eliminate particle N3       
          c1 = (exp(2*zeta(n1))-1)/(exp(2*zeta(n1))+1)
          POUT(4,n1+1) = zkt(n1)/(1-c1**2)**0.5
          POUT(1,n1+1) = zkt(n1) * Cos(zphi(n1))
          POUT(2,n1+1) = zkt(n1) * Sin(zphi(n1))
          POUT(3,n1+1) = zkt(n1)/(1-c1**2)**0.5 * c1
C          
          c2 = (exp(2*zeta(n2))-1)/(exp(2*zeta(n2))+1)
          POUT(4,n2+1) = zkt(n2)/(1-c2**2)**0.5
          POUT(1,n2+1) = zkt(n2) * Cos(zphi(n2))
          POUT(2,n2+1) = zkt(n2) * Sin(zphi(n2))
          POUT(3,n2+1) = zkt(n2)/(1-c2**2)**0.5 * c2
C                    
          POUT(1,n3+1) = 0D0
          POUT(2,n3+1) = 0D0
          POUT(3,n3+1) = 0D0
          POUT(4,n3+1) = 0D0
        ELSEIF(IMERGE.EQ.2) THEN  ! ONE JET, eliminate particles N2,N3       
          c1 = (exp(2*zeta(n1))-1)/(exp(2*zeta(n1))+1)
          POUT(4,n1+1) = zkt(n1)/(1-c1**2)**0.5
          POUT(1,n1+1) = zkt(n1) * Cos(zphi(n1))
          POUT(2,n1+1) = zkt(n1) * Sin(zphi(n1))
          POUT(3,n1+1) = zkt(n1)/(1-c1**2)**0.5 * c1
C          
          POUT(4,n2+1) = 0D0
          POUT(1,n2+1) = 0D0
          POUT(2,n2+1) = 0D0
          POUT(3,n2+1) = 0D0
C                    
          POUT(1,n3+1) = 0D0
          POUT(2,n3+1) = 0D0
          POUT(3,n3+1) = 0D0
          POUT(4,n3+1) = 0D0
        ELSE
           WRITE(6,*)"Wrong imerge in JETALG" ,imerge 
        ENDIF

C     SET JET SCALES IN THE USER ROUTINE
      CALL SCALES(SCALEREN,SCALEFACT,Q2)
      end
      

      function rdist(eta1,phi1,eta2,phi2)
c
c Square root of eq. (2.22)
c
      implicit DOUBLE PRECISION (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
c
      dteta=eta1-eta2
      dtphi=phi1-phi2
      dtphi=dacos(dcos(dtphi))
      rdist=dsqrt(dteta**2+dtphi**2)
      return
      end


      function ptdist(pt1,pt2)
c
c (pt1+pt2)*rcone/max(pt1,pt2), for the cone algorithm
c
      implicit DOUBLE PRECISION (a-h,o-z)
      common/rcone/rcone
c
      tmp=(pt1+pt2)/max(pt1,pt2)
      ptdist=tmp*rcone
      return
      end
c
c End of the cone algorithm
c
c 
c  KT algorithm
c
      subroutine kt(ykt,yeta,yphi,zkt,zeta,zphi,d_jet,imerge)
c This subroutine returns jet variables (zkt, zeta, zphi) defined
c in terms of partonic variables (ykt, yeta, yphi) in the KT algorithm
c algorithm. The parameter D must be set in the function ddist (in this file)
      implicit none
      common/djet/djet
      common/imtt/imtt
      common/imtt2/imtt2
      double PRECISION ykt(1:3),yeta(1:3),yphi(1:3)
      double PRECISION zkt(1:3),zeta(1:3),zphi(1:3)
      double PRECISION xd(1:3),xd2(1:2,2:3)
      double PRECISION xmin,tmpphi,tmpeta,tmpkt,d_jet,djet,pi,ddist
      integer imtt(3:4,3:5)
      integer imtt2(3:5,3:5)
      integer imerge,imerget,i,j,n1,n2,n3,i2bd
      parameter (pi=3.14159265358979312D0)
c
      djet=d_jet
      i2bd=1
      do i=1,3
        if(ykt(i).eq.0.d0) then
          i2bd=0
          n1=i
        endif
      enddo
    
      if (i2bd.eq.0) then
        n2=imtt(3,n1+2)-2
        n3=imtt(4,n1+2)-2
c  1-body kinematics
        if (ykt(n2).eq.0 .or. ykt(n3).eq.0) then
          do i=1,3
            zkt(i)=ykt(i)
            zeta(i)=yeta(i)
            zphi(i)=yphi(i)
          enddo
          imerge=0
        else
c  2-body kinematics
          xd(n2)=ykt(n2)**2
          xd(n3)=ykt(n3)**2
          xd2(n2,n3)=ddist(ykt(n2),yeta(n2),yphi(n2),
     $                     ykt(n3),yeta(n3),yphi(n3))
c         merging case 
          if (xd2(n2,n3).lt.min(xd(n2),xd(n3))) then
              zkt(n1)=ykt(n1)
              zeta(n1)=yeta(n1)
              zphi(n1)=yphi(n1)
              call xmerge(ykt(n2),ykt(n3),tmpkt,
     #                    yeta(n2),yeta(n3),tmpeta,
     #                    yphi(n2),yphi(n3),tmpphi)
              zkt(n2)=tmpkt
              zeta(n2)=tmpeta
              zphi(n2)=tmpphi
              zkt(n3)=0.d0
              zeta(n3)=1.d8
              zphi(n3)=0.d0
              imerge=1
c              write(6,*)'esto es lab frame (1)'
c         no-merging       
          else
            do i=1,3
              zkt(i)=ykt(i)
              zeta(i)=yeta(i)
              zphi(i)=yphi(i)
            enddo
            imerge=0
          endif
        endif
c  3-body kinematics
      elseif (i2bd.eq.1) then
        xd(1)=ykt(1)**2
        xd(2)=ykt(2)**2
        xd(3)=ykt(3)**2
        xd2(1,2)=ddist(ykt(1),yeta(1),yphi(1),ykt(2),yeta(2),yphi(2))
        xd2(1,3)=ddist(ykt(1),yeta(1),yphi(1),ykt(3),yeta(3),yphi(3))
        xd2(2,3)=ddist(ykt(2),yeta(2),yphi(2),ykt(3),yeta(3),yphi(3))
        xmin=min(xd(1),xd(2),xd(3),xd2(1,2),xd2(1,3),xd2(2,3))
        imerget=-1 !temporal value

c       check no merging case (0)
        do i=1,3
          if (xd(i).eq.xmin) then
            imerget=0
            n1=i
          endif
        enddo
c       check merging case (1)
        if (imerget.eq.-1) then
          do i=1,2
            do j=i+1,3
              if (xd2(i,j).eq.xmin) then
                imerget=1
                n1=i
                n2=j
              endif
            enddo
          enddo
        endif

c       no merging (0)     
        if (imerget.eq.0) then
          n2=imtt(3,n1+2)-2
          n3=imtt(4,n1+2)-2
c          xd(n2)=ykt(n2)**2
c          xd(n3)=ykt(n3)**2
c          xd2(n2,n3)=ddist(ykt(n2),yeta(n2),yphi(n2),
c     $                     ykt(n3),yeta(n3),yphi(n3))
c         no-merging to merging (01) - one jet
          if (xd2(n2,n3).lt.min(xd(n2),xd(n3))) then
            zkt(n1)=ykt(n1)
            zeta(n1)=yeta(n1)
            zphi(n1)=yphi(n1)
            call xmerge(ykt(n2),ykt(n3),tmpkt,
     #                  yeta(n2),yeta(n3),tmpeta,
     #                  yphi(n2),yphi(n3),tmpphi)
            zkt(n2)=tmpkt
            zeta(n2)=tmpeta
            zphi(n2)=tmpphi
            zkt(n3)=0.d0
            zeta(n3)=1.d8
            zphi(n3)=0.d0
            imerge=1
c            write(6,*)'esto es lab frame yeah'
c         no-merging to no-merging (00)          
          else
            do i=1,3
              zkt(i)=ykt(i)
              zeta(i)=yeta(i)
              zphi(i)=yphi(i)
            enddo
            imerge=0
          endif
c       merging case (1)   
        elseif (imerget.eq.1) then
          n3=imtt2(n1+2,n2+2)-2
          zkt(n3)=ykt(n3)
          zeta(n3)=yeta(n3)
          zphi(n3)=yphi(n3)
          call xmerge(ykt(n1),ykt(n2),tmpkt,
     #                yeta(n1),yeta(n2),tmpeta,
     #                yphi(n1),yphi(n2),tmpphi)
          zkt(n1)=tmpkt
          zeta(n1)=tmpeta
          zphi(n1)=tmpphi
          zkt(n2)=0.d0
          zeta(n2)=1.d8
          zphi(n2)=0.d0

c         new iteration of the algorithm
          xd(n1)=zkt(n1)**2
          xd(n3)=zkt(n3)**2
          xd2(n1,n3)=ddist(zkt(n1),zeta(n1),zphi(n1),
     $        	           zkt(n3),zeta(n3),zphi(n3))
c         merging to merging (11) - one jet
          if (xd2(n1,n3).lt.min(xd(n1),xd(n3))) then
            call xmerge(zkt(n1),zkt(n3),tmpkt,
     #                  zeta(n1),zeta(n3),tmpeta,
     #                  zphi(n1),zphi(n3),tmpphi)
            zkt(n1)=tmpkt
            zeta(n1)=tmpeta
            zphi(n1)=tmpphi
            zkt(n3)=0.d0
            zeta(n3)=1.d8
            zphi(n3)=0.d0
            imerge=2
c            write(6,*)'esto es lab frame'
c            write(6,*)'zkt1',  zkt(n1), zeta(n1), zphi(n1)
c            Stop
c         merging to no-merging (10) - two jet
          else
            imerge=1
          endif
        else
          write(6,*)'Fatal error in subroutine esjet'
          write(6,*)'xmin', xmin
          write(6,*)'xd(1)',xd(1)
          write(6,*)'xd(2)',xd(2)
          write(6,*)'xd(3)',xd(3)
          write(6,*)'xd2(1,2)',xd2(1,2)
          write(6,*)'xd2(1,3)',xd2(1,3)
          write(6,*)'xd2(2,3)',xd2(2,3)
          stop
        endif
      endif  
      return
      end

      function rdist2(eta1,phi1,eta2,phi2)
      implicit DOUBLE PRECISION (a-h,o-z)
      Double Precision pi,eta1,eta2,phi1,phi2,dteta,dtphi,rdist2
      parameter (pi=3.14159265358979312D0)
c
      dteta=eta1-eta2
      dtphi=phi1-phi2
      dtphi=dacos(dcos(dtphi))
      rdist2=dteta**2+dtphi**2
      return
      end

      function ddist(pt1,eta1,phi1,pt2,eta2,phi2)
      implicit DOUBLE PRECISION (a-h,o-z)
      Double Precision ddist,djet,pt1,eta1,phi1,pt2,eta2,phi2
      common/djet/djet
      tmp=min(pt1**2,pt2**2)*rdist2(eta1,phi1,eta2,phi2)
      ddist=tmp/djet**2
      return
      end
c
c End of the KT algorithm
c
      subroutine xmerge
     #    (xkt1,xkt2,xkt3,xeta1,xeta2,xeta3,xphi1,xphi2,xphi3)
c
c Defines the jet variables (xkt3,xeta3,xphi3) in terms of the
c partonic variables when two partons are merged
c     ET weighted recombination scheme
c
      implicit DOUBLE PRECISION (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
c
      xkt3=xkt1+xkt2
      xeta3=(xkt1*xeta1+xkt2*xeta2)/xkt3
      if(abs(xphi1-xphi2).lt.pi)then
        xphi3=(xkt1*xphi1+xkt2*xphi2)/xkt3
      else
        if(xphi2.lt.xphi1)then
          xphi3=(xkt1*xphi1+xkt2*(xphi2+2*pi))/xkt3
        else
          xphi3=(xkt1*xphi1+xkt2*(xphi2-2*pi))/xkt3
        endif
      endif
      xphi3=atan2(sin(xphi3),cos(xphi3))
      return
      end

c

      subroutine xmerge2
     #    (xkt1,xkt2,xkt3,xeta1,xeta2,xeta3,xphi1,xphi2,xphi3,xy3,tmpm3)
c
c Defines the jet variables (xkt3,xeta3,xphi3) in terms of the
c partonic variables when two partons are merged
c  Recombination in E-scheme for merging
c  xeta3= true rapidity (not pseudorapidity)
c
      implicit DOUBLE PRECISION (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
c
c energy and momentum of parton 1
      c1 = (exp(2*xeta1)-1)/(exp(2*xeta1)+1)
      e1 = xkt1/(1-c1**2)**0.5
      px1 = xkt1 * Cos(xphi1)
      py1 = xkt1 * sin(xphi1)
      pz1 =  e1 * c1
c energy and momentum of parton 2
      c2 = (exp(2*xeta2)-1)/(exp(2*xeta2)+1)
      e2 = xkt2/(1-c2**2)**0.5
      px2 = xkt2 * Cos(xphi2)
      py2 = xkt2 * sin(xphi2)
      pz2 =  e2 * c2
c energy and momentum of parton 1+2
      e3 = e1 + e2
      px3 = px1 + px2
      py3 = py1 + py2
      pz3 = pz1 + pz2
c phi,eta and Et for merged jet
      xphi3 = atan2(py3,px3)
      xy3 = 0.5d0 * log( (e3+pz3)/(e3-pz3)) !y = rapidity
      pmod=(pz3**2 + px3**2 + py3**2)**0.5
      tmpm3=(e3**2 - pz3**2 - px3**2 - py3**2)**0.5
      xeta3= 0.5d0 * log( (pmod+pz3)/(pmod-pz3)) ! eta pseudorapidity
      xkt3=(px3**2+py3**2)**0.5
c
      return
      end

      subroutine xmerge3
     #    (xkt1,xkt2,xkt3,xeta1,xeta2,xeta3,xphi1,xphi2,xphi3,xy3,tmpm3)
c
c Defines the jet variables (xkt3,xeta3,xphi3) in terms of the
c partonic variables when one parton (2) and one massive (tmpm3) jet (1) 
c are merged
c  Recombination in E-scheme for merging
c  xeta3= true rapidity (not pseudorapidity)
c
      implicit DOUBLE PRECISION (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
c
c energy and momentum of parton 1
      c1 = (exp(2*xeta1)-1)/(exp(2*xeta1)+1)
      e1 = (xkt1**2+tmp3**2)**0.5/(1-c1**2)**0.5
      px1 = xkt1 * Cos(xphi1)
      py1 = xkt1 * sin(xphi1)
      pz1 =  e1 * c1
c energy and momentum of parton 2
      c2 = (exp(2*xeta2)-1)/(exp(2*xeta2)+1)
      e2 = xkt2/(1-c2**2)**0.5
      px2 = xkt2 * Cos(xphi2)
      py2 = xkt2 * sin(xphi2)
      pz2 =  e2 * c2
c energy and momentum of parton 1+2
      e3 = e1 + e2
      px3 = px1 + px2
      py3 = py1 + py2
      pz3 = pz1 + pz2
c phi,eta and Et for merged jet
      xphi3 = atan2(py3,px3)
      xy3 = 0.5d0 * log( (e3+pz3)/(e3-pz3)) !y = rapidity
      pmod=(pz3**2 + px3**2 + py3**2)**0.5
      xeta3= 0.5d0 * log( (pmod+pz3)/(pmod-pz3)) ! eta pseudorapidity
      xkt3=(px3**2+py3**2)**0.5
      tmpm3=(e3**2 - pz3**2 - px3**2 - py3**2)**0.5 !mass
c
      return
      end

c 
c anti-kt algorithm
c
      subroutine antikt(ykt,yeta,yphi,zkt,zeta,zphi,d_jet,iimerge,ischm)
c This subroutine returns jet variables (zkt, zeta, zphi) defined
c in terms of partonic variables (ykt, yeta, yphi) in the anti-kt algorithm
      implicit DOUBLE PRECISION (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      common/djet/djet
      common/imtt/imtt
      common/imtt2/imtt2
      DOUBLE PRECISION ykt(1:3),yeta(1:3),yphi(1:3)
      DOUBLE PRECISION zkt(1:3),zeta(1:3),zphi(1:3)
      DOUBLE PRECISION xd(1:3),xd2(1:2,2:3)
      integer imtt(3:4,3:5),imtt2(3:5,3:5),ischm
c
      imerge=-14
      iimerge=0
      djet=d_jet
      i2bd=1
      do i=1,3
        if (ykt(i).eq.0.d0) then
        i2bd=0
        n1=i
        endif
      enddo

      if(i2bd.eq.0)then
        n2=imtt(3,n1+2)-2
        n3=imtt(4,n1+2)-2
c  1-body kinematics
        if (ykt(n2).eq.0 .or. ykt(n3).eq.0) then
          do i=1,3
            zkt(i)=ykt(i)
            zeta(i)=yeta(i)
            zphi(i)=yphi(i)
          enddo
          iimerge=0
        else
c 2-body kinematics
          xd(n2)=1d0/ykt(n2)**2
          xd(n3)=1d0/ykt(n3)**2
          xd2(n2,n3)=adist(ykt(n2),yeta(n2),yphi(n2),
     $                     ykt(n3),yeta(n3),yphi(n3))
c         merging case 
          if (xd2(n2,n3).lt.min(xd(n2),xd(n3))) then
              zkt(n1)=ykt(n1)
              zeta(n1)=yeta(n1)
              zphi(n1)=yphi(n1)
              if (ischm.eq.0) then
               call xmerge2(ykt(n2),ykt(n3),tmpkt,
     #                     yeta(n2),yeta(n3),tmpeta,
     #                     yphi(n2),yphi(n3),tmpphi,tmpy,tmptmpm3)
              elseif (ischm.eq.1) then
               call xmerge(ykt(n2),ykt(n3),tmpkt,
     #                     yeta(n2),yeta(n3),tmpeta,
     #                     yphi(n2),yphi(n3),tmpphi)
              endif
              zkt(n2)=tmpkt
              zeta(n2)=tmpeta
              zphi(n2)=tmpphi
              zkt(n3)=0.d0
              zeta(n3)=1.d8
              zphi(n3)=0.d0
              iimerge=1
c         no-merging       
          else
            do i=1,3
              zkt(i)=ykt(i)
              zeta(i)=yeta(i)
              zphi(i)=yphi(i)
            enddo
            iimerge=0
          endif
        endif
c 3-body kinematics
      elseif (i2bd.eq.1) then
        xd(1)=1d0/ykt(1)**2
        xd(2)=1d0/ykt(2)**2
        xd(3)=1d0/ykt(3)**2
        xd2(1,2)=adist(ykt(1),yeta(1),yphi(1),ykt(2),yeta(2),yphi(2))
        xd2(1,3)=adist(ykt(1),yeta(1),yphi(1),ykt(3),yeta(3),yphi(3))
        xd2(2,3)=adist(ykt(2),yeta(2),yphi(2),ykt(3),yeta(3),yphi(3))
        xmin=min(xd(1),xd(2),xd(3),xd2(1,2),xd2(1,3),xd2(2,3))
        imerge=-1

c       check no merging case (0)        
        do i=1,3
          if (xd(i).eq.xmin) then
            imerge=0
            n1=i
          endif
        enddo
c       check merging case (1)        
        if (imerge.eq.-1) then
          do i=1,2
            do j=i+1,3
              if(xd2(i,j).eq.xmin)then
                imerge=1
                n1=i
                n2=j
              endif
            enddo
          enddo
        endif

c       no merging (0)        
        if(imerge.eq.0)then
c no merging yet eliminate n1 from the list (n1 = jet)
          n2=imtt(3,n1+2)-2
          n3=imtt(4,n1+2)-2

          zkt(n1)=ykt(n1)
          zeta(n1)=yeta(n1)
          zphi(n1)=yphi(n1)
c         no-merging to merging (01)          
          if (xd2(n2,n3).lt.min(xd(n2),xd(n3))) then
          if (ischm.eq.0) then 
           call xmerge2(ykt(n3),ykt(n2),tmpkt,
     #                 yeta(n3),yeta(n2),tmpeta,
     #                 yphi(n3),yphi(n2),tmpphi,tmpy,tmpm3)
          elseif (ischm.eq.1) then
           call xmerge(ykt(n3),ykt(n2),tmpkt,
     #                 yeta(n3),yeta(n2),tmpeta,
     #                 yphi(n3),yphi(n2),tmpphi)
          endif
          zkt(n2)=tmpkt
          zeta(n2)=tmpeta
          zphi(n2)=tmpphi
          zkt(n3)=0.d0
          zeta(n3)=1.d8
          zphi(n3)=0.d0
          iimerge=1          

c         no-merging to no-merging (00) 
          else
            zkt(n2)=ykt(n2)
            zeta(n2)=yeta(n2)
            zphi(n2)=yphi(n2)
            zkt(n3)=ykt(n3)
            zeta(n3)=yeta(n3)
            zphi(n3)=yphi(n3)
          endif
                              
        elseif (imerge.eq.1) then
c merge n1 and n2 (into the new n1), n3 is the other jet        
          n3=imtt2(n1+2,n2+2)-2
          zkt(n3)=ykt(n3)
          zeta(n3)=yeta(n3)
          zphi(n3)=yphi(n3)
          if (ischm.eq.0) then
           call xmerge2(ykt(n1),ykt(n2),tmpkt,
     #                 yeta(n1),yeta(n2),tmpeta,
     #                 yphi(n1),yphi(n2),tmpphi,tmpy,tmpm3)
          elseif (ischm.eq.1) then
           call xmerge(ykt(n1),ykt(n2),tmpkt,
     #                 yeta(n1),yeta(n2),tmpeta,
     #                 yphi(n1),yphi(n2),tmpphi)
          endif
          zkt(n1)=tmpkt
          zeta(n1)=tmpeta
          zphi(n1)=tmpphi
          zkt(n2)=0.d0
          zeta(n2)=1.d8
          zphi(n2)=0.d0

c         new iteration of the algorithm
          xd(n1)=1d0/zkt(n1)**2
          xd(n3)=1d0/zkt(n3)**2
          if (ischm.eq.0) then
           xd2(n1,n3)=adist(zkt(n1),tmpy,zphi(n1),
     $                      zkt(n3),zeta(n3),zphi(n3))
          elseif (ischm.eq.1) then
           xd2(n1,n3)=adist(zkt(n1),zeta(n1),zphi(n1),
     $                      zkt(n3),zeta(n3),zphi(n3))
          endif
c         merging to merging (11) - one jet
          if (xd2(n1,n3).lt.min(xd(n1),xd(n3))) then
            if (ischm.eq.0) then
             call xmerge3(zkt(n1),zkt(n3),tmpkt,
     #                   tmpy,zeta(n3),tmpeta,
     #                   zphi(n1),zphi(n3),tmpphi,tmpy,tmpm3)
            elseif (ischm.eq.1) then
             call xmerge(zkt(n1),zkt(n3),tmpkt,
     #                   zeta(n1),zeta(n3),tmpeta,
     #                   zphi(n1),zphi(n3),tmpphi)
            endif
            zkt(n1)=tmpkt
            zeta(n1)=tmpeta
            zphi(n1)=tmpphi
            zkt(n3)=0.d0
            zeta(n3)=1.d8
            zphi(n3)=0.d0
            iimerge=2

c         merging to no-merging (10) - two jet
          else
            iimerge=1
          endif 
        else
          write(6,*)'Fatal error in subroutine antikt'
          write(6,*)'imerge=',imerge
          stop
        endif
      endif
      return
      end


      function adist(pt1,eta1,phi1,pt2,eta2,phi2)
      implicit DOUBLE PRECISION (a-h,o-z)
      common/djet/djet
c
      tmp=min(1d0/pt1**2,1d0/pt2**2)*rdist2(eta1,phi1,eta2,phi2)
      adist=tmp/djet**2
      return
      end
c
c End of anti-KT algorithm
c
