c
c     CIRCLES approximates possible crossing points of two different 
c     distance circles around two different points on the Earth's 
c     surface
c
c     INPUT:
c
c           dlat1, dlon1   - geographic coordinates of middle point 1
c
c           dlat2, dlon2   - geographic coordinates of middle point 2
c
c           d1, d2         - radius of circles 1 & 2 in [km]
c
c           dd1, dd2       - uncertainty of d1 & d2 in [km]
c
c           nc             - number of already found crossings
c
c     OUTPUT:
c
c           clat, clon     - geographic coordinates of crossing poins
c           
c           dcros          - uncertainty of crossing point(s) [km]
c
c           nc             - number of crossings including eventual 
c                            new ones
c
c     calls: depi, delazd, coomean
c     uses : q2
c
c     Johannes Schweitzer, NORSAR, April 2026
c

      subroutine circles(dlat1, dlon1, dlat2, dlon2, d1, d2, 
     +                   dd1, dd2, clat, clon, dcros, dist, nhit, nc)

c

      IMPLICIT real*8 (A-H,O-Z)
      IMPLICIT integer (I-N)

      dimension nhit(1000), clat(1000), clon(1000), dcros(1000),
     +          dist(1000), xlat(2,1001), xlon(2,1001)

      real*8    q2, fcheck

      call depi(dlat1,dlon1,dlat2,dlon2,dum,delstat,azid,bazd,dkm)

      ddd = dd1 + dd2
      if(delstat.le.ddd) go to 500

      ddel1 = d1 + d2
      ddel2 = ddel1 + ddd

      if(delstat.gt.ddel2) go to 500

      ddel3 = dabs((d1-dd1)-(d2-dd2))
      if(delstat.lt.ddel3) go to 500

      dcros0 = dpythag(dd1,dd2)
      r2d = 45.d0/datan(1.d0)
      ind = 2

      if(delstat.ge. ddel3 .and. 
     +   delstat.le. dabs((d1+dd1)-(d2+dd2))) then

        if(d1.le.d2) then
          call delazd(dlat1,dlon1,azid,d1,ind,xlat1,xlon1)
          call delazd(dlat2,dlon2,bazd+180.d0,d2,ind,xlat2,xlon2)
        else
          call delazd(dlat1,dlon1,azid+180.d0,d1,ind,xlat1,xlon1)
          call delazd(dlat2,dlon2,bazd,d2,ind,xlat2,xlon2)
        endif

        call depi(xlat1,xlon1,xlat2,xlon2,dum,delk,dum,dum,dkm)
        call coomean(xlat1,xlon1,xlat2,xlon2,glat,glon)

        ic = nc 

        if(nc.eq.0) then

           ic = 1
           nhit(ic) = 1
           clat(ic) = glat
           clon(ic) = glon
           dcros(ic) = dcros0 / dkm
           dist(ic) = delk

        else

           do 50 k = 1,ic

           call depi(clat(k),clon(k),glat,glon,dum,delk,dum,dum,dkm)

           if(delk.le.ddd) then
              if(delk.lt.dist(k)) then
                 clat(k) = glat
                 clon(k) = glon
                 dist(k) = delk
                 nhit(k) = nhit(k) + 1
              endif
              goto 50
           else if (k.eq.nc) then
              nc = nc + 1
              clat(nc) = glat
              clon(nc) = glon
              dcros(nc)  = dcros0 / dkm
              dist(nc) = delk
              nhit(nc) = 1
           endif

50         continue

        endif

        goto 400

      endif

      dphi1 = r2d*fcheck(1.d0 - q2(dd1*0.75d0)/2.d0/q2(d1),2)
      dphi2 = r2d*fcheck(1.d0 - q2(dd2*0.75d0)/2.d0/q2(d2),2)

      dphi = dmin1(dphi1,dphi2)
      nphi = int(360.d0/dphi)

      if(nphi.gt.1001) then
         nphi = 1001
         dphi = 0.36d0
      endif

      do 100 i=1,nphi

      azi = i * dphi

      call delazd(dlat1,dlon1,azi,d1,ind,xlat(1,i),xlon(1,i))
      call delazd(dlat2,dlon2,azi,d2,ind,xlat(2,i),xlon(2,i))

100   continue

      ic = nc

      do 300 i=1,nphi-1

      xlat1 = xlat(1,i)
      xlon1 = xlon(1,i)

      do 200 j=i+1,nphi

         xlat2 = xlat(2,i)
         xlon2 = xlon(2,i)

         call depi(xlat1,xlon1,xlat2,xlon2,dum,delk,
     +             dum,dum,dkm)

         if(delk.gt.ddd) goto 200

         if(ic.eq.0) then

            call coomean(xlat1,xlon1,xlat2,xlon2,glat,glon)
            clat(1) = glat
            clon(1) = glon
            dcros(1)  = (delk/2.d0 +dcros0) / dkm
            dist(1) = delk
            nhit(1) = 1
            ic = 1
            goto 200

         endif

         call coomean(xlat1,xlon1,xlat2,xlon2,glat,glon)

         do 150 k = 1,ic

            call depi(clat(k),clon(k),glat,glon,dum,delk,dum,dum,dkm)

            if(delk.le.ddd) then
               if(delk.lt.dist(k)) then
                  clat(k) = glat
                  clon(k) = glon
                  dist(k) = delk
                  dcros(k) =  (delk/2.d0 +dcros0) / dkm
               endif
               if(k.le.nc) nhit(k) = nhit(k) + 1
               goto 150
            endif

            if (k.eq.ic) then
               ij = ic + 1
               clat(ij) = glat
               clon(ij) = glon
               dcros(ij) = (delk/2.d0 +dcros0) / dkm
               dist(ij) = delk
               nhit(ij) = 1
               ic = ij
               goto 200
            endif

150      continue

200   continue

300   continue

400   continue

      nc = ic

500   continue
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine coomean(xlat,xlon,ylat,ylon,zlat,zlon)
c
      subroutine coomean(xlat,xlon,ylat,ylon,zlat,zlon)
      real*8 xlat,xlon,ylat,ylon,zlat,zlon
      real*8 x1, x2
      real*8 alpha1

      zlat = (xlat + ylat) / 2.d0

      x1 = xlon
      x2 = ylon
      if(x1*x2.lt.0.d0) then
        if (x1.lt.0.d0) x1 = x1 + 360.d0
        if (x2.lt.0.d0) x2 = x2 + 360.d0
      endif
      zlon = alpha1 ((x1 + x2) / 2.d0)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     driver subroutine for subroutine circles
c
      subroutine get_circ(sla,slo,nstat,del,dels,ndel,
     *                   blat,blats,blon,blons,nepi,typctl)

      IMPLICIT real*8 (A-H,O-Z)
      IMPLICIT integer (I-N)

      dimension sla(*),slo(*),del(*),dels(*), ndel(*)

      integer   typctl

      dimension nhit(1000), clat(1000), clon(1000), dcros(1000),
     +          dist(1000)

      nepi = 0

      nc = 0
      nc0 = 0

      do 200 i = 1,nstat-1

      if(ndel(i).eq.0) go to 200

      if(nepi.gt.5) go to 250

      sla1 = sla(i)
      slo1 = slo(i)
      del1 = del(i)
      delt  = del(i)*0.1d0
      dels1 = dmax1(delt,dels(i))

      do 100 j = i+1,nstat

         if(ndel(j).eq.0) go to 100

         delt  = del(j)*0.1d0
         dels2 = dmax1(delt,dels(j))

         call circles(sla1,slo1,sla(j),slo(j),del1,del(j),
     +                dels1,dels2,clat,clon,dcros,dist,nhit,nc)

         if(nc.gt.nc0) then
            nepi = nepi + 1
            nc0 = nc
         endif

100   continue

200   continue

250   if(nc.gt.0 .and. nepi.gt.2) then
         nhit0 = 0
         do 300 i = 1,nc
         if(nhit(i).gt.nhit0) then
            nx = i
            nhit0 = nhit(i)
         endif
300      continue

         blat = clat(nx)
         blon = clon(nx)
         blats = dcros(nx)
         if(dabs(blat).lt. 90.d0) then
            blons = dabs(blats) / dcos(blat*datan(1.d0)/45.d0)
            if(blons.gt.180.d0) blons = 180.d0
         else
            blons = 180.d0
         endif

         if(typctl.gt.5) print *,'GET_CIRCLE:',nx,blat,blats,blon,blons

      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
