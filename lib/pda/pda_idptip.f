      SUBROUTINE PDA_IDPTIP(xd,yd,zd,nt,ipt,nl,ipl,pdd,iti,xii, 
     :                      yii,zii,istat)
 
c this subroutine performs punctual interpolation or extrapola-
c tion, i.e., determines the z value at a point.
 
c the input parameters are
c     xd,yd,zd = arrays of dimension ndp containing the x,
c           y, and z coordinates of the data points, where
c           ndp is the number of the data points,
c     nt  = number of triangles,
c     ipt = integer array of dimension 3*nt containing the
c           point numbers of the vertices of the triangles,
c     nl  = number of border line segments,
c     ipl = integer array of dimension 3*nl containing the
c           point numbers of the end points of the border
c           line segments and their respective triangle
c           numbers,
c     pdd = array of dimension 5*ndp containing the partial
c           derivatives at the data points,
c     iti = triangle number of the triangle in which lies
c           the point for which interpolation is to be
c           performed,
c     xii,yii = x and y coordinates of the point for which
c           interpolation is to be performed.
 
c the output parameter is
c     zii = interpolated z value.
*     istat = starlink error message.
 
 
c declaration statements
c      dimension xd(100), yd(100), zd(100), ipt(585), ipl(300), pdd(500)
      dimension xd(*), yd(*), zd(*), ipt(*), ipl(*), pdd(*)
      common /idpi  / itpv
      dimension x(3), y(3), z(3), pd(15), zu(3), zv(3), zuu(3), zuv(3), 
     :          zvv(3)
      real lu, lv
      equivalence (p5, p50)
      integer istat
 
*   check the inherited error status.
      if ( istat.ne.0 ) return
 
c preliminary processing
      it0 = iti
      ntl = nt + nl
      if ( it0.le.ntl ) then
 
c calculation of zii by interpolation.
c checks if the necessary coefficients have been calculated.
         if ( it0.ne.itpv ) then
 
c loads coordinate and partial derivative values at the
c vertices.
            jipt = 3*(it0-1)
            jpd = 0
            do 20 i = 1, 3
               jipt = jipt + 1
               idp = ipt(jipt)
               x(i) = xd(idp)
               y(i) = yd(idp)
               z(i) = zd(idp)
               jpdd = 5*(idp-1)
               do 10 kpd = 1, 5
                  jpd = jpd + 1
                  jpdd = jpdd + 1
                  pd(jpd) = pdd(jpdd)
 10            continue
 20         continue
 
c determines the coefficients for the coordinate system
c transformation from the x-y system to the u-v system
c and vice versa.
            x0 = x(1)
            y0 = y(1)
            a = x(2) - x0
            b = x(3) - x0
            c = y(2) - y0
            d = y(3) - y0
            ad = a*d
            bc = b*c
            dlt = ad - bc
            ap = d/dlt
            bp = -b/dlt
            cp = -c/dlt
            dp = a/dlt
 
c converts the partial derivatives at the vertices of the
c triangle for the u-v coordinate system.
            aa = a*a
            act2 = 2.0*a*c
            cc = c*c
            ab = a*b
            adbc = ad + bc
            cd = c*d
            bb = b*b
            bdt2 = 2.0*b*d
            dd = d*d
            do 40 i = 1, 3
               jpd = 5*i
               zu(i) = a*pd(jpd-4) + c*pd(jpd-3)
               zv(i) = b*pd(jpd-4) + d*pd(jpd-3)
               zuu(i) = aa*pd(jpd-2) + act2*pd(jpd-1) + cc*pd(jpd)
               zuv(i) = ab*pd(jpd-2) + adbc*pd(jpd-1) + cd*pd(jpd)
               zvv(i) = bb*pd(jpd-2) + bdt2*pd(jpd-1) + dd*pd(jpd)
 40         continue
 
c calculates the coefficients of the polynomial.
            p00 = z(1)
            p10 = zu(1)
            p01 = zv(1)
            p20 = 0.5*zuu(1)
            p11 = zuv(1)
            p02 = 0.5*zvv(1)
            h1 = z(2) - p00 - p10 - p20
            h2 = zu(2) - p10 - zuu(1)
            h3 = zuu(2) - zuu(1)
            p30 = 10.0*h1 - 4.0*h2 + 0.5*h3
            p40 = -15.0*h1 + 7.0*h2 - h3
            p50 = 6.0*h1 - 3.0*h2 + 0.5*h3
            h1 = z(3) - p00 - p01 - p02
            h2 = zv(3) - p01 - zvv(1)
            h3 = zvv(3) - zvv(1)
            p03 = 10.0*h1 - 4.0*h2 + 0.5*h3
            p04 = -15.0*h1 + 7.0*h2 - h3
            p05 = 6.0*h1 - 3.0*h2 + 0.5*h3
            lu = sqrt(aa+cc)
            lv = sqrt(bb+dd)
            thxu = atan2(c, a)
            thuv = atan2(d, b) - thxu
            csuv = cos(thuv)
            p41 = 5.0*lv*csuv/lu*p50
            p14 = 5.0*lu*csuv/lv*p05
            h1 = zv(2) - p01 - p11 - p41
            h2 = zuv(2) - p11 - 4.0*p41
            p21 = 3.0*h1 - h2
            p31 = -2.0*h1 + h2
            h1 = zu(3) - p10 - p11 - p14
            h2 = zuv(3) - p11 - 4.0*p14
            p12 = 3.0*h1 - h2
            p13 = -2.0*h1 + h2
            thus = atan2(d-c, b-a) - thxu
            thsv = thuv - thus
            aa = sin(thsv)/lu
            bb = -cos(thsv)/lu
            cc = sin(thus)/lv
            dd = cos(thus)/lv
            ac = aa*cc
            ad = aa*dd
            bc = bb*cc
            g1 = aa*ac*(3.0*bc+2.0*ad)
            g2 = cc*ac*(3.0*ad+2.0*bc)
            h1 = -aa*aa*aa*(5.0*aa*bb*p50+(4.0*bc+ad)*p41)
     :           - cc*cc*cc*(5.0*cc*dd*p05+(4.0*ad+bc)*p14)
            h2 = 0.5*zvv(2) - p02 - p12
            h3 = 0.5*zuu(3) - p20 - p21
            p22 = (g1*h2+g2*h3-h1)/(g1+g2)
            p32 = h2 - p22
            p23 = h3 - p22
            itpv = it0
         end if
 
c converts xii and yii to u-v system.
         dx = xii - x0
         dy = yii - y0
         u = ap*dx + bp*dy
         v = cp*dx + dp*dy
 
c evaluates the polynomial.
         p0 = p00 + v*(p01+v*(p02+v*(p03+v*(p04+v*p05))))
         p1 = p10 + v*(p11+v*(p12+v*(p13+v*p14)))
         p2 = p20 + v*(p21+v*(p22+v*p23))
         p3 = p30 + v*(p31+v*p32)
         p4 = p40 + v*p41
         zii = p0 + u*(p1+u*(p2+u*(p3+u*(p4+u*p5))))
         return
      else
 
c calculation of zii by extrapolation was removed by James, May 00,
c because we don't want to extrapolate for the IFU images!
         zii = 0.0
         return

      end if
      end
