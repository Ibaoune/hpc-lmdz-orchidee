      subroutine minmaxsource(zq,qmin,qmax,comment)

      USE dimphy
      USE infotrac

#include "dimensions.h"

!      character*20 comment
      character*(*) comment
      real qmin,qmax
      real zq(klon,nbtr)

      integer ijmin, lmin, ijlmin
      integer ijmax, lmax, ijlmax

      integer ismin,ismax

      ijlmin=ismin(klon*nbtr,zq,1)
      lmin=(ijlmin-1)/klon+1
      ijmin=ijlmin-(lmin-1)*klon
      zqmin=zq(ijmin,lmin)

      ijlmax=ismax(klon*nbtr,zq,1)
      lmax=(ijlmax-1)/klon+1
      ijmax=ijlmax-(lmax-1)*klon
      zqmax=zq(ijmax,lmax)

      if(zqmin.lt.qmin.or.zqmax.gt.qmax) 
     s     write(*,9999) comment,
     s     ijmin,lmin,zqmin,ijmax,lmax,zqmax

      return
9999  format(a20,2('  q(',i4,',',i2,')=',e12.5))
      end
