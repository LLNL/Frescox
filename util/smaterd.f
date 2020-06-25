!***********************************************************************
! 
!    Copyright 2018, I.J. Thompson
!
!    This file is part of FRESCOX.
!
!    FRESCOX is free software; you can redistribute it and/or
!    modify it under the terms of the GNU General Public License
!    as published by the Free Software Foundation; either version 2
!    of the License, or (at your option) any later version.
!    
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!    
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
!    MA  02110-1301, USA.
!
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!
!    The precise terms and conditions for copying, distribution and
!    modification are contained in the file COPYING.
!
!***********************************************************************
      implicit real*8(a-h,o-z)
      parameter(MXP=1, MXX=1, LMAX=1010, MIJ=4)
      complex*16 s,smat(MXP,MXX,0:LMAX,0:MIJ,0:MIJ,0:MIJ,0:MIJ),sl
      real*8 jtot,jin,j,jp(MXP,MXX),jt(MXP,MXX),jpin,jtin
      integer p,e,pmax,emax(MXP),jttmax(MXP,MXX)
      logical some
      nint(x) = int(x+0.5)

      pmax=0     
      do 3 p=1,MXP
      do 2 e=1,MXX
      do 1 jtt=0,LMAX
      do 1 ijin=0,MIJ
      do 1 ilin=0,MIJ
      do 1 ij=0,MIJ
      do 1 il=0,MIJ
1     smat(p,e,jtt,ijin,ilin,ij,il) = 0.0
2     jttmax(p,e) = 0
3     emax(p) = 0

4     read(5,5,end=20) s,l,j,jtot,p,e,lin,jin,
     *                 jp(p,e),jt(p,e),jpin,jtin
5     format(2f15.10,i6,2f6.1,2i3,i6,f6.1,4f4.1)
      jtt = int(jtot)
      ijin = jin-(jtot-jtin)
      ilin = lin-(jin-jpin)
      ij   = j-(jtot-jt(p,e))
      il   = l-(j-jp(p,e))
      if(jtt.gt.LMAX.or.p.gt.MXP.or.e.gt.MXX.or.
     *                  max(ijin,ilin,ij,il).gt.MIJ) go to 10
      smat(p,e,jtt,ijin,ilin,ij,il) = s
      pmax = max(p,pmax)
      emax(p) = max(e,emax(p))
      jttmax(p,e) = max(jtt,jttmax(p,e))
      go to 4
10    write(0,*) 'Insufficient arrays!'
      write(0,*) jtt,'>',LMAX,' or ',p,'>',MXP,' or ',
     *            e,'>',MXX,' or ',max(ijin,ilin,ij,il),'>',MIJ
      stop
20    do 40 p=1,pmax
      do 40 e=1,emax(p)
      write(6,*) '# From partition ',p,', excitation pair ',e
      do 35 ijin=0,nint(2.*jtin)
      do 35 ilin=0,nint(2.*jpin)
      do 35 ij=0,nint(2.*jt(p,e))
      do 35 il=0,nint(2.*jp(p,e))    
      write(6,*) '# ',jttmax(p,e)+1,' S diffs for IJIN,ILIN,IJ,IL =',
     *           ijin,ilin,it,il
      some=.false.
      sl = smat(p,e,0,ijin,ilin,ij,il) 
      do 30 jtt=1,jttmax(p,e)
      s = smat(p,e,jtt,ijin,ilin,ij,il) 
      if(abs(s).gt.1e-20.and.jtt.ge.2) then
	  write(6,32) jtt,s-sl
          some=.true.
	  endif
       sl = s
30    continue
      if(some) write(6,*) '&'
32    format(I6,2F15.10)
35    continue
40    continue
      stop
      end