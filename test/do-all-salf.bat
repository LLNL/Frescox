#***********************************************************************
# 
#    Copyright 2018, I.J. Thompson
#
#    This file is part of FRESCOX.
#
#    FRESCOX is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    FRESCOX is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with FRESCOX. If not, see <http://www.gnu.org/licenses/>.
#
#    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
#    LICENSE
#
#    The precise terms and conditions for copying, distribution and
#    modification are contained in the file COPYING.
#
#***********************************************************************
del frescox.in

frescox < be11.in > be11.out

frescox < cadprsc.in > capdrsc.out

copy cadprsc-n.in frescox.in
frescox  > capdrsc-n.out
del frescox.in

frescox < f19xfs.in > f19xfs.out

copy f19xfs-n.in frescox.in
frescox > f19xfs-n.out
del frescox.in

frescox < lane20.in > lane20.out

frescox < on2.in > on2.out

frescox < xeta.in > xeta.out

frescox < e80f49b.in > e80f49b.out

copy e80f49b-n.in frescox.in
frescox > e80f49b-n.out
del frescox.in

frescox < b8ex.in > b8ex.out

frescox < dex.in > dex.out
frescox < dexcc.in > dexcc.out

frescox < li7p-EM.in > li7p-EM.out

sfrescox <  na.min > na.out
sfrescox <  ss.min > ss.out
\rm for*.*
del for*.*
