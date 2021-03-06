        parameter (maxve=200000)
        parameter (maxre=35)
        parameter (maxep=50)
        parameter (maxepr=15)
        parameter (maxdat=25)
        parameter (maxplv=1500)
        parameter (maxplp=700)
        parameter(pi=3.1415926535)
        parameter(rtod=180./pi)
        parameter(dtor=pi/180.)
        parameter (nulumax=100)
	parameter (nuseumax=100)
        parameter (maxq=35)
        parameter (maxx=201)
        parameter (maxy=201)
        parameter (maxz=201)
c	include 'snake-incl.f' version snake900 
        real sind, cosd, tand
        real asind, acosd, atand, atan2d
        external sind, cosd, tand
        external asind, acosd, atand, atan2d
        real*8 dsind, dcosd
        external dsind,dcosd
