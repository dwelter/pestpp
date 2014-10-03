	program stor

	implicit none

	real recharge, conductance, storage, inithead, time,
     +  coeff, factor

	open(unit=10,file='input.dat',status='old')
	open(unit=20,file='output.dat')
	read(10,*) recharge, conductance, storage
	read(10,*) inithead
	read(10,*)

	coeff=(recharge/conductance-inithead)
	factor=conductance/storage
	write(20,10)
10	format('   Time           Water_Level')
	do 
	  read(10,*,end=100) time
	  write(20,20) time,inithead+coeff*(1.0-exp(-factor*time))
20	  format(1x,1pg14.7,2x,1pg14.7)
	end do

100	end
