with open("long_in.dat",'r') as f:
	v = float(f.readline().strip())

new_val = v**2

f = open('long_out.dat','w')
f.write(' {0:20.8E} \n'.format(new_val) )
f.close()