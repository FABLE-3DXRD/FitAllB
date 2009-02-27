import sys
import ImageD11.columnfile as ic

#usage python omega_rotate.py in.flt out.flt

infile = sys.argv[1]
outfile = sys.argv[2]
flt = ic.columnfile(infile)
omega = flt.getcolumn('omega')
Min_o = flt.getcolumn('Min_o')
Max_o = flt.getcolumn('Max_o')
IMax_o = flt.getcolumn('IMax_o')
offset = -360.
flt.setcolumn(omega+offset,'omega')
flt.setcolumn(Min_o+offset,'Min_o')
flt.setcolumn(Max_o+offset,'Max_o')
flt.setcolumn(IMax_o+offset,'IMax_o')
flt.writefile(outfile)

