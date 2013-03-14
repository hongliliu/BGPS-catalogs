import plfit
import numpy as np
from agpy import readcol
import atpy

def bgps_to_reg(filename,outfile):
    bolocat = atpy.Table(filename)

    outfile = open(outfile,'w')
    print >>outfile,"galactic"

    for ii in xrange(len(bolocat)):
        flux = float(bolocat.flux_40[ii])
        if flux < 0.2: color='darkblue'
        elif flux < 0.5: color='blue'
        elif flux < 1.0: color='cyan'
        elif flux < 2.0: color='turquoise'
        elif flux < 3.0: color='green'
        elif flux < 5.0: color='yellow'
        elif flux < 10.0: color='orange'
        elif flux < 20.0: color='red'
        elif flux < 50.0: color='maroon'
        else: color='purple'

        try:
            major = bolocat.maj[ii]
            minor = bolocat.min[ii]
            posang = bolocat.pos_ang[ii] + 90.
        except AttributeError:
            major = bolocat.mommaj_as[ii]
            minor = bolocat.mommin_as[ii]
            posang = bolocat.posang[ii] + 90.

        print >>outfile,"ellipse(%s,%s,%s\",%s\",%s) # text={%s} width=%i color=%s" % (bolocat.glon[ii],bolocat.glat[ii],major,minor,posang,bolocat.name[ii],np.round(np.log(flux+1)+1),color)

    outfile.close()

if __name__ == "__main__":
    import optparse

    parser=optparse.OptionParser()
    #parser.add_option("--filename","-f",help="Input Filename",default="/Users/adam/work/catalogs/bolocam_gps_v1_0_1.reg")
    #parser.add_option("--outname","-o",help="Output Filename",default="/Users/adam/work/catalogs/bgps_ellipses.reg")
    options,args = parser.parse_args()

    bgps_to_reg(*args)#,options.filename,options.outname)
