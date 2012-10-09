"""
Find the closest SCUBA GPS point
"""
import numpy as np

def nearest_source(x1,y1,x2,y2):
    d = np.sqrt((x1-x2)**2 + (y1-y2)**2)
    return d.argmin(),d.min()

import atpy
from agpy import readcol

bgps = atpy.Table('bgps_iras_langston_match.tbl',type='ascii')
magpis = readcol('/Users/adam/work/catalogs/merge6_20.cat',fixedformat=[ 7,  7, 13, 13,  6,  9, 10,  8,  7,  7,  6,  5,  9, 10,  8,  9,  2],asStruct=True,skipline=10,nullval='          ')

scuba = atpy.Table('J_ApJS_175_277_table2.dat.fits')

dbest,bestmatch = np.zeros(len(bgps)),np.zeros(len(bgps))
for ii in xrange(len(bgps)):
    bestmatch[ii],dbest[ii] = nearest_source(bgps.glon_peak[ii],bgps.glat_peak[ii],scuba.GLON,scuba.GLAT)

bgps.add_column('F850',scuba['F850'][bestmatch.astype('int')],unit="Jy/beam")
bgps.add_column('P850',scuba['P850'][bestmatch.astype('int')],unit="Jy")
bgps.add_column('F450',scuba['F450'][bestmatch.astype('int')],unit="Jy/beam")
bgps.add_column('P450',scuba['P450'][bestmatch.astype('int')],unit="Jy")
bgps.add_column('BGPS-Scubadist',dbest,unit='degrees')

bgps.write("/Users/adam/work/catalogs/bgps_iras_langston_scuba_match.tbl",overwrite=True)
