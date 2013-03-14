"""
Find the closest magpis GPS point
"""
import numpy as np

def nearest_source(x1,y1,x2,y2):
    d = np.sqrt((x1-x2)**2 + (y1-y2)**2)
    return d.argmin(),d.min()

import atpy
from agpy import readcol

bgps = atpy.Table('bgps_iras_langston_scuba_match.tbl',type='ascii')
magpis = readcol('/Users/adam/work/catalogs/merge6_20.cat',fixedformat=[ 7,  7, 13, 13,  6,  9, 10,  8,  7,  7,  6,  5,  9, 10,  8,  9,  2],asStruct=True,skipline=10,nullval='          ')


dbest,bestmatch = np.zeros(len(bgps)),np.zeros(len(bgps))
for ii in xrange(len(bgps)):
    bestmatch[ii],dbest[ii] = nearest_source(bgps.glon_peak[ii],bgps.glat_peak[ii],magpis.Long,magpis.Lat)

bgps.add_column('Fint6',magpis.Fint6[bestmatch.astype('int')],unit="Jy/beam")
bgps.add_column('Fpeak6',magpis.Fpeak6[bestmatch.astype('int')],unit="Jy")
bgps.add_column('Fint20',magpis.Fint20[bestmatch.astype('int')],unit="Jy/beam")
bgps.add_column('Fpeak20',magpis.Fpeak20[bestmatch.astype('int')],unit="Jy")
bgps.add_column('BGPS-magpisdist',dbest,unit='degrees')

bgps.write("/Users/adam/work/catalogs/bgps_iras_langston_scuba_magpis_match.tbl",overwrite=True)
