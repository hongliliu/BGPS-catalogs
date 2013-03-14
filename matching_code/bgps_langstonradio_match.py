"""
Find the closest Langston 2000 8/14 GHz point, report flux
"""
import numpy as np

def nearest_source(x1,y1,x2,y2):
    d = np.sqrt((x1-x2)**2 + (y1-y2)**2)
    return d.argmin(),d.min()

import atpy
from agpy import readcol

bgps = atpy.Table('bgps_iras_match.tbl',type='ascii')

for frequency in (8,14):
    langston = atpy.Table('J_AJ_119_2801_s%i.dat.gz.fits' % frequency)
    glon_lang,glat_lang = zip(*[(float(coord[:6]),float(coord[6:])) for coord in langston['[LMD2000b]']])

    dbest,bestmatch = np.zeros(len(bgps)),np.zeros(len(bgps))
    for ii in xrange(len(bgps)):
        bestmatch[ii],dbest[ii] = nearest_source(bgps.glon_peak[ii],bgps.glat_peak[ii],glon_lang,glat_lang)

    bgps.add_column('Speak%iGHz' % frequency,langston['Sp'][bestmatch.astype('int')],unit="Jy")
    bgps.add_column('Sint%iGHz' % frequency,langston['Sint'][bestmatch.astype('int')],unit="Jy")
    bgps.add_column('MAJD%iGHz' % frequency,langston['MajDiam'][bestmatch.astype('int')],unit="Arcmin")
    bgps.add_column('MIND%iGHz' % frequency,langston['MinDiam'][bestmatch.astype('int')],unit="Arcmin")
    bgps.add_column('BGPS-Lang%idist' % frequency,dbest,unit='degrees')
    try:
        bgps.add_column('SpInd',langston['Sp+Index'][bestmatch.astype('int')])
    except ValueError:
        pass

bgps.write("/Users/adam/work/catalogs/bgps_iras_langston_match.tbl",overwrite=True)
