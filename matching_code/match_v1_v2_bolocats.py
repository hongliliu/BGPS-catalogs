import atpy # version 0.9.7
import numpy as np
import pyfits

import sys
import pyspherematch

if 'bolocatv1' not in locals():
    bolocatv1 = atpy.Table('../v1.0/bolocam_gps_v1_0_1.tbl', type='ipac')


if 'bgpsv2_nonew' not in locals():
    bgpsv2_nonew = atpy.Table('../v2.0/bgps_v2.1.tbl', type='ipac')

def merge_on_dist(cat1, cat2, cat1x='glon_peak', cat1y='glat_peak',
        cat2x='glon_peak', cat2y='glat_peak', cat1prefix='', cat2prefix=''):
    """
    Merge cat2 onto cat1 based on nearest match to cat1 source
    """

    # set up dtypes for new recordarray
    # cat1.dtypes = not valid (why?).  cat1[0] instead...
    dtypes = (
            [(cat1prefix+fn,dt[0]) for fn,dt in cat1[0].dtype.fields.iteritems()] + 
            [(cat2prefix+fn,dt[0]) for fn,dt in cat2[0].dtype.fields.iteritems()] +
            [('v1v2dist','float'),('v1v2dx','float'),('v1v2dy','float')])
    newrec = np.empty(len(cat1), dtype=dtypes)

    # copy over all fields from original "cat1"
    for field in cat1[0].dtype.fields:
        newrec[field] = cat1[field]

    for ii,(x1,y1) in enumerate(zip(cat1[cat1x],cat1[cat1y])):

        dx = x1 - cat2[cat2x]
        dy = y1 - cat2[cat2y]
        dd = ( dx**2 + dy**2 ) **0.5
        closest = dd.argmin()
        
        newrec['v1v2dist'][ii] = dd[closest]
        newrec['v1v2dx'][ii] = dx[closest]
        newrec['v1v2dy'][ii] = dy[closest]
        for field in cat2[0].dtype.fields:
            newrec[cat2prefix+field][ii] = cat2[field][closest]

    return newrec

v1v2_onv1 = merge_on_dist(bolocatv1,bgpsv2_nonew,cat2prefix='v2',cat2x='glon_peak',cat2y='glat_peak')
v1v2_onv2 = merge_on_dist(bgpsv2_nonew,bolocatv1,cat2prefix='v1',cat1x='glon_peak',cat1y='glat_peak')
v1v2_onv1_cen = merge_on_dist(bolocatv1,bgpsv2_nonew,cat2prefix='v2',cat2x='glon_cen',cat2y='glat_cen',cat1x='glon_cen',cat1y='glat_cen')
v1v2_onv2_cen = merge_on_dist(bgpsv2_nonew,bolocatv1,cat2prefix='v1',cat1x='glon_cen',cat1y='glat_cen',cat2x='glon_cen',cat2y='glat_cen')

matchnames_peak_v1 = []
matchdist_peak_v1 = []
matchxdiff_peak_v1 = []
matchydiff_peak_v1 = []
matchcatnum_peak_v1 = []
v1matches_peak = bolocatv1.flux*0

glonv1_peak = bolocatv1.glon_peak - (360*(bolocatv1.glon_peak > 345))
glonv2_peak = bgpsv2_nonew.glon_peak - (360*(bgpsv2_nonew.glon_peak > 345))

# find nearest v1 source to each v2 source
indv2v1, matchindv2v1, distv2v1 = pyspherematch.spherematch(bgpsv2_nonew['glon_peak'],bgpsv2_nonew['glat_peak'],bolocatv1['glon_peak'],bolocatv1['glat_peak'])

for name, catnum, glon_peak, glat_peak, ii, mi, dd in \
        zip(bgpsv2_nonew['name'], bgpsv2_nonew['cnum'], glonv2_peak,
                bgpsv2_nonew['glat_peak'], indv2v1, matchindv2v1, distv2v1):
    # finding nearest v1 to v2

    dx = (glonv1_peak[mi]-glon_peak)
    dy = (bolocatv1['glat_peak'][mi]-glat_peak)
    #dd = ( dx**2 + dy**2 ) **0.5
    #closest = dd.argmin()

    matchnames_peak_v1.append(bolocatv1.name[mi])
    matchdist_peak_v1.append(dd)
    matchxdiff_peak_v1.append(dx)
    matchydiff_peak_v1.append(dy)
    matchcatnum_peak_v1.append(bolocatv1['cnum'][mi])
    v1matches_peak[mi]+=1

print "There are %i v2 sources in the v1/v2 overlap region" % (len(bgpsv2_nonew))
print "Of these, %i have a neighbor within 100 arcsec (%i do not)" % ((distv2v1 < 100./3600.).sum(),len(bgpsv2_nonew)-(distv2v1 < 100./3600.).sum())
print "Out of %i v1 sources, %i have a unique match in v2" % (len(bolocatv1), len(np.unique(matchindv2v1)))
print "There are %i v1 sources with a match in v2" % ((v1matches_peak>0).sum())
    
#duplicates_peak_v1 = np.array([(cn,matchcatnum_peak_v1.count(cn)) for cn in set(matchcatnum_peak_v1) if matchcatnum_peak_v1.count(cn) > 1])


# find v2 matches to v1 now
#
matchnames_peak_v2 = []
matchdist_peak_v2 = []
matchxdiff_peak_v2 = []
matchydiff_peak_v2 = []
matchcatnum_peak_v2 = []
v2matches_peak = bgpsv2_nonew.flux*0


# find nearest v2 source to each v1 source
indv1v2, matchindv1v2, distv1v2 = pyspherematch.spherematch(bolocatv1['glon_peak'],bolocatv1['glat_peak'],bgpsv2_nonew['glon_peak'],bgpsv2_nonew['glat_peak'])

for name, catnum, glon_peak, glat_peak, ii, mi, dd in \
        zip(bolocatv1.name, bolocatv1['cnum'], glonv1_peak,
                bolocatv1['glat_peak'], indv1v2, matchindv1v2, distv1v2):
    # finding nearest v2 to v1

    dx = (glonv2_peak[mi]-glon_peak)
    dy = (bgpsv2_nonew['glat_peak'][mi]-glat_peak)
    #dd = ( dx**2 + dy**2 ) **0.5
    #closest = dd.argmin()

    matchnames_peak_v2.append(bgpsv2_nonew.name[mi])
    matchdist_peak_v2.append(dd)
    matchxdiff_peak_v2.append(dx)
    matchydiff_peak_v2.append(dy)
    matchcatnum_peak_v2.append(bgpsv2_nonew['cnum'][mi])
    v2matches_peak[mi]+=1

print "There are %i v1 sources in the v2/v1 overlap region" % (len(bolocatv1))
print "Of these, %i have a neighbor within 100 arcsec (%i do not)" % ((distv2v1 < 100./3600.).sum(),len(bolocatv1)-(distv1v2 < 100./3600.).sum())
print "Out of %i v2 sources, %i have a unique match in v1" % (len(bgpsv2_nonew), len(np.unique(matchindv1v2)))
print "There are %i v2 sources with a match in v1" % ((v2matches_peak>0).sum())
    
#duplicates_peak_v2 = np.array([(cn,matchcatnum_peak_v2.count(cn)) for cn in set(matchcatnum_peak_v2) if matchcatnum_peak_v2.count(cn) > 1])

# where v1 and v2 point to each other
double_match = [(v1ind,v2ind) for v2ind,v1ind in zip(indv2v1,matchindv2v1) if matchindv1v2[v1ind]==v2ind]
print "There are %i sources with unique matches between v1 and v2" % (len(double_match))





glonv1_cen = bolocatv1.glon_cen - (360*(bolocatv1.glon_cen > 345))
glonv2_cen = bgpsv2_nonew.glon_cen - (360*(bgpsv2_nonew.glon_cen > 345))

matchnames_cen_v1 = []
matchdist_cen_v1 = []
matchxdiff_cen_v1 = []
matchydiff_cen_v1 = []
matchcatnum_cen_v1 = []
v1matches_cen = bolocatv1.flux*0

for name,catnum,glon_cen,glat_cen in zip(bgpsv2_nonew.name,bgpsv2_nonew.cnum,glonv2_cen,bgpsv2_nonew.glat_cen):

    dx = (glonv1_cen-glon_cen)
    dy = (bolocatv1.glat_cen-glat_cen)
    dd = ( dx**2 + dy**2 ) **0.5
    closest = dd.argmin()

    matchnames_cen_v1.append(bolocatv1.name[closest])
    matchdist_cen_v1.append(dd.min())
    matchxdiff_cen_v1.append(dx[dd.argmin()])
    matchydiff_cen_v1.append(dy[dd.argmin()])
    matchcatnum_cen_v1.append(bolocatv1.cnum[closest])
    v1matches_cen[closest]+=1

duplicates_cen_v1 = np.array([(cn,matchcatnum_cen_v1.count(cn)) for cn in set(matchcatnum_cen_v1) if matchcatnum_cen_v1.count(cn) > 1])



for name,catnum,glon_peak,glat_peak in zip(bolocatv1.name,bolocatv1.cnum,glonv1_peak,bolocatv1.glat_peak):
    # finding nearest v2 source to a v1 source

    dx = (glonv2_peak-glon_peak)
    dy = (bgpsv2_nonew.glat_peak-glat_peak)
    dd = ( dx**2 + dy**2 ) **0.5
    closest = dd.argmin()

    matchnames_peak_v2.append(bgpsv2_nonew.name[closest])
    matchdist_peak_v2.append(dd.min())
    matchxdiff_peak_v2.append(dx[dd.argmin()])
    matchydiff_peak_v2.append(dy[dd.argmin()])
    matchcatnum_peak_v2.append(bgpsv2_nonew.cnum[closest])
    v2matches_peak[closest]+=1

    
duplicates_peak_v2 = np.array([(cn,matchcatnum_peak_v2.count(cn)) for cn in set(matchcatnum_peak_v2) if matchcatnum_peak_v2.count(cn) > 1])


matchnames_cen_v2 = []
matchdist_cen_v2 = []
matchxdiff_cen_v2 = []
matchydiff_cen_v2 = []
matchcatnum_cen_v2 = []
v2matches_cen = bgpsv2_nonew['flux']*0

for name,catnum,glon_cen,glat_cen in zip(bolocatv1.name,bolocatv1.cnum,glonv1_cen,bolocatv1.glat_cen):

    dx = (glonv2_cen-glon_cen)
    dy = (bgpsv2_nonew.glat_cen-glat_cen)
    dd = ( dx**2 + dy**2 ) **0.5
    closest = dd.argmin()

    matchnames_cen_v2.append(bgpsv2_nonew.name[closest])
    matchdist_cen_v2.append(dd.min())
    matchxdiff_cen_v2.append(dx[dd.argmin()])
    matchydiff_cen_v2.append(dy[dd.argmin()])
    matchcatnum_cen_v2.append(bgpsv2_nonew.cnum[closest])
    v2matches_cen[closest]+=1

duplicates_cen_v2 = np.array([(cn,matchcatnum_cen_v2.count(cn)) for cn in set(matchcatnum_cen_v2) if matchcatnum_cen_v2.count(cn) > 1])

from pylab import *
rcParams['font.size']=24
figure(1)
clf()
title("v1 sources with no v2 match")
plot(glonv1_peak[v1matches_peak==0], bolocatv1.glat_peak[v1matches_peak==0], 'b.')
plot(glonv1_peak, bolocatv1.glat_peak, 'b.')

for ii,v2entry in enumerate(bgpsv2_nonew):
    if v2entry['name'] not in matchnames_peak_v2:
        v1match = bolocatv1[matchnames_peak_v1[ii]==bolocatv1.name][0]
        plot([v1match['glon_peak']-(360*(v1match['glon_peak']>345)),v2entry['glon_peak']-(360*(v2entry['glon_peak']>345))],[v1match['glat_peak'],v2entry['glat_peak']], color='r')
gca().axis([-10,50,-0.5,0.5])

figure(2)
clf()
title("v2 sources with no v1 match")
plot(glonv2_peak[v2matches_peak==0], bgpsv2_nonew.glat_peak[v2matches_peak==0], 'r.')
plot(glonv2_peak, bgpsv2_nonew.glat_peak, 'r.')
for ii,v1entry in enumerate(bolocatv1):
    if v1entry['name'] not in matchnames_peak_v1:
        v2match = bgpsv2_nonew[matchnames_peak_v2[ii]==bgpsv2_nonew.name][0]
        plot([v2match['glon_peak']-(360*(v2match['glon_peak']>345)),v1entry['glon_peak']-(360*(v1entry['glon_peak']>345))],[v2match['glat_peak'],v1entry['glat_peak']], color='b')
gca().axis([-10,50,-0.5,0.5])


figure(3)
clf()
bins = np.logspace(np.log10(0.01),np.log10(100))
hist(bgpsv2_nonew.flux_40[np.array(double_match)[:,1]], bins=bins, label="v1-v2 pairs",
        log=True, facecolor=(1,0.8,0.8,0.1), linewidth=0,
        edgecolor=(1,0,0,0.5), histtype='stepfilled')
hist(bgpsv2_nonew.flux_40[v2matches_peak>0], bins=bins, label="v2 with v1 match",
        log=True, facecolor=(1,0.8,0.8,0.1), linewidth=4,
        edgecolor=(1,0,0,0.5), histtype='step')
hist(bgpsv2_nonew.flux_40[v2matches_peak==0], bins=bins, label="v2 (no v1 match)",
        alpha=0.5,log=True,facecolor='none', linewidth=4, edgecolor='b', linestyle="dotted",
        histtype='step')
hist(bolocatv1.flux_40[v1matches_peak==0]*1.5,  bins=bins,  label="v1 (no v2 match)", alpha=0.5, log=True, facecolor='none',  
    linewidth=4, edgecolor='g', linestyle='dashed',
    histtype='step')
gca().set_xscale('log')
gca().set_yscale('log')
gca().set_ylim(0.5,gca().get_ylim()[1])
xlabel('Flux Density in 40\" apertures (Jy)',fontsize='24')
ylabel('Number of sources',fontsize='24')
legend(loc='best',prop={'size':24})

savefig('/Users/adam/work/bolocam/bolocat/bolocat_match_v1v2_histogram.png',bbox_inches='tight')
savefig('/Users/adam/work/bolocam/bolocat/bolocat_match_v1v2_histogram.pdf',bbox_inches='tight')
draw()

# create a merged catalog...
v2withv1matchesbypeak = [(v2entry,bolocatv1[bolocatv1.cnum==v1catnum],v1dist,dx,dy)
        for (v2entry,v1name,v1dist,dx,dy,v1catnum) in zip(
                                bgpsv2_nonew,
                                matchnames_peak_v1,
                                matchdist_peak_v1,
                                matchxdiff_peak_v1,
                                matchydiff_peak_v1,
                                matchcatnum_peak_v1)]
dtypes = (
        [(fn,dt[0]) for fn,dt in (v2entry.dtype.fields.iteritems())] + 
        [('v1'+fn,dt[0]) for fn,dt in bolocatv1[0].dtype.fields.iteritems()] +
        [('v1v2dist','float'),('v1v2dx','float'),('v1v2dy','float')])
newrec = np.empty(len(v2withv1matchesbypeak), dtype=dtypes)
for field in bgpsv2_nonew[0].dtype.fields:
    newrec[field] = bgpsv2_nonew[field]
for ii,line in enumerate(v2withv1matchesbypeak):
    v1entry = line[1]
    for field in v1entry.dtype.fields:
        newrec['v1'+field][ii] = v1entry[field][0]
    newrec['v1v2dist'][ii] = line[2]
    newrec['v1v2dx'][ii] = line[3]
    newrec['v1v2dy'][ii] = line[4]
v2base_v1v2merge = newrec

figure(4)
clf()
hist(v1v2_onv2['v1v2dist']*3600,logspace(-4,1)*3600)
vlines(40,gca().get_ylim()[0],gca().get_ylim()[1],linestyle='--',color='r')
gca().set_xscale('log')
xlabel("V1-V2 distance (arcsec)")
title("V2 catalog with V1 nearest-matches")

figure(5)
clf()
hist(v1v2_onv1['v1v2dist']*3600,logspace(-4,1)*3600)
vlines(40,gca().get_ylim()[0],gca().get_ylim()[1],linestyle='--',color='r')
gca().set_xscale('log')
xlabel("V1-V2 distance (arcsec)")
title("V1 catalog with V2 nearest-matches")


def filenames_to_fieldnames(filenames):
    import os
    return [os.path.basename(fn).replace("v2.0_ds2_","").replace("_13pca_map20.fits","")
            for fn in filenames]

def compare_distances(cat=v1v2_onv2,fig1=6,
        titlestr="V2 catalog with V1 nearest-matches"):
    figure(fig1)
    clf()
    plot(cat['v1v2dx']*3600,cat['v1v2dy']*3600,'o',mec='none',mfc='r',alpha=0.05)
    plot(cat['v1v2dx'].mean()*3600,cat['v1v2dy'].mean()*3600,'gx',mew=2)
    xlabel('v1-v2 x distance (arcsec)')
    ylabel('v1-v2 y distance (arcsec)')
    title(titlestr)

    lt40dist = cat['v1v2dist']*3600 < 40
    cat_lt40 = cat[lt40dist]
    plot(cat_lt40['v1v2dx'].mean()*3600,cat_lt40['v1v2dy'].mean()*3600,'b+',mew=2)
    axis([-1000,1000,-1000,1000])

    figure(fig1+1)
    clf()
    plot(cat_lt40['v1v2dx']*3600,cat_lt40['v1v2dy']*3600,'o',mec='none',mfc='b',alpha=0.1)
    plot(cat_lt40['v1v2dx'].mean()*3600,cat_lt40['v1v2dy'].mean()*3600,'rx',mew=2)
    xlabel('v1-v2 x distance (arcsec)')
    ylabel('v1-v2 y distance (arcsec)')
    title(titlestr+"(only d<40)")
    print titlestr
    print "Mean offset: dx=%f dy=%f" % (cat_lt40['v1v2dx'].mean()*3600,cat_lt40['v1v2dy'].mean()*3600)

atlasgal_overlap = (v1v2_onv2['glon_peak'] < 20) * (v1v2_onv2['glon_peak'] > -10)

def recarray_to_table(recarray):
    T = atpy.Table()
    for field,(dtype,throwaway) in recarray.dtype.fields.iteritems():
        T.add_column(field, recarray[field], dtype=dtype)
    return T

recarray_to_table(v1v2_onv2).write('../cross/v1v2_onv2_peak_table.fits',overwrite=True)
recarray_to_table(v1v2_onv1).write('../cross/v1v2_onv1_peak_table.fits',overwrite=True)
recarray_to_table(v1v2_onv2_cen).write('../cross/v1v2_onv2_cen_table.fits',overwrite=True)
recarray_to_table(v1v2_onv1_cen).write('../cross/v1v2_onv1_cen_table.fits',overwrite=True)
