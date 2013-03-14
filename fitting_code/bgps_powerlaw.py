from agpy import plfit
from agpy import readcol

if not globals().has_key('dotests'): dotests = raw_input('Do tests? ') not in ['','n','N']

rcParams['font.size'] = 36

bcnames,bgps = readcol('/Users/adam/work/catalogs/bolocam_gps_v1_0.tbl',names=True,skipline=45)
bgpsd = readcol('/Users/adam/work/catalogs/bolocam_gps_v1_0.tbl',names=True,skipline=45,asdict=True)

bgflux = bgps[:,18].astype('float')
bgflux40 = bgps[:,12].astype('float')
bgrad = bgps[:,11]
bgrad = bgrad[bgrad!='null'].astype('float')

print "plBGFLUX"
plBGFLUX = plfit(bgflux,nosmall=True)
if dotests: plBGFLUXp,plBGFLUXksarr = plBGFLUX.test_pl()
# xmin: 4.877  n(>xmin): 369  alpha: 2.60503 +/- 0.0835548  Likelihood: -1009  ks: 0.943031
# p(2500) = .52
figure(1); clf()
xlabel('Flux Density (Jy)'); ylabel('N(S)')
plBGFLUX.plotpdf(nbins=40,dolog=True)
savefig('/Users/adam/work/massfunc/BGPS_bolocat_flux_Nhist.png',papertype='b3')
figure(2); clf();
xlabel('Flux Density (Jy)'); ylabel('$\Delta$N/$\Delta$S')
plBGFLUX.plotpdf(nbins=40,dolog=True,dnds=True)
savefig('/Users/adam/work/massfunc/BGPS_bolocat_flux_dndshist.png',papertype='b3')
clf(); plBGFLUX.alphavsks()
savefig('/Users/adam/work/massfunc/BGPS_bolocat_flux_alphaks.png',papertype='b3')


print "plBGFLUX40"
plBGFLUX40 = plfit(bgflux40,nosmall=True)
if dotests: plBGFLUX40p,plBGFLUX40ksarr = plBGFLUX40.test_pl()
# xmin: 0.293  n(>xmin): 1582  alpha: 2.40205 +/- 0.03525  Likelihood: -233.702  ks: 0.985827
# p(2500) = 0.218
figure(1); clf()
xlabel('Flux Density (Jy)'); ylabel('N(S)')
plBGFLUX40.plotpdf(nbins=40,dolog=True)
savefig('/Users/adam/work/massfunc/BGPS_bolocat_flux40_Nhist.png',papertype='b3')
figure(2); clf();
xlabel('Flux Density (Jy)'); ylabel('$\Delta$N/$\Delta$S')
plBGFLUX40.plotpdf(nbins=40,dolog=True,dnds=True)
savefig('/Users/adam/work/massfunc/BGPS_bolocat_flux40_dndshist.png',papertype='b3')
clf(); plBGFLUX40.alphavsks()
savefig('/Users/adam/work/massfunc/BGPS_bolocat_flux40_alphaks.png',papertype='b3')


gctf=readcol('/Users/adam/work/catalogs/gc/Fig4_Mass_r40_Tflux.txt')
gctr=readcol('/Users/adam/work/catalogs/gc/Fig4_Mass_r40_Tradial.txt')
gc20=readcol('/Users/adam/work/catalogs/gc/Fig4_Mass_r40_T20.txt')
gcr120t20=readcol('/Users/adam/work/catalogs/gc/Fig4_Mass_r120_T20.txt')

print "GCTflux"
GCTflux = plfit(gctf[:,1],nosmall=True)
if dotests: GCTfluxp,GCTfluxksarr = GCTflux.test_pl()
# xmin: 4.877  alpha: 2.60503 +/- 0.0835548  Likelihood: -1009  ks: 0.943031
# 0.36359999999999998
print "GCTrad"
GCTrad = plfit(gctr[:,1],nosmall=True)
if dotests: GCTradp,GCTradksarr = GCTrad.test_pl()
# xmin: 198.1  alpha: 2.6885 +/- 0.102006  Likelihood: -1741.87  ks: 0.953082
# 0.35320000000000001
print "GCT20"
GCT20 = plfit(gc20[:,1],nosmall=True)
if dotests: GCT20p,GCT20ksarr = GCT20.test_pl()
# xmin: 368  alpha: 2.57735 +/- 0.10098  Likelihood: -1729.06  ks: 0.959864
# 0.4052
print "GCr120T20"
GCr120T20 = plfit(gcr120t20[:,1],nosmall=True)
if dotests: GCr120T20p,GCr120T20ksarr = GCr120T20.test_pl()
# xmin: 2294.8  alpha: 2.60296 +/- 0.113917  Likelihood: -1760.3  ks: 0.954278
# 0.388

figure(1); clf(); xlabel('Mass (M$_\odot$)'); ylabel('N(M)')
GCTflux.plotpdf(dolog=True)
savefig('/Users/adam/work/massfunc/GC_mass_Tflux_Nhist.png',papertype='b3')
figure(1); clf(); xlabel('Mass (M$_\odot$)'); ylabel('N(M)')
GCTrad.plotpdf(dolog=True)
savefig('/Users/adam/work/massfunc/GC_mass_Trad_Nhist.png',papertype='b3')
figure(1); clf(); xlabel('Mass (M$_\odot$)'); ylabel('N(M)')
GCT20.plotpdf(dolog=True)
savefig('/Users/adam/work/massfunc/GC_mass_T20_Nhist.png',papertype='b3')
figure(1); clf(); xlabel('Mass (M$_\odot$)'); ylabel('N(M)')
GCr120T20.plotpdf(dolog=True)
savefig('/Users/adam/work/massfunc/GC_mass_r120_T20_Nhist.png',papertype='b3')

figure(2); clf(); xlabel('Mass (M$_\odot$)'); ylabel('N($>$M)')
GCTflux.plotcdf()
savefig('/Users/adam/work/massfunc/GC_mass_Tflux_CDF.png',papertype='b3')
figure(2); clf(); xlabel('Mass (M$_\odot$)'); ylabel('N($>$M)')
GCTrad.plotcdf()
savefig('/Users/adam/work/massfunc/GC_mass_Trad_CDF.png',papertype='b3')
figure(2); clf(); xlabel('Mass (M$_\odot$)'); ylabel('N($>$M)')
GCT20.plotcdf()
savefig('/Users/adam/work/massfunc/GC_mass_T20_CDF.png',papertype='b3')
figure(2); clf(); xlabel('Mass (M$_\odot$)'); ylabel('N($>$M)')
GCr120T20.plotcdf()
savefig('/Users/adam/work/massfunc/GC_mass_r120_T20_CDF.png',papertype='b3')


figure(2); clf();
GCTflux.alphavsks()
savefig('/Users/adam/work/massfunc/GC_mass_Tflux_alphaks.png',papertype='b3')
figure(2); clf();
GCTrad.alphavsks()
savefig('/Users/adam/work/massfunc/GC_mass_Trad_alphaks.png',papertype='b3')
figure(2); clf();
GCT20.alphavsks();
savefig('/Users/adam/work/massfunc/GC_mass_T20_alphaks.png',papertype='b3')
figure(2); clf();
GCr120T20.alphavsks()
savefig('/Users/adam/work/massfunc/GC_mass_r120_T20_alphaks.png',papertype='b3')








