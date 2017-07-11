
# This makefile creates the 3D rbelt code
# setenv machine <machine name> (for c-shell)
# set machine=<machine name> (for BASH & other shells)
# also modify file make.$(machine) as needed

#include make.$(machine)
include make.atrus

.SUFFIXIES: .F .c .o .f
.c.o:
	$(cc) $(CFLAGS) -c $<
.F.o:
	$(COMPILER) $(FFLAGS) $(MPIINCDIR) $(HDFINCDIR) -c $<
.f.o:
	$(COMPILER) $(FFLAGS) $(MPIINCDIR) $(HDFINCDIR) -c $<

#Source code files

1p_dip_src = rbelt-main1t.f \
      rbelt-dipole.f \
      rbelt-nogrid.f \
      rbelt-dist.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-boris.f \
      rbelt-gcenter.f \
      rbelt-io.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      rbelt-coortrans.f \
      Geopack-2008.F \

1p_cdip_src = rbelt-main1t.f \
      rbelt-cdipole.f \
      rbelt-nogrid.f \
      rbelt-dist.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-boris.f \
      rbelt-gcenter.f \
      rbelt-io.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      rbelt-coortrans.f \
      Geopack-2008.F \

1p_ts05analytic_src = rbelt-main1t.f \
      rbelt-ts05.f \
      TS04c.f \
      rbelt-nogrid.f \
      rbelt-dist.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-io.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      Geopack-2008.F \

1p_hdf_src = rbelt-main.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-dist.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-boris.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp4d.f \
      rbelt-io.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      rbelt-coortrans.f \
      Geopack-2008.f \

1p_hdf1_src = rbelt-main1t.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-dist.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-io.f  \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      Geopack-2008.f \

1p_ts05_src = rbelt-main1t.f \
      rbelt-grid.f \
      rbelt-ts05_samp.f \
      TS04c.f \
      rbelt-dist.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-io.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      Geopack-2008.f \

rbelt_dip_src = rbelt-main1t.f \
      rbelt-dipole.f \
      rbelt-nogrid.f \
      rbelt-dist_disk.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-io.f  \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-coortrans.f \
      rbelt-lstar.f \
      Geopack-2008.F \

rbelt_cdip_src = rbelt-main1t.f \
      rbelt-cdipole.f \
      rbelt-nogrid.f \
      rbelt-dist_disk.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-io.f  \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-coortrans.f \
      rbelt-lstar.f \
      Geopack-2008.F \

rbelt_dst_src = rbelt-main1t.f \
      rbelt-dstfld.f \
      rbelt-nogrid.f \
      rbelt-dist_disk.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-io.f  \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-coortrans.f \
      rbelt-lstar.f \
      Geopack-2008.F \

rbelt_ulf_src = rbelt-main1t.f \
      rbelt-ulffld.f \
      rbelt-nogrid.f \
      rbelt-dist_disk.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-io.f  \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-coortrans.f \
      rbelt-lstar.f \
      Geopack-2008.F \

rbelt_hdf_src = rbelt-main.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-dist_disk.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-boris.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp4d.f \
      rbelt-io.f   \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-coortrans.f \
      rbelt-lstar.f \
      Geopack-2008.F \

rbelt_lfm_src = rbelt-main.f \
      rbelt-grid.f \
      rbelt-dist_tail.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp4d.f \
      rbelt-io.f   \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      Geopack-2008.F \
      rbelt-lfm_samp.f \
      hdftake_mpi.F \
      ionset.F \
      grid.F \
      cord.F \
      therest.F \
      jcalc2.F \
      flux.F \
      interp-all.F \
      intcoef1d.F \
      cart-int.F \

rbelt_hdf1_src = rbelt-main1t.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-dist_disk.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-io.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      Geopack-2008.F \

rbelt_ts05_src = rbelt-main.f \
      rbelt-grid.f \
      rbelt-ts05_samp.f \
      TS04c.f \
      rbelt-dist_tail.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp4d.f \
      rbelt-io.f   \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      Geopack-2008.F \

isopt_dip_src = rbelt-main1t.f \
      rbelt-dipole.f \
      rbelt-nogrid.f \
      rbelt-dist_isopt.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-io.f \
      rbelt-string.f \

isopt_hdf1_src = rbelt-main1t.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-dist_isopt.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-io.f \
      rbelt-string.f \
      rbelt-time.f \

isosphr_dip_src = rbelt-main1t.f \
      rbelt-dipole.f \
      rbelt-nogrid.f \
      rbelt-dist_sphr.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-io.f \
      rbelt-string.f \

isosphr_hdf_src = rbelt-main.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-dist_sphr.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp4d.f \
      rbelt-io.f \
      rbelt-string.f \
      rbelt-time.f \

isosphr_hdf1_src = rbelt-main1t.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-dist_sphr.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-io.f \
      rbelt-string.f \
      rbelt-time.f \

isosphr_ts05_src = rbelt-main1t.f \
      rbelt-grid.f \
      rbelt-ts05_samp.f \
      TS04c.f \
      rbelt-dist_sphr.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-io.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      Geopack-2008.F \

isosphr_nob_src = rbelt-main1t.f \
      rbelt-zero_fields.f \
      rbelt-nogrid.f \
      rbelt-dist_sphr.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-io.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      Geopack-2008.F \

isosphr_odip_src = rbelt-main1t.f \
      rbelt-offset_dip.f \
      rbelt-nogrid.f \
      rbelt-dist_sphr.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-io.f \
      rbelt-string.f \

isosphr_t96_src = rbelt-main1t.f \
      rbelt-t96.f \
      rbelt-nogrid.f \
      rbelt-dist_sphr.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-io.f \
      rbelt-string.f \
      T96.f \

gcr1pt_dip_src = rbelt-main1t.f \
      rbelt-gcr1pt_ana.f \
      rbelt-dipole.f \
      rbelt-nogrid.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-string.f \
      Geopack-2008.f \
      
gcr1pt_igrf_src = rbelt-main1t.f \
      rbelt-gcr1pt_ana.f \
      rbelt-igrf.f \
      rbelt-nogrid.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-string.f \
      Geopack-2008.f \

gcr1pt_hdf_src = rbelt-main1t.f \
      rbelt-gcr1pt_hdf.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-string.f \
      rbelt-time.f \
      Geopack-2008.f \

acs_dip_src = rbelt-main1t.f \
      rbelt-access.f \
      rbelt-dipole.f \
      rbelt-nogrid.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-string.f \
      Geopack-2008.f \
      
acs_igrf_src = rbelt-main1t.f \
      rbelt-access.f \
      rbelt-igrf.f \
      rbelt-nogrid.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-string.f \
      Geopack-2008.f \

acs_hdf_src = rbelt-main1t.f \
      rbelt-access.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-string.f \
      rbelt-time.f \
      Geopack-2008.f \

pen_dip_src = rbelt-main1t.f \
      rbelt-penumbra.f \
      rbelt-dipole.f \
      rbelt-nogrid.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-igrf_Lm.f \
      rbelt-lstar.f \
      Geopack-2008.f \
      
pen_igrf_src = rbelt-main1t.f \
      rbelt-penumbra.f \
      rbelt-igrf.f \
      rbelt-nogrid.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-igrf_Lm.f \
      rbelt-lstar.f \
      Geopack-2008.f \

pen_hdf_src = rbelt-main1t.f \
      rbelt-penumbra.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-igrf_Lm.f \
      rbelt-lstar.f \
      Geopack-2008.f \

gs_dip_src = rbelt-main1t.f \
      rbelt-nogrid.f \
      rbelt-dipole.f \
      rbelt-gcrsurf.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-string.f \
      rbelt-time.f \
      Geopack-2008.F \

gs_igrf_src = rbelt-main1t.f \
      rbelt-nogrid.f \
      rbelt-igrf.f \
      rbelt-gcrsurf.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      Geopack-2008.F \

gcrsurf_src = rbelt-main1t.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-gcrsurf.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      Geopack-2008.F \

gcrts05_src = rbelt-main1t.f \
      rbelt-grid.f \
      rbelt-ts05_samp.f \
      TS04c.f \
      rbelt-gcrsurf.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      Geopack-2008.F \

monticar_src = rbelt-main1t.f \
      rbelt-grid.f \
      rbelt-ts05_samp.f \
      TS04c.f \
      rbelt-monticar.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      Geopack-2008.F \

gcr1t_src = rbelt-main1t.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-gcr1t.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-string.f \
      rbelt-time.f \
      Geopack-2008.F \

gcr1ts05_src = rbelt-main1t.f \
      rbelt-grid.f \
      rbelt-ts05_samp2.f \
      TS04c.f \
      rbelt-gcr1t.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      Geopack-2008.F \

clat_src = rbelt-main1t.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-clat1t.ff \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-string.f \
      rbelt-time.f \
      Geopack-2008.F \

gcr3t_src = rbelt-main3t.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-gcr3t.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-string.f \
      rbelt-time.f \
      Geopack-2008.F \

gcr3ts05_src = rbelt-main3t.f \
      rbelt-grid.f \
      rbelt-ts05_samp2.f \
      TS04c.f \
      rbelt-gcr3t.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-string.f \
      rbelt-time.f \
      rbelt-lstar.f \
      Geopack-2008.F \

flight1t_obj = rbelt-main1t.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-flight1t.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-string.f \
      rbelt-time.f \
      Geopack-2008.f \
      rbelt-flt.f \

flight3t_obj = rbelt-main3t.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-flight3t.f \
      rbelt-geocut.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-string.f \
      rbelt-time.f \
      Geopack-2008.f \
      rbelt-flt.f \

flt_src = rbelt-flt.f \
      rbelt-grid.f \
      rbelt-hdffileread.f \
      rbelt-dist.f \
      rbelt-particle_loop.f \
      rbelt-time_loop.f \
      rbelt-rk.f \
      rbelt-lorentz.f \
      rbelt-gcenter.f \
      rbelt-linterp3d.f \
      rbelt-io.f \
      rbelt-string.f \
      second.c

pstprc_src = rbelt-pstprc.f \
      rbelt-ppgrid.f \
      rbelt-string.f \
      rbelt-time.f \
      Geopack-2008.F

setread_src = setread.f \
      rbelt-string.f \
      rbelt-time.f \
      Geopack-2008.F

##########################################################################
###########################   Object files   #############################
##########################################################################


1p_dip_obj = rbelt-main1t.o \
      rbelt-dipole.o \
      rbelt-nogrid.o \
      rbelt-dist.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-boris.o \
      rbelt-gcenter.o \
      rbelt-io.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      rbelt-coortrans.o \
      Geopack-2008.o \

1p_cdip_obj = rbelt-main1t.o \
      rbelt-cdipole.o \
      rbelt-nogrid.o \
      rbelt-dist.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-boris.o \
      rbelt-gcenter.o \
      rbelt-io.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      rbelt-coortrans.o \
      Geopack-2008.o \

1p_ts05analytic_obj = rbelt-main1t.o \
      rbelt-ts05.o \
      TS04c.o \
      rbelt-nogrid.o \
      rbelt-dist.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-io.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      Geopack-2008.o \

1p_hdf_obj = rbelt-main.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-dist.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-boris.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp4d.o \
      rbelt-io.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      rbelt-coortrans.o \
      Geopack-2008.o \

1p_hdf1_obj = rbelt-main1t.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-dist.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-io.o  \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      Geopack-2008.o \

1p_ts05_obj = rbelt-main1t.o \
      rbelt-grid.o \
      rbelt-ts05_samp.o \
      TS04c.o \
      rbelt-dist.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-io.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      Geopack-2008.o \

rbelt_dip_obj = rbelt-main1t.o \
      rbelt-dipole.o \
      rbelt-nogrid.o \
      rbelt-dist_disk.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-io.o  \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-coortrans.o \
      rbelt-lstar.o \
      Geopack-2008.o \

rbelt_cdip_obj = rbelt-main1t.o \
      rbelt-cdipole.o \
      rbelt-nogrid.o \
      rbelt-dist_disk.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-io.o  \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-coortrans.o \
      rbelt-lstar.o \
      Geopack-2008.o \

rbelt_dst_obj = rbelt-main1t.o \
      rbelt-dstfld.o \
      rbelt-nogrid.o \
      rbelt-dist_disk.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-io.o  \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-coortrans.o \
      rbelt-lstar.o \
      Geopack-2008.o \

rbelt_ulf_obj = rbelt-main1t.o \
      rbelt-ulffld.o \
      rbelt-nogrid.o \
      rbelt-dist_disk.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-io.o  \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-coortrans.o \
      rbelt-lstar.o \
      Geopack-2008.o \

rbelt_hdf_obj = rbelt-main.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-dist_disk.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-boris.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp4d.o \
      rbelt-io.o   \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-coortrans.o \
      rbelt-lstar.o \
      Geopack-2008.o \

rbelt_lfm_obj = rbelt-main.o \
      rbelt-grid.o \
      rbelt-dist_tail.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp4d.o \
      rbelt-io.o   \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      Geopack-2008.o \
      rbelt-lfm_samp.o \
      hdftake_mpi.o \
      ionset.o \
      grid.o \
      cord.o \
      therest.o \
      jcalc2.o \
      flux.o \
      interp-all.o \
      intcoef1d.o \
      cart-int.o \

rbelt_hdf1_obj = rbelt-main1t.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-dist_disk.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-io.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      Geopack-2008.o \

rbelt_ts05_obj = rbelt-main.o \
      rbelt-grid.o \
      rbelt-ts05_samp.o \
      TS04c.o \
      rbelt-dist_tail.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp4d.o \
      rbelt-io.o   \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      Geopack-2008.o \

isopt_dip_obj = rbelt-main1t.o \
      rbelt-dipole.o \
      rbelt-nogrid.o \
      rbelt-dist_isopt.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-io.o \
      rbelt-string.o \

isopt_hdf1_obj = rbelt-main1t.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-dist_isopt.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-io.o \
      rbelt-string.o \
      rbelt-time.o \

isosphr_dip_obj = rbelt-main1t.o \
      rbelt-dipole.o \
      rbelt-nogrid.o \
      rbelt-dist_sphr.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-io.o \
      rbelt-string.o \

isosphr_hdf_obj = rbelt-main.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-dist_sphr.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp4d.o \
      rbelt-io.o \
      rbelt-string.o \
      rbelt-time.o \

isosphr_hdf1_obj = rbelt-main1t.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-dist_sphr.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-io.o \
      rbelt-string.o \
      rbelt-time.o \

isosphr_ts05_obj = rbelt-main1t.o \
      rbelt-grid.o \
      rbelt-ts05_samp.o \
      TS04c.o \
      rbelt-dist_sphr.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-io.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      Geopack-2008.o \

isosphr_nob_obj = rbelt-main1t.o \
      rbelt-zero_fields.o \
      rbelt-nogrid.o \
      rbelt-dist_sphr.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-io.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      Geopack-2008.o \

isosphr_odip_obj = rbelt-main1t.o \
      rbelt-offset_dip.o \
      rbelt-nogrid.o \
      rbelt-dist_sphr.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-io.o \
      rbelt-string.o \

isosphr_t96_obj = rbelt-main1t.o \
      rbelt-t96.o \
      rbelt-nogrid.o \
      rbelt-dist_sphr.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-io.o \
      rbelt-string.o \
      T96.o \

gcr1pt_dip_obj = rbelt-main1t.o \
      rbelt-gcr1pt_ana.o \
      rbelt-dipole.o \
      rbelt-nogrid.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-string.o \
      Geopack-2008.o \
      
gcr1pt_igrf_obj = rbelt-main1t.o \
      rbelt-gcr1pt_ana.o \
      rbelt-igrf.o \
      rbelt-nogrid.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-string.o \
      Geopack-2008.o \

gcr1pt_hdf_obj = rbelt-main1t.o \
      rbelt-gcr1pt_hdf.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-string.o \
      rbelt-time.o \
      Geopack-2008.o \

acs_dip_obj = rbelt-main1t.o \
      rbelt-access.o \
      rbelt-dipole.o \
      rbelt-nogrid.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-string.o \
      Geopack-2008.o \
   
acs_igrf_obj = rbelt-main1t.o \
      rbelt-access.o \
      rbelt-igrf.o \
      rbelt-nogrid.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-string.o \
      Geopack-2008.o \

acs_hdf_obj = rbelt-main1t.o \
      rbelt-access.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-string.o \
      rbelt-time.o \
      Geopack-2008.o \

pen_dip_obj = rbelt-main1t.o \
      rbelt-penumbra.o \
      rbelt-dipole.o \
      rbelt-nogrid.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-igrf_Lm.o \
      rbelt-lstar.o \
      Geopack-2008.o \
   
pen_igrf_obj = rbelt-main1t.o \
      rbelt-penumbra.o \
      rbelt-igrf.o \
      rbelt-nogrid.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-igrf_Lm.o \
      rbelt-lstar.o \
      Geopack-2008.o \
      
pen_hdf_obj = rbelt-main1t.o \
      rbelt-penumbra.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-igrf_Lm.o \
      rbelt-lstar.o \
      Geopack-2008.o \

gs_dip_obj = rbelt-main1t.o \
      rbelt-nogrid.o \
      rbelt-dipole.o \
      rbelt-gcrsurf.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-string.o \
      rbelt-time.o \
      Geopack-2008.o \

gs_igrf_obj = rbelt-main1t.o \
      rbelt-nogrid.o \
      rbelt-igrf.o \
      rbelt-gcrsurf.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      Geopack-2008.o \

gcrsurf_obj = rbelt-main1t.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-gcrsurf.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      Geopack-2008.o \

gcrts05_obj = rbelt-main1t.o \
      rbelt-grid.o \
      rbelt-ts05_samp.o \
      TS04c.o \
      rbelt-gcrsurf.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      Geopack-2008.o \

monticar_obj = rbelt-main1t.o \
      rbelt-grid.o \
      rbelt-ts05_samp.o \
      TS04c.o \
      rbelt-monticar.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      Geopack-2008.o \

gcr1t_obj = rbelt-main1t.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-gcr1t.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-string.o \
      rbelt-time.o \
      Geopack-2008.o \

gcr1ts05_obj = rbelt-main1t.o \
      rbelt-grid.o \
      rbelt-ts05_samp2.o \
      TS04c.o \
      rbelt-gcr1t.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      Geopack-2008.o \

clat_obj = rbelt-main1t.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-clat1t.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-string.o \
      rbelt-time.o \
      Geopack-2008.o \

gcr3t_obj = rbelt-main3t.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-gcr3t.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-string.o \
      rbelt-time.o \
      Geopack-2008.o \

gcr3ts05_obj = rbelt-main3t.o \
      rbelt-grid.o \
      rbelt-ts05_samp2.o \
      TS04c.o \
      rbelt-gcr3t.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-string.o \
      rbelt-time.o \
      rbelt-lstar.o \
      Geopack-2008.o \

flight1t_obj = rbelt-main1t.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-flight1t.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-string.o \
      rbelt-time.o \
      Geopack-2008.o \
      rbelt-flt.o \

flight3t_obj = rbelt-main3t.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-flight3t.o \
      rbelt-geocut.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-string.o \
      rbelt-time.o \
      Geopack-2008.o \
      rbelt-flt.o \

flt_obj = rbelt-flt.o \
      rbelt-grid.o \
      rbelt-hdffileread.o \
      rbelt-dist.o \
      rbelt-particle_loop.o \
      rbelt-time_loop.o \
      rbelt-rk.o \
      rbelt-lorentz.o \
      rbelt-gcenter.o \
      rbelt-linterp3d.o \
      rbelt-io.o \
      rbelt-string.o \
      second.o

pstprc_obj = rbelt-pstprc.o \
      rbelt-ppgrid.o \
      rbelt-string.o \
      rbelt-time.o \
      Geopack-2008.o \


setread_obj = setread.o \
      rbelt-string.o \
      rbelt-time.o \
      Geopack-2008.o

#executables


# FOR COMPUTING SINGLE PARTICLE TRAJECTORIES

# to compute a single particle trajectory in a dipole
1p_dip: $(1p_dip_obj)
	$(COMPILER) -o 1p_dip $(1p_dip_obj) \
	$(FFLAGS)

# to compute a single particle trajectory in a compressed dipole
1p_cdip: $(1p_cdip_obj)
	$(COMPILER) -o 1p_cdip $(1p_cdip_obj) \
	$(FFLAGS)

# to compute a single particle trajectory in analytic TS05 field
1p_ts05analytic: $(1p_ts05analytic_obj)
	$(COMPILER) -o 1p_ts05analytic $(1p_ts05analytic_obj) \
	$(FFLAGS)

# to compute a single particle trajectory in time dependent fields
1p_hdf: $(1p_hdf_obj)
	$(COMPILER) -o 1p_hdf $(1p_hdf_obj) \
	$(FFLAGS)
	
# to compute a single particle trajectory in static field snapshots
1p_hdf1: $(1p_hdf1_obj)
	$(COMPILER) -o 1p_hdf1 $(1p_hdf1_obj) \
	$(FFLAGS)

# to compute a single particle trajectory in gridded TS05 field
1p_ts05: $(1p_ts05_obj)
	$(COMPILER) -o 1p_ts05 $(1p_ts05_obj) \
	$(FFLAGS)

# FOR TRAPPED RADIATION BELT DISTRIBUTION

# to launch & compute rbelt dist. in a dipole,
# launched from an annulus in the equatorial plane
rbelt_dip: $(rbelt_dip_obj)
	$(COMPILER) -o rbelt_dip $(rbelt_dip_obj) \
	$(FFLAGS)

# to launch & compute rbelt dist. in a compressed dipole,
# launched from an annulus in the equatorial plane
rbelt_cdip: $(rbelt_cdip_obj)
	$(COMPILER) -o rbelt_cdip $(rbelt_cdip_obj) \
	$(FFLAGS)

# to launch & compute rbelt dist. in a compressed dipole,
# launched from an annulus in the equatorial plane
rbelt_dst: $(rbelt_dst_obj)
	$(COMPILER) -o rbelt_dst $(rbelt_dst_obj) \
	$(FFLAGS)

# to launch & compute rbelt dist. in a ULF model,
# launched from an annulus in the equatorial plane
rbelt_ulf: $(rbelt_ulf_obj)
	$(COMPILER) -o rbelt_ulf $(rbelt_ulf_obj) \
	$(FFLAGS)


# to launch & compute rbelt dist. in time dependent fields,
# launched from an annulus in the equatorial plane
rbelt_hdf: $(rbelt_hdf_obj)
	$(COMPILER) -o rbelt_hdf $(rbelt_hdf_obj) \
	$(FFLAGS)

# to launch & compute rbelt dist. in time dependent fields 
# read directly from LFM output files onto rbelt grid,
# launched from an annulus in the equatorial plane
rbelt_lfm: $(rbelt_lfm_obj)
	$(COMPILER) -o rbelt_lfm $(rbelt_lfm_obj) \
	$(FFLAGS)

# to launch & compute rbelt dist. in static field snapshots,
# launched from an annulus in the equatorial plane
rbelt_hdf1: $(rbelt_hdf1_obj)
	$(COMPILER) -o rbelt_hdf1 $(rbelt_hdf1_obj) \
	$(FFLAGS)

# to launch & compute rbelt dist. in time-dep. fields,
# launched from magnetotail region
rbelt_ts05: $(rbelt_ts05_obj)
	$(COMPILER) -o rbelt_ts05 $(rbelt_ts05_obj) \
	$(FFLAGS)


# FOR DISTRIBUTION ISOTROPICALLY INCIDENT FROM OUTSIDE MAIN DIPOLE FIELD

# to launch & compute isotropic dist. into a dipole,
# launched into spherical shell near outer boundary
isosphr_dip:$(isosphr_dip_obj)
	$(COMPILER) -o isosphr_dip $(isosphr_dip_obj) \
	$(FFLAGS)
	
# to launch & compute isotropic dist. into time dependent fields,
# launched into spherical shell near outer boundary
isosphr_hdf:$(isosphr_hdf_obj)
	$(COMPILER) -o isosphr_ $(isosphr_hdf_obj) \
	$(FFLAGS)

# to launch & compute isotropic dist. into static field snapshots,
# launched into spherical shell near outer boundary
isosphr_hdf1:$(isosphr_hdf1_obj)
	$(COMPILER) -o isosphr_hdf1 $(isosphr_hdf1_obj) \
	$(FFLAGS)

# to launch & compute isotropic dist. into static field snapshots,
# launched into spherical shell near outer boundary
# fields put directly on rbelt grid without HDF field files
isosphr_ts05:$(isosphr_ts05_obj)
	$(COMPILER) -o isosphr_ts05 $(isosphr_ts05_obj) \
	$(FFLAGS)

# to launch & compute isotropic dist. into a zero fields,
# launched into spherical shell near outer boundary
isosphr_nob:$(isosphr_nob_obj)
	$(COMPILER) -o isosphr_nob $(isosphr_nob_obj) \
	$(FFLAGS)

# to launch & compute isotropic dist. into offset dipole field,
# launched into spherical shell near outer boundary
isosphr_odip:$(isosphr_odip_obj)
	$(COMPILER) -o isosphr_odip $(isosphr_odip_obj) \
	$(FFLAGS)

# to launch & compute isotropic dist. into t96 fields,
# launched into spherical shell near outer boundary
isosphr_t96:$(isosphr_t96_obj)
	$(COMPILER) -o isosphr_t96 $(isosphr_t96_obj) \
	$(FFLAGS)  


# FOR COMPUTING PARTICLE TRAJECTORIES LAUNCHED ISOTROPOCALLY FROM A POINT 

# to launch and compute tajectories isotorpically from a single point 
# in a dipole
isopt_dip:$(isopt_dip_obj)
	$(COMPILER) -o isopt_dip $(isopt_dip_obj) \
	$(FFLAGS)

# to launch and compute tajectories isotorpically from a single point 
# in static field snapshots
isopt_hdf1:$(isopt_hdf1_obj)
	$(COMPILER) -o isopt_hdf1 $(isopt_hdf1_obj) \
	$(FFLAGS)


# GEOMAGNITIC CUTOFF RIGIDITY SEARCH ROUTINES

# to find cutoff rigidity at a single point in an analytic dipole field
gcr1pt_dip:$(gcr1pt_dip_obj)
	$(COMPILER) -o gcr1pt_dip $(gcr1pt_dip_obj) \
	$(FFLAGS)
	
# to find cutoff rigidity at a single point in analytic IGRF field
gcr1pt_igrf:$(gcr1pt_igrf_obj)
	$(COMPILER) -o gcr1pt_igrf $(gcr1pt_igrf_obj) \
	$(FFLAGS)

# to find  cutoff rigidities at a single point in static field snapshot(s)
gcr1pt_hdf:$(gcr1pt_hdf_obj)
	$(COMPILER) -o gcr1pt_hdf $(gcr1pt_hdf_obj) \
	$(FFLAGS)


# penumbra sample routines
# launch single particle, see if it escapes to outer boundary
acs_dip:$(acs_dip_obj)
	$(COMPILER) -g -o acs_dip $(acs_dip_obj) \
	$(FFLAGS)
	
# launch single particle, see if it escapes to outer boundary
acs_igrf:$(acs_igrf_obj)
	$(COMPILER) -o acs_igrf $(acs_igrf_obj) \
	$(FFLAGS)

# launch single particle, see if it escapes to outer boundary
acs_hdf:$(acs_hdf_obj)
	$(COMPILER) -o acs_hdf $(acs_hdf_obj) \
	$(FFLAGS)

# launch particles with increasing rigidity to map penumbra region
pen_dip:$(pen_dip_obj)
	$(COMPILER) -g -o pen_dip $(pen_dip_obj) \
	$(FFLAGS)
	
# launch particles with increasing rigidity to map penumbra region
pen_igrf:$(pen_igrf_obj)
	$(COMPILER) -o pen_igrf $(pen_igrf_obj) \
	$(FFLAGS)

# launch particles with increasing rigidity to map penumbra region
pen_hdf:$(pen_hdf_obj)
	$(COMPILER) -o pen_hdf $(pen_hdf_obj) \
	$(FFLAGS)


# additional routines

# to find cutoff rigidities at points on a geographic grid
gs_dip: $(gs_dip_obj)
	$(COMPILER) -o gs_dip $(gs_dip_obj) \
	$(FFLAGS)

# to find cutoff rigidities at points on a geographic grid
gs_igrf: $(gs_igrf_obj)
	$(COMPILER) -o gs_igrf $(gs_igrf_obj) \
	$(FFLAGS)

# to find cutoff rigidities at points on a geographic grid
gcrsurf: $(gcrsurf_obj)
	$(COMPILER) -o gcrsurf $(gcrsurf_obj) \
	$(FFLAGS)

# to find cutoff rigidities at points on a geographic grid
# fields put directly on rbelt grid without HDF field files
gcrts05: $(gcrts05_obj)
	$(COMPILER) -o gcrts05 $(gcrts05_obj) \
	$(FFLAGS)

# to find cutoff rigidities at points on a geographic grid
# fields put directly on rbelt grid without HDF field files
monticar: $(monticar_obj)
	$(COMPILER) -o monticar $(monticar_obj) \
	$(FFLAGS)

# to find cutoff rigidities along a spacecraft trajectory 
# in static field snapshots
gcr1t: $(gcr1t_obj)
	$(COMPILER) -o gcr1t $(gcr1t_obj) \
	$(FFLAGS)

# to find cutoff rigidities along a spacecraft trajectory 
# in static field snapshots read directly to grid from TS05
gcr1ts05: $(gcr1ts05_obj)
	$(COMPILER) -o gcr1ts05 $(gcr1ts05_obj) \
	$(FFLAGS)

# searches for cutoff latitudes along a spacecraft trajectory
clat: $(clat_obj)
	$(COMPILER) -o clat $(clat_obj) \
	$(FFLAGS)

# to find cutoff rigidities along a spacecraft trajectory 
# in field snapshot in rbelt time step grid position 3 
# interpolated from positions 1 & 2.
gcr3t: $(gcr3t_obj)
	$(COMPILER) -o gcr3t $(gcr3t_obj) \
	$(FFLAGS)
	
# to find cutoff rigidities along a spacecraft trajectory 
# in static field snapshots in rbelt time step grid position 3 
# interpolated from positions 1 & 2, and read directly to grid from TS05
gcr3ts05: $(gcr3ts05_obj)
	$(COMPILER) -o gcr3ts05 $(gcr3ts05_obj) \
	$(FFLAGS)

# NAIRAS PROJ. GEOMAGNITIC CUTOFF RIGIDITY SEARCH ROUTINES
	
# to find cutoff rigidities along flight path, specified in NAIRAS 
# project input file, in static field snapshots (e.g. IGRF)
flight1t: $(flight1t_obj)
	$(COMPILER) -o flight1t $(flight1t_obj) \
	$(FFLAGS)
		
# to find cutoff rigidities along flight path, specified in NAIRAS 
# project input file, in rbelt time step grid position 3 
# interpolated from positions 1 & 2.
flight3t: $(flight3t_obj)
	$(COMPILER) -o flight3t $(flight3t_obj) \
	$(FFLAGS)
	
	
# MISC.

# radiation belt post proscessing code
pstprc: $(pstprc_obj)
	$(COMPILER) -o pstprc $(pstprc_obj) \
	$(FFLAGS)

# field line tracing algorythm
flt: $(flt_obj)
	$(COMPILER) -o flt $(flt_obj) \
	$(FFLAGS)
	
setread: $(setread_obj)
	$(COMPILER) -o setread $(setread_obj) \
	$(FFLAGS)


clean:
	- rm *.o 

#-DEBUG:subscript_check=ON:div_check=3:verbose_runtime=ON:fullwarn=ON:conform_check=ON
#setenv TRAP_FPE "ALL=COUNT; UNDERFL=ZERO; OVERFL=IEEE,TRACE(2),ABORT(100); DIVZERO=TRACE(2),ABORT(3)"

