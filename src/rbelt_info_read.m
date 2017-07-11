function f = rbelt_info_read(basename,fnumber)

filename=[basename,'-info-',num2str(fnumber),'.txt']

fp = fopen(filename,'r');

%
%     from rbelt-grid.inc
[s1,c] = fscanf(fp, '%s', [1,1]);
[nx,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[ny,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[nz,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[nt,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[nstep,c] = fscanf(fp, '%u\n', [1,1]);
%
%     from rbelt-const.inc
[s1,c] = fscanf(fp, '%s', [1,1]);
[charge,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[charge_sign,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[m0,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[b0,c] = fscanf(fp, '%f\n', [1,1]);
%
%     from rbelt-y0.inc
[s1,c] = fscanf(fp, '%s', [1,1]);
[num_particles,c] = fscanf(fp, '%u\n', [1,1]);
%
%     from rbelt-io.inc
[s1,c] = fscanf(fp, '%s', [1,1]);
[max_yout,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[num_flts,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[num_ints,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[max_wsteps,c] = fscanf(fp, '%u\n', [1,1]);
%
%     rbelt namelist
[s1,c] = fscanf(fp, '%s', [1,1]);
[basename,c] = fscanf(fp, '%s\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[firstfilenum,c] = fscanf(fp, '%u\n', [1,1]);
%
%     bounds namelist
[s1,c] = fscanf(fp, '%s', [1,1]);
[rmin,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[rmax,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[tmax,c] = fscanf(fp, '%f\n', [1,1]);
%
%     io namelist
[s1,c] = fscanf(fp, '%s', [1,1]);
[dthalt,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[flux_dt,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[init_twrite,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[dtwrite,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[binio,c] = fscanf(fp, '%s\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[flux_out,c] = fscanf(fp, '%s\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[dist_out,c] = fscanf(fp, '%s\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[prcp_out,c] = fscanf(fp, '%s\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[init_out,c] = fscanf(fp, '%s\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[cone_out,c] = fscanf(fp, '%s\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[lcalc,c] = fscanf(fp, '%u\n', [1,1]);
%
%     Lorentz namelist
[s1,c] = fscanf(fp, '%s', [1,1]);
[tstep_lrntz,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[dx_max,c] = fscanf(fp, '%f\n', [1,1]);
%
%     guiding center namelist
[s1,c] = fscanf(fp, '%s', [1,1]);
[tstep_gc,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[dx_max_gc,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[go2lrntz,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[go2gc,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[dtgo2gc,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[seed,c] = fscanf(fp, '%u\n', [1,1]);
%
%     dist namelist
[s1,c] = fscanf(fp, '%s', [1,1]);
[dist_seed,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[dt_dist,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[init_t,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[emin,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[emax,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[epa_min,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[epa_max,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[lmin,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[lmax,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[radius,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[exp,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[factor,c] = fscanf(fp, '%f\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[flag0,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[flag_switch,c] = fscanf(fp, '%s\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[initdist,c] = fscanf(fp, '%u\n', [1,1]);
%
%
[s1,c] = fscanf(fp, '%s', [1,1]);
[year0,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[doy0,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[hour0,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[min0,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[sec0,c] = fscanf(fp, '%u\n', [1,1]);
%
% output file data
[s1,c] = fscanf(fp, '%s', [1,1]);
[num_wsteps,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[wlines_dist,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[wlines_flux,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[wlines_prcp,c] = fscanf(fp, '%u\n', [1,1]);
[s1,c] = fscanf(fp, '%s', [1,1]);
[wlines_cone,c] = fscanf(fp, '%u\n', [1,1]);

f=wlines_dist

fclose(fp);
