#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 12:59:24 2020

@author: ssli

Predict theoretical cosmic shear signal with input cosmological parameters

Modified from Likelihood for the KiDS+VIKING-450 correlation functions 
    https://github.com/fkoehlin/kv450_cf_likelihood_public

"""

from __future__ import print_function

from scipy import interpolate as itp
from scipy import special
import os
import numpy as np
import math

import io_cs


def CSsignalFunc(data, cosmo):
    """
    Modeling the cosmic shear signal from theory 
    """

    # Number of bins
    nzbins = len(data.const['z_bins_min'])
    # Number of correlation
    nzcorrs = int(nzbins * (nzbins + 1) / 2)


    #####################################################################
    # read in data vector
    #####################################################################
    # the theta list is the one actually used
    path_tmp = os.path.join(data.paths['data'], 'data_vector/for_cosmo/' + data.conf['prefix_data_vector'] + data.conf['sample'] + '.dat')
    theta, _, _ = np.loadtxt(path_tmp, unpack=True)
    print('Loaded data vectors from: \n', path_tmp)
    # we assume theta is same for xi_p and xi_m
    theta_bins = np.concatenate((theta, theta))

    # Read angular cut values
    if data.conf['use_cut_theta']:
        mask, mask_indices, mask_suffix = io_cs.ReadCutValueFunc(data, nzbins, nzcorrs, theta_bins)
    else:
        mask = np.ones(2 * nzcorrs * data.const['ntheta'])
        mask_indices = np.where(mask == 1)[0]
        mask_suffix = data.conf['cutvalues_file'][:-4]



    #####################################################################
    # read redshift distribution
    #####################################################################
    z_samples, hist_samples, nz_samples = io_cs.ReadZdistribution(data, nzbins)

    # prevent undersampling of histograms!
    if data.const['nzmax'] < nz_samples:
        print("You are trying to integrate at lower resolution than supplied by the n(z) histograms. \n Increase nzmax! Aborting run now... \n")
        exit()
    # if that's the case, we want to integrate at histogram resolution and need to account for
    # the extra zero entry added
    elif data.const['nzmax'] == nz_samples:
        nzmax = z_samples.shape[1]
        # requires that z-spacing is always the same for all bins...
        z_p = z_samples[0, :]
        print('Redshift integrations performed at resolution of redshift distribution histograms! \n')
    # if we interpolate anyway at arbitrary resolution the extra 0 doesn't matter
    else:
        nzmax = data.const['nzmax'] + 1
        z_p = np.linspace(z_samples.min(), z_samples.max(), nzmax)
        print('Redshift integrations performed at set "nzmax" resolution! \n')

    pz = np.zeros((nzmax, nzbins))
    pz_norm = np.zeros(nzbins)
    for zbin in range(nzbins):
        # we assume that the histograms loaded are given as left-border histograms
        # and that the z-spacing is the same for each histogram
        spline_pz = itp.interp1d(z_samples[zbin, :], hist_samples[zbin, :], kind=data.conf['type_redshift_interp'])
        mask_min = z_p >= z_samples[zbin, :].min()
        mask_max = z_p <= z_samples[zbin, :].max()
        mask_z = mask_min & mask_max
        # points outside the z-range of the histograms are set to 0!
        pz[mask_z, zbin] = spline_pz(z_p[mask_z])
        # Normalize selection functions
        dz = z_p[1:] - z_p[:-1]
        pz_norm[zbin] = np.sum(0.5 * (pz[1:, zbin] + pz[:-1, zbin]) * dz)

    if ('D_z' in data.nuisance_parameters):

        # Create labels for loading of dn/dz-files:
        zbin_labels = []
        for i in range(nzbins):
            zbin_labels += ['{:.1f}t{:.1f}'.format(data.const['z_bins_min'][i], data.const['z_bins_max'][i])]

        pz = np.zeros((nzmax, nzbins))
        pz_norm = np.zeros(nzbins)
        for zbin in range(nzbins):

            param_name = 'D_z{:}'.format(zbin + 1)
            z_mod = z_p + data.nuisance_parameters['D_z'][param_name]

            # Load n(z) again:
            fname = os.path.join(
            data.paths['data'], 'redshift/' + data.conf['sample'] + '/Nz_{0:}/Nz_{0:}_Mean/Nz_{0:}_z{1:}.asc'.format(data.conf['nz_method'], zbin_labels[zbin]))
            zptemp, hist_pz = np.loadtxt(fname, usecols=(0, 1), unpack=True)
            shift_to_midpoint = np.diff(zptemp)[0] / 2.
            
            spline_pz = itp.interp1d(zptemp + shift_to_midpoint, hist_pz, kind=data.conf['type_redshift_interp'])
            
            mask_min = z_mod >= zptemp.min() + shift_to_midpoint
            mask_max = z_mod <= zptemp.max() + shift_to_midpoint
            mask_tmp = mask_min & mask_max
            # points outside the z-range of the histograms are set to 0!
            pz[mask_tmp, zbin] = spline_pz(z_mod[mask_tmp])
            # Normalize selection functions
            dz = z_p[1:] - z_p[:-1]
            pz_norm[zbin] = np.sum(0.5 * (pz[1:, zbin] + pz[:-1, zbin]) * dz)


    #####################################################################
    # CLASS cosmo initialization
    #####################################################################
    # the last cosmo arguments  
    data.cosmo_arguments['z_max_pk'] = z_p.max()
    print('Intitial cosmological parameters passed to CLASS code:')
    print(data.cosmo_arguments)

    # Prepare the cosmological module with the input parameters
    cosmo.set(data.cosmo_arguments)
    cosmo.compute(["lensing"])
    print('sigma8 =', cosmo.sigma8())

    # Omega_m contains all species!
    Omega_m = cosmo.Omega_m()
    print('Omega_m =', Omega_m)
    small_h = cosmo.h()
    print('h =', small_h)
    # One wants to obtain here the relation between z and r, this is done
    # by asking the cosmological module with the function z_of_r
    r, dzdr = cosmo.z_of_r(z_p)


    ################################################
    # discrete theta values (to convert C_l to xi's)
    ################################################
    if data.conf['use_theory_binning']:
        thetamin = np.min(data.const['theta_bin_min_val']) * 0.8
        thetamax = np.max(data.const['theta_bin_max_val']) * 1.2
    else:
        thetamin = np.min(theta_bins) * 0.8
        thetamax = np.max(theta_bins) * 1.2

    nthetatot = np.ceil(math.log(thetamax / thetamin) / data.const['dlntheta']) + 1
    nthetatot = np.int32(nthetatot)
    theta = np.zeros(nthetatot)
    a2r = math.pi / (180. * 60.)

    # define an array of thetas
    for it in range(nthetatot):
        theta[it] = thetamin * math.exp(data.const['dlntheta'] * it)



    ################################################################
    # discrete l values used in the integral to convert C_l to xi's)
    ################################################################
    # wavenumber l for Cl-integration
    # It is a logspace
    # find nlmax in order to reach lmax with logarithmic steps dlnl
    nlmax = np.int(np.log(data.const['lmax']) / data.const['dlnl']) + 1
    # redefine slightly dlnl so that the last point is always exactly lmax
    dlnl = np.log(data.const['lmax']) / (nlmax - 1)
    l = np.exp(dlnl * np.arange(nlmax))

    if data.conf['integrate_Bessel_with'] in ['brute_force', 'cut_off']:
        # l = x / theta / a2r
        # x = l * theta * a2r

        # We start by considering the largest theta, theta[-1], and for that value we infer
        # a list of l's from the requirement that corresponding x values are spaced linearly with a given stepsize, until xmax.
        # Then we loop over smaller theta values, in decreasing order, and for each of them we complete the previous list of l's,
        # always requiuring the same dx stepsize (so that dl does vary) up to xmax.
        #
        # We first apply this to a running value ll, in order to count the total numbner of ll's, called nl.
        # Then we create the array lll[nl] and we fill it with the same values.
        #
        # we also compute on the fly the critical index il_max[it] such that ll[il_max[it]]*theta[it]*a2r
        # is the first value of x above xmax

        ll=1.
        il=0
        while (ll*theta[-1]*a2r < data.const['dx_threshold']):
            ll += data.const['dx_below_threshold']/theta[-1]/a2r
            il += 1
        for it  in range(nthetatot):
            while (ll*theta[nthetatot-1-it]*a2r < data.const['xmax']) and (ll+data.const['dx_above_threshold']/theta[nthetatot-1-it]/a2r < data.const['lmax']):
                ll += data.const['dx_above_threshold']/theta[nthetatot-1-it]/a2r
                il += 1
        nl = il+1

        lll = np.zeros(nl)
        il_max = np.zeros(nthetatot)
        il=0
        lll[il]=1.
        while (lll[il]*theta[-1]*a2r < data.const['dx_threshold']):
            il += 1
            lll[il] = lll[il-1] + data.const['dx_below_threshold']/theta[-1]/a2r
        for it in range(nthetatot):
            while (lll[il]*theta[nthetatot-1-it]*a2r < data.const['xmax']) and (lll[il] + data.const['dx_above_threshold']/theta[nthetatot-1-it]/a2r < data.const['lmax']):
                il += 1
                lll[il] = lll[il-1] + data.const['dx_above_threshold']/theta[nthetatot-1-it]/a2r
            il_max[nthetatot-1-it] = il

        # finally we compute the array l*dl that will be used in the trapezoidal integration
        # (l is a factor in the integrand [l * C_l * Bessel], and dl is like a weight)
        ldl = np.zeros(nl)
        ldl[0]=lll[0]*0.5*(lll[1]-lll[0])
        for il in range(1,nl-1):
            ldl[il]=lll[il]*0.5*(lll[il+1]-lll[il-1])
        ldl[-1]=lll[-1]*0.5*(lll[-1]-lll[-2])
    else:
        try:
            import pycl2xi.fftlog as fftlog
        except:
            print('FFTLog was requested as integration method for the Bessel functions but is not installed. \n Download it from "https://github.com/tilmantroester/pycl2xi" and follow the installation instructions there (also requires the fftw3 library). \n Aborting run now... \n')
            exit()

        # this has to be declared a self, otherwise fftlog won't be available
        Cl2xi = fftlog.Cl2xi

        # this is sufficient (FFTLog only uses 5k points internally anyways...)
        ell_lin = np.arange(1., 501., 1)
        ell_log = np.logspace(np.log10(501.), np.log10(data.const['lmax']), 5000 - len(ell_lin))
        lll = np.concatenate((ell_lin, ell_log))
        # linspace --> overkill and too slow!
        #lll = np.arange(1., lmax + 1., 1)
        nl = lll.size


    #####################################################################
    # arrays and some integrations necessary for the theory binning:
    #####################################################################
    if data.conf['use_theory_binning']:
        if data.conf['read_weight_func_for_binning']:
            fname = os.path.join(data.paths['data'], data.conf['theory_weight_func_file'])
            thetas, weights = np.loadtxt(fname, unpack=True)
            theory_weight_func = itp.splrep(thetas, weights)
        else:
            thetas = np.linspace(data.const['theta_bin_min_val'], data.const['theta_bin_max_val'], data.const['ntheta'] * int(data.const['theta_nodes_theory']))
            weights = a2r * thetas * data.const['theory_binning_const']
            theory_weight_func = itp.splrep(thetas, weights)

        # first get the theta-bin borders based on ntheta and absolute min and absolute max values
        a = np.linspace(np.log10(data.const['theta_bin_min_val']), np.log10(data.const['theta_bin_max_val']), data.const['ntheta'] + 1)
        theta_bins_tmp = 10.**a
        theta_bin_min = theta_bins_tmp[:-1]
        theta_bin_max = theta_bins_tmp[1:]

        int_weight_func = np.zeros(data.const['ntheta'])
        thetas_for_theory_binning = np.zeros((data.const['ntheta'], int(data.const['theta_nodes_theory'])))
        for idx_theta in range(data.const['ntheta']):
            theta_tmp = np.linspace(theta_bin_min[idx_theta], theta_bin_max[idx_theta], int(data.const['theta_nodes_theory']))
            dtheta = (theta_tmp[1:] - theta_tmp[:-1]) * a2r

            weight_func_integrand = itp.splev(theta_tmp, theory_weight_func)

            int_weight_func[idx_theta] = np.sum(0.5 * (weight_func_integrand[1:] + weight_func_integrand[:-1]) * dtheta)
            # for convenience:
            thetas_for_theory_binning[idx_theta, :] = theta_tmp


    #####################################################################
    # Allocation of various arrays filled and used later
    #####################################################################
    g = np.zeros((nzmax, nzbins))
    pk = np.zeros((nlmax, nzmax))
    pk_lin = np.zeros((nlmax, nzmax))
    Cl_integrand = np.zeros((nzmax, nzcorrs))
    Cl = np.zeros((nlmax, nzcorrs))
    spline_Cl = np.empty(nzcorrs, dtype=(list, 3))
    Cll = np.zeros((nzcorrs, nl))
    BBessel0 = np.zeros(nl)
    BBessel4 = np.zeros(nl)
    xi1 = np.zeros((nthetatot, nzcorrs))
    xi2 = np.zeros((nthetatot, nzcorrs))
    xi1_theta = np.empty(nzcorrs, dtype=(list, 3))
    xi2_theta = np.empty(nzcorrs, dtype=(list, 3))


    ################################################
    # intrinsic alignment
    ################################################
    # needed for IA modelling:
    if ('A_IA' in data.nuisance_parameters) and ('exp_IA' in data.nuisance_parameters):
        amp_IA = data.nuisance_parameters['A_IA']
        exp_IA = data.nuisance_parameters['exp_IA']
        intrinsic_alignment = True
    elif ('A_IA' in data.nuisance_parameters) and ('exp_IA' not in data.nuisance_parameters):
        amp_IA = data.nuisance_parameters['A_IA']
        # redshift-scaling is turned off:
        exp_IA = 0.

        intrinsic_alignment = True
    else:
        intrinsic_alignment = False

    # get linear growth rate if IA are modelled:
    if intrinsic_alignment:
        rho_crit = RhoCriticalFunc(small_h)
        # derive the linear growth factor D(z)
        linear_growth_rate = np.zeros_like(z_p)
        #print(redshifts)
        for index_z, z in enumerate(z_p):
            try:
                # for CLASS ver >= 2.6:
                linear_growth_rate[index_z] = cosmo.scale_independent_growth_factor(z)
            except:
                # my own function from private CLASS modification:
                linear_growth_rate[index_z] = cosmo.growth_factor_at_z(z)
        # normalize to unity at z=0:
        try:
            # for CLASS ver >= 2.6:
            linear_growth_rate /= cosmo.scale_independent_growth_factor(0.)
        except:
            # my own function from private CLASS modification:
            linear_growth_rate /= cosmo.growth_factor_at_z(0.)


    ################################################
    # the lens efficiency
    ################################################
    # that depends on r and the bin
    # g_i(r) = 2r(1+z(r)) int_r^+\infty drs p_r(rs) (rs-r)/rs

    # Compute now the selection function p(r) = p(z) dz/dr normalized
    # to one. The np.newaxis helps to broadcast the one-dimensional array
    # dzdr to the proper shape. Note that p_norm is also broadcasted as
    # an array of the same shape as p_z
    pr = pz * (dzdr[:, np.newaxis] / pz_norm)

    for Bin in range(nzbins):
        # shift from first entry only useful if first enrty is 0!
        for nr in range(1, nzmax-1):
            fun = pr[nr:, Bin] * (r[nr:] - r[nr]) / r[nr:]
            g[nr, Bin] = np.sum(0.5 * (fun[1:] + fun[:-1]) * (r[nr + 1:] - r[nr:-1]))
            g[nr, Bin] *= 2. * r[nr] * (1. + z_p[nr])



    ################################################
    # matter power spectrum
    ################################################    
    # P(k=l/r,z(r)) from cosmological module
    kmax_in_inv_Mpc = data.const['k_max_h_by_Mpc'] * small_h
    for index_l in range(nlmax):
        for index_z in range(1, nzmax):

            k_in_inv_Mpc = (l[index_l] + 0.5) / r[index_z]
            if (k_in_inv_Mpc > kmax_in_inv_Mpc):
                pk_dm = 0.
                pk_lin_dm = 0.
            else:
                pk_dm = cosmo.pk(k_in_inv_Mpc, z_p[index_z])
                pk_lin_dm = cosmo.pk_lin(k_in_inv_Mpc, z_p[index_z])

            pk[index_l, index_z] = pk_dm
            pk_lin[index_l, index_z] = pk_lin_dm


    ################################################
    # convergence power spectrum
    ################################################    
    Cl_GG_integrand = np.zeros_like(Cl_integrand)
    Cl_GG = np.zeros_like(Cl)

    if intrinsic_alignment:
        Cl_II_integrand = np.zeros_like(Cl_integrand)
        Cl_II = np.zeros_like(Cl)

        Cl_GI_integrand = np.zeros_like(Cl_integrand)
        Cl_GI = np.zeros_like(Cl)

    dr = r[1:] - r[:-1]
    # Start loop over l for computation of C_l^shear
    # Start loop over l for computation of E_l
    for il in range(nlmax):
        # find Cl_integrand = (g(r) / r)**2 * P(l/r,z(r))
        for Bin1 in range(nzbins):
            for Bin2 in range(Bin1, nzbins):
                Cl_GG_integrand[1:, one_dim_index(Bin1,Bin2, nzbins)] = g[1:, Bin1] * g[1:, Bin2] / r[1:]**2 * pk[il, 1:]
                if intrinsic_alignment:
                    factor_IA = IAFactorFunc(small_h, Omega_m, rho_crit, \
                                                z_p, linear_growth_rate, amp_IA, exp_IA)
                        
                    if data.conf['use_linear_pk_for_IA']:
                        # this term (II) uses the linear matter power spectrum P_lin(k, z)
                        Cl_II_integrand[1:, one_dim_index(Bin1,Bin2, nzbins)] = pr[1:, Bin1] * pr[1:, Bin2] * factor_IA[1:]**2 / r[1:]**2 * pk_lin[il, 1:]
                        # this term (GI) uses sqrt(P_lin(k, z) * P_nl(k, z))
                        Cl_GI_integrand[1:, one_dim_index(Bin1,Bin2, nzbins)] = (g[1:, Bin1] * pr[1:, Bin2] + g[1:, Bin2] * pr[1:, Bin1]) * factor_IA[1:] / r[1:]**2 * np.sqrt(pk_lin[il, 1:] * pk[il, 1:])
                    else:
                        # both II and GI terms use the non-linear matter power spectrum P_nl(k, z)
                        Cl_II_integrand[1:, one_dim_index(Bin1,Bin2, nzbins)] = pr[1:, Bin1] * pr[1:, Bin2] * factor_IA[1:]**2 / r[1:]**2 * pk[il, 1:]
                        Cl_GI_integrand[1:, one_dim_index(Bin1,Bin2, nzbins)] = (g[1:, Bin1] * pr[1:, Bin2] + g[1:, Bin2] * pr[1:, Bin1]) * factor_IA[1:] / r[1:]**2 * pk[il, 1:]

        # Integrate over r to get C_l^shear_ij = P_ij(l)
        # C_l^shear_ij = 9/16 Omega0_m^2 H_0^4 \sum_0^rmax dr (g_i(r)
        # g_j(r) /r**2) P(k=l/r,z(r)) dr
        # It is then multiplied by 9/16*Omega_m**2
        # and then by (h/2997.9)**4 to be dimensionless
        # (since P(k)*dr is in units of Mpc**4)
        for Bin in range(nzcorrs):
            Cl_GG[il, Bin] = np.sum(0.5 * (Cl_GG_integrand[1:, Bin] + Cl_GG_integrand[:-1, Bin]) * dr)
            Cl_GG[il, Bin] *= 9. / 16. * Omega_m**2
            Cl_GG[il, Bin] *= (small_h / 2997.9)**4

            if intrinsic_alignment:
                Cl_II[il, Bin] = np.sum(0.5 * (Cl_II_integrand[1:, Bin] + Cl_II_integrand[:-1, Bin]) * dr)

                Cl_GI[il, Bin] = np.sum(0.5 * (Cl_GI_integrand[1:, Bin] + Cl_GI_integrand[:-1, Bin]) * dr)
                # here we divide by 4, because we get a 2 from g(r)!
                Cl_GI[il, Bin] *= 3. / 4. * Omega_m
                Cl_GI[il, Bin] *= (small_h / 2997.9)**2

    if intrinsic_alignment:
        Cl = Cl_GG + Cl_GI + Cl_II
    else:
        Cl = Cl_GG

    # Spline Cl[il,Bin1,Bin2] along l
    for Bin in range(nzcorrs):
        spline_Cl[Bin] = list(itp.splrep(l, Cl[:, Bin]))

    # Interpolate Cl at values lll and store results in Cll
    for Bin in range(nzcorrs):
        Cll[Bin,:] = itp.splev(lll[:], spline_Cl[Bin])



    ################################################
    # shear correlation function
    ################################################    
    if data.conf['integrate_Bessel_with'] == 'brute_force':
        # this seems to produce closest match in comparison with CCL
        # I still don't like the approach of just integrating the Bessel
        # functions over some predefined multipole range...
        #t0 = timer()
        # Start loop over theta values
        for it in range(nthetatot):
            BBessel0[:] = special.j0(lll[:] * theta[it] * a2r)
            BBessel4[:] = special.jv(4, lll[:] * theta[it] * a2r)

            # Here is the actual trapezoidal integral giving the xi's:
            # - in more explicit style:
            # for Bin in range(nzbin_pairs):
            #     for il in range(ilmax):
            #         xi1[it, Bin] = np.sum(ldl[il]*Cll[Bin,il]*BBessel0[il])
            #         xi2[it, Bin] = np.sum(ldl[il]*Cll[Bin,il]*BBessel4[il])
            # - in more compact and vectorizable style:
            xi1[it, :] = np.sum(ldl[:] * Cll[:, :] * BBessel0[:], axis=1)
            xi2[it, :] = np.sum(ldl[:] * Cll[:, :] * BBessel4[:], axis=1)
        #dt = timer() - t0
        #print('dt = {:.6f}'.format(dt))
        #print(lll.min(), lll.max(), lll.shape)
        #exit()

        # normalize xis
        xi1 = xi1 / (2. * math.pi)
        xi2 = xi2 / (2. * math.pi)

    elif data.conf['integrate_Bessel_with'] == 'fftlog':
        #t0 = timer()
        for zcorr in range(nzcorrs):

            # convert theta from arcmin to deg; xis are already normalized!
            xi1[:, zcorr] = Cl2xi(Cll[zcorr, :], lll[:], theta[:] / 60., bessel_order=0) #, ell_min_fftlog=lll.min(), ell_max_fftlog=lll.max() + 1e4)
            xi2[:, zcorr] = Cl2xi(Cll[zcorr, :], lll[:], theta[:] / 60., bessel_order=4) #, ell_min_fftlog=lll.min(), ell_max_fftlog=lll.max() + 1e4)
        #dt = timer() - t0
        #print('dt = {:.6f}'.format(dt))
        #print(lll.min(), lll.max(), lll.shape)
        #exit()

    else:
        #t0 = timer()
        for it in range(nthetatot):
            ilmax = il_max[it]

            BBessel0[:ilmax] = special.j0(lll[:ilmax] * theta[it] * a2r)
            BBessel4[:ilmax] = special.jv(4, lll[:ilmax] * theta[it] * a2r)

            # Here is the actual trapezoidal integral giving the xi's:
            # - in more explicit style:
            # for Bin in range(nzcorrs):
            #     for il in range(ilmax):
            #         xi1[it, Bin] = np.sum(ldl[il]*Cll[Bin,il]*BBessel0[il])
            #         xi2[it, Bin] = np.sum(ldl[il]*Cll[Bin,il]*BBessel4[il])
            # - in more compact and vectorizable style:
            xi1[it, :] = np.sum(ldl[:ilmax] * Cll[:, :ilmax] * BBessel0[:ilmax], axis=1)
            xi2[it, :] = np.sum(ldl[:ilmax] * Cll[:, :ilmax] * BBessel4[:ilmax], axis=1)
        #dt = timer() - t0
        #print('dt = {:.6f}'.format(dt))
        #print(lll.min(), lll.max(), lll.shape)
        #exit()
    
        # normalize xis
        xi1 = xi1 / (2. * math.pi)
        xi2 = xi2 / (2. * math.pi)

    # Spline the xi's
    for Bin in range(nzcorrs):
        xi1_theta[Bin] = list(itp.splrep(theta, xi1[:,Bin]))
        xi2_theta[Bin] = list(itp.splrep(theta, xi2[:,Bin]))

    xi_p = np.zeros((data.const['ntheta'], nzcorrs))
    xi_m = np.zeros((data.const['ntheta'], nzcorrs))
    if data.conf['use_theory_binning']:

        #t0 = timer()
        # roughly 0.01s to 0.02s extra...
        for idx_theta in range(data.const['ntheta']):
            theta = thetas_for_theory_binning[idx_theta, :]
            dtheta = (theta[1:] - theta[:-1]) * a2r

            for idx_bin in range(nzcorrs):

                xi_p_integrand = itp.splev(theta, xi1_theta[idx_bin]) * itp.splev(theta, theory_weight_func)
                xi_m_integrand = itp.splev(theta, xi2_theta[idx_bin]) * itp.splev(theta, theory_weight_func)

                xi_p[idx_theta, idx_bin] = np.sum(0.5 * (xi_p_integrand[1:] + xi_p_integrand[:-1]) * dtheta) / int_weight_func[idx_theta]
                xi_m[idx_theta, idx_bin] = np.sum(0.5 * (xi_m_integrand[1:] + xi_m_integrand[:-1]) * dtheta) / int_weight_func[idx_theta]

        # now mix xi_p and xi_m back into xi_obs:
        temp = np.concatenate((xi_p, xi_m))
        xi = OrderXiFunc(temp, data.const['ntheta'], nzcorrs)
        #dt = timer() - t0
        #print(dt)

    else:
        # Get xi's in same column vector format as the data
        #iz = 0
        #for Bin in range(nzcorrs):
        #    iz = iz + 1  # this counts the bin combinations
        #    for i in range(ntheta):
        #        j = (iz-1)*2*ntheta
        #        xi[j+i] = itp.splev(
        #            theta_bins[i], xi1_theta[Bin])
        #        xi[ntheta + j+i] = itp.splev(
        #            theta_bins[i], xi2_theta[Bin])
        # or in more compact/vectorizable form:
        iz = 0
        for Bin in range(nzcorrs):
            iz = iz + 1  # this counts the bin combinations
            j = (iz - 1) * 2 * ntheta
            xi[j:j + data.const['ntheta']] = itp.splev(theta_bins[:data.const['ntheta']], xi1_theta[Bin])
            xi[j + data.const['ntheta']:j + 2 * data.const['ntheta']] = itp.splev(theta_bins[:data.const['ntheta']], xi2_theta[Bin])






    ################################################
    # Shear bias
    ################################################
    # nuisance parameter for m-correction:
    dm_per_zbin = np.zeros((data.const['ntheta'], nzbins))
    if ('dm' in data.nuisance_parameters):
        for zbin in range(nzbins):
            dm_per_zbin[:, zbin] = np.ones(data.const['ntheta']) * data.nuisance_parameters['dm'][zbin]

    # nuisance parameters for constant c-correction:
    dc1_per_zbin = np.zeros((data.const['ntheta'], nzbins))
    dc2_per_zbin = np.zeros((data.const['ntheta'], nzbins))
    if ('dc1' in data.nuisance_parameters):
        for zbin in range(nzbins):
            dc1_per_zbin[:, zbin] = np.ones(data.const['ntheta']) * data.nuisance_parameters['dc1'][zbin]
    if ('dc2' in data.nuisance_parameters):
        for zbin in range(nzbins):
            dc2_per_zbin[:, zbin] = np.ones(data.const['ntheta']) * data.nuisance_parameters['dc2'][zbin]

    # correlate dc1/2_per_zbin in tomographic order of xi1/2:
    dc1_sqr = np.zeros((data.const['ntheta'], nzcorrs))
    dc2_sqr = np.zeros((data.const['ntheta'], nzcorrs))
    # correlate dm_per_zbin in tomographic order of xi1/2:
    dm_plus_one_sqr = np.zeros((data.const['ntheta'], nzcorrs))
    index_corr = 0
    for zbin1 in range(nzbins):
        for zbin2 in range(zbin1, nzbins):

            # c-correction:
            dc1_sqr[:, index_corr] = dc1_per_zbin[:, zbin1] * dc1_per_zbin[:, zbin2]
            dc2_sqr[:, index_corr] = dc2_per_zbin[:, zbin1] * dc2_per_zbin[:, zbin2]

            # m-correction:
            dm_plus_one_sqr[:, index_corr] = (1. + dm_per_zbin[:, zbin1]) * (1. + dm_per_zbin[:, zbin2])

            index_corr += 1

    # get c-correction into form of xi_obs
    temp = np.concatenate((dc1_sqr, dc2_sqr))
    dc_sqr = OrderXiFunc(temp, data.const['ntheta'], nzcorrs)

    # get m-correction into form of xi_obs
    temp = np.concatenate((dm_plus_one_sqr, dm_plus_one_sqr))
    dm_plus_one_sqr_obs = OrderXiFunc(temp, data.const['ntheta'], nzcorrs)

    # Below we construct a theta-dependent c-correction function from
    # measured data (for one z-bin) and scale it with an amplitude per z-bin
    # which is to be fitted
    # this is all independent of the constant c-correction calculated above
    xip_c = np.zeros((data.const['ntheta'], nzcorrs))
    xim_c = np.zeros((data.const['ntheta'], nzcorrs))
    # load theta-dependent c-term function if requested
    # file is assumed to contain values for the same theta values as used
    # for xi_pm!
    if data.conf['use_cterm_function']:
        fname = os.path.join(data.paths['data'], 'SUPPLEMENTARY_FILES/KV450_xi_pm_c_term.dat')
        # function is measured over same theta scales as xip, xim
        xip_c_per_zbin, xim_c_per_zbin = np.loadtxt(fname, usecols=(3, 4), unpack=True)
        print('Loaded (angular) scale-dependent c-term function from: \n', fname, '\n')

        amps_cfunc = np.ones(nzbins)
        for zbin in range(nzbins):
            if ('Ac'in data.nuisance_parameters):
                amps_cfunc[zbin] = data.nuisance_parameters['Ac']

        index_corr = 0
        for zbin1 in range(nzbins):
            for zbin2 in range(zbin1, nzbins):
                xip_c[:, index_corr] = amps_cfunc[zbin1] * amps_cfunc[zbin2] * xip_c_per_zbin
                # TODO: we leave xim_c set to 0 for now!
                #xim_c[:, index_corr] = amps_cfunc[zbin1] * amps_cfunc[zbin2] * xim_c_per_zbin
                index_corr += 1
    # get it into order of xi_obs
    # contains only zeros if function is not requested
    temp = np.concatenate((xip_c, xim_c))
    xipm_c = OrderXiFunc(temp, data.const['ntheta'], nzcorrs)

    xi = xi * dm_plus_one_sqr_obs + xipm_c + dc_sqr

    # write out masked theory vector in list format:    
    io_cs.WriteVectorFunc(data, nzbins, theta_bins, mask_indices, mask_suffix, xi, data.conf['sample'])
    print('Predicted vector saved.')
    print('All Done.')


def OrderXiFunc(temp, ntheta, nzcorrs):
        """
        This function takes xi_pm and constructs
        the xi_pm vector in its observed ordering:
         xi_obs = {xi_p(theta1, z1xz1)... xi_p(thetaK, z1xz1), xi_m(theta1, z1xz1)...
                   xi_m(thetaK, z1xz1);... xi_p(theta1, zNxzN)... xi_p(thetaK, zNxzN),
                   xi_m(theta1, zNxzN)... xi_m(thetaK, zNxN)}
        """

        xi_obs = np.zeros(ntheta * nzcorrs * 2)

        # create the data-vector:
        k = 0
        for j in range(nzcorrs):
            for i in range(2 * ntheta):
                xi_obs[k] = temp[i, j]
                k += 1

        return xi_obs


def IAFactorFunc(small_h, Omega_m, rho_crit,
                    z, linear_growth_rate, amplitude, exponent):

    const = 5e-14 / small_h**2 # Mpc^3 / M_sol

    # arbitrary convention
    z0 = 0.3
    factor = -1. * amplitude * const * rho_crit * Omega_m / linear_growth_rate * ((1. + z) / (1. + z0))**exponent

    return factor


def RhoCriticalFunc(small_h):
    """
    The critical density of the Universe at redshift 0.

    Returns
    -------
    rho_crit in solar masses per cubic Megaparsec.

    """

    # yay, constants...
    Mpc_cm = 3.08568025e24 # cm
    M_sun_g = 1.98892e33 # g
    G_const_Mpc_Msun_s = M_sun_g * (6.673e-8) / Mpc_cm**3.
    H100_s = 100. / (Mpc_cm * 1.0e-5) # s^-1

    rho_crit_0 = 3. * (small_h * H100_s)**2. / (8. * np.pi * G_const_Mpc_Msun_s)

    return rho_crit_0


#######################################################################################################
# This function is used to convert 2D sums over the two indices (Bin1, Bin2) of an N*N symmetric matrix
# into 1D sums over one index with N(N+1)/2 possible values
def one_dim_index(Bin1, Bin2, nzbins):
    if Bin1 <= Bin2:
        return int(Bin2 + nzbins * Bin1 - (Bin1 * (Bin1 + 1)) / 2) 
    else:
        return int(Bin1 + nzbins * Bin2 - (Bin2 * (Bin2 + 1)) / 2) 
