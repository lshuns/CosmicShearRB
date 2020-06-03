#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 11:19:58 2017

@author: fkoehlin

@modified: Shun-Sheng Li
"""

import os
import sys
import glob
import numpy as np

# Bayesian way of defining confidence intervals:
# What's the difference to percentiles?
def minimum_credible_intervals(values, central_value, weights, bins=40):
    """
    Extract minimum credible intervals (method from Jan Haman) FIXME

    copy & paste from Monte Python (2.1.2) with own modifications
    --> checked that this function returns same output as Monte Python; modifications are all okay!!!

    """
    #histogram = info.hist
    #bincenters = info.bincenters
    #levels = info.levels

    histogram, bin_edges = np.histogram(values, bins=bins, weights=weights, normed=False)
    bincenters = 0.5*(bin_edges[1:]+bin_edges[:-1])

    # Defining the sigma contours (1, 2 and 3-sigma)
    levels = np.array([68.27, 95.45, 99.73])/100.

    bounds = np.zeros((len(levels), 2))
    j = 0
    delta = bincenters[1]-bincenters[0]
    left_edge = np.max(int(histogram[0] - 0.5*(histogram[1]-histogram[0])), 0)
    right_edge = np.max(int(histogram[-1] + 0.5*(histogram[-1]-histogram[-2])), 0)
    failed = False
    for level in levels:
        norm = float(
            (np.sum(histogram)-0.5*(histogram[0]+histogram[-1]))*delta)
        norm += 0.25*(left_edge+histogram[0])*delta
        norm += 0.25*(right_edge+histogram[-1])*delta
        water_level_up = np.max(histogram)*1.0
        water_level_down = np.min(histogram)*1.0
        top = 0.

        iterations = 0
        while (abs((top/norm)-level) > 0.0001) and not failed:
            top = 0.
            water_level = (water_level_up + water_level_down)/2.
            #ontop = [elem for elem in histogram if elem > water_level]
            indices = [i for i in range(len(histogram))
                       if histogram[i] > water_level]
            # check for multimodal posteriors
            '''
            if ((indices[-1]-indices[0]+1) != len(indices)):
                print('Could not derive minimum credible intervals for this multimodal posterior!')
                failed = True
                break
            '''
            top = (np.sum(histogram[indices]) -
                   0.5*(histogram[indices[0]]+histogram[indices[-1]]))*(delta)

            # left
            if indices[0] > 0:
                top += (0.5*(water_level+histogram[indices[0]]) *
                        delta*(histogram[indices[0]]-water_level) /
                        (histogram[indices[0]]-histogram[indices[0]-1]))
            else:
                if (left_edge > water_level):
                    top += 0.25*(left_edge+histogram[indices[0]])*delta
                else:
                    top += (0.25*(water_level + histogram[indices[0]]) *
                            delta*(histogram[indices[0]]-water_level) /
                            (histogram[indices[0]]-left_edge))

            # right
            if indices[-1] < (len(histogram)-1):
                top += (0.5*(water_level + histogram[indices[-1]]) *
                        delta*(histogram[indices[-1]]-water_level) /
                        (histogram[indices[-1]]-histogram[indices[-1]+1]))
            else:
                if (right_edge > water_level):
                    top += 0.25*(right_edge+histogram[indices[-1]])*delta
                else:
                    top += (0.25*(water_level + histogram[indices[-1]]) *
                            delta * (histogram[indices[-1]]-water_level) /
                            (histogram[indices[-1]]-right_edge))

            if top/norm >= level:
                water_level_down = water_level
            else:
                water_level_up = water_level
            # safeguard, just in case
            iterations += 1
            if (iterations > 1e4):
                print('The loop to check for sigma deviations was taking too long to converge.')
                break

        # min
        if indices[0] > 0:
            bounds[j][0] = bincenters[indices[0]] - delta*(histogram[indices[0]]-water_level)/(histogram[indices[0]]-histogram[indices[0]-1])
        else:
            if (left_edge > water_level):
                bounds[j][0] = bincenters[0]-0.5*delta
            else:
                bounds[j][0] = bincenters[indices[0]] - 0.5*delta*(histogram[indices[0]]-water_level)/(histogram[indices[0]]-left_edge)

        # max
        if indices[-1] < (len(histogram)-1):
            bounds[j][1] = bincenters[indices[-1]] + delta*(histogram[indices[-1]]-water_level)/(histogram[indices[-1]]-histogram[indices[-1]+1])
        else:
            if (right_edge > water_level):
                bounds[j][1] = bincenters[-1]+0.5*delta
            else:
                bounds[j][1] = bincenters[indices[-1]] + \
                    0.5*delta*(histogram[indices[-1]]-water_level) / \
                    (histogram[indices[-1]]-right_edge)

        j += 1

    for elem in bounds:
        for j in (0, 1):
            elem[j] -= central_value

    return bounds

def weighted_mean(values, weights=None):

    if weights is None:
        weights = np.ones_like(values)

    return np.sum(weights*values)/np.sum(weights)

def quantile(x, q, weights=None):
    """
    Like numpy.percentile, but:
    * Values of q are quantiles [0., 1.] rather than percentiles [0., 100.]
    * scalar q not supported (q must be iterable)
    * optional weights on x
    """
    if weights is None:
        return np.percentile(x, [100. * qi for qi in q])
    else:
        idx = np.argsort(x)
        xsorted = x[idx]
        cdf = np.add.accumulate(weights[idx])
        cdf /= cdf[-1]

        return np.interp(q, cdf, xsorted).tolist()

def get_values_and_intervals(parameters, weights, use_median=False):

    param_values = np.zeros((len(parameters), 7))
    confidence_values = np.zeros((len(parameters), 6))

    for idx, param in enumerate(parameters):

        if use_median:
            central_value = quantile(param, [0.5], weights=weights)[0]
        else:
            central_value = weighted_mean(param, weights=weights)

        # bounds returns [[-1sigma, +1sigma],[-2sigma, +2sigma], [-3sigma, +3sigma]]
        bounds = minimum_credible_intervals(param, central_value, weights, bins=50)

        param_values[idx, :] = np.concatenate(([central_value], bounds[:,0], bounds[:,1]))
        confidence_values[idx, :] = central_value + bounds.flatten()

    return param_values, confidence_values


def write_parameters_to_file(fname, best_fit_params, fit_statistics, param_values_mean, confidence_values_mean, param_values_median, confidence_values_median, labels, labels_tex):

    with open(fname, 'w') as f:
        f.write('# Best fitting values: \n')
        f.write('\chi^2 = {:.4f}, \chi^2_red = {:.4f} ({:} d.o.f.), index in chain = {:.0f} \n'.format(fit_statistics[0], fit_statistics[1], int(fit_statistics[2]), fit_statistics[3]))
        for index, label in enumerate(labels):
            name = label +':'
            f.write(name.ljust(20, ' ')+'{:.4f} \n'.format(best_fit_params[index]))
        ### (weighted) MEAN ###
        f.write('\n'+'# parameter, MEAN, err_minus (68%), err_plus (68%), MEAN, err_minus (95%), err_plus (95%), MEAN, err_minus (99%), err_plus (99%) \n')
        for index, label in enumerate(labels):
            name = label +':'
            f.write(name.ljust(20, ' ') + '{0:.4f} {1:.4f} +{2:.4f}, {0:.4f} {3:.4f} +{4:.4f}, {0:.4f} {5:.4f} +{6:.4f} \n'.format(param_values_mean[index, 0], param_values_mean[index, 1], param_values_mean[index, 4], param_values_mean[index, 2], param_values_mean[index, 5], param_values_mean[index, 3], param_values_mean[index, 6]))
        f.write('\n'+'# parameter, lower bound (68%), upper bound (68%), lower bound (95%), upper bound (95%), lower bound (99%), upper bound (99%) \n')
        for index, label in enumerate(labels):
            name = label +':'
            f.write(name.ljust(20, ' ')+'1sigma >{:.4f}, 1sigma <{:.4f}, 2sigma >{:.4f}, 2sigma <{:.4f}, 3sigma >{:.4f}, 3sigma <{:.4f} \n'.format(confidence_values_mean[index, 0], confidence_values_mean[index, 1], confidence_values_mean[index, 2], confidence_values_mean[index, 3], confidence_values_mean[index, 4], confidence_values_mean[index, 5]))
        ### (weighted) MEDIAN ###
        f.write('\n'+'# parameter, MEDIAN, err_minus (68%), err_plus (68%), MEDIAN, err_minus (95%), err_plus (95%), MEDIAN, err_minus (99%), err_plus (99%) \n')
        for index, label in enumerate(labels):
            name = label +':'
            f.write(name.ljust(20, ' ') + '{0:.4f} {1:.4f} +{2:.4f}, {0:.4f} {3:.4f} +{4:.4f}, {0:.4f} {5:.4f} +{6:.4f} \n'.format(param_values_median[index, 0], param_values_median[index, 1], param_values_median[index, 4], param_values_median[index, 2], param_values_median[index, 5], param_values_median[index, 3], param_values_median[index, 6]))
        f.write('\n'+'# parameter, lower bound (68%), upper bound (68%), lower bound (95%), upper bound (95%), lower bound (99%), upper bound (99%) \n')
        for index, label in enumerate(labels):
            name = label +':'
            f.write(name.ljust(20, ' ')+'1sigma >{:.4f}, 1sigma <{:.4f}, 2sigma >{:.4f}, 2sigma <{:.4f}, 3sigma >{:.4f}, 3sigma <{:.4f} \n'.format(confidence_values_median[index, 0], confidence_values_median[index, 1], confidence_values_median[index, 2], confidence_values_median[index, 3], confidence_values_median[index, 4], confidence_values_median[index, 5]))
        ### (weighted) MEAN (TeX) ###
        f.write('\n'+'\n'+'\n'+'### TeX ###'+'\n'+'# parameter, MEAN, err_minus (68%), err_plus (68%), MEAN, err_minus (95%), err_plus (95%), MEAN, err_minus (99%), err_plus (99%) \n')
        for index, label in enumerate(labels_tex):
            name = label +':'
            f.write(name.ljust(20, ' ')+'{0:.3f}_{{{1:.3f}}}^{{+{2:.3f}}}, {0:.3f}_{{{3:.3f}}}^{{+{4:.3f}}}, {0:.3f}_{{{5:.3f}}}^{{+{6:.3f}}} \n'.format(param_values_mean[index, 0], param_values_mean[index, 1], param_values_mean[index, 4], param_values_mean[index, 2], param_values_mean[index, 5], param_values_mean[index, 3], param_values_mean[index, 6]))
        ### (weighted) MEDIAN (TeX) ###
        f.write('\n'+'\n'+'\n'+'### TeX ###'+'\n'+'# parameter, MEDIAN, err_minus (68%), err_plus (68%), MEDIAN, err_minus (95%), err_plus (95%), MEDIAN, err_minus (99%), err_plus (99%) \n')
        for index, label in enumerate(labels_tex):
            name = label +':'
            f.write(name.ljust(20, ' ')+'{0:.3f}_{{{1:.3f}}}^{{+{2:.3f}}}, {0:.3f}_{{{3:.3f}}}^{{+{4:.3f}}}, {0:.3f}_{{{5:.3f}}}^{{+{6:.3f}}} \n'.format(param_values_median[index, 0], param_values_median[index, 1], param_values_median[index, 4], param_values_median[index, 2], param_values_median[index, 5], param_values_median[index, 3], param_values_median[index, 6]))

    print 'File saved to: \n', fname

    return


if __name__ == '__main__':

    # path_to_chain = './KV450_H0_Dz_IA/'
    path_to_chain = sys.argv[1]

    fname =  glob.glob(path_to_chain + '*.txt')[0]
    data = np.loadtxt(fname)
    weights = data[:, 0]
    #print data
    #print data[:, -1]
    # glob can expand names with *-operator!
    fname =  glob.glob(path_to_chain + '*.paramnames')[0]
    print fname
    names = np.loadtxt(fname, dtype=str, delimiter='\t')
    print names, names.shape
    labels = names[:, 0]
    labels_tex = names[:, 1]

    chi2 = 2. * data[:, 1]
    min_chi2 = chi2.min()
    best_fit_index = np.where(data[:, 1] == data[:, 1].min())
    print best_fit_index
    best_fit_params = data[best_fit_index]
    fit_statistics = np.array([min_chi2, 0., 0., int(best_fit_index[0])])

    params_mean, conf_mean = get_values_and_intervals(data[:, 2:].T, weights, use_median=False)
    params_median, conf_median = get_values_and_intervals(data[:, 2:].T, weights, use_median=True)

    fname = os.path.join(path_to_chain, 'parameter_table.txt')
    write_parameters_to_file(fname, best_fit_params[0, 2:], fit_statistics, params_mean, conf_mean, params_median, conf_median, labels, labels_tex)