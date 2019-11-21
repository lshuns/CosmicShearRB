import getdist
import getdist.plots
import sys
import string
import matplotlib.pyplot as plt
from getdist import loadMCSamples
import numpy as np

#plt.rcParams["mathtext.fontset"] = "cm"

chain_dir = sys.argv[1]
chains = string.split(sys.argv[2])
labels = string.split(sys.argv[3])
color_choice = string.split(sys.argv[4])
burn_remove = string.split(sys.argv[5])
plot_file = sys.argv[6]

### define colours ###

# Colorbrewer2 8-class Accent
#colors_all = ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f', '#bf5b17', '#666666']

# Colorbrewer2 8-class Set1 + gold
colors_all = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#FFD700']

colors = []
for i in range(len(chains)):
    colors.append(colors_all[int(color_choice[i])])    

###  load the chains and add S8 ###
samples = []

for i in range(len(chains)):
    if int(burn_remove[i]) == 1:
        #samples.append(loadMCSamples(chain_dir+chains[i], settings={'ignore_rows':0.3, "smooth_scale_2D":1.0})) # for Marika's emcee chain
        #samples.append(loadMCSamples(chain_dir+chains[i], settings={'ignore_rows':0.54, "smooth_scale_2D":0.2})) # for Tilman's emcee chain
        #samples.append(loadMCSamples(chain_dir+chains[i], settings={'ignore_rows':0.3, "smooth_scale_2D":0.4})) # for Shahab's chain
        samples.append(loadMCSamples(chain_dir+chains[i], settings={'ignore_rows':0.3}))
    else:
        samples.append(loadMCSamples(chain_dir+chains[i]))
    p = samples[i].getParams()
    samples[i].addDerived(p.sigma8*np.sqrt(p.omegam/0.3), name='S_8', label='S_8')
        
### plot ###
g = getdist.plots.getSinglePlotter(chain_dir = chain_dir)
g.settings.rcSizes(12.0,15.0,12.0)
params = g.get_param_array(samples[0], ['omegam', 'sigma8'])
g.plot_2d(samples,'omegam','sigma8',filled=True,colors=colors)
#g.plot_2d(samples[0],'omegam','sigma8',filled=True,colors=(colors[0],))
#g.plot_2d(samples[4],'omegam','sigma8',filled=True,colors=(colors[4],))
#g.add_legend(labels, legend_loc='upper right', colored_text=True, fontsize=12, label_order=(4,1,2,0,3))
g.add_legend(labels, legend_loc='upper right', colored_text=True, fontsize=12)
g.setAxes(params, lims=[0.1, 0.52, 0.51, 1.28])
g.export(plot_file)
