## STAN DIAGNOSTICS

import matplotlib.pyplot as plt
from matplotlib.pyplot import rcParams
import numpy as np
import itertools
import configparser

config = configparser.ConfigParser()
config.read('../../stan_plots_colors.config')
colors = config

def autolabel(rects, ax):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()*0.5, height,
                '%.3g' % height,
                ha='center', va='bottom')

        
def parameter_plot(parameter_samples, chains, parameter_name, ax = None, show = True, true_value = None):
    ind = np.arange(chains)  # the x locations for the groups
    width = 0.75             # the width of the bars

    if ax is None:
        fig, ax = plt.subplots()
    parameter_means = []
    parameter_stds  = []
    parameter_5perc = []
    parameter_95perc = []
    l = int(len(parameter_samples)/chains)
    for c in range(chains):
        parameter_means.append(np.mean(parameter_samples[c*(l):(c+1)*l]))
        parameter_stds.append(2*np.std(parameter_samples[c*(l):(c+1)*l]))
        parameter_5perc.append(np.percentile(parameter_samples[c*(l):(c+1)*l],5))
        parameter_95perc.append(np.percentile(parameter_samples[c*(l):(c+1)*l],95))
    rects = ax.bar(ind, parameter_means, width, color='#'+colors['PARAMETER_BAR_PLOT']['bar'],
                   zorder = 2) #, yerr=parameter_stds)
    if true_value is not None:
        ax.axhline(y=true_value, xmin=0, xmax=chains, linewidth=2, color = '#'+colors['PARAMETER_BAR_PLOT']['true_value'])
    # add some text for labels, title and axes ticks
    ax.set_title(parameter_name)
    autolabel(rects, ax)
    # add percentiles
    for i, rect in enumerate(rects):
        height = rect.get_height()
        ax.scatter([rect.get_x() + rect.get_width()*0.5,
                    rect.get_x() + rect.get_width()*0.5], 
                   [parameter_5perc[i], parameter_95perc[i]],
                   c = '#'+colors['PARAMETER_BAR_PLOT']['error'],
                   zorder = 3)
    if show: 
        plt.show()

def chain_plot(parameter_samples, chains, parameter_name, ax = None, show = True, true_value = None):
    if ax is None:
        fig, ax = plt.subplots()
    L = len(parameter_samples)
    for c in range(chains):
        ax.plot(parameter_samples[(c*int(L/chains)):((c+1)*int(L/chains))], color = '#'+colors['CHAIN_PLOT']['chain'+str(c)])
    ax.set_title(parameter_name)
    if true_value is not None:
        ax.axhline(y=true_value, xmin=0, xmax=len(parameter_samples), linewidth=2, color = '#'+colors['CHAIN_PLOT']['true_value'])
    if show:
        plt.show()
    
def scatterplot_matrix(data, names, **kwargs):
    """Plots a scatterplot matrix of subplots.  Each row of "data" is plotted
    against other rows, resulting in a nrows by nrows grid of subplots with the
    diagonal subplots labeled with "names".  Additional keyword arguments are
    passed on to matplotlib's "plot" command. Returns the matplotlib figure
    object containg the subplot grid."""
    numvars, numdata = data.shape
    fig, axes = plt.subplots(nrows=numvars, ncols=numvars, figsize=(numvars*3,numvars*3))
    fig.subplots_adjust(hspace=0.05, wspace=0.05)

    for ax in axes.flat:
        # Hide all ticks and labels
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        # Set up ticks only on one side for the "edge" subplots...
        if ax.is_first_col():
            ax.yaxis.set_ticks_position('left')
        if ax.is_last_col():
            ax.yaxis.set_ticks_position('right')
        if ax.is_first_row():
            ax.xaxis.set_ticks_position('top')
        if ax.is_last_row():
            ax.xaxis.set_ticks_position('bottom')

    # Plot the data.
    for i, j in zip(*np.triu_indices_from(axes, k=1)):
        for x, y in [(i,j), (j,i)]:
            axes[x,y].plot(data[x], data[y], **kwargs)

    # Label the diagonal subplots...
    for i, label in enumerate(names):
        axes[i,i].annotate(label, (0.5, 0.5), xycoords='axes fraction',
                ha='center', va='center')

    # Turn on the proper x or y axes ticks.
    for i, j in zip(range(numvars), itertools.cycle((-1, 0))):
        axes[j,i].xaxis.set_visible(True)
        axes[i,j].yaxis.set_visible(True)

    return fig

def sample_plots(samples, chains, true_values = dict(), plots = ['bar','chains','scatter_matrix']):
    if 'bar' in plots:
        rcParams['figure.figsize'] = 3*len(samples), 3-0.22
        fig, axes = plt.subplots(1, len(samples))
        fig.subplots_adjust(hspace=0.22, wspace=0.22)
        p = 0
        for key, value in samples.items():
            if key in true_values:
                true_value = true_values[key]
            else:
                true_value = None
            parameter_plot(value, chains, key, axes[p], False, true_value)
            p += 1
        plt.show()
    
    if 'chains' in plots:
        rcParams['figure.figsize'] = 3*len(samples), 3*len(samples)
        fig, axes = plt.subplots(len(samples),1)
        fig.subplots_adjust(hspace=0.22, wspace=0.22)
        p = 0
        for key, value in samples.items():
            if key in true_values:
                true_value = true_values[key]
            else:
                true_value = None
            chain_plot(value, chains, key, axes[p], False, true_value)
            p += 1
        plt.show()
    
    if 'scatter_matrix' in plots:
        l = []
        names = samples.keys()
        for n in names:
            l.append(samples[n])
        data= np.array(l)
        scatterplot_matrix(data, names, markersize = colors['SCATTER_MATRIX']['markersize'], linewidth = 0, 
                           marker = '.',
                           markerfacecolor = '#'+colors['SCATTER_MATRIX']['markerfacecolor'], 
                           markeredgecolor =  '#'+colors['SCATTER_MATRIX']['markeredgecolor'],
                           )
        plt.show()