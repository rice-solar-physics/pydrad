"""
Plotting methods to easily visualize HYDRAD results
"""
import copy

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.colors

__all__ = ['plot_strand',
           'plot_profile',
           'plot_time_distance',
           'plot_histogram',
           'plot_time_mesh',]


def plot_histogram(vals, bins, ax=None, **kwargs):
    """
    Given a set of bin edges and the values in each bin, plot
    the histogram.
    
    # Parameters
    ax (`matplotlib.pyplot): Matplotlib axis instance
    vals (array-like): value in each bin
    bins (array-like): Bin edges, including the rightmost edge
    kwargs : Plotting keyword arguments
    """
    if ax is None:
        fig = plt.figure()
        ax = plt.gca()
    ymin = ax.get_ylim()[0]
    ax.step(bins[:-1], vals, where='post', **kwargs)
    if 'label' in kwargs:
        del kwargs['label']
    ax.step(bins[-2:], [vals[-1], vals[-1]], where='pre', **kwargs)
    ax.vlines(bins[0], ymin, vals[0], **kwargs)
    ax.vlines(bins[-1], ymin, vals[-1], **kwargs)
    return ax


@u.quantity_input
def plot_time_distance(strand, quantities, delta_s: u.cm, **kwargs):
    """
    Make a time-distance plot for a particular quantity or
    quantities on a uniform grid.

    # Parameters
    strand (`#hydrad_tools.parse.Strand`):
    quantities (`str`, `list`, `tuple`): Name of quantity or quantities to plot.
    Optionally, you can also pass in a tuple of `(str,array-like)`, where `str` is
    the label and the second entry is the quantity to plot, already interpolated
    onto a common grid.
    norm (`dict`, optional): Dictionary of colormap normalizations; one per
    quantity.
    """
    grid = strand.get_uniform_grid(delta_s).to(kwargs.pop('space_unit', 'cm'))
    # Interpolate quantities to constant grid as needed
    quantities = copy.deepcopy(quantities)
    for i, q in enumerate(quantities):
        if type(q) is str:
            quantities[i] = (q, strand.to_constant_grid(q, grid))
    plot_time_mesh(strand, quantities, grid, r'$s$', **kwargs)


def plot_time_mesh(strand, quantities, y_grid, y_label, **kwargs):
    """
    Plot a given quantity as a function of some variable and time
    for a given strand.
    """
    y_mesh, t_mesh = np.meshgrid(y_grid.value, strand.time.value,)
    t_mesh = (t_mesh * strand.time.unit).to(kwargs.pop('time_unit', 's'))
    y_mesh = y_mesh * y_grid.unit
    quantities = copy.deepcopy(quantities)
    if type(quantities) is not list:
        quantities = [quantities]
    fig, ax = plt.subplots(
        len(quantities), 1,
        figsize=kwargs.pop('figsize', (10, 2.5*len(quantities))),
        sharex=True,
        sharey=True,
    )
    if len(quantities) == 1:
        ax = [ax]
    # NOTE: remove these here so we can send the rest to pcolormesh
    norm = kwargs.pop('norm', {})
    cmap = kwargs.pop('cmap', {})
    yscale = kwargs.pop('yscale', 'linear')
    for i, q in enumerate(quantities):
        label, data = q  # If q is a tuple of label, array
        im = ax[i].pcolormesh(
            t_mesh.value,
            y_mesh.value,
            data.value,
            cmap=cmap.get(label, 'viridis'),
            norm=norm.get(label, None),
            **kwargs,
        )
        cbar = fig.colorbar(im, ax=ax[i])
        if data.unit is u.dimensionless_unscaled:
            cbar.ax.set_ylabel(f'{label}')
        else:
            cbar.ax.set_ylabel(f'{label} [{data.unit}]')
    ax[-1].set_xlabel(f'$t$ [{t_mesh.unit}]')
    ax[-1].set_ylabel(f'{y_label} [{y_mesh.unit}]')
    ax[-1].set_yscale(yscale)


def plot_strand(strand, limits=None, cmap='viridis', **kwargs):
    """
    Plot hydrodynamic quantities at multiple timesteps

    # Parameters
    strand (#pydrad.parse.Strand): Loop strand object
    limits (`dict`): Set axes limits for hydrodynamic quantities, optional
    cmap (`str`): The colormap to map the timestep index to
    plot_kwargs (`dict`): Any keyword arguments used matplotlib.plot, optional
    figsize (`tuple`): Width and height of figure, optional
    """
    limits = {} if limits is None else limits
    plot_kwargs = kwargs.get('plot_kwargs', {})
    fig, axes = _setup_figure(strand[0], limits, **kwargs)
    colors = matplotlib.colors.LinearSegmentedColormap.from_list(
        '', plt.get_cmap(cmap).colors, N=len(strand))
    for i, p in enumerate(strand):
        plot_kwargs['color'] = colors(i)
        _ = _plot_profile(p, axes, **plot_kwargs)
    plt.show()


def plot_profile(profile, **kwargs):
    """
    Plot hydrodynamic quantites at a single timestep

    # Parameters
    profile (#pydrad.parse.Strand): Loop profile object
    limits (`dict`): Set axes limits for hydrodynamic quantities, optional
    plot_kwargs (`dict`): Any keyword arguments used matplotlib.plot, optional
    figsize (`tuple`): Width and height of figure, optional
    """
    limits = kwargs.pop('limits', {})
    plot_kwargs = kwargs.get('plot_kwargs', {})
    fig, axes = _setup_figure(profile, limits, **kwargs)
    _plot_profile(profile, axes, **plot_kwargs)
    plt.show()


def _setup_figure(profile, limits, **kwargs):
    # Setup frame
    fig, axes = plt.subplots(
        2, 2, figsize=kwargs.get('figsize', (10, 10)), sharex=True)
    # Limits
    axes[0, 0].set_ylim(limits.get('temperature', (0, 15)))
    axes[0, 1].set_ylim(limits.get('density', (1e8, 1e13)))
    axes[0, 1].set_yscale('log')
    axes[1, 0].set_ylim(limits.get('pressure', (0.1, 1e2)))
    axes[1, 0].set_yscale('log')
    axes[1, 1].set_ylim(limits.get('velocity', (-1e2, 1e2)))
    axes[1, 1].set_xlim(profile.coordinate[[0, -1]].to(u.Mm).value)
    # Labels
    axes[0, 0].set_ylabel(r'$T$ [MK]')
    axes[0, 1].set_ylabel(r'$n$ [cm$^{-3}$]')
    axes[1, 0].set_ylabel(r'$P$ [dyne cm$^{-2}$]')
    axes[1, 1].set_ylabel(r'$v$ [km s$^{-1}$]')
    axes[1, 0].set_xlabel(r'$s$ [Mm]')
    axes[1, 1].set_xlabel(r'$s$ [Mm]')

    return fig, axes


def _plot_profile(profile, axes, **kwargs):
    line1a, = axes[0, 0].plot(
        profile.coordinate.to(u.Mm),
        profile.electron_temperature.to(u.MK),
        **kwargs,
        ls='-'
    )
    line1b, = axes[0, 0].plot(
        profile.coordinate.to(u.Mm),
        profile.ion_temperature.to(u.MK),
        **kwargs,
        ls='--'
    )
    line2a, = axes[0, 1].plot(
        profile.coordinate.to(u.Mm),
        profile.electron_density.to(u.cm**(-3)),
        **kwargs,
        ls='-'
    )
    line2b, = axes[0, 1].plot(
        profile.coordinate.to(u.Mm),
        profile.ion_density.to(u.cm**(-3)),
        **kwargs,
        ls='--'
    )
    line3a, = axes[1, 0].plot(
        profile.coordinate.to(u.Mm),
        profile.electron_pressure.to(u.dyne / (u.cm**2)),
        **kwargs,
        ls='-'
    )
    line3b, = axes[1, 0].plot(
        profile.coordinate.to(u.Mm),
        profile.ion_pressure.to(u.dyne / (u.cm**2)),
        **kwargs,
        ls='--'
    )
    line4, = axes[1, 1].plot(
        profile.coordinate.to(u.Mm),
        profile.velocity.to(u.km/u.s),
        **kwargs
    )
    return line1a, line1b, line2a, line2b, line3a, line3b, line4
