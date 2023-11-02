"""
Plotting methods to easily visualize HYDRAD results
"""
import copy

import astropy.units as u
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

__all__ = ['plot_strand',
           'plot_profile',
           'plot_time_distance',
           'plot_histogram',
           'plot_time_mesh']


def plot_histogram(vals, bins, ax=None, **kwargs):
    """
    Given a set of bin edges and the values in each bin, plot
    the histogram.

    Parameters
    ----------
    vals : array-like
        value in each bin
    bins : array-like
        Bin edges, including the rightmost edge
    ax : `matplotlib.pyplot`, optional
        Matplotlib axis instance
    kwargs : `dict`
        Plotting keyword arguments

    Returns
    -------
    ax : `matplotlib.pyplot`
        Axes instance with histogram plot attached
    """
    if ax is None:
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

    Parameters
    ----------
    strand : `~pydrad.parse.Strand`
    quantities : `str`, `list`, `tuple`
        Name of quantity or quantities to plot. Optionally, you can also pass
        in a tuple of `(str, array-like)`, where `str` is the label and the
        second entry is the quantity to plot, already interpolated onto a
        common grid.
    delta_s : `~astropy.units.Quantity`
        Spacing of the uniform spatial grid to interpolate the quantity onto
    space_unit : `str` or `astropy.quantity.Unit`, optional
        Unit for the spatial axis

    Optional Parameters
    -------------------
    All other parameters are passed to `~pydrad.visualize.plot_time_mesh`.
    """
    # NOTE: I'm setting this as the default here so that the quadrilaterals
    # (i.e. the cells represented in the mesh) are centered on grid points.
    # This is because the time and spatial grids here represent the centers
    # of the grid points and thus should be centered on their faces.
    shading = kwargs.pop('shading', 'nearest')
    grid = strand.get_uniform_grid(delta_s).to(kwargs.pop('space_unit', 'cm'))
    # Interpolate quantities to constant grid as needed
    quantities = copy.deepcopy(quantities)
    if type(quantities) is not list:
        quantities = [quantities]
    for i, q in enumerate(quantities):
        if type(q) is str:
            quantities[i] = (q, strand.to_constant_grid(q, grid))
    fig, ax = plot_time_mesh(strand, quantities, grid, r'$s$', shading=shading, **kwargs)
    return fig, ax


def plot_time_mesh(strand, quantities, y_grid, y_label, **kwargs):
    """
    Plot a given quantity as a function of some variable and time
    for a given strand.

    Parameters
    ----------
    strand : `~pydrad.parse.Strand`
    quantities : `list`, `tuple`
        List of or single tuple of `(str, array-like)`, where `str`
        is the label and the second entry is the quantity to plot,
        already interpolated onto a common grid.
    y_grid : `~astropy.units.Quantity`
        Grid other than time to interpolate the quantity onto.
    y_label : `str`
        Axis label for the other dimension against which to plot the quantity
        or quantities.
    norm : normalization or `dict`, optional
        Colormap normalization or dictionary of normalizations with keys
        corresponding to the quantity names.
    cmap : `str` or colormap instance or `dict`, optional
        Colormap to use for all quantities except velocity which will always
        use the ``RdBu_r`` diverging colormap. If `dict`, mapping of colormaps
        to different quantity names.
    figsize : `tuple`, optional
        Dimensions of the resulting figure
    time_unit : `str` or `astropy.quantity.Unit`, optional
        Unit for the time axis

    Optional Parameters
    -------------------
    All other keyword arguments are passed to `matplotlib.pcolormesh`.
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
    units = kwargs.pop('units', {})
    labels = kwargs.pop('labels', {})
    yscale = kwargs.pop('yscale', 'linear')
    cbar_size = kwargs.pop('cbar_size', '5%')
    cbar_pad = kwargs.pop('cbar_pad', '1%')
    for i, (label, data) in enumerate(quantities):
        data = data.to(units.get(label, data.unit))
        im = ax[i].pcolormesh(
            t_mesh.value,
            y_mesh.value,
            data.value,
            cmap=cmap.get(label, 'viridis') if isinstance(cmap, dict) else cmap,
            norm=norm.get(label, None) if isinstance(norm, dict) else norm,
            **kwargs,
        )
        cax = make_axes_locatable(ax[i]).append_axes(
            "right",
            size=cbar_size,
            pad=cbar_pad,
        )
        cbar = fig.colorbar(im, cax=cax)
        cbar_label = labels.get(label, label)
        if data.unit != u.dimensionless_unscaled:
            cbar_label += f' [{data.unit}]'
        cbar.set_label(cbar_label)
        ax[i].set_ylabel(f'{y_label} [{y_mesh.unit}]')
    ax[-1].set_xlabel(f'$t$ [{t_mesh.unit}]')
    ax[-1].set_yscale(yscale)

    return fig, ax


def plot_strand(strand, limits=None, cmap='viridis', **kwargs):
    """
    Plot hydrodynamic quantities at multiple timesteps

    This function takes all of the same keyword arguments as
    `~pydrad.visualize.plot_profile`.

    Parameters
    ----------
    strand : `~pydrad.parse.Strand`
    limits : `dict`, optional
        Axes limits for hydrodynamic quantities
    cmap : `str` or colormap instance
        The colormap to map the timestep index to

    See Also
    --------
    plot_profile
    """
    limits = {} if limits is None else limits
    plot_kwargs = kwargs.get('plot_kwargs', {})
    fig, axes = _setup_figure(strand[0], limits, **kwargs)
    colors = matplotlib.colors.LinearSegmentedColormap.from_list(
        '', plt.get_cmap(cmap).colors, N=len(strand))
    for i, p in enumerate(strand):
        plot_kwargs['color'] = colors(i)
        _ = _plot_profile(p, axes, **plot_kwargs)


def plot_profile(profile, **kwargs):
    """
    Plot hydrodynamic quantites at a single timestep

    Parameters
    ----------
    profile : `~pydrad.parse.Profile`
    limits : `dict`, optional
        Axes limits for hydrodynamic quantities
    plot_kwargs : `dict`, optional
        Any keyword arguments used `~matplotlib.pyplot.plot`
    figsize : `tuple`, optional
        Width and height of figure
    """
    limits = kwargs.pop('limits', {})
    plot_kwargs = kwargs.get('plot_kwargs', {})
    fig, axes = _setup_figure(profile, limits, **kwargs)
    _plot_profile(profile, axes, **plot_kwargs)


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
