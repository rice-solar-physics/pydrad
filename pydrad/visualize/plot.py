"""
Plotting methods to easily visualize HYDRAD results
"""
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.colors

__all__ = ['plot_strand', 'plot_profile', 'plot_time_distance']


@u.quantity_input
def plot_time_distance(strand, quantities, delta_s: u.cm, **kwargs):
    """
    Make a time-distance plot for a particular quantity or
    quantities on a uniform grid.

    # Parameters
    strand (`#hydrad_tools.parse.Strand`):
    quantities (`str`, `list`): Name of quantity or quantities to plot
    """
    grid = strand.get_uniform_grid(delta_s)
    s_mesh, t_mesh = np.meshgrid(grid.value, strand.time.value,)
    t_mesh = (t_mesh*strand.time.unit).to(kwargs.get('time_unit', 's'))
    s_mesh = (s_mesh*grid.unit).to(kwargs.get('space_unit', 'cm'))
    if type(quantities) is str:
        quantities = [quantities]
    fig, ax = plt.subplots(
        len(quantities), 1,
        figsize=kwargs.get('figsize', (10, 2.5*len(quantities))),
        sharex=True,
        sharey=True,
    )
    if len(quantities) == 1:
        ax = [ax]
    for i, q in enumerate(quantities):
        q_uni = strand.to_constant_grid(q, grid)
        im = ax[i].pcolormesh(
            t_mesh.value,
            s_mesh.value,
            q_uni.value,
            cmap=kwargs.get('cmap', 'plasma'),
            norm=kwargs.get('norm', {}).get(q, None),
        )
        cbar = fig.colorbar(im, ax=ax[i])
        cbar.ax.set_ylabel(f'{q} [{q_uni.unit}]')
    ax[i].set_xlabel(f'$t$ [{t_mesh.unit}]')
    ax[i].set_ylabel(f'$s$ [{s_mesh.unit}]')


def plot_strand(strand, limits=None, cmap='viridis', **kwargs):
    """
    Plot hydrodynamic quantities at multiple timesteps

    # Parameters
    strand (#pydrad.parse.Strand): Loop strand object
    limits (`dict`): Set axes limits for hydrodynamic quantities, optional
    plot_kwargs (`dict`): Any keyword arguments used matplotlib.plot, optional
    figsize (`tuple`): Width and height of figure, optional
    """
    limits = {} if limits is None else limits
    plot_kwargs = kwargs.get('plot_kwargs', {})
    fig, axes = _setup_figure(strand[0], limits, **kwargs)
    colors = matplotlib.colors.LinearSegmentedColormap.from_list(
        '', plt.get_cmap(cmap).colors, N=len(strand))
    # NOTE: once strand indexing is fixed, we can index it directly
    for i, p in enumerate(strand):
        plot_kwargs['color'] = colors(i)
        _ = _plot_profile(p, axes, **plot_kwargs)
    plt.show()


def plot_profile(profile, **kwargs):
    limits = kwargs.get('limits', {})
    if 'limits' in kwargs:
        del kwargs['limits']
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
    axes[1, 1].set_ylim(limits.get('velocity', (-5e7, 5e7)))
    axes[1, 1].set_xlim(0, profile.coordinate[-1].to(u.cm).value)
    # Labels
    axes[0, 0].set_ylabel(r'$T$ [MK]')
    axes[0, 1].set_ylabel(r'$n$ [cm$^{-3}$]')
    axes[1, 0].set_ylabel(r'$P$ [dyne cm$^{-2}$ s$^{-1}$]')
    axes[1, 1].set_ylabel(r'$v$ [cm s$^{-1}$]')
    axes[1, 0].set_xlabel(r'$s$ [Mm]')
    axes[1, 1].set_xlabel(r'$s$ [Mm]')

    return fig, axes


def _plot_profile(profile, axes, **kwargs):
    line1a, = axes[0, 0].plot(
        profile.coordinate.to(u.cm),
        profile.electron_temperature.to(u.MK),
        **kwargs,
        ls='-'
    )
    line1b, = axes[0, 0].plot(
        profile.coordinate.to(u.cm),
        profile.ion_temperature.to(u.MK),
        **kwargs,
        ls='--'
    )
    line2a, = axes[0, 1].plot(
        profile.coordinate.to(u.cm),
        profile.electron_density,
        **kwargs,
        ls='-'
    )
    line2b, = axes[0, 1].plot(
        profile.coordinate.to(u.cm),
        profile.ion_density,
        **kwargs,
        ls='--'
    )
    line3a, = axes[1, 0].plot(
        profile.coordinate.to(u.cm),
        profile.electron_pressure,
        **kwargs,
        ls='-'
    )
    line3b, = axes[1, 0].plot(
        profile.coordinate.to(u.cm),
        profile.ion_pressure,
        **kwargs,
        ls='--'
    )
    line4, = axes[1, 1].plot(
        profile.coordinate.to(u.cm),
        profile.velocity,
        **kwargs
    )
    return line1a, line1b, line2a, line2b, line3a, line3b, line4
