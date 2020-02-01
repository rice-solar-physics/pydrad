"""
Plotting methods to easily visualize HYDRAD results
"""
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.colors

__all__ = ['plot_strand', 'plot_profile']


def plot_time_distance(strand, name, **kwargs):
    """
    Make a time-distance plot for a particular quantity

    # Parameters
    strand (`#hydrad_tools.parse.Strand`):
    name (`str`): Name of quantity
    """
    ...


def _setup_time_distance_figure():
    ...


def _setup_time_distance_axis():
    ...


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
