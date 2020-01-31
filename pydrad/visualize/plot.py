"""
Plotting methods to easily visualize HYDRAD results
"""
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.colors

__all__ = ['plot_strand']


def plot_strand(strand, start=0, stop=None, step=1, **kwargs):
    """
    Plot hydrodynamic quantities at multiple timesteps

    # Parameters
    strand (#pydrad.parse.Strand): Loop strand object
    start (`int`): Starting time index, optional
    stop (`int`): Final time index, optional
    step (`int`): Number of steps between successive timesteps
    limits (`dict`): Set axes limits for hydrodynamic quantities, optional
    plot_kwargs (`dict`): Any keyword arguments used matplotlib.plot, optional
    figsize (`tuple`): Width and height of figure, optional
    """
    if stop is None:
        stop = strand.time.shape[0] + 1
    limits = kwargs.get('limits', {})
    if 'limits' in kwargs:
        del kwargs['limits']
    plot_kwargs = kwargs.get('plot_kwargs', {})
    fig, axes = _setup_figure(strand[0], limits, **kwargs)
    colors = matplotlib.colors.LinearSegmentedColormap.from_list(
        '', plt.get_cmap(kwargs.get('cmap', 'viridis')).colors,
        N=strand.time[start:stop:step].shape[0])
    # NOTE: once strand indexing is fixed, we can index it directly
    for i, t in enumerate(strand.time[start:stop:step]):
        j = np.where(t == strand.time)[0][0]
        plot_kwargs['color'] = colors(i)
        _ = _plot_profile(strand[j], axes, **plot_kwargs)
    plt.show()


def _setup_figure(profile, limits, **kwargs):
    # Setup frame
    fig, axes = plt.subplots(2, 2, figsize=kwargs.get('figsize', (10, 10)), sharex=True)
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
    line1a, = axes[0, 0].plot(profile.coordinate.to(u.cm), profile.electron_temperature.to(u.MK),
                              **kwargs, ls='-')
    line1b, = axes[0, 0].plot(profile.coordinate.to(u.cm), profile.ion_temperature.to(u.MK),
                              **kwargs, ls='--')
    line2a, = axes[0, 1].plot(profile.coordinate.to(u.cm), profile.electron_density,
                              **kwargs, ls='-')
    line2b, = axes[0, 1].plot(profile.coordinate.to(u.cm), profile.ion_density, **kwargs, ls='--')
    line3a, = axes[1, 0].plot(profile.coordinate.to(u.cm), profile.electron_pressure,
                              **kwargs, ls='-')
    line3b, = axes[1, 0].plot(profile.coordinate.to(u.cm), profile.ion_pressure, **kwargs, ls='--')
    line4, = axes[1, 1].plot(profile.coordinate.to(u.cm), profile.velocity, **kwargs)
    return line1a, line1b, line2a, line2b, line3a, line3b, line4
