"""
Animate evolution of hydrodynamic quantities
"""
import numpy as np
from matplotlib.animation import FuncAnimation
import astropy.units as u

from .plot import _plot_profile, _setup_figure

__all__ = ['animate_strand']


def animate_strand(strand, start=0, stop=None, step=1, **kwargs):
    """
    Return a matplotlib animation of time-dependent hydrodynamic quantities of a
    strand. Takes the same arguments as #pydrad.visualize.plot_strand(). See the
    [matplotlib animation docs](https://matplotlib.org/api/animation_api.html) 
    for examples on how to write the movie to a file.
    """
    if stop is None:
        stop = strand.time.shape[0] + 1
    plot_kwargs = kwargs.get('plot_kwargs', {})
    limits = kwargs.get('limits', {})
    if 'limits' in kwargs:
        del kwargs['limits']
    # Initialize figure
    fig, axes = _setup_figure(strand[start], limits, **kwargs)
    if 'color' not in plot_kwargs:
        plot_kwargs['color'] = 'C0'
    l1a, l1b, l2a, l2b, l3a, l3b, l4 = _plot_profile(strand[start], axes, **plot_kwargs)
    # Define updater
    def update_plot(i):
        profile = strand[i]
        l1a.set_data(profile.coordinate.to(u.cm), profile.electron_temperature.to(u.MK))
        l1b.set_data(profile.coordinate.to(u.cm), profile.ion_temperature.to(u.MK))
        l2a.set_data(profile.coordinate.to(u.cm), profile.electron_density)
        l2b.set_data(profile.coordinate.to(u.cm), profile.ion_density)
        l3a.set_data(profile.coordinate.to(u.cm), profile.electron_pressure)
        l3b.set_data(profile.coordinate.to(u.cm), profile.ion_pressure)
        l4.set_data(profile.coordinate.to(u.cm), profile.velocity)
        fig.suptitle(r'$t={:4.0f}$ {}'.format(strand.time[i].value, strand.time[i].unit), y=0.905)
        return l1a, l1b, l2a, l2b, l3a, l3b, l4
    frames = [np.where(strand.time == t)[0][0] for t in strand.time[start:stop:step]]
    return FuncAnimation(
        fig, update_plot, blit=kwargs.get('blit', True), frames=frames,
        interval=kwargs.get('interval', 10), repeat=kwargs.get('repeat', True))
