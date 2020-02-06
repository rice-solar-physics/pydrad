"""
Animate evolution of hydrodynamic quantities
"""
from matplotlib.animation import FuncAnimation
import astropy.units as u

from .plot import _plot_profile, _setup_figure

__all__ = ['animate_strand']


def animate_strand(strand, **kwargs):
    """
    Return a matplotlib animation of time-dependent hydrodynamic quantities of a
    strand. Takes the same arguments as #pydrad.visualize.plot_strand(). See the
    [matplotlib animation docs](https://matplotlib.org/api/animation_api.html) 
    for examples on how to write the movie to a file.
    """
    plot_kwargs = kwargs.get('plot_kwargs', {})
    limits = kwargs.get('limits', {})
    if 'limits' in kwargs:
        del kwargs['limits']
    fig, axes = _setup_figure(strand[0], limits, **kwargs)
    # Make sure the lines stay the same color
    if 'color' not in plot_kwargs:
        plot_kwargs['color'] = 'C0'
    l1a, l1b, l2a, l2b, l3a, l3b, l4 = _plot_profile(strand[0],
                                                     axes,
                                                     **plot_kwargs)

    def update_plot(i):
        p = strand[i]
        l1a.set_data(p.coordinate.to(u.Mm), p.electron_temperature.to(u.MK))
        l1b.set_data(p.coordinate.to(u.Mm), p.ion_temperature.to(u.MK))
        l2a.set_data(p.coordinate.to(u.Mm), p.electron_density.to(u.cm**(-3)))
        l2b.set_data(p.coordinate.to(u.Mm), p.ion_density.to(u.cm**(-3)))
        l3a.set_data(p.coordinate.to(u.Mm), p.electron_pressure.to(u.dyne/(u.cm**2)))
        l3b.set_data(p.coordinate.to(u.Mm), p.ion_pressure.to(u.dyne/(u.cm**2)))
        l4.set_data(p.coordinate.to(u.Mm), p.velocity.to(u.km/u.s))
        fig.suptitle(r'$t={:.0f}$ {}'.format(
            strand.time[i].value, strand.time[i].unit), y=0.905)
        return l1a, l1b, l2a, l2b, l3a, l3b, l4
    
    return FuncAnimation(
        fig, update_plot, blit=kwargs.get('blit', True), frames=len(strand),
        interval=kwargs.get('interval', 10), repeat=kwargs.get('repeat', True))
