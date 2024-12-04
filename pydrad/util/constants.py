"""
Constants used throughout `pydrad`
"""
from astropy.constants import Constant

__all__ = [
    'gamma',
    'm_avg_ion',
]


gamma = Constant(
    abbrev='gamma',
    name='heat capacity ratio of a monatomic ideal gas',
    value=5/3,
    unit='',
    uncertainty=0,
    reference='https://en.wikipedia.org/wiki/Heat_capacity_ratio',
    system=None,
)


# NOTE: This is pulled directly from the HYDRAD source code
# https://github.com/rice-solar-physics/HYDRAD/blob/master/Resources/source/constants.h
# and is the average ion mass for a H-He plasma computed by weighting the H and He ion
# masses with a particular set of abundances.
# m_avg_ion = m_H + \frac{n_{He}}{n_{H}} m_{He}
m_avg_ion = Constant(
    abbrev='m_avg_ion',
    name='average ion mass for a H-He plasma',
    value=2.171e-24,
    unit='g',
    uncertainty=0,
    reference='https://github.com/rice-solar-physics/HYDRAD/blob/master/Resources/source/constants.h',
    system='cgs',
)
