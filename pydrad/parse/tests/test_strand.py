"""
Test parsing HYDRAD results
"""
from pydrad.parse import Strand


def test_parse_initial_conditions(hydrad):
    s = Strand(hydrad)
    assert hasattr(s, 'initial_conditions')
