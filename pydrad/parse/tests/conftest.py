"""
Common fixtures for testing pydrad.parse
"""
import pytest

from pydrad.parse import Strand


@pytest.fixture
def strand(hydrad):
    return Strand(hydrad)
