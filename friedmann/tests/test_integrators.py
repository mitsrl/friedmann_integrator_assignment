"""Tests for Friedmann equation integrators."""

import pytest
import numpy as np

from friedmann import integrators

# Used google unit conversions.
HUBBLE_DISTANCE = 2997.9246     # Mpc/h


def test_matter_dominated():
    # Matter dominated, using units with little h in them.
    integ = integrators.CumulativeSumIntegrator(0, 1, 100)

    # From analytic calulation.
    correct_comoving_distance = 2 * HUBBLE_DISTANCE * (1 - np.sqrt(integ.a))
    assert np.allclose(integ.comoving_distance, correct_comoving_distance, rtol=1e-2)

    assert integ.hubble_distance == pytest.approx(HUBBLE_DISTANCE, 0.1)

def test_lambda_dominated():
    integ = integrators.CumulativeSumIntegrator(0, 0, 100)
    # For lambda domination Hubble parameter is constant.
    assert np.allclose(integ.hubble_parameter, 100)

def test_radiation_dominated():
    integ = integrators.CumulativeSumIntegrator(1, 0, 100)

