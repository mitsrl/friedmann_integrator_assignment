# Friedmann Equation Integrator Assignment
Starter code for 8.942

This repo already contains the shell of a Friedmann Integrator, that actually
works for a matter dominated universe. Take a look at `integrators.py` to see
the code and explanations therein. Also look at `test_integrators.py`, which
in addition to testing the code, shows how it is meant to be used. Currently
the tests fail. Improve the code such that the tests pass.

You should be able to make the current tests pass without actually coding up
the numerical integral. Once the tests pass, write new tests based on analytic
solutions to the Friedmann Equations. Once you have new tests, write the
code that makes these tests pass.

This workflow - of writing tests first then writting the code that passes the
tests - is called test driven devellopment. For this assignment it might seem
like overkill, but for more involved projects it can save you days (or even
weeks) of frustration.
