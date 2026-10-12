"""Aerosonde UAV parameter set -- the single source of truth for all numbers.

Source: Beard & McLain, *Small Unmanned Aircraft: Theory and Practice*
(Princeton, 2012), Appendix E (Aerosonde parameters).

Rules for this file:
  * Every value carries SI units and a comment citing where it came from.
  * Values not taken directly from the book are marked ``# NOT FROM BOOK``
    with the reason, so they can be audited.

Parameters are added phase by phase (mass/inertia in Phase 1, aerodynamics,
propulsion and actuators in Phase 2). Intentionally empty in Phase 0.
"""
