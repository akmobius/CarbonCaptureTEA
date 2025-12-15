# Carbon Capture Technoeconomic Analysis (System LCOE vs. Renewables Share (India, 2023 base))
Models the LCOE for India's power system % as the share of renewables increases, comparing a baseline coal-only % system with various CCS scenarios.

# Edit History (pre-github)
Initial commit: 12/10/25 Emma Rutherford, Rev 1: 12/15/25 Ariel Mobius

# Outputs
FIGURE 1:
  • No CCS vs Coal+CCS (low/mid/high)
  • Current mix-weighted system LCOE as reference
  • RE leg includes integration + storage cost that grows with RE share
  Compares the overall system LCOE for No CCS vs. CCS (low, mid, high cost) against the current mix-weighted system LCOE. The RE component includes integration and battery storage costs that scale with RE share


FIGURE 2:
  • Mid CCS with adjustments for:
      - Health + mortality costs (combined),
      - Crop-loss damages,
      - Net-of-SOx/NOx (remove 25% of CCS capex),
      - Carbon credits (low / mid / high, in ₹/tCO2),
      - Composite "ALL reductions" (configurable via toggles)

All editable inputs and toggles are at the top.
