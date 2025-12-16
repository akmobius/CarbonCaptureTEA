# Carbon Capture Technoeconomic Analysis (System LCOE vs. Renewables Share (India, 2023 base))
Models the LCOE for India's power system % as the share of renewables increases, comparing a baseline coal-only % system with various CCS scenarios. It calculates and visualizes the system-wide LCOE (in ₹/kWh) as a function of the Renewables Share (ranging from 0% to 100% of the coal+RE generating capacity)
# Outputs
FIGURE 1: Compares the system LCOE across a base case and three main CCS cost scenarios against the current average system LCOE (the reference line). (1) No CCS vs Coal+CCS (low/mid/high), (2) Current mix-weighted system LCOE as reference, (3) RE leg includes integration + storage cost that grows with RE share. Compares the overall system LCOE for No CCS vs. CCS (low, mid, high cost) against the current mix-weighted system LCOE. The RE component includes integration and battery storage costs that scale with RE share. 

Crossover Point: Demonstrates the RE share at which the Coal + CCS system LCOE falls below the current system LCOE, achieving cost parity.


FIGURE 2: Focuses on the Mid-Cost CCS scenario and applies specific economic adjustments to model the impact of different policy and valuation choices. Avoided Externalities: Plots the LCOE reduction gained by monetizing the avoided damages from coal pollution, specifically: Health and Mortality costs, Crop-Loss damages. 

Carbon Credit Revenue: Plots the LCOE reduction achieved by incorporating a revenue stream from selling captured $\text{CO}_2$ (carbon credits) at low, mid, and high price points (₹/t$\text{CO}_2$).Co-Benefit Adjustment: Includes an adjustment for the capital cost overlap between CCS and required $\text{SO}_x$/$\text{NO}_x$ scrubbers, reducing the effective CCS capital cost.

Composite Scenario: Allows for the toggling and visualization of a custom scenario combining selected externalities and credit revenues.

All editable inputs and toggles are at the top.

# Edit History
Rev 1: 12/15/25 Ariel Mobius

Initial commit: 12/10/25 Emma Rutherford
