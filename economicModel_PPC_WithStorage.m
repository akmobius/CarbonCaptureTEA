% -------------------------------------------------------
% System LCOE vs. Renewables Share (India, 2023 base)
% 
% This script models the Levelized Cost of Energy (LCOE) for India's power system
% as the share of Renewables (RE) increases, comparing a baseline Coal-only
% system with various Carbon Capture and Storage (CCS) scenarios.
% 
% Original Script: 12/10/25 Emma Rutherford
% Rev 1: 12/15/25 Ariel Mobius
% 
% FIGURE 1:
%   • No CCS vs Coal+CCS (low/mid/high)
%   • Current mix-weighted system LCOE as reference
%   • RE leg includes integration + storage cost that grows with RE share
%   Compares the overall system LCOE for No CCS vs. CCS (low, mid, high cost)
%           against the current mix-weighted system LCOE. The RE component
%           includes integration and battery storage costs that scale with RE share.
% 
% FIGURE 2:
%   • Mid CCS with adjustments for:
%       - Health + mortality costs (combined)
%       - Crop-loss damages
%       - Net-of-SOx/NOx (remove 25% of CCS capex)
%       - Carbon credits (low / mid / high, in ₹/tCO2)
%       - Composite "ALL reductions" (configurable via toggles)
% 
% All editable inputs and toggles are at the top.
% -------------------------------------------------------
close all; clear; clc;

%% ================== PLOT TOGGLES (TURN LINES ON/OFF) ====================

% ---- Figure 1 line visibility ----
showFig1_noCCS      = true;   % No CCS baseline LCOE curve
showFig1_CCS_low    = true;   % Coal+CCS LCOE curve (low cost scenario)
showFig1_CCS_mid    = true;   % Coal+CCS LCOE curve (mid cost scenario)
showFig1_CCS_high   = true;   % Coal+CCS LCOE curve (high cost scenario)
showFig1_refLine    = true;   % Current average system LCOE (horizontal reference)

% ---- Figure 2 line visibility ----
showFig2_noCCS      = true;   % No CCS baseline LCOE curve (for context)
showFig2_CCS_mid    = true;   % Mid-Cost Coal+CCS LCOE curve (base for adjustments)
showFig2_healthMort = false;  % Mid CCS with health/mortality externality subtracted
showFig2_crops      = false;  % Mid CCS with crop-loss externality subtracted
showFig2_netSOx     = false;  % Mid CCS with CAPEX reduced for SOx/NOx co-benefit
showFig2_creditLow  = true;   % Mid CCS with low carbon credit revenue added
showFig2_creditMid  = true;   % Mid CCS with mid carbon credit revenue added
showFig2_creditHigh = true;   % Mid CCS with high carbon credit revenue added
showFig2_allRed     = false;  % Composite line including selected reductions/credits
showFig2_refLine    = true;   % Current average system LCOE (horizontal reference)

%% ====== COMPOSITE ("ALL REDUCTIONS") CONTENT TOGGLES ====================
% These control what is subtracted / credited in the composite mid-CCS line

% These control which economic factors are included in the 'showFig2_allRed' composite line.
use_netSOxNOx_inComposite    = false;   % Factor in CAPEX reduction due to SOx/NOx co-benefit
use_healthMort_inComposite   = true;    % Factor in avoided health + mortality costs (externality)
use_crops_inComposite        = true;    % Factor in avoided crop-loss damages (externality)
use_carbonCredit_inComposite = false;   % Factor in carbon credit revenue

% Choose which carbon credit level the composite uses: 'low' | 'mid' | 'high'
composite_credit_level = 'mid';

%% ================== USER INPUTS (EDIT THESE ONLY) =======================

% --- Observed current average system PPC (₹/kWh) ---
% Set this to your empirical current system PPC (from govt / DISCOM data)
% Subsidy value calculated based on (Observed tariff) - (Calculated plant-gate LCOE)
subsidy_Rs = 7.3-4.2;     % <<< EDIT THIS VALUE Subsidized cost to consumer (₹/kWh)
subsidize = false;        % Logical: Apply subsidy to all calculated LCOE curves?
% --- Currency conversion ---
usdToRupee = 89.96;   % Rs per USD

% --- Technology LCOEs (USD/kWh) ---
% Order: [Coal, Oil, Gas, Nuclear, Hydro, Solar, Wind, Bio, Waste]
LCOE_source_USD = [ ...
    0.0879;   % Coal (IRENA 2023, India)
    0.3568;     % Oil  (IRENA 2023, World)
    0.1157;     % Gas 
    5.8/usdToRupee;    % Nuclear
    0.047;    % Hydropower (large) (IRENA 2023, India)
    0.048;    % Solar PV utility-scale (IRENA 2023, India)
    0.046;    % Onshore wind (IRENA 2023, India)
    0.062;    % Bio (IRENA 2023, India)
    0.062     % Waste
];
% LCOE_source_USD = [ ...
%     4.34/usdToRupee;   % Coal (IRENA 2023, India)
%     5.35/usdToRupee;     % Oil  (IRENA 2023, World)
%     5.35/usdToRupee;     % Gas 
%     3.6/usdToRupee;    % Nuclear
%     2.45/usdToRupee;    % Hydropower (large) (IRENA 2023, India)
%     3.97/usdToRupee;    % Solar PV utility-scale (IRENA 2023, India)
%     4.04/usdToRupee;    % Onshore wind (IRENA 2023, India)
%     5.65/usdToRupee;    % Bio
%     0.062     % Waste
% ];

% Index map for readability (do not edit)
idxCoal   = 1;
idxOil    = 2;
idxGas    = 3;
idxNucl   = 4;
idxHydro  = 5;
idxSolar  = 6;
idxWind   = 7;
idxBio    = 8;
idxWaste  = 9;

% --- CURRENT GENERATION MIX — INDIA 2023 (GWh/year) ---
% Electricity generation sources, India, 2023 (GWh)
% Coal        1478683
% Oil         4153
% Natural gas 58764
% Biofuels    37373
% Waste       3981
% Nuclear     47937
% Hydropower  143890
% Solar PV    118528
% Wind        93465
%
% Map to model order [Coal, Oil, Gas, Nuclear, Hydro, Solar, Wind, Bio, Waste]
gen_TWh = [ ...
    1478683;   % Coal
    4153;      % Oil
    58764;     % Natural gas
    47937;     % Nuclear
    143890;    % Hydropower
    118528;    % Solar PV
    93465;     % Wind
    37373;     % Biofuels
    3981       % Waste
] / 1000;   % convert GWh -> TWh

% --- Renewables mix for future (solar + wind only) ---
solar_share_RE = 0.75;    % fraction of RE generation that is solar (0–1)

% --- Integration cost (system-level) for VRE (USD/kWh per unit RE share) ---
% RE share r ranges from 0 to 1. Set =0 to ignore integration costs.
% Accounts for costs like grid upgrades due to variable RE (VRE).
a_integ_USD = 0.000;      % Linear cost increase (USD/kWh per unit of RE share r)

% --- Battery storage economics (for firming RE) ---
storage_capex_USD_per_kWh      = 152;    % $/kWh of battery capacity
storage_lifetime_years         = 15;     % years
storage_cycles_per_year        = 365;    % full cycles per year
storage_roundtrip_efficiency   = 0.95;    % fraction
storage_opex_USD_per_kWh_thru  = 0.002;  % variable O&M per kWh throughput

% Storage requirement vs RE share r
r_storage_start       = 0.3;   % below this RE share, ~0 storage
r_storage_full        = 1;   % at or above this, storage requirement = max
storage_share_RE_max  = 0.5;   % at high RE, fraction of RE energy through storage

% --- CCS cost model based on NITI Aayog CCUS (2022), coal-based power ---
% Table 6-4: Coal-based power (800 MW, 5 mtpa CCU capacity)
% Capital charges: 700–1,000 ₹/tCO2
% Cash (operating) costs: 2,100–2,500 ₹/tCO2
% Total capture cost: 2,800–3,500 ₹/tCO2
% (Excludes transport, sequestration, monitoring, which add ~US$10–15/tCO2)

% Capture costs are defined in ₹/tCO2 for low, mid, and high scenarios.
% Capital cost (Rs/tCO2)
CCS_capex_low_Rs_per_t  = 700;
CCS_capex_mid_Rs_per_t  = (700 + 1000)/2;   % 850 (midpoint assumption)
CCS_capex_high_Rs_per_t = 1000;

% Cash / O&M cost (Rs/tCO2)
CCS_opex_low_Rs_per_t   = 2100;
CCS_opex_mid_Rs_per_t   = (2100 + 2500)/2;  % 2300
CCS_opex_high_Rs_per_t  = 2500;

% Convert capture costs to USD/tCO2 using usdToRupee
CCS_capex_low_USD_per_t  = CCS_capex_low_Rs_per_t  / usdToRupee;
CCS_capex_mid_USD_per_t  = CCS_capex_mid_Rs_per_t  / usdToRupee;
CCS_capex_high_USD_per_t = CCS_capex_high_Rs_per_t / usdToRupee;

CCS_opex_low_USD_per_t   = CCS_opex_low_Rs_per_t   / usdToRupee;
CCS_opex_mid_USD_per_t   = CCS_opex_mid_Rs_per_t   / usdToRupee;
CCS_opex_high_USD_per_t  = CCS_opex_high_Rs_per_t  / usdToRupee;

% Transport & storage (and monitoring) costs, from NITI Aayog:
% Additional US$10–15/tCO2 on top of capture cash costs
CCS_TS_low_USD_per_t  = 5;
CCS_TS_mid_USD_per_t  = 15;
CCS_TS_high_USD_per_t = 45;

% --- CO2 emissions intensity of coal generation ---
CO2_intensity_t_per_MWh = 0.95;       % tCO2/MWh
CO2_intensity_t_per_kWh = CO2_intensity_t_per_MWh / 1000;  % tCO2/kWh

% --- Fraction of CO2 captured by CCS ---
capture_rate = 0.95;     % e.g., 90% capture

% --- Externalities & damages as TOTAL annual costs (USD/year, attributed to coal) ---
healthMort_cost_total_USD_per_year = 4.5e9;   % combined health + mortality cost
crop_loss_total_USD_per_year       = 2.0e9;   % crop-loss damages

% --- Carbon credit / CO2 price assumptions (for CCS cases, in ₹/tCO2) ---
carbon_price_low_Rs_per_t  = 250;    % low voluntary credit
carbon_price_mid_Rs_per_t  = 800;    % mid compliance credit
carbon_price_high_Rs_per_t = 1500;   % high credit case

% --- Fraction of CCS CAPEX that overlaps with SOx/NOx scrubber requirement ---
soxNOx_capex_fraction = 0.15;        % X% of CCS CAPEX assumed needed anyway

%% ================== DERIVED PARAMETERS (DO NOT EDIT) ====================

% Coal, solar, wind LCOEs in USD
LCOE_coal_noCCS_USD = LCOE_source_USD(idxCoal);
LCOE_solar_USD      = LCOE_source_USD(idxSolar);
LCOE_wind_USD       = LCOE_source_USD(idxWind); % Weighted RE LCOE

% Weighted RE LCOE in USD
LCOE_RE_USD = solar_share_RE * LCOE_solar_USD + ...
              (1 - solar_share_RE) * LCOE_wind_USD;

% Current mix-weighted system LCOE (USD & Rs/kWh)
shares       = gen_TWh ./ sum(gen_TWh);
LCOE_ref_USD = sum(shares .* LCOE_source_USD); %Current mix average LCOE
LCOE_ref_Rs  = LCOE_ref_USD * usdToRupee;

if subsidize
    LCOE_ref_Rs = LCOE_ref_Rs - subsidy_Rs; %Adjust for susbidy if toggled
end
fprintf('Current mix-weighted plant-gate LCOE ≈ $%.3f/kWh (≈ ₹%.2f/kWh)\n', ...
    LCOE_ref_USD, LCOE_ref_Rs);

% Battery storage LCOE (USD/kWh delivered)
lifetime_throughput_kWh = storage_cycles_per_year * storage_lifetime_years * storage_roundtrip_efficiency;
LCOE_storage_USD = (storage_capex_USD_per_kWh / lifetime_throughput_kWh) + ...
                   storage_opex_USD_per_kWh_thru;
LCOE_storage_Rs  = LCOE_storage_USD * usdToRupee;

% --- CCS Uplift Calculation (The increase in LCOE due to CCS) ---
% Uplift (dLCOE) is calculated by converting the per-tonne CCS costs (CAPEX, OPEX, T&S)
% into a per-kWh cost using CO2 intensity and capture rate.

% Capture CAPEX contribution
dLCOE_CCS_low_capex_USD  = CO2_intensity_t_per_kWh * capture_rate * CCS_capex_low_USD_per_t;
dLCOE_CCS_mid_capex_USD  = CO2_intensity_t_per_kWh * capture_rate * CCS_capex_mid_USD_per_t;
dLCOE_CCS_high_capex_USD = CO2_intensity_t_per_kWh * capture_rate * CCS_capex_high_USD_per_t;

% Capture OPEX (cash cost) contribution
dLCOE_CCS_low_opex_USD   = CO2_intensity_t_per_kWh * capture_rate * CCS_opex_low_USD_per_t;
dLCOE_CCS_mid_opex_USD   = CO2_intensity_t_per_kWh * capture_rate * CCS_opex_mid_USD_per_t;
dLCOE_CCS_high_opex_USD  = CO2_intensity_t_per_kWh * capture_rate * CCS_opex_high_USD_per_t;

% Transport & storage contribution
dLCOE_CCS_low_TS_USD     = CO2_intensity_t_per_kWh * capture_rate * CCS_TS_low_USD_per_t;
dLCOE_CCS_mid_TS_USD     = CO2_intensity_t_per_kWh * capture_rate * CCS_TS_mid_USD_per_t;
dLCOE_CCS_high_TS_USD    = CO2_intensity_t_per_kWh * capture_rate * CCS_TS_high_USD_per_t;

% Total CCS uplift (capture + T&S)
dLCOE_CCS_low_USD  = dLCOE_CCS_low_capex_USD  + dLCOE_CCS_low_opex_USD  + dLCOE_CCS_low_TS_USD;
dLCOE_CCS_mid_USD  = dLCOE_CCS_mid_capex_USD  + dLCOE_CCS_mid_opex_USD  + dLCOE_CCS_mid_TS_USD;
dLCOE_CCS_high_USD = dLCOE_CCS_high_capex_USD + dLCOE_CCS_high_opex_USD + dLCOE_CCS_high_TS_USD;

% Convert costs to Rs/kWh
LCOE_coal_noCCS = LCOE_coal_noCCS_USD * usdToRupee;
LCOE_RE         = LCOE_RE_USD         * usdToRupee;
a_integ         = a_integ_USD         * usdToRupee;

dLCOE_CCS_low   = dLCOE_CCS_low_USD   * usdToRupee;
dLCOE_CCS_mid   = dLCOE_CCS_mid_USD   * usdToRupee;
dLCOE_CCS_high  = dLCOE_CCS_high_USD  * usdToRupee;

% Final LCOE of Coal+CCS (Rs/kWh)
LCOE_coal_CCS_low  = LCOE_coal_noCCS + dLCOE_CCS_low;
LCOE_coal_CCS_mid  = LCOE_coal_noCCS + dLCOE_CCS_mid;
LCOE_coal_CCS_high = LCOE_coal_noCCS + dLCOE_CCS_high;

% Mid CCS net-of-SOx/NOx (subtract soxNOx_capex_fraction of CAPEX in Rs)
dLCOE_CCS_mid_capex_Rs  = dLCOE_CCS_mid_capex_USD  * usdToRupee;
dLCOE_CCS_mid_opex_Rs   = dLCOE_CCS_mid_opex_USD   * usdToRupee;
dLCOE_CCS_mid_capex_net = (1 - soxNOx_capex_fraction) * dLCOE_CCS_mid_capex_Rs;
LCOE_coal_CCS_mid_netSOx = LCOE_coal_noCCS + dLCOE_CCS_mid_opex_Rs + dLCOE_CCS_mid_capex_net;

% --- Externalities Conversion (Rs/kWh) ---
% Convert total annual damage costs (USD/year) into a per-kWh cost based on baseline coal generation.
gen_TWh_coal_baseline = gen_TWh(idxCoal);            % coal TWh/year
coal_kWh_baseline     = gen_TWh_coal_baseline * 1e9; % 1 TWh = 1e9 kWh

healthMort_cost_USD_per_kWh = healthMort_cost_total_USD_per_year / coal_kWh_baseline;
crop_loss_cost_USD_per_kWh  = crop_loss_total_USD_per_year       / coal_kWh_baseline;

healthMort_cost_Rs_per_kWh  = healthMort_cost_USD_per_kWh * usdToRupee;
crop_loss_cost_Rs_per_kWh   = crop_loss_cost_USD_per_kWh  * usdToRupee;

% --- Carbon Credits Conversion (Rs/kWh) ---
% Convert per-tonne revenue into a per-kWh revenue stream for CCS coal.
% Carbon credits (Rs/kWh) for coal+CCS generation (low/mid/high)
credit_low_Rs_per_kWh  = CO2_intensity_t_per_kWh * capture_rate * carbon_price_low_Rs_per_t;
credit_mid_Rs_per_kWh  = CO2_intensity_t_per_kWh * capture_rate * carbon_price_mid_Rs_per_t;
credit_high_Rs_per_kWh = CO2_intensity_t_per_kWh * capture_rate * carbon_price_high_Rs_per_t;

%% ================== STORAGE REQUIREMENT VS RE SHARE =====================
% Calculates the storage cost adder for the RE LCOE.
% Renewables share axis (0–100% of coal+RE system)
r = linspace(0,1,101);

% Logic to determine the fraction of RE energy that must pass through storage (firming cost)
% Fraction of RE energy that passes through storage as function of r
storage_share_RE = zeros(size(r));
idx_mid = (r > r_storage_start) & (r < r_storage_full);
% Linear ramp-up of storage requirement between r_storage_start (30%) and r_storage_full (100%)
storage_share_RE(idx_mid) = storage_share_RE_max * ...
    (r(idx_mid) - r_storage_start) ./ (r_storage_full - r_storage_start);
storage_share_RE(r >= r_storage_full) = storage_share_RE_max;

% Storage cost adder on RE leg (Rs/kWh of RE)
storageAdder = storage_share_RE * LCOE_storage_Rs;

% Integration cost adder
integ = a_integ * r;

%% ================== SYSTEM LCOE CALCULATIONS ===========================
% The core LCOE calculation: System LCOE = r*(RE LCOE + Adders) + (1 - r)*Coal LCOE
% where r is the renewables share.

% No-CCS baseline
LCOE_sys_noCCS    = r .* (LCOE_RE + integ + storageAdder) + (1 - r) .* LCOE_coal_noCCS;

% CCS cases (low / mid / high)
LCOE_sys_CCS_low  = r .* (LCOE_RE + integ + storageAdder) + (1 - r) .* LCOE_coal_CCS_low;
LCOE_sys_CCS_mid  = r .* (LCOE_RE + integ + storageAdder) + (1 - r) .* LCOE_coal_CCS_mid;
LCOE_sys_CCS_high = r .* (LCOE_RE + integ + storageAdder) + (1 - r) .* LCOE_coal_CCS_high;

% Mid CCS net-of-SOx/NOx (applies only to the coal portion (1-r))
LCOE_sys_CCS_mid_netSOx = r .* (LCOE_RE + integ + storageAdder) + (1 - r) .* LCOE_coal_CCS_mid_netSOx;

% Externalities Reductions (applied as a subtraction from the mid-CCS LCOE)
% Note: The health/crop costs are *avoided* only in the coal share (1-r).
LCOE_sys_mid_healthMort = LCOE_sys_CCS_mid - (1 - r) .* healthMort_cost_Rs_per_kWh;
LCOE_sys_mid_crops = LCOE_sys_CCS_mid - (1 - r) .* crop_loss_cost_Rs_per_kWh;

% Carbon Credit Reductions (applied only to the coal portion (1-r))
% Revenue is subtracted from the coal+CCS LCOE, with a floor of zero.
LCOE_coal_CCS_mid_creditLow = max(LCOE_coal_CCS_mid - credit_low_Rs_per_kWh, 0);
LCOE_sys_mid_creditLow = r .* (LCOE_RE + integ + storageAdder) + ...
                         (1 - r) .* LCOE_coal_CCS_mid_creditLow;

% Carbon credit (mid) applied to coal+CCS leg
LCOE_coal_CCS_mid_creditMid = max(LCOE_coal_CCS_mid - credit_mid_Rs_per_kWh, 0);
LCOE_sys_mid_creditMid = r .* (LCOE_RE + integ + storageAdder) + ...
                         (1 - r) .* LCOE_coal_CCS_mid_creditMid;

% Carbon credit (high) applied to coal+CCS leg
LCOE_coal_CCS_mid_creditHigh = max(LCOE_coal_CCS_mid - credit_high_Rs_per_kWh, 0);
LCOE_sys_mid_creditHigh = r .* (LCOE_RE + integ + storageAdder) + ...
                          (1 - r) .* LCOE_coal_CCS_mid_creditHigh;

% Subsidy adjustment (applied if the 'subsidize' toggle is true)
if subsidize
    LCOE_sys_noCCS = LCOE_sys_noCCS - subsidy_Rs;
    LCOE_sys_CCS_low = LCOE_sys_CCS_low - subsidy_Rs;
    LCOE_sys_CCS_mid = LCOE_sys_CCS_mid - subsidy_Rs;
    LCOE_sys_CCS_high = LCOE_sys_CCS_high - subsidy_Rs;
    LCOE_sys_CCS_mid_netSOx = LCOE_sys_CCS_mid_netSOx - subsidy_Rs;
    LCOE_sys_mid_healthMort = LCOE_sys_mid_healthMort - subsidy_Rs;
    LCOE_sys_mid_crops = LCOE_sys_mid_crops - subsidy_Rs;
    LCOE_sys_mid_creditLow = LCOE_sys_mid_creditLow - subsidy_Rs;
    LCOE_sys_mid_creditMid = LCOE_sys_mid_creditMid - subsidy_Rs;
    LCOE_sys_mid_creditHigh = LCOE_sys_mid_creditHigh - subsidy_Rs;
end
%% ===== BUILD COMPOSITE MID-CCS COAL LCOE FOR "ALL REDUCTIONS" =========
% Logic to build the 'ALL REDUCTIONS' composite line based on the user's composite toggles.

% Start from plain mid-CCS coal (Rs/kWh)
LCOE_coal_CCS_mid_composite = LCOE_coal_CCS_mid;
externalityAdder = zeros(size(r));             % Accumulator for externality benefits

% Apply Net SOx/NOx adjustment if toggled
if use_netSOxNOx_inComposite
    LCOE_coal_CCS_mid_composite = LCOE_coal_CCS_mid_netSOx;
end

if use_carbonCredit_inComposite
    % Selects the correct credit amount based on 'composite_credit_level'
    switch lower(composite_credit_level)
        case 'low'
            credit_to_use = credit_low_Rs_per_kWh;
        case 'mid'
            credit_to_use = credit_mid_Rs_per_kWh;
        case 'high'
            credit_to_use = credit_high_Rs_per_kWh;
        otherwise
            credit_to_use = 0;
    end
    % Subtracts credit revenue (with a floor of 0)
    LCOE_coal_CCS_mid_composite = max( ...
        LCOE_coal_CCS_mid_composite - credit_to_use, 0);
end

compParts = {};  % for legend text

% Accumulate externality benefits (subtracted at the end)
if use_netSOxNOx_inComposite
    compParts{end+1} = '- SOx/NOx';
end
if use_carbonCredit_inComposite
    switch lower(composite_credit_level)
        case 'low'
            compParts{end+1} = 'low credit';
        case 'mid'
            compParts{end+1} = 'mid credit';
        case 'high'
            compParts{end+1} = 'high credit';
    end
end
if use_healthMort_inComposite
    externalityAdder = externalityAdder + healthMort_cost_Rs_per_kWh;
    compParts{end+1} = 'healthcare, mortality,';
end
if use_crops_inComposite
    externalityAdder = externalityAdder + crop_loss_cost_Rs_per_kWh;
    compParts{end+1} = 'crops';
end

% Final calculation for the composite LCOE curve
LCOE_sys_mid_allRed = r .* (LCOE_RE + integ + storageAdder) + ...
                      (1 - r) .* LCOE_coal_CCS_mid_composite - ...
                      (1 - r) .* externalityAdder;

% Logic to create a detailed label for the composite line
if isempty(compParts)
    compositeLabel = 'Mid CCS composite (no reductions)';
else
    if numel(compParts) == 1
        compositeLabel = ['Mid CCS composite (' compParts{1} ')'];
    else
        compositeLabel = ['Mid CCS - (' strjoin(compParts(1:end-1), ', ') ...
                          ' ' compParts{end} ')'];
    end
end

%% ================== CROSSOVER POINTS VS CURRENT LCOE ====================
% Calculates the RE share (r) at which the system LCOE curve intersects the current LCOE reference line.
% This is the "cost parity" point.

[~, idxLow]       = min(abs(LCOE_sys_CCS_low        - LCOE_ref_Rs));
[~, idxMid]       = min(abs(LCOE_sys_CCS_mid        - LCOE_ref_Rs));
[~, idxHigh]      = min(abs(LCOE_sys_CCS_high       - LCOE_ref_Rs));
[~, idxMidNetSOx] = min(abs(LCOE_sys_CCS_mid_netSOx - LCOE_ref_Rs));

r_cross_low       = r(idxLow);
r_cross_mid       = r(idxMid);
r_cross_high      = r(idxHigh);
r_cross_midNetSOx = r(idxMidNetSOx);

% Print crossover points to the command window
fprintf('Low  CCS crossover vs current LCOE at ~%.1f%% RE share.\n',  r_cross_low*100);
fprintf('Mid  CCS crossover vs current LCOE at ~%.1f%% RE share.\n',  r_cross_mid*100);
fprintf('High CCS crossover vs current LCOE at ~%.1f%% RE share.\n', r_cross_high*100);
fprintf('Mid CCS (net SOx/NOx) crossover at ~%.1f%% RE share.\n',     r_cross_midNetSOx*100);

%% ================== FIGURE 1 — LCOE-ONLY COMPARISON ====================

figure(1); clf; hold on; box on; grid on;

% Shaded CCS band (low–high) to visualize the uncertainty range
xBand = [r*100, fliplr(r*100)];
yBand = [LCOE_sys_CCS_low, fliplr(LCOE_sys_CCS_high)];
fill(xBand, yBand, [0.8 0.9 1.0], ...
    'EdgeColor','none', 'FaceAlpha',0.5);

h1 = []; leg1 = {};

% Plotting individual CCS and No-CCS lines based on toggles
if showFig1_noCCS
    p_noCCS = plot(r*100, LCOE_sys_noCCS, ':', 'LineWidth', 3);
    h1(end+1)   = p_noCCS;
    leg1{end+1} = 'No CCS (coal + RE)';
end

if showFig1_CCS_low
    p_low = plot(r*100, LCOE_sys_CCS_low, '-', 'LineWidth', 3);
    h1(end+1)   = p_low;
    leg1{end+1} = 'Coal + CCS (low cost)';
end

if showFig1_CCS_mid
    p_mid = plot(r*100, LCOE_sys_CCS_mid, '-', 'LineWidth', 3);
    h1(end+1)   = p_mid;
    leg1{end+1} = 'Coal + CCS (mid cost)';
end

if showFig1_CCS_high
    p_high = plot(r*100, LCOE_sys_CCS_high, '-', 'LineWidth', 3);
    h1(end+1)   = p_high;
    leg1{end+1} = 'Coal + CCS (high cost)';
end

% Plotting the Current System LCOE reference line
if showFig1_refLine
    refLabel = sprintf('Current System LCOE (₹%.1f/kWh)', LCOE_ref_Rs);
    yline(LCOE_ref_Rs, 'k--', 'LineWidth', 2, ...
        'Label', refLabel, ...
        'LabelHorizontalAlignment','left', ...
        'LabelVerticalAlignment','bottom', ...
        'LabelOrientation','aligned', ...
        'FontSize', 14);
end

% Set labels, title, limits, and legend for Figure 1
xlabel('Renewables (with storage) share of generation in system (%)');
%ylabel('System Power Purchase Cost (PPC) (₹/kWh)');
ylabel('Levelized Cost of Energy (₹/kWh)');
title('System LCOE vs Renewables Share — Coal Plants with and without CCS (India)');

if ~isempty(h1)
    legend(h1, leg1, 'Location','best');
end
ylim([2,15])
grid on;
if exist('improvePlot','file'), improvePlot(); end

%% = FIGURE 2 — EXTERNALITIES, NET SOx/NOx, CREDITS, COMPOSITE =========

figure(2); clf; hold on; box on; grid on;

h2 = []; leg2 = {};

% Plotting lines based on toggles, focusing on adjustments to the Mid-CCS line.
if showFig2_noCCS
    p2_noCCS = plot(r*100, LCOE_sys_noCCS, ':', 'LineWidth', 3);
    h2(end+1)   = p2_noCCS;
    leg2{end+1} = 'No CCS baseline';
end

if showFig2_CCS_mid
    p2_mid = plot(r*100, LCOE_sys_CCS_mid, '-', 'LineWidth', 3);
    h2(end+1)   = p2_mid;
    leg2{end+1} = 'Coal + CCS (mid cost)';
end

if showFig2_healthMort
    p2_healthMort = plot(r*100, LCOE_sys_mid_healthMort, '--', 'LineWidth', 3);
    h2(end+1)   = p2_healthMort;
    leg2{end+1} = 'Mid CCS - health + mortality';
end

if showFig2_crops
    p2_crops = plot(r*100, LCOE_sys_mid_crops, '-.', 'LineWidth', 3);
    h2(end+1)   = p2_crops;
    leg2{end+1} = 'Mid CCS - crop-loss damages';
end

if showFig2_netSOx
    p2_netSOx = plot(r*100, LCOE_sys_CCS_mid_netSOx, '-', 'LineWidth', 3);
    h2(end+1)   = p2_netSOx;
    leg2{end+1} = 'Mid CCS - SOx/NOx (15% CC CAPEX removed)';
end

if showFig2_creditLow
    p2_creditLow = plot(r*100, LCOE_sys_mid_creditLow, '-', 'LineWidth', 3);
    h2(end+1)   = p2_creditLow;
    leg2{end+1} = 'Mid CCS - low carbon credit (₹/250)';
end

if showFig2_creditMid
    p2_creditMid = plot(r*100, LCOE_sys_mid_creditMid, '-', 'LineWidth', 3);
    h2(end+1)   = p2_creditMid;
    leg2{end+1} = 'Mid CCS - mid carbon credit (₹/800)';
end

if showFig2_creditHigh
    p2_creditHigh = plot(r*100, LCOE_sys_mid_creditHigh, '-', 'LineWidth', 3);
    h2(end+1)   = p2_creditHigh;
    leg2{end+1} = 'Mid CCS - high carbon credit (₹/1500)';
end

if showFig2_allRed
    p2_all = plot(r*100, LCOE_sys_mid_allRed, '-', 'LineWidth', 3);
    h2(end+1)   = p2_all;
    leg2{end+1} = compositeLabel;
end

% Plotting the Current System LCOE reference line
if showFig2_refLine
    refLabel2 = sprintf('Current avg system LCOE (₹%.1f/kWh)', LCOE_ref_Rs);
    yline(LCOE_ref_Rs, 'k--', 'LineWidth', 2, ...
        'Label', refLabel2, ...
        'LabelHorizontalAlignment','left', ...
        'LabelVerticalAlignment','bottom', ...
        'LabelOrientation','aligned', ...
        'FontSize', 14);
end

% Set labels, title, limits, and legend for Figure 2
xlabel('Renewables (with storage) share of generation in system (%)');
ylabel('System LCOE (₹/kWh)');
title('System LCOE with CCS — Carbon Credits');

if ~isempty(h2)
    legend(h2, leg2, 'Location','best');
end
ylim([5,13]);
grid on;
if exist('improvePlot','file'), improvePlot(); end


