#pragma once

// =============================================================================
// CC1e0pi Selection Configuration
// =============================================================================

// -----------------------------------------------------------------------------
// Reconstruction mode toggle
// Comment out USE_NUGRAPH to use legacy pandora track score + chi2 PID
// Uncomment to use NuGraph2 semantic classification
// -----------------------------------------------------------------------------
#define USE_NUGRAPH

// -----------------------------------------------------------------------------
// Energy thresholds (nominal values)
// These are used in the truth signal definition and reco cuts
// For threshold scans, use the ntuple output and vary in Python
// -----------------------------------------------------------------------------
const double THRESHOLD_E          = 0.2;    // electron KE threshold [GeV]
const double VISIBILTY_THRESHOLD_P  = 0.05;  // proton visibility threshold [GeV]
const double VISIBILTY_THRESHOLD_PI = 0.025; // pion visibility threshold [GeV]

// -----------------------------------------------------------------------------
// Flash match window
// -----------------------------------------------------------------------------
const double FM_DELTA_Z_MIN  = 0.0;   // [cm]
const double FM_DELTA_Z_MAX  = 100.0; // [cm]
const double FM_FLASH_T_MIN  = -1.0;  // [us]
const double FM_FLASH_T_MAX  = 11.0;  // [us]

// -----------------------------------------------------------------------------
// Shower cuts (nominal values)
// -----------------------------------------------------------------------------
const double SHOWER_ENERGY_MIN  = 0.200; // [GeV]
const double SHOWER_DEDX_MAX    = 3.5;   // [MeV/cm]
const double SHOWER_ANGLE_MAX   = 10.0;  // [deg]
const double SHOWER_CONVGAP_MAX = 5.0;   // [cm]

// -----------------------------------------------------------------------------
// Muon veto (legacy mode only)
// -----------------------------------------------------------------------------
const double MUON_CHI2_MUON_MAX   = 30.0;
const double MUON_CHI2_PROTON_MIN = 60.0;

// -----------------------------------------------------------------------------
// Proton/pion PID (legacy mode only)
// -----------------------------------------------------------------------------
const double PROTON_CHI2_MAX = 100.0;
const double PION_CHI2_MIN   = 100.0;

// -----------------------------------------------------------------------------
// Fiducial and containment margins [cm]
// -----------------------------------------------------------------------------
const double FV_MARGIN_X    = 25.0;
const double FV_MARGIN_Y    = 25.0;
const double FV_MARGIN_Z_UP = 30.0;
const double FV_MARGIN_Z_DN = 50.0;

const double CONT_MARGIN = 5.0;
