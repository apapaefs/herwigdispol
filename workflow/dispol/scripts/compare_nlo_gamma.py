#!/usr/bin/env python3
"""
Compare Herwig gamma-only NLO POWHEG cross sections with POLDIS predictions.

Usage:
    python compare_nlo_gamma.py

Reads the .out files from the POSNLO-GAMMA and NEGNLO-GAMMA runs and extracts
the total cross sections. Computes Total = PositiveNLO - NegativeNLO (the true
NLO cross section), then compares PP-PM and K_pol with POLDIS.

Key formula: sigma_NLO = sigma_PositiveNLO - sigma_NegativeNLO
(because NLOWeight returns max(0, wgt) for Pos and max(0, -wgt) for Neg)
"""
import re
import sys
import os
import math

def extract_xsec(filename):
    """Extract total cross section from Herwig .out file.
    Returns (value_nb, error_nb) in nanobarns, or (None, None) if not found."""
    if not os.path.exists(filename):
        return None, None
    with open(filename) as f:
        text = f.read()
    # Herwig format (nb units from header, not on data line):
    # "Total (from generated events):  1e+06  1e+06  6.346(4)e+00"
    # The value is X.XXX(Y)e+ZZ where (Y) is the error on the last digit(s)
    m = re.search(r'Total\s+\(from generated events\).*?([\d.]+)\((\d+)\)e([+-]?\d+)', text)
    if m:
        mantissa = float(m.group(1))
        err_digits = int(m.group(2))
        exponent = int(m.group(3))
        val = mantissa * 10**exponent
        # Error: digits in () = uncertainty on last digits of mantissa
        # e.g., 6.346(4)e+00 means 6.346 +- 0.004
        n_decimals = len(m.group(1).split('.')[-1]) if '.' in m.group(1) else 0
        err = err_digits * 10**(exponent - n_decimals)
        return val, err
    return None, None


def pb(val_nb, err_nb=None):
    """Convert nb to pb."""
    if err_nb is not None:
        return val_nb * 1e3, err_nb * 1e3
    return val_nb * 1e3


def print_comparison(label, herwig_val, herwig_err, poldis_val, poldis_err, unit="pb"):
    """Print a comparison line with ratio and significance."""
    ratio = herwig_val / poldis_val
    delta_pct = (ratio - 1.0) * 100
    diff = herwig_val - poldis_val
    diff_err = math.sqrt(herwig_err**2 + poldis_err**2)
    nsigma = abs(diff) / diff_err if diff_err > 0 else float('inf')
    print(f"    {label}")
    print(f"      Herwig = {herwig_val:.1f} +- {herwig_err:.1f} {unit}")
    print(f"      POLDIS = {poldis_val:.1f} +- {poldis_err:.1f} {unit}")
    print(f"      Ratio  = {ratio:.4f}  ({delta_pct:+.2f}%,  {nsigma:.1f}sigma)")
    return ratio, delta_pct, nsigma


def main():
    basedir = os.path.dirname(os.path.abspath(__file__))

    # ================================================================
    # POLDIS reference values (gamma-only, 20M points)
    # ================================================================
    poldis = {
        'lo_unpol':    6600.787,   # pb
        'nlo_unpol':   6197.993,   # pb
        'nnlo_unpol':  5923.464,   # pb
        'lo_pol':        52.185,   # pb (sigma_pol = (PP-PM)/2)
        'nlo_pol':       53.235,   # pb
        'nnlo_pol':      53.523,   # pb
        # Errors
        'lo_unpol_e':    1.483,
        'nlo_unpol_e':   1.856,
        'nnlo_unpol_e':  2.723,
        'lo_pol_e':      0.011,
        'nlo_pol_e':     0.013,
        'nnlo_pol_e':    0.017,
    }
    poldis['k_unpol'] = poldis['nlo_unpol'] / poldis['lo_unpol']
    poldis['k_pol']   = poldis['nlo_pol']   / poldis['lo_pol']

    print("=" * 72)
    print("  Herwig vs POLDIS: gamma-only NLO polarized DIS")
    print("  e-(18 GeV) + p+(275 GeV), Q2=[25,2500], y=[0.2,0.6]")
    print("  PDF: PDF4LHC15_nnlo_100_pdfas, alpha_EM=1/137.036")
    print("  POLDIS: 20M points  |  Herwig: 100M events")
    print("=" * 72)
    print()
    print("  POLDIS reference (gamma-only, 20M):")
    print(f"    Unpol:  LO = {poldis['lo_unpol']:.1f} +- {poldis['lo_unpol_e']:.1f}")
    print(f"            NLO = {poldis['nlo_unpol']:.1f} +- {poldis['nlo_unpol_e']:.1f},  K = {poldis['k_unpol']:.4f}")
    print(f"            NNLO = {poldis['nnlo_unpol']:.1f} +- {poldis['nnlo_unpol_e']:.1f}")
    print(f"    Pol:    LO = {poldis['lo_pol']:.3f} +- {poldis['lo_pol_e']:.3f}")
    print(f"            NLO = {poldis['nlo_pol']:.3f} +- {poldis['nlo_pol_e']:.3f},  K = {poldis['k_pol']:.4f}")
    print(f"            NNLO = {poldis['nnlo_pol']:.3f} +- {poldis['nnlo_pol_e']:.3f}")
    print(f"    PP-PM:  LO = {2*poldis['lo_pol']:.2f} pb,  NLO = {2*poldis['nlo_pol']:.2f} pb")
    print()

    # ================================================================
    # Read Herwig output files
    # ================================================================
    pos_files = {
        'PP': 'DIS-POL-POWHEG_PP-POSNLO-GAMMA.out',
        'PM': 'DIS-POL-POWHEG_PM-POSNLO-GAMMA.out',
        '00': 'DIS-POL-POWHEG_00-POSNLO-GAMMA.out',
    }
    neg_files = {
        'PP': 'DIS-POL-POWHEG_PP-NEGNLO-GAMMA.out',
        'PM': 'DIS-POL-POWHEG_PM-NEGNLO-GAMMA.out',
        '00': 'DIS-POL-POWHEG_00-NEGNLO-GAMMA.out',
    }

    pos = {}  # PositiveNLO cross sections in nb
    neg = {}  # NegativeNLO cross sections in nb
    for tag in ['PP', 'PM', '00']:
        pos[tag] = extract_xsec(os.path.join(basedir, pos_files[tag]))
        neg[tag] = extract_xsec(os.path.join(basedir, neg_files[tag]))

    has_pos = {t: pos[t][0] is not None for t in ['PP', 'PM', '00']}
    has_neg = {t: neg[t][0] is not None for t in ['PP', 'PM', '00']}

    # ================================================================
    # Print raw cross sections
    # ================================================================
    print("-" * 72)
    print("  Raw Herwig cross sections (pb):")
    print("-" * 72)
    print(f"  {'':>4s}  {'PositiveNLO':>20s}  {'NegativeNLO':>20s}  {'Total (Pos-Neg)':>20s}")
    for tag in ['PP', 'PM', '00']:
        p_str = f"{pb(pos[tag][0]):.1f} +- {pb(0, pos[tag][1])[1]:.1f}" if has_pos[tag] else "NOT RUN"
        n_str = f"{pb(neg[tag][0]):.1f} +- {pb(0, neg[tag][1])[1]:.1f}" if has_neg[tag] else "NOT RUN"
        if has_pos[tag] and has_neg[tag]:
            tot = pb(pos[tag][0]) - pb(neg[tag][0])
            tot_e = math.sqrt(pb(0, pos[tag][1])[1]**2 + pb(0, neg[tag][1])[1]**2)
            t_str = f"{tot:.1f} +- {tot_e:.1f}"
        elif has_pos[tag]:
            t_str = "(need NegNLO)"
        else:
            t_str = "-"
        print(f"  {tag:>4s}  {p_str:>20s}  {n_str:>20s}  {t_str:>20s}")
    print()

    # ================================================================
    # Comparison with POLDIS
    # ================================================================

    # --- Unpolarized ---
    print("-" * 72)
    print("  UNPOLARIZED NLO COMPARISON:")
    print("-" * 72)

    if has_pos['00']:
        pos_00, pos_00_e = pb(pos['00'][0], pos['00'][1])
        print()
        print_comparison("PositiveNLO only (00) vs POLDIS NLO:",
                        pos_00, pos_00_e,
                        poldis['nlo_unpol'], poldis['nlo_unpol_e'])

        if has_neg['00']:
            neg_00, neg_00_e = pb(neg['00'][0], neg['00'][1])
            tot_00 = pos_00 - neg_00
            tot_00_e = math.sqrt(pos_00_e**2 + neg_00_e**2)
            print()
            print_comparison("Total NLO (Pos-Neg, 00) vs POLDIS NLO:",
                            tot_00, tot_00_e,
                            poldis['nlo_unpol'], poldis['nlo_unpol_e'])
            print()
            print(f"    NegativeNLO contribution: {neg_00:.1f} +- {neg_00_e:.1f} pb")
            print(f"    (= {neg_00/pos_00*100:.2f}% of PositiveNLO)")
    else:
        print("    00 PositiveNLO: NOT YET RUN")
    print()

    # --- Polarized ---
    print("-" * 72)
    print("  POLARIZED NLO COMPARISON:")
    print("-" * 72)

    if has_pos['PP'] and has_pos['PM']:
        pp_pos, pp_pos_e = pb(pos['PP'][0], pos['PP'][1])
        pm_pos, pm_pos_e = pb(pos['PM'][0], pos['PM'][1])
        diff_pos = pp_pos - pm_pos
        diff_pos_e = math.sqrt(pp_pos_e**2 + pm_pos_e**2)
        sigma_pol_pos = diff_pos / 2.0
        sigma_pol_pos_e = diff_pos_e / 2.0

        print()
        print("  From PositiveNLO only:")
        print_comparison("sigma_pol = (PP-PM)/2 vs POLDIS:",
                        sigma_pol_pos, sigma_pol_pos_e,
                        poldis['nlo_pol'], poldis['nlo_pol_e'])

        if has_neg['PP'] and has_neg['PM']:
            pp_neg, pp_neg_e = pb(neg['PP'][0], neg['PP'][1])
            pm_neg, pm_neg_e = pb(neg['PM'][0], neg['PM'][1])

            # Total = Pos - Neg for each polarization
            pp_tot = pp_pos - pp_neg
            pp_tot_e = math.sqrt(pp_pos_e**2 + pp_neg_e**2)
            pm_tot = pm_pos - pm_neg
            pm_tot_e = math.sqrt(pm_pos_e**2 + pm_neg_e**2)

            diff_tot = pp_tot - pm_tot
            diff_tot_e = math.sqrt(pp_tot_e**2 + pm_tot_e**2)
            sigma_pol_tot = diff_tot / 2.0
            sigma_pol_tot_e = diff_tot_e / 2.0

            print()
            print("  From Total NLO (Pos - Neg):")
            print_comparison("sigma_pol = (PP-PM)/2 vs POLDIS:",
                            sigma_pol_tot, sigma_pol_tot_e,
                            poldis['nlo_pol'], poldis['nlo_pol_e'])

            # Check NegativeNLO polarized component
            diff_neg = pp_neg - pm_neg
            diff_neg_e = math.sqrt(pp_neg_e**2 + pm_neg_e**2)
            print()
            print(f"    NegNLO PP-PM = {diff_neg:.2f} +- {diff_neg_e:.2f} pb")
            print(f"    (non-zero means NegNLO has polarized component)")

        print()
        # K_pol using POLDIS LO baseline
        kpol_pos = sigma_pol_pos / poldis['lo_pol']
        print(f"    K_pol (Herwig, PosNLO only)  = {kpol_pos:.4f}")
        if has_neg['PP'] and has_neg['PM']:
            kpol_tot = sigma_pol_tot / poldis['lo_pol']
            print(f"    K_pol (Herwig, Total Pos-Neg) = {kpol_tot:.4f}")
        print(f"    K_pol (POLDIS)                = {poldis['k_pol']:.4f}")
    else:
        print("    PP/PM PositiveNLO: NOT YET RUN")
    print()

    # --- Consistency checks ---
    print("-" * 72)
    print("  CONSISTENCY CHECKS:")
    print("-" * 72)
    if has_pos['PP'] and has_pos['PM'] and has_pos['00']:
        pp_v = pb(pos['PP'][0])
        pm_v = pb(pos['PM'][0])
        z_v  = pb(pos['00'][0])
        avg = (pp_v + pm_v) / 2.0
        print(f"    PositiveNLO: (PP+PM)/2 = {avg:.1f},  00 = {z_v:.1f},  diff = {avg-z_v:.1f} pb ({(avg/z_v-1)*100:.2f}%)")

    if all(has_neg[t] for t in ['PP', 'PM', '00']):
        pp_n = pb(neg['PP'][0])
        pm_n = pb(neg['PM'][0])
        z_n  = pb(neg['00'][0])
        avg_n = (pp_n + pm_n) / 2.0
        print(f"    NegativeNLO: (PP+PM)/2 = {avg_n:.1f},  00 = {z_n:.1f},  diff = {avg_n-z_n:.1f} pb ({(avg_n/z_n-1)*100:.2f}%)")

    if has_pos['PP'] and has_pos['PM'] and has_pos['00'] and all(has_neg[t] for t in ['PP', 'PM', '00']):
        pp_t = pb(pos['PP'][0]) - pb(neg['PP'][0])
        pm_t = pb(pos['PM'][0]) - pb(neg['PM'][0])
        z_t  = pb(pos['00'][0]) - pb(neg['00'][0])
        avg_t = (pp_t + pm_t) / 2.0
        print(f"    Total NLO:   (PP+PM)/2 = {avg_t:.1f},  00 = {z_t:.1f},  diff = {avg_t-z_t:.1f} pb ({(avg_t/z_t-1)*100:.2f}%)")

    print()
    print("=" * 72)


if __name__ == '__main__':
    main()
