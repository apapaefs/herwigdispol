#!/usr/bin/env python3
"""
Parse NLO_CUM diagnostic output from Herwig PP, PM, 00 runs.

The C++ diagnostic uses Born-weighted ratio estimators:
  K_X = Σ(X_i × w_Born_i) / Σ(w_Born_i)
where w_Born = eq² × f_q(xB) × (1 + a_born×l + l²) ∝ Born cross section.
This gives correct K-factors regardless of VEGAS importance sampling.

Extracts the final cumulative K-values and computes:
  1. Component-by-component K comparison between PP, PM, 00
  2. Polarized K-factor decomposition: K_pol for each NLO component

Usage:
    python parse_nlo_cum.py PP.log PM.log [00.log]
"""
import re
import sys


def parse_nlo_cum(filename):
    """Parse NLO_CUM lines from a Herwig log file.
    Returns the last (most converged) set of K-values."""
    data = {}
    data_real = {}
    # Pattern for a single float: optional sign, digits/dot, optional exponent
    _fp = r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?'
    with open(filename) as f:
        for line in f:
            if line.startswith('NLO_CUM_REAL '):
                for m in re.finditer(r'(\w+)=(' + _fp + r')', line):
                    key, val = m.group(1), float(m.group(2))
                    data_real[key] = val
                m2 = re.match(r'NLO_CUM_REAL\s+(\d+)', line)
                if m2:
                    data_real['N'] = int(m2.group(1))
            elif line.startswith('NLO_CUM '):
                # Handle K_wgt=VALUE+-ERROR format
                for m in re.finditer(r'(\w+)=(' + _fp + r')(?:\+-(' + _fp + r'))?', line):
                    key = m.group(1)
                    data[key] = float(m.group(2))
                    if m.group(3):
                        data[key + '_err'] = float(m.group(3))
                m2 = re.match(r'NLO_CUM\s+(\d+)', line)
                if m2:
                    data['N'] = int(m2.group(1))
    data.update(data_real)
    return data


def main():
    if len(sys.argv) < 3:
        print("Usage: python parse_nlo_cum.py PP.log PM.log [00.log]")
        sys.exit(1)

    pp_file = sys.argv[1]
    pm_file = sys.argv[2]
    oo_file = sys.argv[3] if len(sys.argv) > 3 else None

    pp = parse_nlo_cum(pp_file)
    pm = parse_nlo_cum(pm_file)
    oo = parse_nlo_cum(oo_file) if oo_file else None

    if not pp or not pm:
        print("ERROR: Could not parse NLO_CUM data from log files")
        sys.exit(1)

    print("=" * 80)
    print("  Born-Weighted K-Factor Decomposition: PP vs PM" +
          (" vs 00" if oo else ""))
    print("=" * 80)
    print()
    print(f"  Events: PP={pp.get('N','?')}, PM={pm.get('N','?')}" +
          (f", 00={oo.get('N','?')}" if oo else ""))
    print()
    print("  K_X = Born-weighted average of NLO component X")
    print("  K_wgt = K_virt + K_cqu + K_cqp + K_cgu + K_cgp + K_rq + K_rg")
    print("  K_wgt ≈ σ_NLO / σ_Born (true K-factor)")
    print()

    # Components to analyze
    components = [
        ('K_virt',  'Virtual + soft'),
        ('K_cqu',   'Coll. quark (unpol)'),
        ('K_cqp',   'Coll. quark (pol)'),
        ('K_cgu',   'Coll. gluon (unpol)'),
        ('K_cgp',   'Coll. gluon (pol)'),
        ('K_rq',    'Real quark (QCDC)'),
        ('K_rg',    'Real gluon (BGF)'),
        ('K_wgt',   'TOTAL K'),
        ('K_wold',  'Old F2 K (check)'),
    ]

    # --- Table of K-factors ---
    print("-" * 80)
    header = f"  {'Component':<22s}  {'K_PP':>14s}  {'K_PM':>14s}"
    if oo:
        header += f"  {'K_00':>14s}"
    header += f"  {'K_PP - K_PM':>14s}"
    print(header)
    print("-" * 80)

    for key, label in components:
        pp_val = pp.get(key, 0.0)
        pm_val = pm.get(key, 0.0)
        diff = pp_val - pm_val
        line = f"  {label:<22s}  {pp_val:>14.8f}  {pm_val:>14.8f}"
        if oo:
            oo_val = oo.get(key, 0.0)
            line += f"  {oo_val:>14.8f}"
        line += f"  {diff:>14.8f}"
        print(line)

    print("-" * 80)

    # --- Real emission decomposition ---
    real_comps = [
        ('K_rqKu', 'realq unpol kernel'),
        ('K_rqKp', 'realq pol kernel'),
        ('K_rgKu', 'realg unpol kernel'),
        ('K_rgKp', 'realg pol kernel'),
    ]

    print()
    print("-" * 80)
    print("  Real emission kernel decomposition:")
    print("-" * 80)
    header2 = f"  {'Component':<22s}  {'K_PP':>14s}  {'K_PM':>14s}"
    if oo:
        header2 += f"  {'K_00':>14s}"
    header2 += f"  {'K_PP - K_PM':>14s}"
    print(header2)
    print("-" * 80)
    for key, label in real_comps:
        pp_val = pp.get(key, 0.0)
        pm_val = pm.get(key, 0.0)
        diff = pp_val - pm_val
        line = f"  {label:<22s}  {pp_val:>14.8f}  {pm_val:>14.8f}"
        if oo:
            oo_val = oo.get(key, 0.0)
            line += f"  {oo_val:>14.8f}"
        line += f"  {diff:>14.8f}"
        print(line)
    print("-" * 80)

    # --- Polarized K-factor ---
    print()
    print("=" * 80)
    print("  Polarized K-factor (K_pol)")
    print("=" * 80)
    print()
    print("  K_pol = [K_PP × σ_Born_PP − K_PM × σ_Born_PM] / [σ_Born_PP − σ_Born_PM]")
    print()

    # POLDIS reference Born cross sections
    sigma_unpol = 6600.787  # pb (POLDIS LO)
    sigma_pol   = 52.185    # pb (POLDIS LO)
    sigma_pp = sigma_unpol + sigma_pol   # 6652.972
    sigma_pm = sigma_unpol - sigma_pol   # 6548.602
    delta_sigma = sigma_pp - sigma_pm    # = 2 × σ_pol = 104.370

    print(f"  σ_Born (POLDIS LO, pb):")
    print(f"    PP  = {sigma_pp:.1f}  (σ_unpol + σ_pol)")
    print(f"    PM  = {sigma_pm:.1f}  (σ_unpol − σ_pol)")
    print(f"    PP−PM = {delta_sigma:.1f}  (= 2×σ_pol)")
    print()

    nlo_components = [
        ('K_virt',  'Virtual + soft'),
        ('K_cqu',   'Coll. quark (unpol)'),
        ('K_cqp',   'Coll. quark (pol)'),
        ('K_cgu',   'Coll. gluon (unpol)'),
        ('K_cgp',   'Coll. gluon (pol)'),
        ('K_rq',    'Real quark (QCDC)'),
        ('K_rg',    'Real gluon (BGF)'),
        ('K_wgt',   'TOTAL'),
    ]

    print(f"  {'Component':<22s}  {'K_PP×σ_PP':>12s}  {'K_PM×σ_PM':>12s}  {'PP−PM (pb)':>12s}  {'K_pol_comp':>12s}")
    print("-" * 80)

    for key, label in nlo_components:
        kpp = pp.get(key, 0.0)
        kpm = pm.get(key, 0.0)
        s_pp = kpp * sigma_pp
        s_pm = kpm * sigma_pm
        pp_m_pm = s_pp - s_pm
        k_pol_comp = pp_m_pm / delta_sigma if delta_sigma != 0 else 0
        print(f"  {label:<22s}  {s_pp:>12.2f}  {s_pm:>12.2f}  {pp_m_pm:>12.4f}  {k_pol_comp:>12.6f}")

    print("-" * 80)
    print()

    # Total K_pol
    K_pp = pp.get('K_wgt', 0.0)
    K_pm = pm.get('K_wgt', 0.0)
    K_pp_err = pp.get('K_wgt_err', 0.0)
    K_pm_err = pm.get('K_wgt_err', 0.0)
    K_pol = (K_pp * sigma_pp - K_pm * sigma_pm) / delta_sigma
    # Error propagation: δK_pol = sqrt((σ_PP × δK_PP)² + (σ_PM × δK_PM)²) / Δσ
    K_pol_err = ((sigma_pp * K_pp_err)**2 + (sigma_pm * K_pm_err)**2)**0.5 / delta_sigma

    K_unpol = K_pp if oo is None else oo.get('K_wgt', K_pp)

    print(f"  K_pol (Herwig)  = {K_pol:.6f} ± {K_pol_err:.6f}")
    print(f"  K_pol (POLDIS)  = 1.020100")
    print(f"  K_unpol (Herwig) = {K_unpol:.6f}")
    print(f"  K_unpol (POLDIS) = 0.938900")
    print()
    if K_pol_err > 0:
        tension = abs(K_pol - 1.0201) / K_pol_err
        print(f"  Tension: {tension:.1f}σ")
    print()

    # Consistency checks
    print("=" * 80)
    print("  Consistency checks:")
    print("=" * 80)

    # Components sum to total
    for tag, d in [('PP', pp), ('PM', pm)] + ([('00', oo)] if oo else []):
        total = sum(d.get(k, 0) for k in ['K_virt', 'K_cqu', 'K_cqp', 'K_cgu', 'K_cgp', 'K_rq', 'K_rg'])
        wgt = d.get('K_wgt', 0)
        print(f"  {tag}: Σ(K_components) = {total:.8f}, K_wgt = {wgt:.8f}, diff = {total-wgt:.2e}")

    # Real emission decomposition
    for tag, d in [('PP', pp), ('PM', pm)] + ([('00', oo)] if oo else []):
        rq = d.get('K_rq', 0)
        rq_sum = d.get('K_rqKu', 0) + d.get('K_rqKp', 0)
        rg = d.get('K_rg', 0)
        rg_sum = d.get('K_rgKu', 0) + d.get('K_rgKp', 0)
        print(f"  {tag}: K_rq={rq:.8f} vs Ku+Kp={rq_sum:.8f}; "
              f"K_rg={rg:.8f} vs Ku+Kp={rg_sum:.8f}")

    if oo:
        # For 00: K_wgt includes F_L at NLO; K_wold is F₂-only.
        # K_wold should match POLDIS K_unpol = 0.939 (F₂-only reference).
        # K_wgt - K_wold = NLO F_L contribution (~1-2% effect).
        k_wgt = oo.get('K_wgt', 0)
        k_wold = oo.get('K_wold', 0)
        print(f"  00: K_wgt = {k_wgt:.8f} (F₂+F_L), K_wold = {k_wold:.8f} (F₂ only)")
        print(f"      diff = {k_wgt - k_wold:.2e} (= NLO F_L contribution)")
        print(f"      K_wold should match POLDIS K_unpol = 0.9389")
        # Polarized components should be zero for 00
        for key in ['K_cqp', 'K_cgp', 'K_rqKp', 'K_rgKp']:
            val = oo.get(key, 0)
            print(f"  00: {key} = {val:.2e} (should be ~0)")

    print()


if __name__ == '__main__':
    main()
