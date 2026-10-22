#!/usr/bin/env python3
"""Compare actual Herwig HardVertex amplitudes with an independent MadLoop export.

No analytic Herwig amplitude helper is imported. Polarization overlaps determine
all helicity-label/phase conversions. No amplitude-dependent phase is fitted.
"""
import argparse
import cmath
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import subprocess

HELICITIES = list(itertools.product((0, 2), repeat=4))
COLOUR = 5 / 24


def dot(a, b):
    return a[0]*b[0] - sum(x*y for x, y in zip(a[1:], b[1:]))


def pairs(fields):
    return [complex(float(a), float(b)) for a, b in zip(fields[::2], fields[1::2])]


def require_finite(values, context):
    if any(not (math.isfinite(complex(z).real) and math.isfinite(complex(z).imag))
           for z in values):
        raise ValueError(f'Non-finite {context}')


def read_native(path):
    events = []
    event = None
    for line in path.read_text().splitlines():
        f = line.split()
        if not f:
            continue
        tag = f[0]
        if tag == 'EVENT':
            if event is not None:
                raise ValueError('Unterminated native event')
            event = dict(number=int(f[1]), momenta={}, pol={}, amp={}, rho={}, outrho={})
        elif tag == 'END':
            if (len(event['momenta']), len(event['pol']), len(event['amp']),
                    len(event['rho']), len(event['outrho'])) != (4, 8, 81, 18, 18):
                raise ValueError('Incomplete native helicity record')
            require_finite(event['param'] + event['me2'] +
                           [z for key in ('momenta', 'pol') for row in event[key].values() for z in row] +
                           [z for key in ('amp', 'rho', 'outrho') for z in event[key].values()],
                           'native record')
            events.append(event)
            event = None
        elif event is None:
            raise ValueError('Record outside native event')
        elif tag in ('PARAM', 'ME2'):
            event[tag.lower()] = list(map(float, f[1:]))
        elif tag == 'MOM':
            event['momenta'][int(f[1])] = list(map(float, f[2:]))
        elif tag == 'POL':
            event['pol'][tuple(map(int, f[1:3]))] = pairs(f[3:])
        elif tag in ('RHO', 'OUTRHO'):
            event[tag.lower()][tuple(map(int, f[1:4]))] = pairs(f[4:])[0]
        elif tag == 'AMP':
            event['amp'][tuple(map(int, f[1:5]))] = pairs(f[5:])[0]
        else:
            raise ValueError(f'Unknown record {tag}')
    if event is not None or not events:
        raise ValueError('Empty or truncated native export')
    return events


def read_madloop(text):
    result = dict(amp={}, pol={}, precision=[], hel_me2={}, momenta_used={})
    hel = None
    for line in text.splitlines():
        f = line.split()
        if not f:
            continue
        if f[0] == 'COUPLING':
            result['g'], result['e'] = map(float, f[1:])
        elif f[0] == 'TOTAL':
            result['me2'], result['total_precision'] = map(float, f[1:3])
            result['total_return_code'] = int(f[3])
        elif f[0] == 'HEL':
            hel = tuple(0 if int(x) == -1 else 2 for x in f[2:6])
            result['precision'].append(float(f[-2]))
            result['hel_me2'][hel] = float(f[6])
        elif f[0] == 'MOM_USED':
            leg = int(f[1])-1
            value = list(map(float, f[2:]))
            if leg in result['momenta_used'] and value != result['momenta_used'][leg]:
                raise ValueError('MadLoop used different kinematics across helicities')
            result['momenta_used'][leg] = value
        elif f[0] in ('AMP', 'POL'):
            vals = [complex(float(a), float(b)) for a, b in
                    re.findall(r'\(([^,]+),([^\)]+)\)', line)]
            if f[0] == 'AMP':
                result['amp'][hel] = vals
            else:
                leg = int(f[1]) - 1
                key = (leg, hel[leg])
                if key in result['pol'] and max(abs(a-b) for a, b in
                                               zip(result['pol'][key], vals)) > 1e-12:
                    raise ValueError('Inconsistent MadLoop polarization vectors')
                result['pol'][key] = vals
    if set(result['amp']) != set(HELICITIES) or len(result['pol']) != 8 or len(result['momenta_used']) != 4:
        raise ValueError('Missing MadLoop helicities or polarization vectors')
    if any(len(row) != 2 for row in result['amp'].values()) or any(len(row) != 4 for row in result['pol'].values()):
        raise ValueError('Unexpected MadLoop colour or vector dimensions')
    require_finite([result['g'], result['e'], result['me2']] +
                   [z for key in ('amp', 'pol', 'momenta_used') for row in result[key].values() for z in row],
                   'MadLoop record')
    return result


def convert_basis(mg, native):
    """Map MadLoop's colour-flow coefficient to Herwig's d^{abc} coefficient."""
    overlaps = {}
    residual = 0.
    unitary_error = 0.
    for leg in range(4):
        for h in (0, 2):
            hpol = native['pol'][leg, h]
            for m in (0, 2):
                overlaps[leg, h, m] = -dot(hpol, [z.conjugate() for z in mg['pol'][leg, m]])
            remainder = [hpol[k] - sum(overlaps[leg, h, m] * mg['pol'][leg, m][k]
                                      for m in (0, 2)) for k in range(4)]
            p = native['momenta'][leg]
            pivot = max(range(4), key=lambda k: abs(p[k]))
            gauge = remainder[pivot] / p[pivot]
            residual = max(residual, max(abs(remainder[k] - gauge*p[k]) for k in range(4)))
        for h in (0, 2):
            for k in (0, 2):
                gram = sum(overlaps[leg, h, m] * overlaps[leg, k, m].conjugate() for m in (0, 2))
                unitary_error = max(unitary_error, abs(gram - (h == k)))
    # Tr(Ta Tb Tc) + Tr(Tc Tb Ta) = d^{abc}/2. Both generated
    # colour-flow coefficients are equal; each is twice Herwig's coefficient.
    scale = max(abs(z) for row in mg['amp'].values() for z in row)
    colour_error = max(abs(row[0]-row[1]) for row in mg['amp'].values()) / scale
    amplitudes = {}
    for h in HELICITIES:
        amplitudes[h] = sum(mg['amp'][m][0] / 2 * math.prod(
            overlaps[leg, h[leg], m[leg]] for leg in range(4)) for m in HELICITIES)
    momentum_scale = max(p[0] for p in native['momenta'].values())
    momentum_error = max(abs(mg['momenta_used'][i][k]-native['momenta'][i][k])
                         for i in range(4) for k in range(4))/momentum_scale
    return amplitudes, dict(polarization_gauge_residual=residual,
                            polarization_unitarity_error=unitary_error,
                            momentum_repair_relative=momentum_error,
                            colour_flow_error=colour_error)


def weighted(amp, rho):
    total = 0j
    for a, b, ap, bp, c, d in itertools.product((0, 2), repeat=6):
        total += rho[0, a, ap] * rho[1, b, bp] * amp[a, b, c, d] * amp[ap, bp, c, d].conjugate()
    return COLOUR * total.real


def density(amp, rho, leg):
    result = {(i, j): 0j for i in (0, 2) for j in (0, 2)}
    for a, b, ap, bp, i, j, other in itertools.product((0, 2), repeat=7):
        h = (a, b, i, other) if leg == 0 else (a, b, other, i)
        hp = (ap, bp, j, other) if leg == 0 else (ap, bp, other, j)
        result[i, j] += rho[0, a, ap] * rho[1, b, bp] * amp[h] * amp[hp].conjugate()
    trace = result[0, 0] + result[2, 2]
    return {key: value/trace for key, value in result.items()}


def metrics(native, mg, mapped, phase, basis):
    actual = {h: native['amp'][h] for h in HELICITIES}
    scale = max(map(abs, mapped.values()))
    aligned = {h: phase*z for h, z in mapped.items()}
    errors = {h: abs(actual[h]-aligned[h]) for h in HELICITIES}
    longitudinal = max(abs(z) for h, z in native['amp'].items() if 1 in h)
    unpolarized = COLOUR / 4 * sum(abs(z)**2 for z in actual.values())
    mg_weighted = weighted(mapped, native['rho'])
    relative_error = max(errors[h]/max(abs(mapped[h]), 1e-6*scale) for h in HELICITIES)
    density_error = 0.
    for leg in (0, 1):
        expected = density(mapped, native['rho'], leg)
        for i in range(3):
            for j in range(3):
                density_error = max(density_error, abs(native['outrho'][leg, i, j] - expected.get((i, j), 0j)))
    result = dict(event=native['number'], **basis,
                  native_me2=native['me2'][0], native_cached_me2=native['me2'][1],
                  native_vertex_me2=native['me2'][2], madloop_weighted_me2=mg_weighted,
                  unpolarized_native_me2=unpolarized, unpolarized_madloop_me2=mg['me2'],
                  amplitude_error_scaled=max(errors.values())/scale,
                  amplitude_relative_error=relative_error,
                  interference_error_scaled=max(abs(actual[h]*actual[k].conjugate() -
                        mapped[h]*mapped[k].conjugate()) for h in HELICITIES for k in HELICITIES)/scale**2,
                  longitudinal_amplitude_max=longitudinal,
                  unpolarized_me2_relative_error=abs(unpolarized/mg['me2']-1),
                  weighted_me2_relative_error=abs(native['me2'][0]/mg_weighted-1),
                  weight_vertex_relative_error=abs(native['me2'][0]/native['me2'][2]-1),
                  cached_weight_relative_error=(abs(native['me2'][0]/native['me2'][1]-1)
                      if native['me2'][1] else None),
                  outgoing_density_max_error=density_error,
                  madloop_precision=max(mg['precision']),
                  madloop_total_return_code=mg['total_return_code'],
                  mu2=native['param'][0], alpha_s=native['param'][1], alpha_em=native['param'][2])
    result['helicities'] = [dict(helicity=list(h), herwig=[actual[h].real, actual[h].imag],
        madloop_aligned=[aligned[h].real, aligned[h].imag],
        abs_error=errors[h], relative_error=errors[h]/max(abs(mapped[h]), 1e-6*scale)) for h in HELICITIES]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('native', type=Path)
    parser.add_argument('--madloop', type=Path, required=True, help='Diagnostic SubProcesses/PV... directory')
    parser.add_argument('--output', type=Path, required=True, help='New comparison directory')
    parser.add_argument('--mode', choices=('Massive', 'Massless5'), required=True)
    parser.add_argument('--reuse', action='store_true', help='Reuse check outputs only if input and executable hashes match')
    parser.add_argument('--tolerance', type=float, default=1e-6)
    args = parser.parse_args()
    args.madloop = args.madloop.resolve()
    args.output.mkdir(parents=True, exist_ok=args.reuse)
    events = read_native(args.native)
    proc_card = (args.madloop/'../../Cards/proc_card_mg5.dat').resolve().read_text().replace('\\\n', '')
    if args.mode == 'Massless5' and not re.search(r'generate g g > a g / t ', proc_card):
        raise ValueError('Massless5 reference must explicitly exclude the top loop')
    if args.mode == 'Massive' and 'import model loop_sm-c_mass' not in proc_card:
        raise ValueError('Massive reference must use loop_sm-c_mass')
    mass_block = False
    masses = {}
    param_path = args.madloop/'param_card.dat'
    for line in param_path.read_text().splitlines():
        fields = line.split('#')[0].split()
        if not fields:
            continue
        if fields[0].lower() in ('block', 'decay'):
            mass_block = fields[0].lower() == 'block' and fields[1].lower() == 'mass'
        elif mass_block:
            masses[int(fields[0])] = float(fields[1].lower().replace('d', 'e'))
    for native in events:
        expected = native['param'][3:] if args.mode == 'Massive' else (0., 0.)
        for flavour, mass in enumerate(expected, 4):
            if abs(masses.get(flavour, 0.) - mass) > 1e-12*max(1., abs(mass)):
                raise ValueError(f'Reference quark mass mismatch for PDG {flavour}')
        if any(masses.get(q, 0.) != 0. for q in (1, 2, 3)):
            raise ValueError('Reference u,d,s must be massless')
    binary = args.madloop/'check'
    binary_hash = hashlib.sha256(binary.read_bytes()).hexdigest()
    param_hash = hashlib.sha256(param_path.read_bytes()).hexdigest()
    settings_hash = hashlib.sha256((args.madloop/'../MadLoop5_resources/MadLoopParams.dat').read_bytes()).hexdigest()
    prepared = []
    for index, native in enumerate(events):
        payload = ' '.join(map(str, native['param'][:3]))+'\n'
        payload += '\n'.join(' '.join(map(str, native['momenta'][i])) for i in range(4))+'\n'
        fingerprint = hashlib.sha256((binary_hash+param_hash+settings_hash+payload).encode()).hexdigest()
        output = args.output/f'point-{index:04d}.txt'
        stamp = output.with_suffix('.sha256')
        if not (args.reuse and output.exists() and stamp.exists() and stamp.read_text() == fingerprint):
            process = subprocess.run([str(binary)], cwd=args.madloop, input=payload,
                                     text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
            output.write_text(process.stdout)
            process.check_returncode()
            stamp.write_text(fingerprint)
        mg = read_madloop(output.read_text())
        if abs(mg['g']**2/(4*math.pi)/native['param'][1]-1) > 1e-12 or abs(mg['e']**2/(4*math.pi)/native['param'][2]-1) > 1e-12:
            raise ValueError('MadLoop couplings do not match native couplings')
        mapped, basis = convert_basis(mg, native)
        prepared.append((native, mg, mapped, basis))
        print(f'MadLoop {index+1}/{len(events)} evaluated', flush=True)
    # The stored-vector comparison fixes the convention to +1. Report the
    # measured common phase as a diagnostic, but never fit it away.
    overlap = sum(sum(native['amp'][h] * mapped[h].conjugate() for h in HELICITIES) /
                  sum(abs(z)**2 for z in mapped.values()) for native, _, mapped, _ in prepared)
    phase = complex(1.)
    results = [metrics(native, mg, mapped, phase, basis) for native, mg, mapped, basis in prepared]
    keys = ['amplitude_error_scaled', 'amplitude_relative_error', 'interference_error_scaled',
            'unpolarized_me2_relative_error', 'weighted_me2_relative_error',
            'weight_vertex_relative_error',
            'outgoing_density_max_error', 'longitudinal_amplitude_max',
            'polarization_gauge_residual', 'polarization_unitarity_error',
            'momentum_repair_relative', 'colour_flow_error']
    maxima = {key: max(r[key] for r in results) for key in keys}
    report = dict(mode=args.mode, events=len(results), helicities=16*len(results),
                  interference_products=256*len(results), common_phase_radians=cmath.phase(phase),
                  observed_common_phase_radians=cmath.phase(overlap),
                  source=str(args.native.resolve()), source_sha256=hashlib.sha256(args.native.read_bytes()).hexdigest(),
                  comparator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  madloop_check_sha256=binary_hash, madloop_param_card_sha256=param_hash,
                  madloop_settings_sha256=settings_hash, maxima=maxima, results=results,
                  tolerance=args.tolerance, passed=all(v < args.tolerance for v in maxima.values()))
    (args.output/'comparison.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'results'}, indent=2))
    if not report['passed']:
        raise SystemExit('Actual native helicity comparison failed; inspect comparison.json')


if __name__ == '__main__':
    main()
