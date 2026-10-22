#!/usr/bin/env python3
"""Compare complex amplitudes, interference products and absolute ME2.

The independently generated MG colour-flow coefficient equals twice the
d^{abc} coefficient. Only ONE common phase is aligned, never a helicity phase.
"""
import argparse
import cmath
import json
import math
from pathlib import Path
import re
import subprocess


def compare(path, probe, masses):
    momenta, polarizations, amplitudes, precisions = [], [], [], []
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == 'COUPLING':
            g, e = map(float, fields[1:])
        elif fields[0] == 'TOTAL':
            total = float(fields[1])
        elif fields[0] == 'MOMENTUM':
            momenta.extend(map(float, fields[2:]))
        elif fields[0] == 'HEL':
            precisions.append(float(fields[-2]))
        elif fields[0] in ('AMP', 'POL'):
            values = [complex(float(a), float(b)) for a, b in
                      re.findall(r'\(([^,]+),([^\)]+)\)', line)]
            if fields[0] == 'AMP':
                amplitudes.append(values)
            else:
                for z in values:
                    polarizations.extend((z.real, z.imag))
    assert len(amplitudes) == 16 and len(polarizations) == 512
    inputs = ' '.join(map(str, [*masses, *momenta, *polarizations]))
    output = subprocess.check_output([str(probe)], input=inputs, text=True)
    reduced = [complex(*map(float, line.split())) for line in output.splitlines()]
    coupling = 2 * e * g**3 / (4 * math.pi)**2
    calculated = [2 * coupling * z for z in reduced]
    reference = [z[0] for z in amplitudes]
    overlap = sum(a * b.conjugate() for a, b in zip(reference, calculated))
    phase = overlap / abs(overlap)
    scale = max(map(abs, reference))
    amplitude_error = max(abs(a-phase*b) for a, b in zip(reference, calculated))/scale
    interference_error = max(abs(a*c.conjugate()-b*d.conjugate())
                             for a, b in zip(reference, calculated)
                             for c, d in zip(reference, calculated))/scale**2
    me2 = 5/24 * coupling**2 * sum(abs(z)**2 for z in reduced)/4
    result = dict(source=str(path), me2=me2, madloop_me2=total,
                  relative_me2_error=abs(me2/total-1),
                  amplitude_error=amplitude_error,
                  interference_error=interference_error,
                  common_phase=cmath.phase(phase),
                  madloop_precision=max(precisions))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files', type=Path, nargs='+')
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--masses', type=float, nargs=3, default=(1.27, 4.18, 172.69))
    parser.add_argument('--tolerance', type=float, default=1e-6)
    args = parser.parse_args()
    results = [compare(p, args.probe.resolve(), args.masses) for p in args.files]
    print(json.dumps(results, indent=2))
    assert all(max(r[k] for k in ('relative_me2_error', 'amplitude_error',
                                 'interference_error')) < args.tolerance for r in results)
