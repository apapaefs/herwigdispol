# LO NC DIS Calculator Note

## Purpose

`calc_lo_dis_nc.py` is a small LO neutral-current DIS reference calculator for
the current DISPOL workflow. It is meant to provide a POLDIS-style baseline for
`GAMMA`, `Z`, `INT = ALL - GAMMA - Z`, and `ALL`.

V1 reports:

- `sigma_unpol`
- `sigma_LL`

It does not attempt to reconstruct Herwig `00`, `PP`, `PM`, `MP`, or `MM`.

## Default plain39-style inputs

- Beam energies: `E_e = 18 GeV`, `E_p = 275 GeV`
- Phase space: `49 <= Q^2 <= 2500 GeV^2`, `0.2 <= y <= 0.6`
- Electroweak inputs:
  - `alpha_EM(M_Z) = 0.00729927`
  - `sin^2(theta_W) = 0.221639970483740179`
  - `M_Z = 91.1876 GeV`
  - `Gamma_Z = 2.4952 GeV`
- PDFs:
  - unpolarized: `PDF4LHC15_nnlo_100_pdfas`
  - polarized difference: `BDSSV24-NNLO`

The active flavor sum uses `u,d,s,c,b` and antiquarks.

## LO formulas

For each channel, the code constructs the POLDIS-style vector and axial channel
pieces:

- `GAMMA`
- `Z`
- `INT`
- `ALL = GAMMA + Z + INT`

and then evaluates

- `sigma_unpol ~ Y_+ * vec + Y_- * ax`
- `sigma_LL ~ Y_- * vec + Y_+ * ax`

with

- `Y_+ = 1 + (1-y)^2`
- `Y_- = 1 - (1-y)^2`
- `x_B = Q^2 / (s y)`, `s = 4 E_e E_p`

The absolute LO differential normalization matches the POLDIS LO convention:

- `d^2sigma / dQ^2 dy = 2*pi*alpha_EM^2 / (y Q^4) * structure * NRM`
- `NRM = 389385 * 1000`

## Conventions

Default behavior is intentionally POLDIS-style. In particular, the default `Z`
treatment follows the fixed-width DIS convention used in the POLDIS formulas:

- `PROPZ = Q^4 / ((Q^2 + M_Z^2)^2 + (M_Z Gamma_Z)^2)`
- `INTGZ = LEPCHARGE * (Q^2 + M_Z^2) / Q^2`

This is intentionally the POLDIS-style spacelike fixed-width choice.

The calculator also supports `--herwig-lo-mode`, which switches to a closer
Herwig-LO convention:

- spacelike `Z` propagator with zero width
- `sin^2(theta_W) = 1 - (M_W/M_Z)^2`
- Herwig-like `ep` cut enforcement with
  - `x1 = 1`
  - `x2 = x_B`
  - optional `X1Min`, `X2Min`, `W2Min`, `W2Max`

This mode is meant for comparison studies, not as the default physics
reference.

The calculator also supports an optional approximate quark-mass model for
comparison studies:

- `--quark-mass-mode massive-born`
- `--quark-mass-mode slowrescale`
- `--quark-mass-mode slowrescale-jacobian`

`massive-born` is the preferred comparison option when you want to keep PDFs at
`x_B` and modify only the LO partonic kernel/kinematics. It uses a numerical
2->2 Born helicity trace with massive incoming/outgoing quarks to build a
channel-by-channel kernel ratio and applies that to the exact massless
POLDIS-style kernels.

The slow-rescaling options are cruder sensitivity models built on the massless
LO kernels:

- `x_pdf = x_B * (1 + m_q^2 / Q^2)`
- optional extra Jacobian factor `Q^2 / (Q^2 + m_q^2)`

When enabled, the default masses are the Herwig nominal quark masses:

- `m_u = 0.00216 GeV`
- `m_d = 0.00467 GeV`
- `m_s = 0.0934 GeV`
- `m_c = 1.27 GeV`
- `m_b = 4.18 GeV`

They can be overridden with `--mq-u`, `--mq-d`, `--mq-s`, `--mq-c`, and
`--mq-b`.

## Modes

- Default run: one representative pointwise report plus one integrated report
- `point`: fixed `(Q^2, y)` evaluation
- `integrate`: cut-integrated LO cross section with deterministic
  Gauss-Legendre quadrature in `log(Q^2)` and `y`
