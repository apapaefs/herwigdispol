#!/usr/bin/env python3.10
"""LO neutral-current DIS calculator for POLDIS-style gamma/Z/ALL references."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from typing import Iterable, List, Sequence

from lo_dis_nc_common import (
    CHANNEL_ALL,
    CHANNEL_GAMMA,
    CHANNEL_INT,
    CHANNEL_ORDER,
    CHANNEL_Z,
    BeamInputs,
    DEFAULT_ALPHA_EM,
    DEFAULT_BEAM_EA_GEV,
    DEFAULT_BEAM_EB_GEV,
    DEFAULT_DIFF_PDF,
    DEFAULT_GAMMA_Z,
    DEFAULT_LEPTON_CHARGE,
    DEFAULT_MODE_HERWIG_LO,
    DEFAULT_MODE_POLDIS,
    DEFAULT_QUARK_MASS_MODEL,
    DEFAULT_QUARK_MASS_MODEL_MASSIVE_BORN,
    DEFAULT_QUARK_MASS_MODEL_SLOWRESCALE,
    DEFAULT_QUARK_MASS_MODEL_SLOWRESCALE_JACOBIAN,
    DEFAULT_MW,
    DEFAULT_MZ,
    DEFAULT_Q2_MAX,
    DEFAULT_Q2_MIN,
    DEFAULT_Q2_POINT,
    DEFAULT_SUM_PDF,
    DEFAULT_Y_MAX,
    DEFAULT_Y_MIN,
    DEFAULT_Y_POINT,
    EWInputs,
    LOCalculatorMode,
    QuarkMassModel,
    filter_rows_by_channel,
    herwig_default_quark_mass_map,
    import_lhapdf,
    integrate_contributions,
    point_audit_observables,
    point_contributions,
    resolve_sin2_theta_w,
    serialize_dataclass_rows,
    sum_integrated_contributions,
    sum_point_contributions,
    x_bjorken,
)


DEFAULT_MODE = "both"
DEFAULT_NQ2 = 64
DEFAULT_NY = 48


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        nargs="?",
        choices=("point", "integrate", "both"),
        default=DEFAULT_MODE,
        help="Run a fixed point, an integrated cross section, or both (default).",
    )
    parser.add_argument(
        "--channel",
        default="all",
        help="Select one channel (GAMMA, Z, INT, ALL) or use 'all' for all four.",
    )
    parser.add_argument(
        "--flavor-scope",
        choices=("total", "per-flavor", "both"),
        default="both",
        help="Whether to print only totals, only per-flavor rows, or both.",
    )
    parser.add_argument("--q2", type=float, default=DEFAULT_Q2_POINT, help="Point-mode Q^2 in GeV^2.")
    parser.add_argument("--y", type=float, default=DEFAULT_Y_POINT, help="Point-mode y.")
    parser.add_argument("--q2-min", type=float, default=DEFAULT_Q2_MIN, help="Integrated minimum Q^2 in GeV^2.")
    parser.add_argument("--q2-max", type=float, default=DEFAULT_Q2_MAX, help="Integrated maximum Q^2 in GeV^2.")
    parser.add_argument("--y-min", type=float, default=DEFAULT_Y_MIN, help="Integrated minimum y.")
    parser.add_argument("--y-max", type=float, default=DEFAULT_Y_MAX, help="Integrated maximum y.")
    parser.add_argument("--nq2", type=int, default=DEFAULT_NQ2, help="Gauss-Legendre points in log(Q^2).")
    parser.add_argument("--ny", type=int, default=DEFAULT_NY, help="Gauss-Legendre points in y.")
    parser.add_argument("--theta", type=float, default=None, help="Weak mixing angle theta_W in radians.")
    parser.add_argument(
        "--sin2-theta-w",
        type=float,
        default=None,
        help="sin^2(theta_W). Defaults to the plain39 POLDIS value unless --herwig-lo-mode is used.",
    )
    parser.add_argument("--alpha-em", type=float, default=DEFAULT_ALPHA_EM, help="Fixed alpha_EM input.")
    parser.add_argument("--mw", type=float, default=DEFAULT_MW, help="W mass in GeV.")
    parser.add_argument("--mz", type=float, default=DEFAULT_MZ, help="Z mass in GeV.")
    parser.add_argument("--gamma-z", type=float, default=DEFAULT_GAMMA_Z, help="Z width in GeV.")
    parser.add_argument("--beam-ea", type=float, default=DEFAULT_BEAM_EA_GEV, help="Lepton beam energy in GeV.")
    parser.add_argument("--beam-eb", type=float, default=DEFAULT_BEAM_EB_GEV, help="Hadron beam energy in GeV.")
    parser.add_argument("--sum-pdf", default=DEFAULT_SUM_PDF, help="Unpolarized LHAPDF set.")
    parser.add_argument("--diff-pdf", default=DEFAULT_DIFF_PDF, help="Polarized-difference LHAPDF set.")
    parser.add_argument(
        "--audit-qid",
        type=int,
        default=None,
        help="Optional single-flavor PDG id for the LO_GAMMA_POINT-style audit summary in point mode.",
    )
    parser.add_argument(
        "--audit-jacobian",
        type=float,
        default=None,
        help="Optional Herwig jacobian value used to reconstruct sigma_hat in point mode.",
    )
    parser.add_argument(
        "--audit-s-hat",
        type=float,
        default=None,
        help="Optional Herwig sHat in GeV^2 used to reconstruct sigma_hat in point mode.",
    )
    parser.add_argument(
        "--audit-xb",
        type=float,
        default=None,
        help="Optional Herwig Bjorken x_B used for the audit PDF lookup when it differs from Q^2/(y s).",
    )
    parser.add_argument(
        "--herwig-lo-mode",
        action="store_true",
        help="Use a Herwig-LO-like convention: spacelike Z with zero width, sin^2(theta_W)=1-(mW/mZ)^2, and optional ep x/W^2 cuts.",
    )
    parser.add_argument(
        "--quark-mass-mode",
        choices=("none", "slowrescale", "slowrescale-jacobian", "massive-born"),
        default="none",
        help="Optional quark-mass model. 'massive-born' keeps PDFs at x_B and modifies only the LO partonic kernel/kinematics. The slow-rescaling modes use Herwig nominal quark masses unless overridden.",
    )
    parser.add_argument("--mq-u", type=float, default=None, help="Override the u/ubar mass in GeV.")
    parser.add_argument("--mq-d", type=float, default=None, help="Override the d/dbar mass in GeV.")
    parser.add_argument("--mq-s", type=float, default=None, help="Override the s/sbar mass in GeV.")
    parser.add_argument("--mq-c", type=float, default=None, help="Override the c/cbar mass in GeV.")
    parser.add_argument("--mq-b", type=float, default=None, help="Override the b/bbar mass in GeV.")
    parser.add_argument("--x1-min", type=float, default=None, help="Override the Herwig-like X1Min cut used in --herwig-lo-mode.")
    parser.add_argument("--x2-min", type=float, default=None, help="Override the Herwig-like X2Min cut used in --herwig-lo-mode.")
    parser.add_argument("--w2-min", type=float, default=None, help="Override the Herwig-like W^2 minimum cut in GeV^2.")
    parser.add_argument("--w2-max", type=float, default=None, help="Override the Herwig-like W^2 maximum cut in GeV^2.")
    parser.add_argument(
        "--lepton-charge",
        type=float,
        default=DEFAULT_LEPTON_CHARGE,
        help="Incoming lepton charge: -1 for e-, +1 for e+.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Append a machine-readable JSON dump after the text report.",
    )
    return parser.parse_args(argv)


def resolve_channel_selection(selection: str) -> List[str]:
    if selection == "all":
        return list(CHANNEL_ORDER)
    upper = selection.upper()
    if upper in CHANNEL_ORDER:
        return [upper]
    raise SystemExit(
        f"Unknown channel selection '{selection}'. Use GAMMA, Z, INT, ALL, or all."
    )


def build_table(headers: Sequence[str], rows: Iterable[Sequence[str]]) -> str:
    data = [list(row) for row in rows]
    widths = [len(header) for header in headers]
    for row in data:
        for idx, cell in enumerate(row):
            widths[idx] = max(widths[idx], len(cell))

    def fmt(row: Sequence[str]) -> str:
        return "  ".join(cell.ljust(widths[idx]) for idx, cell in enumerate(row))

    lines = [fmt(headers), "  ".join("-" * width for width in widths)]
    lines.extend(fmt(row) for row in data)
    return "\n".join(lines)


def fmt_float(value: float) -> str:
    return f"{value:+.12e}"


def make_calculator_mode(args: argparse.Namespace) -> LOCalculatorMode:
    if not args.herwig_lo_mode:
        return DEFAULT_MODE_POLDIS
    return LOCalculatorMode(
        name=DEFAULT_MODE_HERWIG_LO.name,
        spacelike_z_width=DEFAULT_MODE_HERWIG_LO.spacelike_z_width,
        enforce_ep_x_cuts=DEFAULT_MODE_HERWIG_LO.enforce_ep_x_cuts,
        x1_min=DEFAULT_MODE_HERWIG_LO.x1_min if args.x1_min is None else args.x1_min,
        x2_min=DEFAULT_MODE_HERWIG_LO.x2_min if args.x2_min is None else args.x2_min,
        w2_min=DEFAULT_MODE_HERWIG_LO.w2_min if args.w2_min is None else args.w2_min,
        w2_max=DEFAULT_MODE_HERWIG_LO.w2_max if args.w2_max is None else args.w2_max,
    )


def make_quark_mass_model(args: argparse.Namespace) -> QuarkMassModel:
    if args.quark_mass_mode == "none":
        return DEFAULT_QUARK_MASS_MODEL
    if args.quark_mass_mode == "massive-born":
        return DEFAULT_QUARK_MASS_MODEL_MASSIVE_BORN
    if args.quark_mass_mode == "slowrescale":
        return DEFAULT_QUARK_MASS_MODEL_SLOWRESCALE
    return DEFAULT_QUARK_MASS_MODEL_SLOWRESCALE_JACOBIAN


def make_quark_masses(args: argparse.Namespace) -> dict[str, float]:
    masses = herwig_default_quark_mass_map()
    overrides = {
        "u": args.mq_u,
        "d": args.mq_d,
        "s": args.mq_s,
        "c": args.mq_c,
        "b": args.mq_b,
    }
    for base_name, value in overrides.items():
        if value is None:
            continue
        masses[base_name] = value
        masses[f"{base_name}bar"] = value
    return masses


def make_inputs(
    args: argparse.Namespace,
) -> tuple[EWInputs, BeamInputs, LOCalculatorMode, QuarkMassModel, dict[str, float]]:
    calculator_mode = make_calculator_mode(args)
    quark_mass_model = make_quark_mass_model(args)
    quark_masses = make_quark_masses(args)
    if args.herwig_lo_mode:
        sin2_theta_w = 1.0 - (args.mw / args.mz) ** 2
        if args.theta is not None:
            theta_implied = resolve_sin2_theta_w(args.theta, None)
            if abs(theta_implied - sin2_theta_w) > 1.0e-12:
                raise SystemExit(
                    f"--theta is inconsistent with --herwig-lo-mode masses: "
                    f"sin^2(theta_W)={theta_implied:.15f} vs 1-(mW/mZ)^2={sin2_theta_w:.15f}."
                )
        if args.sin2_theta_w is not None and abs(args.sin2_theta_w - sin2_theta_w) > 1.0e-12:
            raise SystemExit(
                f"--sin2-theta-w={args.sin2_theta_w:.15f} is inconsistent with "
                f"--herwig-lo-mode masses, which imply {sin2_theta_w:.15f}."
            )
    else:
        sin2_theta_w = resolve_sin2_theta_w(args.theta, args.sin2_theta_w)
    ew_inputs = EWInputs(
        sin2_theta_w=sin2_theta_w,
        alpha_em=args.alpha_em,
        mw=args.mw,
        mz=args.mz,
        gamma_z=args.gamma_z,
        lepton_charge=args.lepton_charge,
    )
    beams = BeamInputs(
        beam_ea_gev=args.beam_ea,
        beam_eb_gev=args.beam_eb,
    )
    return ew_inputs, beams, calculator_mode, quark_mass_model, quark_masses


def render_point_report(
    args: argparse.Namespace,
    ew_inputs: EWInputs,
    beams: BeamInputs,
    calculator_mode: LOCalculatorMode,
    quark_mass_model: QuarkMassModel,
    quark_masses: dict[str, float],
    sum_pdf,
    diff_pdf,
) -> tuple[str, dict]:
    rows = point_contributions(
        args.q2,
        args.y,
        ew_inputs,
        beams,
        sum_pdf,
        diff_pdf,
        mode=calculator_mode,
        quark_masses_gev=quark_masses,
        mass_model=quark_mass_model,
    )
    channels = resolve_channel_selection(args.channel)
    filtered_rows = filter_rows_by_channel(rows, channels)
    totals = sum_point_contributions(filtered_rows)
    x_b = x_bjorken(args.q2, args.y, beams.s_hadronic)

    sections = [
        "Pointwise LO NC DIS",
        "-------------------",
        (
            f"Inputs: Q^2={args.q2:.6f} GeV^2  y={args.y:.6f}  x_B={x_b:.12f}  "
            f"s={beams.s_hadronic:.6f} GeV^2"
        ),
        f"LO convention: {calculator_mode.name}",
        f"Quark mass model: {quark_mass_model.name}",
        f"Channel selection: {', '.join(channels)}",
        "",
    ]

    if args.flavor_scope in ("total", "both"):
        total_rows = []
        for channel in channels:
            row = totals[channel]
            total_rows.append(
                [
                    channel,
                    fmt_float(row.structure_unpol),
                    fmt_float(row.structure_ll),
                    fmt_float(row.d2sigma_unpol_pb_per_gev2),
                    fmt_float(row.d2sigma_ll_pb_per_gev2),
                ]
            )
        sections.extend(
            [
                "Total PDF-weighted point contributions",
                build_table(
                    (
                        "channel",
                        "structure_unpol",
                        "structure_ll",
                        "d2sigma_unpol [pb/GeV^2]",
                        "d2sigma_LL [pb/GeV^2]",
                    ),
                    total_rows,
                ),
                "",
            ]
        )

    if args.flavor_scope in ("per-flavor", "both"):
        flavor_rows = []
        for row in filtered_rows:
            flavor_rows.append(
                [
                    row.flavor,
                    row.channel,
                    fmt_float(row.vec),
                    fmt_float(row.ax),
                    fmt_float(row.sigma_unpol_kernel),
                    fmt_float(row.sigma_ll_kernel),
                    fmt_float(row.x_pdf),
                    fmt_float(row.mass_gev),
                    fmt_float(row.mass_weight),
                    fmt_float(row.xf_unpolarized),
                    fmt_float(row.xf_polarized),
                    fmt_float(row.d2sigma_unpol_pb_per_gev2),
                    fmt_float(row.d2sigma_ll_pb_per_gev2),
                ]
            )
        sections.extend(
            [
                "Per-flavor Born kernels and differential contributions",
                build_table(
                    (
                        "flavor",
                        "channel",
                        "vec",
                        "ax",
                        "kernel_unpol",
                        "kernel_LL",
                        "x_pdf",
                        "mass [GeV]",
                        "mass_weight",
                        "xf_unpol",
                        "xf_pol",
                        "d2sigma_unpol [pb/GeV^2]",
                        "d2sigma_LL [pb/GeV^2]",
                    ),
                    flavor_rows,
                ),
                "",
            ]
        )

    audit_payload = None
    if args.audit_qid is not None:
        audit_rows = []
        for channel in channels:
            audit = point_audit_observables(
                args.q2,
                args.y,
                args.audit_qid,
                ew_inputs,
                beams,
                sum_pdf,
                diff_pdf,
                channel=channel,
                mode=calculator_mode,
                pdf_q2=args.q2,
                jacobian=args.audit_jacobian,
                s_hat=args.audit_s_hat,
                x_b=args.audit_xb,
            )
            audit_rows.append(audit)

        headers = ["channel", "qid", "pdf_sum"]
        if args.audit_jacobian is not None:
            headers.append("sigma_hat [nb]")
        headers.append("d2sigma_unpol [pb/GeV^2]")

        audit_table_rows = []
        for audit in audit_rows:
            row = [
                str(audit["channel"]),
                str(audit["qid"]),
                fmt_float(float(audit["pdf_sum"])),
            ]
            if args.audit_jacobian is not None:
                row.append(fmt_float(float(audit["sigma_hat_nb"])))
            row.append(fmt_float(float(audit["d2sigma_unpol_pb_per_gev2"])))
            audit_table_rows.append(row)

        sections.extend(
            [
                "Audit Observables",
                build_table(tuple(headers), audit_table_rows),
                "",
            ]
        )
        audit_payload = audit_rows

    payload = {
        "mode": "point",
        "inputs": {
            "q2": args.q2,
            "y": args.y,
            "x_b": x_b,
            "s_hadronic": beams.s_hadronic,
        },
        "calculator_mode": asdict(calculator_mode),
        "quark_mass_model": asdict(quark_mass_model),
        "quark_masses_gev": quark_masses,
        "audit": audit_payload,
        "totals": [
            {
                "channel": channel,
                "structure_unpol": totals[channel].structure_unpol,
                "structure_ll": totals[channel].structure_ll,
                "d2sigma_unpol_pb_per_gev2": totals[channel].d2sigma_unpol_pb_per_gev2,
                "d2sigma_ll_pb_per_gev2": totals[channel].d2sigma_ll_pb_per_gev2,
            }
            for channel in channels
        ],
        "per_flavor": serialize_dataclass_rows(filtered_rows),
    }
    return "\n".join(sections).rstrip(), payload


def render_integration_report(
    args: argparse.Namespace,
    ew_inputs: EWInputs,
    beams: BeamInputs,
    calculator_mode: LOCalculatorMode,
    quark_mass_model: QuarkMassModel,
    quark_masses: dict[str, float],
    sum_pdf,
    diff_pdf,
) -> tuple[str, dict]:
    rows = integrate_contributions(
        args.q2_min,
        args.q2_max,
        args.y_min,
        args.y_max,
        args.nq2,
        args.ny,
        ew_inputs,
        beams,
        sum_pdf,
        diff_pdf,
        mode=calculator_mode,
        quark_masses_gev=quark_masses,
        mass_model=quark_mass_model,
    )
    channels = resolve_channel_selection(args.channel)
    filtered_rows = filter_rows_by_channel(rows, channels)
    totals = sum_integrated_contributions(filtered_rows)

    sections = [
        "Integrated LO NC DIS",
        "--------------------",
        (
            f"Cuts: Q^2 in [{args.q2_min:.6f}, {args.q2_max:.6f}] GeV^2, "
            f"y in [{args.y_min:.6f}, {args.y_max:.6f}]"
        ),
        f"LO convention: {calculator_mode.name}",
        f"Quark mass model: {quark_mass_model.name}",
        f"Quadrature: nq2={args.nq2}  ny={args.ny}",
        f"Channel selection: {', '.join(channels)}",
        "",
    ]

    if args.flavor_scope in ("total", "both"):
        total_rows = []
        for channel in channels:
            row = totals[channel]
            total_rows.append(
                [
                    channel,
                    fmt_float(row.sigma_unpol_pb),
                    fmt_float(row.sigma_ll_pb),
                ]
            )
        sections.extend(
            [
                "Integrated totals",
                build_table(
                    ("channel", "sigma_unpol [pb]", "sigma_LL [pb]"),
                    total_rows,
                ),
                "",
            ]
        )

    if args.flavor_scope in ("per-flavor", "both"):
        flavor_rows = []
        for row in filtered_rows:
            flavor_rows.append(
                [
                    row.flavor,
                    row.channel,
                    fmt_float(row.sigma_unpol_pb),
                    fmt_float(row.sigma_ll_pb),
                ]
            )
        sections.extend(
            [
                "Integrated per-flavor contributions",
                build_table(
                    ("flavor", "channel", "sigma_unpol [pb]", "sigma_LL [pb]"),
                    flavor_rows,
                ),
                "",
            ]
        )

    payload = {
        "mode": "integrate",
        "cuts": {
            "q2_min": args.q2_min,
            "q2_max": args.q2_max,
            "y_min": args.y_min,
            "y_max": args.y_max,
            "nq2": args.nq2,
            "ny": args.ny,
        },
        "calculator_mode": asdict(calculator_mode),
        "quark_mass_model": asdict(quark_mass_model),
        "quark_masses_gev": quark_masses,
        "totals": [asdict(totals[channel]) for channel in channels],
        "per_flavor": serialize_dataclass_rows(filtered_rows),
    }
    return "\n".join(sections).rstrip(), payload


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if (args.audit_jacobian is None) != (args.audit_s_hat is None):
        raise SystemExit("Use --audit-jacobian and --audit-s-hat together.")
    ew_inputs, beams, calculator_mode, quark_mass_model, quark_masses = make_inputs(args)
    lhapdf = import_lhapdf()
    sum_pdf = lhapdf.mkPDF(args.sum_pdf, 0)
    diff_pdf = lhapdf.mkPDF(args.diff_pdf, 0)

    reports: List[str] = []
    payload = {
        "defaults": {
            "sum_pdf": args.sum_pdf,
            "diff_pdf": args.diff_pdf,
            "ew_inputs": asdict(ew_inputs),
            "beam_inputs": asdict(beams),
            "calculator_mode": asdict(calculator_mode),
            "quark_mass_model": asdict(quark_mass_model),
            "quark_masses_gev": quark_masses,
            "channel_selection": args.channel,
            "flavor_scope": args.flavor_scope,
        }
    }

    if args.mode in ("point", "both"):
        report, point_payload = render_point_report(
            args,
            ew_inputs,
            beams,
            calculator_mode,
            quark_mass_model,
            quark_masses,
            sum_pdf,
            diff_pdf,
        )
        reports.append(report)
        payload["point"] = point_payload

    if args.mode in ("integrate", "both"):
        report, integration_payload = render_integration_report(
            args,
            ew_inputs,
            beams,
            calculator_mode,
            quark_mass_model,
            quark_masses,
            sum_pdf,
            diff_pdf,
        )
        reports.append(report)
        payload["integrate"] = integration_payload

    print("\n\n".join(reports))

    if args.json:
        print("\nJSON")
        print("----")
        print(json.dumps(payload, indent=2, sort_keys=True))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
