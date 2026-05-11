from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from .vectors import FourVector


ME_STATUS_ROLE2 = 990072
ME_STATUS_ROLE3 = 990073


@dataclass(frozen=True)
class EventParticle:
    pid: int
    momentum: FourVector
    status: int
    role: int | None = None
    barcode: int | None = None


@dataclass(frozen=True)
class HepMCEvent:
    event_number: int
    weight: float
    incoming: tuple[EventParticle, ...]
    outgoing: tuple[EventParticle, ...]
    attributes: dict[str, str] = field(default_factory=dict)
    weights: dict[str, float] = field(default_factory=dict)

    @property
    def particles(self) -> tuple[EventParticle, ...]:
        return self.incoming + self.outgoing

    def with_weights(self, weights: dict[str, float]) -> "HepMCEvent":
        return HepMCEvent(
            event_number=self.event_number,
            weight=self.weight,
            incoming=self.incoming,
            outgoing=self.outgoing,
            attributes=dict(self.attributes),
            weights=dict(weights),
        )


def _format_float(value: float) -> str:
    return f"{value:.12e}"


def _particle_line(barcode: int, vertex: int, particle: EventParticle) -> str:
    p = particle.momentum
    return (
        f"P {barcode} {vertex} {particle.pid} "
        f"{_format_float(p.px)} {_format_float(p.py)} {_format_float(p.pz)} "
        f"{_format_float(p.e)} {_format_float(p.mass())} {particle.status}"
    )


def _cross_section_line(cross_section_pb: float) -> str:
    values = [_format_float(cross_section_pb), _format_float(0.0), "-1", "-1"]
    return "A 0 GenCrossSection " + " ".join(values)


def write_hepmc3(path: str | Path, events: Iterable[HepMCEvent], cross_section_pb: float | None = None) -> float:
    event_list = list(events)
    if cross_section_pb is None:
        cross_section_pb = sum(event.weight for event in event_list)

    weight_names: list[str] = []
    for event in event_list:
        for name in event.weights:
            if name not in weight_names:
                weight_names.append(name)

    lines = [
        "HepMC::Version 3.03.00",
        "HepMC::Asciiv3-START_EVENT_LISTING",
        "W " + " ".join(["Weight", *weight_names]),
        "T HerwigFO-Python\\|0.1\\|Standalone fixed-order NC DIS validator",
    ]

    for event in event_list:
        particle_count = len(event.particles)
        lines.append(f"E {event.event_number} 1 {particle_count}")
        lines.append("U GEV MM")
        weight_values = [event.weight]
        weight_values.extend(float(event.weights.get(name, 0.0)) for name in weight_names)
        lines.append("W " + " ".join(_format_float(value) for value in weight_values))
        lines.append(_cross_section_line(cross_section_pb))

        next_barcode = 1
        incoming_barcodes: list[tuple[int, EventParticle]] = []
        for particle in event.incoming:
            barcode = particle.barcode or next_barcode
            next_barcode = max(next_barcode, barcode + 1)
            incoming_barcodes.append((barcode, particle))

        outgoing_barcodes: list[tuple[int, EventParticle]] = []
        for particle in event.outgoing:
            barcode = particle.barcode or next_barcode
            next_barcode = max(next_barcode, barcode + 1)
            outgoing_barcodes.append((barcode, particle))

        for key, value in sorted(event.attributes.items()):
            lines.append(f"A 0 {key} {value}")
        for barcode, particle in outgoing_barcodes:
            if particle.role is not None:
                lines.append(f"A {barcode} herwig_powheg_me_parton {particle.role}")

        for barcode, particle in incoming_barcodes:
            lines.append(_particle_line(barcode, 0, particle))

        incoming_ids = ",".join(str(barcode) for barcode, _ in incoming_barcodes)
        lines.append(f"V -1 0 [{incoming_ids}]")

        for barcode, particle in outgoing_barcodes:
            lines.append(_particle_line(barcode, -1, particle))

    lines.append("HepMC::Asciiv3-END_EVENT_LISTING")
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n")
    return float(cross_section_pb)
