#!/usr/bin/env python3
"""Relocate only a saved test's loader paths; preserve the serialized object payload.

Validation only. Production must regenerate its run files. This avoids testing
the old library accidentally when checking version-0 object compatibility.
"""
import argparse
import hashlib
import json
from pathlib import Path


def relocate(data, old_prefix, new_prefix):
    lines = data.splitlines(keepends=True)
    if lines[:3] != [b"ThePEG version 1 Database\n", b"0\n", b"3\n"]:
        raise ValueError("unrecognized saved-file header; refusing to modify it")
    old, new = old_prefix.encode(), new_prefix.encode()
    cursor = 3
    changed = 0
    for _ in range(2):
        count = int(lines[cursor])
        cursor += 1
        for index in range(cursor,cursor+count):
            if lines[index].startswith(old+b"/"):
                lines[index] = new+lines[index][len(old):]
                changed += 1
        cursor += count
    if not changed:
        raise ValueError("no loader paths matched the requested old prefix")
    payload = b"".join(lines[cursor:])
    original_payload = b"".join(data.splitlines(keepends=True)[cursor:])
    assert payload == original_payload
    return b"".join(lines), {"changed_loader_paths":changed,
        "unchanged_object_payload_sha256":hashlib.sha256(payload).hexdigest()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",type=Path)
    parser.add_argument("destination",type=Path)
    parser.add_argument("--old-prefix",required=True)
    parser.add_argument("--new-prefix",required=True)
    args=parser.parse_args()
    data, report=relocate(args.source.read_bytes(),args.old_prefix,args.new_prefix)
    with args.destination.open("xb") as stream:
        stream.write(data)
    print(json.dumps(report))


if __name__ == "__main__":
    main()
