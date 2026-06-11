#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path


def replace_interface_filename(text: str, class_name: str, member_name: str) -> tuple[str, int]:
    lines = text.splitlines()
    start = end = None
    for i, line in enumerate(lines):
        if 'interfaceFilename' in line and class_name in ''.join(lines[max(0, i-1):min(len(lines), i+2)]):
            start = i
            break
    if start is None:
        for i, line in enumerate(lines):
            if 'interfaceFilename' in line:
                start = i
                break
    if start is None:
        return text, 0

    for j in range(start, min(len(lines), start + 30)):
        if f'&{class_name}::{member_name}, "", true, false);' in lines[j]:
            end = j
            break
    if end is None:
        for j in range(start, min(len(lines), start + 30)):
            if lines[j].strip().endswith('"", true, false);'):
                end = j
                break
    if end is None:
        return text, 0

    prefix = lines[start]
    repl = [
        prefix,
        '    ("Filename",',
        '     "The name of the file where the YODA histograms are put. If empty, "',
        '     "the run name will be used instead. \' .yoda\' will in any case be "',
        '     "appended to the file name.",',
        f'     &{class_name}::{member_name}, "", true, false);',
    ]
    repl[3] = '     "the run name will be used instead. \'\\.yoda\' will in any case be "'
    lines = lines[:start] + repl + lines[end + 1 :]
    return '\n'.join(lines) + '\n', 1


def subn(text: str, pattern: str, repl: str, flags: int = 0) -> tuple[str, int]:
    return re.subn(pattern, repl, text, flags=flags)


def patch_rivet_analysis(text: str) -> tuple[str, int]:
    total = 0

    text, n = subn(
        text,
        r'if \( _rivet \)\{\s*#if ThePEG_RIVET_VERSION > 1\s*try \{\s*_rivet->analyze\(\*hepmc\);\s*\} catch \(const YODA::Exception & e\) \{\s*Throw\(\) << "Warning: Rivet/Yoda got the exception: "<< e\.what\(\)<<"\\n"\s*<< Exception::warning;\s*\}\s*#else\s*#error "Unknown ThePEG_RIVET_VERSION"\s*#endif\s*\}',
        'if ( _rivet ){\n'
        '    try {\n'
        '      _rivet->analyze(*hepmc);\n'
        '    } catch (const YODA::Exception & e) {\n'
        '      Throw() << "Warning: Rivet/Yoda got the exception: "<< e.what()<<"\\n"\n'
        '              << Exception::warning;\n'
        '    }\n'
        '  }',
        flags=re.S,
    )
    total += n

    text2, n = replace_interface_filename(text, 'RivetAnalysis', '_filename')
    text = text2
    total += n

    text, n = subn(
        text,
        r'#if ThePEG_RIVET_VERSION > 2\s*_rivet->setCrossSection\(make_pair\(generator\(\)->integratedXSec\(\)/picobarn,\s*generator\(\)->integratedXSecErr\(\)/picobarn\)\);\s*#else\s*_rivet->setCrossSection\(generator\(\)->integratedXSec\(\)/picobarn\);\s*#endif',
        '_rivet->setCrossSection(double(generator()->integratedXSec()/picobarn),\n'
        '                            double(generator()->integratedXSecErr()/picobarn));',
        flags=re.S,
    )
    total += n

    text, n = subn(
        text,
        r'string fname = _filename;\s*#if ThePEG_RIVET_VERSION > 1\s*if \( fname.empty\(\) \) fname = generator\(\)->path\(\) \+ "/" \+ generator\(\)->runName\(\) \+ "\\.yoda";\s*#else\s*#error "Unknown ThePEG_RIVET_VERSION"\s*#endif',
        'string fname = _filename;\n'
        '    if ( fname.empty() ) fname = generator()->path() + "/" + generator()->runName() + ".yoda";',
        flags=re.S,
    )
    total += n

    text, n = subn(
        text,
        r'#if ThePEG_RIVET_VERSION > 3\s*_rivet->setCheckBeams\(_checkBeams\);\s*#else\s*_rivet->checkBeams\(_checkBeams\);\s*#endif',
        '_rivet->setCheckBeams(_checkBeams);',
        flags=re.S,
    )
    total += n

    text, n = subn(text, r'_rivet->checkBeams\(_checkBeams\);', '_rivet->setCheckBeams(_checkBeams);')
    total += n

    text, n = subn(text, r'\n\s*#error "Unknown ThePEG_RIVET_VERSION"', '')
    total += n

    return text, total


def patch_nlo_rivet_analysis(text: str) -> tuple[str, int]:
    total = 0

    text, n = subn(
        text,
        r'if\(_rivet\)\{\s*#if ThePEG_RIVET_VERSION == 1\s*_rivet->analyze\(\*hepmc\);\s*#elif ThePEG_RIVET_VERSION > 1\s*try \{\s*_rivet->analyze\(\*hepmc\);\s*\} catch \(const YODA::Exception & e\) \{\s*Throw\(\) << "Warning: Rivet/Yoda got the exception: "<< e\.what\(\)<<"\\n"\s*<< Exception::warning;\s*\}\s*#else\s*#error "Unknown ThePEG_RIVET_VERSION"\s*#endif\s*\}',
        'if(_rivet){\n'
        '    try {\n'
        '      _rivet->analyze(*hepmc);\n'
        '    } catch (const YODA::Exception & e) {\n'
        '      Throw() << "Warning: Rivet/Yoda got the exception: "<< e.what()<<"\\n"\n'
        '              << Exception::warning;\n'
        '    }\n'
        '  }',
        flags=re.S,
    )
    total += n

    text, n = subn(
        text,
        r'if \( _rivet \)\{\s*#if ThePEG_RIVET_VERSION == 1\s*_rivet->analyze\(\*hepmc\);\s*#elif ThePEG_RIVET_VERSION > 1\s*try \{\s*_rivet->analyze\(\*hepmc\);\s*\} catch \(const YODA::Exception & e\) \{\s*Throw\(\) << "Warning: Rivet/Yoda got the exception: "<< e\.what\(\)<<"\\n"\s*<< Exception::warning;\s*\}\s*#else\s*#error "Unknown ThePEG_RIVET_VERSION"\s*#endif\s*\}',
        'if ( _rivet ){\n'
        '        try {\n'
        '          _rivet->analyze(*hepmc);\n'
        '        } catch (const YODA::Exception & e) {\n'
        '          Throw() << "Warning: Rivet/Yoda got the exception: "<< e.what()<<"\\n"\n'
        '                  << Exception::warning;\n'
        '        }\n'
        '      }',
        flags=re.S,
    )
    total += n

    text2, n = replace_interface_filename(text, 'NLORivetAnalysis', 'filename')
    text = text2
    total += n

    text, n = subn(
        text,
        r'#if ThePEG_RIVET_VERSION > 2\s*_rivet->setCrossSection\(make_pair\(generator\(\)->integratedXSec\(\)/picobarn,\s*generator\(\)->integratedXSecErr\(\)/picobarn\)\);\s*#else\s*_rivet->setCrossSection\(generator\(\)->integratedXSec\(\)/picobarn\);\s*#endif',
        '_rivet->setCrossSection(double(generator()->integratedXSec()/picobarn),\n'
        '                            double(generator()->integratedXSecErr()/picobarn));',
        flags=re.S,
    )
    total += n

    text, n = subn(
        text,
        r'string fname = filename;\s*#if ThePEG_RIVET_VERSION == 1\s*if \( fname.empty\(\) \) fname = generator\(\)->path\(\) \+ "/" \+ generator\(\)->runName\(\) \+ "\\.aida";\s*#elif ThePEG_RIVET_VERSION > 1\s*if \( fname.empty\(\) \) fname = generator\(\)->path\(\) \+ "/" \+ generator\(\)->runName\(\) \+ "\\.yoda";\s*#else\s*#error "Unknown ThePEG_RIVET_VERSION"\s*#endif',
        'string fname = filename;\n'
        '    if ( fname.empty() ) fname = generator()->path() + "/" + generator()->runName() + ".yoda";',
        flags=re.S,
    )
    total += n

    text, n = subn(text, r'\n\s*#error "Unknown ThePEG_RIVET_VERSION"', '')
    total += n

    return text, total


def main() -> int:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('.')
    ra = root / 'Analysis' / 'RivetAnalysis.cc'
    nra = root / 'Analysis' / 'NLORivetAnalysis.cc'

    if not ra.exists() or not nra.exists():
        print('Usage: python3 thepeg_rivet4_patch.py /path/to/ThePEG-2.3.0', file=sys.stderr)
        print('Could not find Analysis/RivetAnalysis.cc and Analysis/NLORivetAnalysis.cc', file=sys.stderr)
        return 2

    text = ra.read_text()
    text2, n1 = patch_rivet_analysis(text)
    ra.write_text(text2)

    text = nra.read_text()
    text2, n2 = patch_nlo_rivet_analysis(text)
    nra.write_text(text2)

    print(f'Patched {ra} ({n1} edits)')
    print(f'Patched {nra} ({n2} edits)')
    if n1 == 0 and n2 == 0:
        print('No edits were applied. The files may already be patched or may differ from expected ThePEG 2.3.0 sources.', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
