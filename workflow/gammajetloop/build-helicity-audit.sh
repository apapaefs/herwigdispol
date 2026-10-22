#!/usr/bin/env bash
# Only the diagnostic observer is built here; the native ME comes from Herwig.
set -euo pipefail
task_source=${1:?usage: build-helicity-audit.sh HERWIG_SOURCE NEW_BUILD_DIRECTORY}
task_build=${2:?provide a new build directory}
if [[ -e "$task_build" ]]; then
  echo "Refusing to overwrite existing directory: $task_build" >&2
  exit 1
fi
mkdir -p "$task_build"
task_helpers=$(cd "$(dirname "$0")" && pwd)
task_prefix=$(herwig-config --prefix || true)
test -f "$task_prefix/lib/Herwig/Herwig.so"
read -r -a task_cxx <<< "${CXX:-c++}"
read -r -a task_cppflags <<< "$(herwig-config --cppflags)"
read -r -a task_ldflags <<< "$(herwig-config --ldflags)"
read -r -a task_ldlibs <<< "$(herwig-config --ldlibs)"
"${task_cxx[@]}" -std=c++17 -O2 -fPIC -shared \
  "${task_cppflags[@]}" -I"$task_source/MatrixElement/Hadron" \
  "$task_helpers/NativeHelicityAudit.cc" \
  "$task_prefix/lib/Herwig/Herwig.so" "${task_ldflags[@]}" "${task_ldlibs[@]}" \
  -Wl,-rpath,"$task_prefix/lib/Herwig" -Wl,-rpath,"$task_prefix/lib/ThePEG" \
  -o "$task_build/NativeHelicityAudit.so"
