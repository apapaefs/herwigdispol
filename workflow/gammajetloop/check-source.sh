#!/usr/bin/env bash
# Load the compatible Herwig/ThePEG environment before running this script.
set -euo pipefail
task_source=${1:?usage: check-source.sh HERWIG_SOURCE NEW_BUILD_DIRECTORY}
task_build=${2:?provide a new build directory}
if [[ -e "$task_build" ]]; then
  echo "Refusing to overwrite existing directory: $task_build" >&2
  exit 1
fi
mkdir -p "$task_build"
task_build=$(cd "$task_build" && pwd)
task_helpers=$(cd "$(dirname "$0")" && pwd)
task_hadron="$task_source/MatrixElement/Hadron"
read -r -a task_cxx <<< "${CXX:-c++}"
read -r -a task_cppflags <<< "$(herwig-config --cppflags)"
read -r -a task_ldflags <<< "$(herwig-config --ldflags)"
read -r -a task_ldlibs <<< "$(herwig-config --ldlibs)"
task_rpaths=()
for task_flag in "${task_ldflags[@]}"; do
  if [[ "$task_flag" == -L* ]]; then
    task_rpaths+=("-Wl,-rpath,${task_flag#-L}")
  fi
done
"${task_cxx[@]}" -std=c++11 -O2 -UNDEBUG \
  "$task_hadron/GammaJetLoopAmplitudesTest.cc" \
  "$task_hadron/GammaJetLoopAmplitudes.cc" -o "$task_build/amplitudes-test"
"$task_build/amplitudes-test"
"${task_cxx[@]}" -std=c++11 -O2 -I"$task_hadron" \
  "$task_helpers/point_probe.cc" "$task_hadron/GammaJetLoopAmplitudes.cc" \
  -o "$task_build/point-probe"
"${task_cxx[@]}" -std=c++17 -O2 -UNDEBUG "${task_cppflags[@]}" \
  "$task_helpers/test_hard_vertex.cc" \
  "$task_source/MatrixElement/ProductionMatrixElement.cc" \
  "$task_source/MatrixElement/HardVertex.cc" \
  "${task_ldflags[@]}" "${task_ldlibs[@]}" "${task_rpaths[@]}" \
  -o "$task_build/hard-vertex-test"
"$task_build/hard-vertex-test"
