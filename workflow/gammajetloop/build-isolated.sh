#!/usr/bin/env bash
# Source a compatible polarized Herwig environment before running this script.
set -euo pipefail
task_source=${1:?usage: build-isolated.sh HERWIG_SOURCE NEW_BUILD_DIRECTORY}
task_build=${2:?provide a new build directory}
if [[ -e "$task_build" ]]; then
  echo "Refusing to overwrite existing directory: $task_build" >&2
  exit 1
fi
mkdir -p "$task_build/src" "$task_build/lib"
task_build=$(cd "$task_build" && pwd)
# Herwig 7.3's shell helper prints valid output but can exit 1 because its
# final option test is false. Validate the returned installation explicitly.
task_prefix=$(herwig-config --prefix || true)
test -f "$task_prefix/lib/Herwig/Herwig.so"
test -x "$task_prefix/bin/Herwig"
task_hadron="$task_source/MatrixElement/Hadron"
cp "$task_hadron"/{MEPP2GammaJetLoop,GammaJetLoopAmplitudes}.{cc,h} "$task_build/src/"
cp "$task_hadron/GammaJetLoopAmplitudesTest.cc" "$task_build/src/"
read -r -a task_cxx <<< "${CXX:-c++}"
# herwig-config returns compiler arguments rather than a single path.
read -r -a task_cppflags <<< "$(herwig-config --cppflags)"
read -r -a task_ldflags <<< "$(herwig-config --ldflags)"
read -r -a task_ldlibs <<< "$(herwig-config --ldlibs)"
for source in MEPP2GammaJetLoop GammaJetLoopAmplitudes; do
  "${task_cxx[@]}" -std=c++17 -O2 -fPIC "${task_cppflags[@]}" \
    -c "$task_build/src/$source.cc" -o "$task_build/$source.o"
done
"${task_cxx[@]}" -shared -o "$task_build/lib/HwGammaJetLoopValidation.so" \
  "$task_build/MEPP2GammaJetLoop.o" "$task_build/GammaJetLoopAmplitudes.o" \
  "$task_prefix/lib/Herwig/Herwig.so" "${task_ldflags[@]}" "${task_ldlibs[@]}" \
  -Wl,-rpath,"$task_prefix/lib/Herwig" -Wl,-rpath,"$task_prefix/lib/ThePEG"
"${task_cxx[@]}" -std=c++17 -O2 "$task_build/src/GammaJetLoopAmplitudesTest.cc" \
  "$task_build/src/GammaJetLoopAmplitudes.cc" -o "$task_build/amplitudes-test"
"$task_build/amplitudes-test"
cat > "$task_build/herwig-loop" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export LD_PRELOAD="$task_build/lib/HwGammaJetLoopValidation.so\${LD_PRELOAD:+:\$LD_PRELOAD}"
exec "$task_prefix/bin/Herwig" "\$@" -L "$task_build/lib"
EOF
chmod +x "$task_build/herwig-loop"
sha256sum "$task_build/src/"* "$task_build/lib/HwGammaJetLoopValidation.so" \
  > "$task_build/SHA256SUMS"
printf 'Isolated launcher: %s/herwig-loop\n' "$task_build"
