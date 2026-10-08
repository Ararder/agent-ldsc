#!/usr/bin/env bash
# Build the synthetic bundle and run one request through the worker with the given runtime.
# Usage: ci-fixture-run.sh docker|apptainer IMAGE_OR_SIF OUTDIR
set -euo pipefail
runtime=$1 image=$2 out=$3
mkdir -p "$out/run/input"
cp tests/fixtures/requests/fixture-request.json "$out/run/input/request.json"
cp tests/fixtures/requests/peaks.bed tests/fixtures/requests/genes.txt "$out/run/input/"
if [ "$runtime" = docker ]; then
  run() { docker run --rm --network none --user "$(id -u):$(id -g)" -v "$PWD:/src:ro" -v "$out:/out" --entrypoint "$1" "$image" "${@:2}"; }
else
  run() { apptainer exec --cleanenv --containall --no-home --bind "$PWD:/src:ro,$out:/out" "$image" "$@"; }
fi
run /opt/envs/worker/bin/python -I /src/tests/fixtures/make_fixture_bundle.py /out/bundle
run agent-ldsc-worker run --request /out/run/input/request.json --refs /out/bundle --work /out/run --jobs 2
test -s "$out/run/COMPLETE.json"
