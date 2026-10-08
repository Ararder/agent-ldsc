#!/usr/bin/env bash
# Run the full test suite inside an image. Tests that need ldsc/liftOver/ldsR must run here.
set -euo pipefail
image=$1
# The launcher is host-side and not in the image; it is mounted only for its unit tests.
docker run --rm --network none -v "$PWD/tests:/src/tests:ro" -v "$PWD/launcher:/src/launcher:ro" \
  -v "$PWD/profiles:/src/profiles:ro" -e AGENT_LDSC_REQUIRE_ALL_TESTS=1 \
  --entrypoint /bin/sh "$image" -c '
    set -e
    ldsc --help >/dev/null
    agent-ldsc-worker versions
    cd /tmp && python -m pytest -q -p no:cacheprovider -o "pythonpath=/src/tests /src/launcher" --rootdir /src/tests /src/tests'
