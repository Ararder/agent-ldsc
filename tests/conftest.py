import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("AGENT_LDSC_HOME", str(ROOT))
# Local development fallbacks (the image provides these at their default locations).
_dev_ldsc = Path("/tmp/agent-ldsc-tools/ldsc27/bin/python")
if "AGENT_LDSC_LDSC_CMD" not in os.environ and _dev_ldsc.exists():
    os.environ["AGENT_LDSC_LDSC_CMD"] = f"arch -x86_64 {_dev_ldsc} /tmp/agent-ldsc-src/ldsc/ldsc.py"
_dev_liftover = Path("/tmp/agent-ldsc-tools/liftover/bin")
if shutil.which("liftOver") is None and _dev_liftover.exists():
    os.environ["PATH"] = f"{_dev_liftover}:{os.environ['PATH']}"


def _ldsc_available() -> bool:
    cmd = os.environ.get("AGENT_LDSC_LDSC_CMD", "/opt/envs/ldsc/bin/python /opt/ldsc/ldsc.py").split()
    return Path(cmd[-1]).exists()


def _rgate_available() -> bool:
    if shutil.which("Rscript") is None:
        return False
    return subprocess.run(["Rscript", "-e", "stopifnot(requireNamespace('ldsR', quietly=TRUE))"],
                          capture_output=True).returncode == 0


def pytest_collection_modifyitems(config, items):
    checks = {"ldsc": _ldsc_available(), "liftover": shutil.which("liftOver") is not None,
              "rgate": _rgate_available()}
    for item in items:
        for mark, ok in checks.items():
            if mark in item.keywords and not ok:
                if os.environ.get("AGENT_LDSC_REQUIRE_ALL_TESTS"):
                    pytest.fail(f"{mark} unavailable but AGENT_LDSC_REQUIRE_ALL_TESTS is set")
                item.add_marker(pytest.mark.skip(reason=f"{mark} unavailable"))


@pytest.fixture(scope="session")
def fixture_bundle(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("bundle") / "fixture-grch37-v1"
    subprocess.run([sys.executable, str(ROOT / "tests/fixtures/make_fixture_bundle.py"), str(out)], check=True)
    return out


@pytest.fixture
def bundle(fixture_bundle):
    from agent_ldsc_worker.refs import load_bundle
    return load_bundle(fixture_bundle)


@pytest.fixture
def strict_profile():
    from agent_ldsc_worker.run import load_profile
    return load_profile("strict-v1")
