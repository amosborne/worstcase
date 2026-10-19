"""Install each built artifact into an isolated environment with runtime deps only."""

import subprocess
import sys
import tempfile
import venv
from pathlib import Path

root = Path(__file__).resolve().parents[1]
artifacts = sorted((root / "dist").glob("worstcase-*.whl")) + sorted(
    (root / "dist").glob("worstcase-*.tar.gz")
)
if not artifacts:
    raise SystemExit("Run poetry build first")
(root / "build").mkdir(exist_ok=True)
for artifact in artifacts:
    with tempfile.TemporaryDirectory(dir=root / "build", prefix="artifact-") as directory:
        target = Path(directory)
        venv.EnvBuilder(with_pip=True).create(target / "venv")
        scripts = "Scripts" if sys.platform == "win32" else "bin"
        python = target / "venv" / scripts / "python"
        subprocess.run(
            [str(python), "-m", "pip", "install", str(artifact)], check=True
        )
        subprocess.run(
            [str(python), str(root / "scripts" / "check_core.py")],
            cwd=target, check=True,
        )
        print(f"Validated {artifact.name}", flush=True)
