"""Audit exact installed dependencies, excluding only this local project."""

import subprocess
import sys
from importlib.metadata import distributions
from pathlib import Path

root = Path(__file__).resolve().parents[1]
requirements = root / "build" / "audit-requirements.txt"
requirements.parent.mkdir(exist_ok=True)
packages = {
    distribution.metadata["Name"]: distribution.version
    for distribution in distributions()
    if distribution.metadata["Name"].lower() != "worstcase"
}
requirements.write_text(
    "".join(f"{name}=={version}\n" for name, version in sorted(packages.items())),
    encoding="utf-8",
)
subprocess.run(
    [
        sys.executable,
        "-m",
        "pip_audit",
        "--requirement",
        str(requirements),
        "--no-deps",
        "--disable-pip",
        "--strict",
        "--progress-spinner",
        "off",
    ],
    check=True,
)
