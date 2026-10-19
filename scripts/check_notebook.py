"""Execute the notebook with this Poetry environment's interpreter."""

import json
import sys
import tempfile
from pathlib import Path

import nbformat
from nbclient import NotebookClient

root = Path(__file__).resolve().parents[1]
(root / "build").mkdir(exist_ok=True)
with tempfile.TemporaryDirectory(prefix="kernel-", dir=root / "build") as directory:
    kernel = Path(directory) / "python3"
    kernel.mkdir()
    (kernel / "kernel.json").write_text(
        json.dumps(
            {
                "argv": [
                    sys.executable,
                    "-m",
                    "ipykernel_launcher",
                    "-f",
                    "{connection_file}",
                ],
                "display_name": "worstcase validation",
                "language": "python",
            }
        ),
        encoding="utf-8",
    )
    from jupyter_client.kernelspec import KernelSpecManager

    notebook = nbformat.read(root / "example.ipynb", as_version=4)
    client = NotebookClient(
        notebook,
        timeout=180,
        kernel_name="python3",
        resources={"metadata": {"path": str(root)}},
    )
    client.create_kernel_manager()
    client.km.kernel_spec_manager = KernelSpecManager(kernel_dirs=[directory])
    client.execute()
    (root / "build").mkdir(exist_ok=True)
    nbformat.write(notebook, root / "build" / "example.executed.ipynb")
print("Notebook executed successfully")
