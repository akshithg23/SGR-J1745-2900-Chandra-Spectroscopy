"""Run CIAO specextract with a temporary fix for SciServer's Python 3.11 build."""

import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main():
    command = ["specextract", *sys.argv[1:]]
    spec = importlib.util.find_spec("ciao_contrib")
    package = Path(next(iter(spec.submodule_search_locations)))
    relative = Path("_tools/merging.py")
    source = package / relative
    code = source.read_text(encoding="utf-8")
    try:
        compile(code, str(source), "exec")
    except SyntaxError as error:
        lines = code.splitlines(keepends=True)
        index = error.lineno - 1
        original = lines[index]
        if "ObsIds" not in original or '{" ".join(' not in original:
            raise
        lines[index] = original.replace('{" ".join(', "{' '.join(", 1)
        fixed = "".join(lines)
        compile(fixed, str(source), "exec")
    else:
        return subprocess.call(command)

    print("Applying temporary CIAO Python 3.11 quoting fix.", flush=True)
    with tempfile.TemporaryDirectory(prefix="spa660-ciao-") as tmp:
        local = Path(tmp) / "ciao_contrib"
        shutil.copytree(package, local, ignore=shutil.ignore_patterns("__pycache__"))
        (local / relative).write_text(fixed, encoding="utf-8")
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join(filter(None, [tmp, env.get("PYTHONPATH")]))
        return subprocess.call(command, env=env)


if __name__ == "__main__":
    sys.exit(main())
