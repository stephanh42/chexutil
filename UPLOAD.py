from subprocess import check_call
from sys import executable as py_exec
from glob import glob

check_call("git clean -x -f -d".split())

for out_fmt, ext in (("plain", ""), ("rst", ".rst")):
    check_call(("pandoc", "--from=markdown", "--to=" + out_fmt, "--output=README" + ext, "README.md"))

check_call((py_exec, "setup.py", "sdist"))
check_call((py_exec, "setup.py", "bdist_wheel"))

print("""1. Upload to testpypi
2. Upload to pypi
9. Exit""")

repos = {"1": "testpypi", "2": "pypi"}

choice = repos.get(input("Make your choice: ").strip())
if choice is not None:
    args = (py_exec, "-m", "twine", "upload", "-r", choice) + tuple(glob("dist/*"))
    print(" ".join(args))
    check_call(args)
