import sys
print("Python:", sys.version)

import importlib
for pkg in ["scipy", "numpy", "pandas", "matplotlib", "pymc", "baycomp", "statsmodels"]:
    try:
        spec = importlib.util.find_spec(pkg)
        if spec is None:
            print(pkg + ": NOT FOUND")
        else:
            m = importlib.import_module(pkg)
            v = getattr(m, "__version__", "?")
            print(pkg + ": " + str(v))
    except Exception as e:
        print(pkg + ": ERROR " + str(e))

print("DONE")
