import sys
sys.stdout.write("start\n")
sys.stdout.flush()
try:
    import numpy as np
    sys.stdout.write("numpy ok\n")
    sys.stdout.flush()
except Exception as e:
    sys.stdout.write(f"numpy error: {e}\n")
    sys.stdout.flush()
try:
    import pandas as pd
    sys.stdout.write("pandas ok\n")
    sys.stdout.flush()
except Exception as e:
    sys.stdout.write(f"pandas error: {e}\n")
    sys.stdout.flush()
try:
    from scipy import stats
    sys.stdout.write("scipy ok\n")
    sys.stdout.flush()
except Exception as e:
    sys.stdout.write(f"scipy error: {e}\n")
    sys.stdout.flush()
sys.stdout.write("end\n")
sys.stdout.flush()
