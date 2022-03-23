from .test_findorb import *
from .raiden import *
from .wrapper import *
#from .eph2ades import *

import sys
import os.path
import shutil
import subprocess

# Check if 'fo' is on the path; if it isn't assume it's
# in the same directory as python (which is true for conda-based installs)
# and prepend that directory to the PATH.
#
# NOTE: we _must_ prepend to the PATH (however dirty it may be), because
# 'fo' behaves differently if called with a full path vs. if called as
# just 'fo'.
#
if shutil.which("fo") is None:
    fo_dir = os.path.split(sys.executable)[0]
    paths = os.environ["PATH"].split(":")
    if fo_dir not in paths:
        os.environ["PATH"] = fo_dir + ":" + os.environ["PATH"]

if shutil.which("fo") is None:
    raise Exception("FATAL: FindOrb not found. FindOrb's `fo` executable must be on $PATH")

subprocess.call(["fo"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
