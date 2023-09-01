from getpass import getuser
from pathlib import Path
import os, socket
import platform
import glob

# Get system information
USER_NAME = getuser()
HOST_NAME = socket.gethostname()
HOST_DOMAIN = socket.gethostbyname_ex(HOST_NAME)[0]
HOST_IS_GBO = HOST_DOMAIN.endswith("gb.nrao.edu")
HOST_OS = platform.uname()[0]

# Collect paths
BASE_DIR = Path(__file__).resolve().parent.parent
GUI_BASE_DIR = Path(__file__).resolve().parent

if HOST_IS_GBO:
    DYSH_DATA_ROOT="/home/gbtdata"
else:
    DYSH_DATA_ROOT=BASE_DIR

sessions = [os.path.basename(x) for x in glob.glob(DYSH_DATA_ROOT+'/*')]
print(sessions)