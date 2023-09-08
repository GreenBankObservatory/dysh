from bs4 import BeautifulSoup
import os, wget, socket, requests, subprocess
from getpass import getuser
from pathlib import Path
import platform
import glob

# Get system information
USER_NAME = getuser()
HOST_NAME = socket.gethostname()
try:    
    HOST_DOMAIN = socket.gethostbyname_ex(HOST_NAME)[0]
    HOST_IS_GBO = HOST_DOMAIN.endswith("gb.nrao.edu")
except:
    HOST_IS_GBO = False
HOST_OS = platform.uname()[0]
print(USER_NAME, HOST_NAME, HOST_OS)

# Collect paths
BASE_DIR = Path(__file__).resolve().parent.parent
GUI_BASE_DIR = Path(__file__).resolve().parent

if HOST_IS_GBO:
    DYSH_PRIMARY_DATA_ROOT = "/home/gbtdata"
    DYSH_EXAMPLE_DATA_ROOT = "https://www.gb.nrao.edu/dysh/example_data"
    sessions = [os.path.basename(x) for x in glob.glob(DYSH_PRIMARY_DATA_ROOT+'/*')]
    print(sessions)
else:
    DYSH_PRIMARY_DATA_ROOT=BASE_DIR
    DYSH_EXAMPLE_DATA_ROOT = "https://www.gb.nrao.edu/dysh/example_data"


url = "https://www.gb.nrao.edu/dysh/example_data/azel/data/AGBT16B_018_01.raw.vegas"
ext = 'fits'

def listFD(url, ext=''):
    page = requests.get(url).text
    print(page)
    soup = BeautifulSoup(page, 'html.parser')
    return [node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]

for file in listFD(DYSH_EXAMPLE_DATA_ROOT, ext):
    print(file)
print(glob.glob(DYSH_EXAMPLE_DATA_ROOT+'/*'))

subprocess.Popen("mkdir ~/mnt/data_dir mount -t data:/dir/ /mnt/data_dir", shell=True)