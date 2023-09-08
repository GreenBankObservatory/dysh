import os

HOME = os.getenv('HOME')

XDG_CONFIG_HOME = os.getenv('XDG_CONFIG_HOME')
if XDG_CONFIG_HOME is None:
    XDG_CONFIG_HOME = os.path.join(HOME, '.config/')
    os.environ['XDG_CONFIG_HOME'] = XDG_CONFIG_HOME
    print(f"No XDG_CONFIG_HOME found. Setting to {XDG_CONFIG_HOME}")

DYSH_CONFIG = os.getenv('DYSH_CONFIG')
if DYSH_CONFIG is None:
    DYSH_CONFIG = os.path.join(XDG_CONFIG_HOME, 'dysh/')
    os.environ['DYSH_CONFIG'] = DYSH_CONFIG
    print(f"No DYSH_CONFIG found. Setting to {DYSH_CONFIG}")

print(HOME, XDG_CONFIG_HOME, DYSH_CONFIG)