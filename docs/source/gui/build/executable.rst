***********************
Building the Executable
***********************

PyInstaller
===========

To build with PyInstaller

.. code:: bash

    (dysh) $ cd gui
    (dysh) $ pyinstaller app.py

Building From Scratch
#####################

If you somehow lose or clear the ``dysh.spec`` file, you need to run a lot more flags through `PyInstaller`

.. code:: bash

    (dysh) $ cd gui
    (dysh) $ pyinstaller --onefile --noconsole --name "dysh" --icon=./static/favicon.ico  --clean -y --collect-all asdf --collect-all asdf_standard --collect-all asdf_transform_schemas --collect-all packaging --collect-all pkg_resources --collect-all astropy --collect-all lz4 --recursive-copy-metadata asdf --recursive-copy-metadata astropy app.py

Troubleshooting
===============

Windows
#######

If you get the following error:

.. code:: bash

    OSError: [WinError 225] Operation did not complete successfully because the file contains a virus or potentially unwanted software.
    ...
    win32ctypes.pywin32.pywintypes.error: (225, 'BeginUpdateResourceW', 'Operation did not complete successfully because the file contains a virus or potentially unwanted software.')

This is the antivirus program thinking you're getting a virus. To circumvent this:

1. Open the Windows Security app

2. Navigate to "Virus & threat protection"

3. Click "Manage settings" under "Virus & threat protection settings"

4. Turn off "Real-time protection"

5. Run the `PyInstaller` build command

6. Turn "Real-time protection" back on

What You Can Ignore
###################

You can safely ignore the following messages:

* ``WARNING: Library {LIBNAME} required via ctypes not found``

  * See https://github.com/pyinstaller/pyinstaller/issues/1403
