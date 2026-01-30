
.. _usersguide-logging

**********************
Using the dysh logger
**********************

dysh has a `built-in logging mechanism <https://dysh.readthedocs.io/en/latest/reference/modules/dysh.log.html>`_
based on the 
`Python logging module <https://docs.python.org/3/library/logging.html>`_ 
whereby messages are categorized by log level (severity) and logged to
a file and/or the screen.   The command line arguments that control logging are

.. code:: bash

   dysh 

     -v {0,1,2,3}, --verbosity {0,1,2,3} Set logging verbosity
     --log LOG                           Specify log path
     -q, --quiet                         Silence DEBUG- and INFO-level logs to stderr

The log verbosity level is set either on the command line with `-v` or `--verbosity`, or can be set in the Python session itself with `init_logging()`. The numeric values correspond to log levels as follows

.. list-table:: Log Levels
   :header-rows: 1

   * - Value
     - Level
   * - 0
     - DEBUG
   * - 1
     - INFO
   * - 2
     - WARNING
   * - 3
     - CRITICAL


Setting the log level interactively while in dysh:

.. code:: Python

   init_logging(2)

