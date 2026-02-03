
.. _usersguide-logging:

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

The log verbosity level is set either on the command line with `-v`
or `--verbosity` or can be set in the Python session itself with `~dysh.log.init_logging`.  You can think of 
verbosity as the complement of the
traditional `log severity level <https://docs.python.org/3/howto/logging.html#logging-basic-tutorial>`_. 
The values of verbosity and severity levels are as follows:

.. list-table:: Log Levels
   :header-rows: 1

   * - Verbosity
     - Log Severity Level
   * - 3
     - DEBUG
   * - 2
     - INFO
   * - 1
     - ERROR
   * - 0
     - WARNING


Setting the log verbosity level interactively while in dysh:

.. code:: Python

   from dysh.log import init_logging
   init_logging(2)  # Log only ERROR or higher 
