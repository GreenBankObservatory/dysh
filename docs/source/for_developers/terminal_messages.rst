*****************
Terminal Messages
*****************

Dysh Rich Console
=================

The `dysh/config/rich_theme.py` file defines a `rich.console.Console` object called `DyshRichConsole`, which employs a `rich.highlighter.RegexHighlighter` object called `DyshRichHighlighter` to highlight syntax.

To use the `DyshRichConsole` to print pretty messages, import the following at the top of your code:

.. code::

    from dysh.config.rich_theme import DyshRichConsole

Then, to use it:

.. code::

    DyshRichConsole.print(foo)

Friendly Messages
=================

`dysh/util/messages.py` has a class called `FriendlyMessages`, which includes general messages like greetings.

System Messages
===============

`dysh/util/messages.py` has another class called `SystemMessages`, which give information about the system.
