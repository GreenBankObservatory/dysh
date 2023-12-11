************
Git Workflow
************

Background
==========

Dysh uses the `OpenAstronomy model for its branching and release structure <https://packaging-guide.openastronomy.org/en/latest/releasing.html#releasing-from-branches>`_. In practice, this essentially means that we are using a standard "feature branch" workflow, where new features are added via Pull Requests based on ``main``.

If you are an external collaborator, you will need to create Pull Requests from a fork of Dysh, but the instructions below should otherwise be fairly similar.

Developers are welcome to manage their git workflow however they want, but using the `GitHub CLI <https://cli.github.com/>`_ or `GitHub Desktop <https://desktop.github.com/>`_ simplify many operations.

How to Contribute
=================


There are two types of workflow:
1. Contributing to a future release. This is the most common workflow
2. Patching a previous release. This involves changes to a "release branch" -- changes will *not* go on `main`

Below is a diagram demonstrating both of these workflows. They are explored in more depth farther down.

.. mermaid::

    %%{init: { 'logLevel': 'debug', 'theme': 'base', 'gitGraph': {'showBranches': true, 'showCommitLabel':true,'mainBranchOrder': 10}} }%%

    gitGraph

        commit id: "freeze 0.2" tag: "v0.3dev"
        branch v0.2 order: 11

        checkout v0.2
        commit id: "start v0.2 beta test" tag: "v0.2.0b1"

        branch v0.2.0b1_fix order: 14
        commit id: "fix bugs from beta"
        checkout v0.2
        merge v0.2.0b1_fix id: "Release v0.2.0rc1" tag: "v0.2.0rc1"
        commit id: "Release v0.2.0" tag: "v0.2.0"
        commit id: "Bug fix for v0.2" tag: "v0.2.1"
        checkout v0.2

        checkout main
        branch feature_x order:3
        commit id: "start feature x"
        commit id: "finish feature x"
        checkout main
        merge feature_x id: "freeze 0.3" tag: "v0.4dev"

        branch v0.3 order:100
        commit id: "start v0.3 beta test" tag: "v0.3.0b1"
        checkout main

        branch feature_y order:1
        checkout feature_y
        commit id:"implement feature y"
        checkout main

Contribute to ``main`` (Add to a Future Release)
++++++++++++++++++++++++++++++++++++++++++++++++


Primary development is always done on a feature branch off of `main`. The only way for code to get into `main` is via a `Pull Request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests>`_ on the `Dysh GitHub <https://github.com/GreenBankObservatory/dysh>`_. A Pull Request is essentially a group of commits with a common purpose, and it is associated with a discussion thread. It encapsulates both the code changes themselves, and the explanation of why they are being made, how they work, etc.

Here is the lifecycle of a simple Pull Request:

.. code-block:: sh

    git checkout main
    git fetch origin main
    git checkout -b my-cool-new-feature
    # <make some code changes>
    git commit -am "Desc of changes"
    # <make some code changes>
    git commit -am "More changes"
    git push origin my-cool-new-feature

At this point, the developer considers the work complete, and the Dysh GitHub has a new branch, ``my-cool-new-feature``, and the developer's repo will look like:

.. mermaid::

    %%{init: { 'logLevel': 'debug', 'theme': 'base', 'gitGraph': {'showBranches': true, 'showCommitLabel':true}} }%%

    gitGraph
        commit id:" "
        branch my-cool-new-feature
        checkout my-cool-new-feature
        commit id: "Desc of changes"
        commit id: "More changes"

Now they can open a Pull Request in the Dysh repo on GitHub.

At this point, a few things happen:

1. GitHub CI kicks off, which runs code quality checks and unit tests
2. Another developer will review your code contribution
3. Discussion of the new feature will occur
4. Eventually, the PR will be accepted or rejected
    - If it is accepted, it will be merged into ``main`` and be a part of the next release
    - If it is rejected, it will go nowhere


If accepted, the Dysh repo will now look like:

.. mermaid::

    %%{init: { 'logLevel': 'debug', 'theme': 'base', 'gitGraph': {'showBranches': true, 'showCommitLabel':true}} }%%

    gitGraph
        commit id:" "
        branch my-cool-new-feature
        checkout my-cool-new-feature
        commit id: "Desc of changes"
        commit id: "More changes"
        checkout main
        merge my-cool-new-feature




Contribute to a Previous Release (Add to a release branch)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In this example, release ``v0.2.0`` has been released, and development has started on ``v0.3.0`` (on ``main``). Then, a bug is reported in ``v0.2.0``, necessitating a bug fix.


.. mermaid::

    %%{init: { 'logLevel': 'debug', 'theme': 'base', 'gitGraph': {'showBranches': true, 'showCommitLabel':true}} }%%

    gitGraph

        commit id: "freeze 0.2" tag: "v0.3dev"
        branch v0.2

        checkout v0.2

        checkout main
        commit id: "  "
        checkout v0.2

        commit id: "Release 0.2.0" tag: "v0.2.0"
        branch v0.2.0-bug-fix
        commit id: "fix bug"
        checkout v0.2
        merge v0.2.0-bug-fix id: "Release v0.2.1"
        commit id: "Release 0.2.1" tag: "v0.2.1"

        checkout v0.2

        checkout main
        commit id: " "
