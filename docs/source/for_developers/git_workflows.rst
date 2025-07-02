*************
Git Workflows
*************

Branches
========

The development team employs a GitFlow workflow with feature-based branching.
This means that code in the `main` branch should always be in a releasable state.
Developers should maintain a single branch for each feature, bug fix, etc. and submit pull requests to the `main` branch when ready.
Pull Requests must be reviewed by a different team member before being allowed to merge.
When it's time to release, a dedicated team member will merge the `release-x.y.z` branch with the `main` branch and tag it accordingly.

.. mermaid::

    ---
    config:
        logLevel: 'debug'
        theme: 'base'
        gitGraph:
            showCommitLabel: false
            mainBranchOrder: 2
    ---


    gitGraph
        commit id: "1"
        branch enhancement-X order: 3
        checkout enhancement-X
        commit id: "2"
        commit id: "3"
        checkout main
        merge enhancement-X
        branch issue-Y-fix  order: 4
        checkout issue-Y-fix
        commit id: "4" 
        commit id: "5"
        commit id: "6"
        checkout main
        merge issue-Y-fix
        branch release-1.B.C  order: 0 
        checkout release-1.B.C
        commit id: "7"
        checkout main
        branch enhancement-Z order: 5
        checkout enhancement-Z
        commit id: "8"
        checkout main
        merge enhancement-Z
        branch release-2.B.C order:1
        checkout release-2.B.C
        commit id: "9"
        checkout main
        commit id: "10"

Releases
========

Release branches will be locked once work on the next release begins. Only quick fixes during the beta period will be allowed.

Setting up your own development branch
======================================

In the directory you want to work in, set up the repo:

.. code-block:: bash

    $ git clone git@github.com:GreenBankObservatory/dysh.git
    $ cd dysh

To check out a branch called `branch-name`, just do

.. code-block:: bash

    $ git checkout branch-name

Development should be done in a separate dedicated branch to that feature or bug fix.
For example, to develop a fix for Github issue YYY that has to do with a spectrum math error, ensure that your `main` branch is up to date and start a new branch.

.. code-block:: bash

    $ git checkout main
    $ git pull
    $ git checkout -b YYY-spectrum-math-fix

Say someone just made changes to `main`. To add them to your branch, do the following:

.. code-block:: bash

    $ git checkout main
    $ git pull
    $ git checkout YYY-spectrum-math-fix
    $ git merge main


Whenever you do some work, make sure you are in your separate development branch.
When you are ready to merge changes made by other developers into your own branch,

When you are ready to commit changes, review what's been changed with

.. code-block:: bash

    $ git status

and then add the intended files using

.. code-block:: bash

    $ git add paths/to/changed_files.py

Check `dysh/.gitignore` to make sure you are not adding ignored files (virtual environment data, `_build/`, etc.). Then commit and push with

.. code-block:: bash

    $ git commit -m "this is my commit message"
    $ git push

The first time you run this, it will give a command about setting the origin upstream. Simply copy and run that command. Users of GitHub Desktop can also achieve all of these above steps using the app interface. Next, go to the `dysh GitHub page <https://github.com/GreenBankObservatory/dysh/>`_ and submit a pull request to the `main` branch.

Now follow the steps in the next page to set up more integrations.
