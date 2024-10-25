*************
Git Workflows
*************

Branches
========

The development team employs a GitFlow workflow with personal branching. This means that code in the `main` branch should always be in a releasable state. Developers should maintain their own development branches and commit changes to a `release-x.y.z` branch. When it's time to release, a dedicated team member will merge the `release-x.y.z` branch with the `main` branch and tag it accordingly.

.. mermaid::

    %%{init: { 'theme': 'base' } }%%

    gitGraph
        commit id: "1"
        commit id: "2"
        branch release-x.y.z order: 1
        commit id: "3"
        branch personB-devel order: 3
        branch personA-devel order: 2
        commit id: "4"
        checkout personB-devel
        commit id: "5"
        commit id: "6"
        checkout personA-devel
        commit id: "7"
        commit id: "8"
        checkout release-x.y.z
        merge personA-devel id: "9"
        commit id: "10"
        checkout personB-devel
        merge release-x.y.z id: "11"
        commit id: "12"
        commit id: "13"
        checkout release-x.y.z
        merge personB-devel id: "14"
        checkout main
        merge release-x.y.z  id: "15" tag: "release-x.y.z"

Releases
========

Release branches will be locked once work on the next release begins.

Setting up your own development branch
======================================

In the directory you want to work in, set up the repo:

.. code-block:: bash

    $ git clone git@github.com:GreenBankObservatory/dysh.git
    $ cd dysh

To check out a branch called {{branch-name}}, just do

.. code-block:: bash

    $ git checkout {{branch-name}}

Current development is done in the `main` branch. To set up your own development branch called `{{your-name}}-devel`, do the following:

.. code-block:: bash

    $ git checkout main
    $ git checkout -b {{your-name}}-devel

Say someone just made changes to `main`. To add them to your branch, do the following:

.. code-block:: bash

    $ git checkout main
    $ git pull
    $ git checkout {{your-name}}-devel
    $ git merge main

Now make and activate a python virtual environment so your work doesn't impede or break other concurrent projects. Whenever you do some work, make sure you are in your own development branch. When you are ready to merge changes made by other developers into your own branch,

When you are ready to commit changes, review what's been changed with

.. code-block:: bash

    $ git status

and then add the intended files using

.. code-block:: bash

    $ git add path/to/changed_file.py

Check `dysh/.gitignore` to make sure you are not adding ignored files (virtual environment data, `_build/`, etc.). Then commit and push with

.. code-block:: bash

    $ git commit -m "this is my commit message"
    $ git push

The first time you run this, it will give a command about setting the origin upstream. Simply copy and run that command. Users of GitHub Desktop can also achieve all of these above steps using the app interface. Next, go to the `dysh GitHub page <https://github.com/GreenBankObservatory/dysh/>`_ and submit a pull request.

Now follow the steps in the next page to set up more integrations.
