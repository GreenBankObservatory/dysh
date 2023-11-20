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

`$ git clone git@github.com:GreenBankObservatory/dysh.git`

`$ cd dysh`

`$ git checkout [branch-name]`

with the [branch-name] replaced by the current release or other branch. Currently, this is "release-0.2.0". Then setup your own development branch with your first name followed by "-devel":

`$ git checkout [name-devel]`

Now make and activate a python virtual environment so your work doesn't impede or break other concurrent projects. Whenever you do some work, make sure you are in your own development branch. When you are ready to merge changes made by other developers into your own branch,

`$ git checkout release-0.2.0`

`$ git pull`

`$ git checkout [name-devel]`

`$ git merge release-0.2.0`

When you are ready to commit changes, review what's been changed with

`$ git status`

and then add the intended files using

`$ git add [path/to/changed/file]`

check `dysh/.gitignore` to make sure you are not adding ignored files (virtual environment data, `_build/`, etc.). Then commit and push with

`$ git commit -m "[this is my commit message]"`

`$ git push`

The first time you run this, it will give a command about setting the origin upstream. Simply copy and run that command. Users of GitHub Desktop can also achieve all of these above steps using the app interface. Next, go to the `dysh GitHub page <https://github.com/GreenBankObservatory/dysh/>`_ and submit a pull request.
