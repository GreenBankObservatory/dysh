The recommended way to develop dysh is to [fork](https://guides.github.com/activities/forking/) dysh
from the *upstream*, work on your own clone of that fork, and issue a pull request when your new code is ready for review.   Below are two ways to do this. Ideally, you should write regression test(s) for your code that can be run in dysh's CI (i.e., by *pytest*). Be sure to occasionally sync your fork with upstream.

Please see our [For Developers guide](https://dysh.readthedocs.io/en/latest/for_developers/index.html) for more detailed information.

## GitHub

Follow the instructions in [github for creating a fork](https://guides.github.com/activities/forking/) and cloning your fork. Once your work is ready for review, [issue a pull request](https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project#making-a-pull-request).   You can optionally select a member of the dysh team to review it. If not, a reviewer will be assigned.
To sync your fork, click the 'Sync Fork' button on the github repo page.

## gh:   github CLI

 If using the **gh** command
from the [github cli](https://cli.github.com/), you should create your own fork, clone locally, and set the *upstream*:

      gh repo fork https://github.com/GreenbankObservatory/dysh

You will be asked if you want to clone the fork (say yes). If all is well, the following commands should show the correct *origin* and *upstream*:

      git remote -v

      origin git@github.com:YOURNAME/dysh.git
      upstream git@github.com:GreenbankObservatory/dysh.git

for both the (fetch) and (push). None of these should be blank! You are now ready for working
in your own branches and issue a pull request (PR).   To sync your fork:

      gh repo sync

To issue a pull request:

      gh pr create

Once the branch has been merged by the upstream, there is no need to keep it.
It can be removed as follows:

      git branch -d your_branch_name
      git push original --delete your_branch_name


## Helpful git options

*  Show all files modified in a branch AAA

       `git diff main...AAA --name-status --diff-filter=M`

*  When was a branch created

      git show --summary \`git merge-base AAA main\`

* To see which files belonged in a commit, find the SHA (e.g. via "git log" or "git log file" if
   you know the file that belonged to it), then

       `git diff-tree --no-commit-id --name-only -r SHA`

* Difference between two SHA's

       `git diff <commit-id> <commit-id>`
