# Developer Install

The dysh README suggests a route to install **dysh** using
**hatch**. As suggested also, there are several ways how you can
install **dysh** for development. We give a few practical examples, all
based on having "dyshN" directories in my ~/GBT directory. It is imperative
that a developer install takes place in a shielded environment, generally
using a virtual environment.  We list a few:

## dysh1: native Python

Here is an example using native python on an vanilla Ubuntu system:

     # first ensure your native python has at least a way to run pip and allow a venv
     sudo pip install python3 python3-pip ipython3 python3.11-venv jupyter-core

     # setup a venv, for example in a $HOME/venv hierarchy
     mkdir $HOME/venv
     python3 -m venv $HOME/venv/dysh
     
     # activate this venv
     source $HOME/venv/dysh/bin/activate
     
     # install hatch
     pip install hatch notebook

After this dysh can be installed in a virtual environment controlled by hatch

     mkdir ~/GBT/dysh1
     cd ~/GBT/dysh1
     git clone https://github.com/GreenBankObservatory/dysh
     cd dysh

     # setup dysh with hatch
     hatch shell
     pip install -r requirements_dev.txt
     pip install -r docs/requirements.txt     
     hatch build
     pip install -e .
     ipython            # this initially gave a matplotlib error, but it went away
     exit               

Any time development is needed:

     source $HOME/venv/dysh/bin/activate
     cd ~/GBT/dysh1/dysh
     hatch shell

and as always, verify it's there:

     python -c 'import dysh; print(dysh.__version__)'
     echo "git BRANCH: $(git branch --show-current)   HEAD: $(git rev-list --count HEAD)"

and when done, exit the hatch sub-shell

     exit

this will still return to the native virtual environment, so one more exit is needed to kill this shell

     exit
     
## dysh2: anaconda3 python

Here is an example using an anaconda3 python, no virtual environments, no hatch, no nothing.

     mkdir ~/GBT/dysh2
     cd ~/GBT/dysh2

     install_anaconda3                # DM me for a copy
     source python_start.sh

     git clone https://github.com/GreenBankObservatory/dysh
     cd dysh
     pip install -r requirements_dev.txt
     pip install -r docs/requirements.txt
     pip install -e .


any time development is needed:

     source ~/GBT/dysh2/python_start.sh

and verify

     python -c 'import dysh; print(dysh.__version__)'
     echo "git BRANCH: $(git branch --show-current)   HEAD: $(git rev-list --count HEAD)"

and when done, exit will terminate the shell

     exit
   

