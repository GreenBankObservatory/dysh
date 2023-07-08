# Developer Install

The dysh README suggests a route to install **dysh** using
**hatch**. As suggested also, there are several ways how you can
install **dysh** for development. We give a few practical examples, all
based on having "dyshN" directories in my ~/GBT directory. It is imperative
that a developer install takes place in a shielded environment, generally
using a virtual environment.  We list a few, but if you found another approach,
please share.

## dysh1: native Python

Here is an example using native python on a vanilla Ubuntu system (python version 3.11 may be different).
You will need initial admin privilages for this.

     # first ensure your native python has at least a way to run pip and allow a venv
     sudo apt install python3 python3-pip ipython3 python3.11-venv jupyter-core

     # setup a venv, for example in a $HOME/venv hierarchy
     mkdir -p $HOME/venv
     python3 -m venv $HOME/venv/dysh1
     
     # activate this venv
     source $HOME/venv/dysh1/bin/activate
     
     # install hatch
     pip install hatch notebook

After this dysh can be installed in a virtual environment controlled by hatch

     mkdir ~/GBT/dysh1
     cd ~/GBT/dysh1
     git clone https://github.com/GreenBankObservatory/dysh
     cd dysh

     # setup dysh with hatch (be sure to be in the dysh directory)
     hatch shell
     pip install -r requirements_dev.txt
         # some warning about running ipython
     pip install -r docs/requirements.txt     
     hatch build
     pip install -e .
     ipython            # this initially gave a matplotlib error, but it went away
     exit               

Any time development is needed:

     source $HOME/venv/dysh1/bin/activate
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
Simple and self-contained, but with an anaconda3 to maintain.

     mkdir ~/GBT/dysh2
     cd ~/GBT/dysh2

     ../install_anaconda3                # DM me for a copy
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
   




## dysh3: anaconda3 python with virtual environment

Here is an example using an anaconda3 python, but now using hatch 

     mkdir ~/GBT/dysh3
     cd ~/GBT/dysh3

     ../install_anaconda3                # DM me for a copy
     source python_start.sh
     
     pip install hatch

After this dysh can be installed in a virtual environment controlled by hatch,
pretty much following what we did in *dysh1*:

     git clone https://github.com/GreenBankObservatory/dysh
     cd dysh

     # setup dysh with hatch (be sure to be in the dysh directory)
     hatch shell
     pip install -r requirements_dev.txt
     pip install -r docs/requirements.txt     
     hatch build
     pip install -e .

and verify

     python -c 'import dysh; print(dysh.__version__)'
     echo "git BRANCH: $(git branch --show-current)   HEAD: $(git rev-list --count HEAD)"

and when done, exit will terminate the shell

     exit
   

Any time development is needed:

     source $HOME/GBT/dysh3/python_start.sh
     cd ~/GBT/dysh3/dysh
     hatch shell





## Sample workflows

Minor issue:  with hatch, if you're not in the code tree (much like git) you don't know
where your code tree is. Do we need peter's "rc" files. Do we need a module file?


###  Simple dysh commands


     python -c 'import dysh; print(dysh.__version__)'
     python -c 'import dysh; print(dysh.__file__)'

###  Building Documentation

     cd dysh/docs
     make html
     xdg-open _build/html/index.html
