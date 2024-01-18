##Steps for publishing a new dysh release

0. Prerequisites
    *  A release branch *release-M.m.p* ( **M**ajor, **m**inor, **p**atch) should have been created a few weeks before the release date.
    * This branch is to be used for betas and other pre-releases, with versioning managed as described below.
1. Housekeeping

     * Switch to the release branch in your sandbox
     * Ensure all required functionality/bug fixes are merged, no outstanding PRs.
     * Ensure all `CI <https://github.com/GreenBankObservatory/dysh/actions>`_ for this branch is passing.
     * Set up readthedocs
         - `Activate the release branch <https://readthedocs.org/projects/dysh/versions/>`_
         - Ensure docs are building
     * Change version string `dysh/__init__.py`
         - we follow the `major.minor.patch[qualifier]` paradigm,
           where `qualifier` is e.g., `b`, `rc1`, etc.
         - push `dysh/__init__.py`
     * Freeze release branch at least week before release

2. Create the release on github.com
     - Follow `these steps to begin the release and generate release notes <https://docs.github.com/en/repositories/releasing-projects-on-github/automatically-generated-release-notes>`_
         -  Create a new tag that matches your version in step 1, e.g. '0.2.0'
         - The release title should be 'v'+[the tag]
         - Choose **pre-release** or **latest release** as appropriate
         - If you want another pair of eyeballs, click *Save Draft*, otherwise click *Publish Release*.
         - The CI should create the release and upload it to pypi.org.

3. Monitor the `CI <https://github.com/GreenBankObservatory/dysh/actions>`_, `pypi.org <https://pypi.org/manage/project/dysh/releases/>`_, and `readthedocs.org <https://readthedocs.org/projects/dysh/>`_ for progress/problems.
