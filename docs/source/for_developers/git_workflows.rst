*************
Git Workflows 
*************

Branches
========

Make your own development/feature branch and stay synced with the current `release-x.y.z` branch. When it's time to release, merge the `release-x.y.z` branch with the `main` branch. Tag the merge accordingly. 

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