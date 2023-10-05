*********
Structure
*********

.. mermaid::

    erDiagram
        SESSION {
            int session_id PK, FK
            int obsblock_id FK
        }
        OBSBLOCK {
            int obsblock_id PK, FK
            int spectrum_id FK
            int flag
        }
        SPECTRUM {
            int spectrum_id PK, FK
            int flag
        }
        SPECTRUM_DATA {
            int spectrum_id PK, FK
            list flags
            list frequencies
            list intensities
        }

        SESSION ||--|{ OBSBLOCK : "has one or more"
        OBSBLOCK ||--|{ SPECTRUM : "has one or more"
        SPECTRUM ||--|| SPECTRUM_DATA : "has one"