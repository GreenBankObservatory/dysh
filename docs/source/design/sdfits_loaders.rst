**************
SDFITS Loaders
**************

Motivation
==========

Although the initial design will be for GBT data, a goal is for `dysh` to be easily modifiable for any single-dish radio telescope. Thus came the idea of an SDFITS loader which would standardize inputs.

ScanBlock
=========

Here's the class diagram for an `ScanBlock` and its derived classes.

.. mermaid::

    classDiagram
        class ScanBlock{
            metadata
            spectra
            summary()
            calibrate()
            baseline()
            timeaverage()
            polaverage()
            finalspectrum()
            undo()
        }
        class PSScan{
            metadata
            spectra
        }
        class FSScan{
            metadata
            spectra
            fold()
        }
        class SubBeamNodScan{
            metadata
            spectra
        }
        class OTFScan{
            metadata
            spectra
        }
        ScanBlock <|-- PSScan
        ScanBlock <|-- FSScan
        ScanBlock <|-- NodScan
        ScanBlock <|-- OTFScan

That's probably not super accurate. I just copied the diagram in the stakeholder presentation from last May.

SDFITSLoad
==========

.. mermaid::

    classDiagram
        class SDFITSLoad{
            _filename
            _bintable
            _index
            _binheader
            _data
            _hdu
            _header
            __len__()
            __repr__()
            _bintable_from_rows()
            _loadlists()
            _summary()
            info()
            bintable()
            binheader()
            filename()
            index()
            reset()
            create_index()
            load()
            fix_meta()
            velocity_convention()
            udata()
            ushow()
            naxis()
            nintegrations()
            rawspectra()
            rawspectrum()
            getrow()
            getspec()
            nrows()
            nchan()
            npol()
            sources()
            scans()
            summary()
            write()
        }


GBTFITSLoad
===========

.. mermaid::

    classDiagram
        class SDFITSLoad{
            method()
        }
        class GBTFITSLoad{
            method()
        }
        SDFITSLoad <|-- GBTFITSLoad
