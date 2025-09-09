.. _scanblocks:

#########
ScanBlock
#########

A `~dysh.spectra.scan.ScanBlock` is a `list` of `~dysh.spectra.scan.ScanBase` objects.
As illustrated by the diagram below, it serves to hold together a series of scans.
One of its purposes is to facilitate working with groups of scans, as it allows time averaging them together (using :py:meth:`ScanBlock.timeaverage <dysh.spectra.scan.ScanBlock.timeaverage>`), or subtracting a baseline from the integrations in a group of scans (using :py:meth:`ScanBlock.subtract_baseline <dysh.spectra.scan.ScanBlock.subtract_baseline>`).

.. mermaid::
    :caption: The contents of a `~dysh.spectra.scan.ScanBlock` are a series of `~dysh.spectra.scan.ScanBase` objects. `~dysh.spectra.scan.ScanBlock` are the return of the calibration routines (e.g., `~dysh.fits.gbtfitsload.GBTFITSLoad.getps`, `~dysh.fits.gbtfitsload.GBTFITSLoad.getfs` or `~dysh.fits.gbtfitsload.GBTFITSLoad.gettp`)

    flowchart TD

    subgraph newLines[GBTFITSLoad
        50 Scans
        Position Switching
        Dual Linear Polarization
        1 Beam
        4 Frequency Windows
        100 integrations each
        ]




    end


    newLines -- getps( scan=45, plnum=1, ifnum=0 ) --> ScanBlock1
    newLines -- gettp( scan=[17,18,19], intnum=np.r_[50:100], ifnum=2 ) --> ScanBlock2



    subgraph ScanBlock1[ScanBlock]
            psscan["spectra.scan.PSScan<br />scans = 44,45<br />plnum = 1<br />ifnum = 0<br />intnum = (0,100)"]
        end
    subgraph ScanBlock2[ScanBlock]
            tpscan1["spectra.scan.TPScan<br />scan=17<br />plnum = 0<br />ifnum = 2<br />intnum=(50,100)"]
            tpscan2["spectra.scan.TPScan<br />scan=18<br />plnum = 0<br />ifnum = 2<br />intnum=(50,100)"]
            tpscan3["spectra.scan.TPScan<br />scan=19<br />plnum = 0<br />ifnum = 2<br />intnum=(50,100)"]

        end

    ScanBlock1[Scan Block] -- timeaverage() --->spectrum1[Spectrum]
    ScanBlock2[Scan Block] -- timeaverage() --->spectrum2[Spectrum]


The return of most calibration methods (e.g., `~dysh.fits.gbtfitsload.GBTFITSLoad.getps`, `~dysh.fits.gbtfitsload.GBTFITSLoad.getfs` or `~dysh.fits.gbtfitsload.GBTFITSLoad.gettp`) is a `~dysh.spectra.scan.ScanBlock`.
Since a `~dysh.spectra.scan.ScanBlock` is a `list`, it is possible to append to a `~dysh.spectra.scan.ScanBlock`.
For example, you can do the following:

.. code:: python3

    from dysh.spectra import ScanBlock
    scan_block = ScanBlock()
    # Calibrate some data. For example:
    # ps_sb = sdfits.getps(scan=1, ifnum=0, plnum=0, fdnum=0)
    # fs_sb = sdfits.getfs(scan=2, ifnum=0, plnum=0, fdnum=0)
    scan_block.append(ps_sb)
    scan_block.append(fs_sb)

This can be useful when working with different observing procedures or for combining polarizations and/or spectral windows (IFs).
