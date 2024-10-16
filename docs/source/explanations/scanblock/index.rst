.. _scanblocks:

#########
ScanBlock
#########

.. mermaid::

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
