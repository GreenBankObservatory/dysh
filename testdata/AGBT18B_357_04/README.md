# AGBT18B_357_04

Argus SubBeamNod observations.

## Making the test data

Retrieve the engineering FITS files from the archive and fill using:

```bash
sdfits -serial -backends=vegas -noprompt AGBT18B_357_04/ScanLog.fits -scans=1:3
```
