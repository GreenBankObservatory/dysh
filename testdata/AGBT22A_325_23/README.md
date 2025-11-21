# AGBT22A_325_23

Argus Nod observations.

## Making the test data

Retrieve the engineering FITS files from the archive and fill using:

```bash
sdfits -serial -backends=vegas -noprompt AGBT22A_325_23/ScanLog.fits -scans=43:46
```
