# AGBT25A_504_03

W-band frequency switched observations.

## Making the test data

Fill with ``sdfits`` and remove any extra files.

```bash
sdfits AGBT25A_504_03 -scans=18
cd AGBT25A_504_03.raw.vegas
rm AGBT25A_504_03.raw.vegas.A.flag AGBT25A_504_03.raw.vegas.E.fits AGBT25A_504_03.raw.vegas.E.flag AGBT25A_504_03.raw.vegas.E.index AGBT25A_504_03.raw.vegas.index
```
