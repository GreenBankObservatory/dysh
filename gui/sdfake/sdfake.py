import json
from astropy.io import fits
from astroquery.simbad import Simbad
from astroquery.splatalogue import Splatalogue
from astroquery.svo_fps import SvoFps
from astropy import units as u
import matplotlib.pyplot as plt

class Source:
    def __init__(self, src_name):
        self.src_name = src_name
        self.query_simbad()

    def query_simbad(self):
        # http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=GAIA&asttype=
        self.filters = {
            'B':{'catalog':'AAVSO','wavelength':4361*u.Angstrom, 'fwhm':650*u.Angstrom, 'color':'b'}, 
            'V':{'catalog':'AAVSO','wavelength':5448*u.Angstrom, 'fwhm':890*u.Angstrom, 'color':'purple'}, 
            'R':{'catalog':'AAVSO','wavelength':6407*u.Angstrom, 'fwhm':1580*u.Angstrom, 'color':'r'}, 
            'G':{'catalog':'GAIA','wavelength':6735.41*u.Angstrom, 'fwhm':4396.69*u.Angstrom, 'color':'g'}, 
            'J':{'catalog':'2MASS','wavelength':12410.51*u.Angstrom, 'fwhm':2149.14*u.Angstrom, 'color':'orange'}, 
            'H':{'catalog':'2MASS','wavelength':16513.66*u.Angstrom, 'fwhm':2609.65*u.Angstrom, 'color':'cyan'}, 
            'K':{'catalog':'2MASS','wavelength':21656.31*u.Angstrom, 'fwhm':2784.55*u.Angstrom, 'color':'k'}, 
            'u':{'catalog':'SDSS','wavelength':3572.18*u.Angstrom, 'fwhm':565.80*u.Angstrom, 'color':'b'}, 
            'g':{'catalog':'SDSS','wavelength':4750.82*u.Angstrom, 'fwhm':1175.63*u.Angstrom, 'color':'g'}, 
            'r':{'catalog':'SDSS','wavelength':6204.29*u.Angstrom, 'fwhm':1130.56*u.Angstrom, 'color':'r'}, 
            'i':{'catalog':'SDSS','wavelength':7519.27*u.Angstrom, 'fwhm':1253.30*u.Angstrom, 'color':'y'}, 
            'z':{'catalog':'SDSS','wavelength':8992.26*u.Angstrom, 'fwhm':998.50*u.Angstrom, 'color':'k'}
            }
        customSimbad = Simbad()
        for fi in self.filters.keys():
            customSimbad.add_votable_fields(f'fluxdata({fi})')
        self.info_tbl = customSimbad.query_object(self.src_name)

class Observation:
    def __init__(self):
        self.set_freqs()
        self.get_lines()

    def set_freqs(self):
        self.fmin = 1.0 * u.GHz
        self.fmax = 2.0 * u.GHz

    def get_lines(self, only_NRAO_rec=True):
        self.line_tbl = Splatalogue.query_lines(
            self.fmin, 
            self.fmax,
            only_NRAO_recommended=only_NRAO_rec
            )

src = Source("HL Tau")
obs = Observation()
filter_list = SvoFps.get_filter_list(facility='Keck')
J_2MASS_transmission = SvoFps.get_transmission_data('2MASS/2MASS.J')
H_2MASS_transmission = SvoFps.get_transmission_data('2MASS/2MASS.H')
K_2MASS_transmission = SvoFps.get_transmission_data('2MASS/2MASS.Ks')
filter_list = SvoFps.get_filter_list(facility='Sloan')
SDSS_u_transmission = SvoFps.get_transmission_data('SLOAN/SDSS.u')
SDSS_g_transmission = SvoFps.get_transmission_data('SLOAN/SDSS.g')
SDSS_r_transmission = SvoFps.get_transmission_data('SLOAN/SDSS.r')
SDSS_i_transmission = SvoFps.get_transmission_data('SLOAN/SDSS.i')
SDSS_z_transmission = SvoFps.get_transmission_data('SLOAN/SDSS.z')

fix, ax = plt.subplots(nrows=2, ncols=1, sharex=True)

cax = ax[0]
for fi in src.filters.keys():
    print(src.info_tbl[f'FILTER_NAME_{fi}'].value[0], src.info_tbl[f'FLUX_{fi}'].value[0])
    lbl_fi = f"{src.info_tbl[f'FILTER_NAME_{fi}'].value[0]} ({src.filters[fi]['catalog']})"
    cax.errorbar(
        src.filters[fi]['wavelength'],
        src.info_tbl[f'FLUX_{fi}'].value[0], 
        xerr=src.filters[fi]['fwhm']/2,
        yerr=src.info_tbl[f'FLUX_ERROR_{fi}'].value[0],
        label=lbl_fi,
        marker='.',
        markersize=5,
        color=src.filters[fi]['color']
        )
cax.legend(loc='upper right')
cax.set_xlabel('Wavelength (A)')
cax.set_ylabel('Flux (mag)')
cax.set_title('Photometry')

cax = ax[1]
cax.plot(J_2MASS_transmission['Wavelength'], J_2MASS_transmission['Transmission'], label="J (2MASS)", color=src.filters['J']['color'])
cax.plot(H_2MASS_transmission['Wavelength'], H_2MASS_transmission['Transmission'], label="H (2MASS)", color=src.filters['H']['color'])
cax.plot(K_2MASS_transmission['Wavelength'], K_2MASS_transmission['Transmission'], label="K (2MASS)", color=src.filters['K']['color'])

cax.plot(SDSS_u_transmission['Wavelength'], SDSS_u_transmission['Transmission'], label="u (SDSS)", color=src.filters['u']['color'])
cax.plot(SDSS_g_transmission['Wavelength'], SDSS_g_transmission['Transmission'], label="g (SDSS)", color=src.filters['g']['color'])
cax.plot(SDSS_r_transmission['Wavelength'], SDSS_r_transmission['Transmission'], label="r (SDSS)", color=src.filters['r']['color'])
cax.plot(SDSS_i_transmission['Wavelength'], SDSS_i_transmission['Transmission'], label="i (SDSS)", color=src.filters['i']['color'])
cax.plot(SDSS_z_transmission['Wavelength'], SDSS_z_transmission['Transmission'], label="z (SDSS)", color=src.filters['z']['color'])
cax.set_xlabel('Wavelength (Angstroms)')
cax.set_ylabel('Transmission Fraction')
cax.set_title('Filter Curves')
cax.legend(loc='upper right')

plt.show()
'''
print(src.info_tbl.columns)
plt.figure()
for fi in src.filters.keys():
    print(src.info_tbl[f'FILTER_NAME_{fi}'].value[0], src.info_tbl[f'FLUX_{fi}'].value[0])
    lbl_fi = f"{src.info_tbl[f'FILTER_NAME_{fi}'].value[0]} ({src.filters[fi]['catalog']})"
    plt.errorbar(
        src.filters[fi]['wavelength'],
        src.info_tbl[f'FLUX_{fi}'].value[0], 
        xerr=src.filters[fi]['fwhm'],
        yerr=src.info_tbl[f'FLUX_ERROR_{fi}'].value[0],
        label=lbl_fi,
        marker='.',
        markersize=5
        )

plt.legend(loc='upper right')
plt.xlabel('Wavelength (A)')
plt.ylabel('Flux (mag)')
plt.show()
'''