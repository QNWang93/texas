import gspread
from oauth2client.service_account import ServiceAccountCredentials
import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

import pygsheets

from astropy.time import Time 

import warnings
warnings.filterwarnings("ignore")

def Rise_covered(Table):
	"""
	Check to see if the rise of the transient is covered by TESS observations.

	Input:
	------
	Table - pandas dataframe 

	Outputs:
	--------
	ind 	- list, indicies that TESS covers
	sector 	- list, TESS sectors of the transients
	"""
	TESS_times = pd.read_csv('TESS_sector_jd.csv').values
	ind = []
	sector = []
	for i in range(len(Table)):
		maxtime = Time(Table['min_date'].iloc[i].replace(' ', 'T')).jd
		disc = Time(Table['disc_date'].iloc[i].replace(' ', 'T')).jd
		
		diff = maxtime - disc 
		if diff > 5:
			buffer = 20 # assuming max and rise is ~18 days
		else:
			buffer = 10 # arbitrary choice
		sneeze = ((TESS_times[:,1] < (maxtime - buffer)) & 
						  (TESS_times[:,2] > (maxtime - 3)))
		if sneeze.any():
			sector += [TESS_times[sneeze,0][0]]
			ind += [i]
	return ind, sector #Table[ind] 


def Check_gal_lat(Table):
	"""
	Calculates galactic latitude, which is used as a flag for event rates.
	"""
	ind = []
	for i in range(len(Table)):
		b = SkyCoord(ra=float(Table['transient_RA'][i])*u.degree, dec=float(Table['transient_Dec'][i])*u.degree, frame='icrs').galactic.b.degree
		if abs(b) >= 10:
			ind += [i]

	return ind

def Check_extinction(Table):
	ind = []
	for i in range(len(Table)):
		if Table['mw_ebv'][i] <= 0.2:
			ind += [i]
	return ind

def Check_point(Table):
	ind = []
	print(Table['point_source_probability'])
	for i in range(len(Table)):
		if np.isfinite(Table['point_source_probability'][i]) and Table['point_source_probability'][i] is not None:
			if Table['point_source_probability'][i] <= 0.8:
				ind += [i]
		else:
			ind += [i]
	return ind


def Gal_coord(Table):
	l = []
	b = []
	for i in range(len(Table)):
		c = SkyCoord(ra=float(Table['transient_RA'][i])*u.degree, dec=float(Table['transient_Dec'][i])*u.degree, frame='icrs')
		l += [c.galactic.l.degree]
		b += [c.galactic.b.degree]
	return l, b


def Check_type(Table):
	ind = np.where((Table['spec_class'] == 'SN Ia') | (Table['spec_class'] == 'None'))[0]
	return ind 

def Check_z(Table):
        return 0

def YSE_list():
	all_cand = pd.read_csv('https://ziggy.ucolick.org/yse/explorer/54/download?format=csv')
	all_cand = all_cand.drop_duplicates(subset='name')
	all_cand['spec_class'] = all_cand['spec_class'].fillna(value = 'None')
	# Spec tyoe 
	ind = Check_type(all_cand)
	good = all_cand.iloc[ind]
	good = good.reset_index(drop=True)
	# galactic latitude
	ind = Check_gal_lat(good)
	good = good.iloc[ind]
	good = good.reset_index(drop=True)
	# extinction
	ind = Check_extinction(good)
	good = good.iloc[ind]
	good = good.reset_index(drop=True)
	# point source
	ind = Check_point(good)
	good = good.iloc[ind]
	good = good.reset_index(drop=True)
	


	df = pd.DataFrame()
	#df['name'] = good['name'] + '#' + url + good['name'] + '/'
	links = []
	#for i in range(len(good['name'])):
	#	links += ['=HYPERLINK("https://ziggy.ucolick.org/yse/transient_detail/{0}/","{0}")'.format(good['name'].iloc[i])]
	l, b = Gal_coord(good)
	df['Name'] = good['name'] 
	df['RA'] = good['transient_RA']
	df['Dec'] = good['transient_Dec']
	df['l'] = l
	df['b'] = b
	df['Peak mag'] = good['min_mag']
	df['Peak time'] = good['min_date']
	df['Disc date'] = good['disc_date']
	df['Type'] = good['spec_class']
	ind = np.where(df['Type'].isna())[0]
	df['Type'].iloc[ind] = 'Phot ' + good['phot_class'].iloc[ind]
	df['PS prob'] = good['point_source_probability']
	df['Redshift'] = good['transient_z']
	df['MW E(B-V)'] = good['mw_ebv']
	
	return df


def Update_sheet():
	"""
	Updates the sheet.

	"""
	filename = './candidates.csv'
	web = pd.read_csv(filename)
	df = YSE_list()
	print(web.keys())
	for i in range(len(df['Name'])):
		name = df['Name'][i]
		if (web['Name'] == name).any():
			ind = np.where(web['Name'] == name)[0]
			for col in df.columns:
				web[col].iloc[ind] = df[col].iloc[i]
		else:
			print('Added ', name)
			web.loc[-1] = df.iloc[i]
			web.index = web.index + 1
			web = web.sort_index()

	web.iloc[:,15:] = web.iloc[:,10:].replace({pd.np.nan: ''})
	
	web.to_csv(filename,index=False)
	print('Updated')
	return 

if __name__ == '__main__':
	Update_sheet()

