import gspread
from oauth2client.service_account import ServiceAccountCredentials
import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import time
import pygsheets
from astropy.io import ascii
from astropy.time import Time 
import warnings
from config import lc_cfg
import texas
import os
warnings.filterwarnings("ignore")

def Save_space(Save):
    """
    Creates a pathm if it doesn't already exist.
    """
    try:
        if not os.path.exists(Save):
            os.makedirs(Save)
    except FileExistsError:
        pass

def get_gglsheet():
# use creds to create a client to interact with the Google Drive API
    scope = ['https://spreadsheets.google.com/feeds',
		'https://www.googleapis.com/auth/drive']
    creds = ServiceAccountCredentials.from_json_keyfile_name('client_secret.json', scope)
    client = gspread.authorize(creds)
    sheet = client.open("texas").sheet1
    return sheet

def TESS_cover(ra, dec, disc):
	before_leeway = 10	 # Days of leeway before date
	after_leeway = 10	 # Days of leeway after date
	tess_date = lc_cfg['tess_date']

	url = 'https://heasarc.gsfc.nasa.gov/cgi-bin/tess/webtess/'
	url += 'wtv.py?Entry={ra}%2C{dec}'
	r = requests.get(url.format(ra=str(ra), dec=str(dec)))
	if r.status_code!=200:
		print('status message:',r.text)
		error = 'ERROR: could not get {url}, status code {code}'
		raise RuntimeError(error.format(url=url, code=r.status_code))
		return(None)

	reg = r"observed in camera \w+.\nSector \w+"
	info = re.findall(reg, r.content.decode())
	sectors=[]
	for k in info:
	    sectors.append(int(re.split(r'\s', k)[5])-1)

	if len(sectors)>0:
		for sector in sectors:
			if (discovery_jd > tess_date[int(sector)-1]-before_leeway and
				discovery_jd < tess_date[int(sector)]+after_leeway):
				return(True)
	return(False)


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
		if Table['mw_ebv'][i] <= 0.3:
			ind += [i]
	return ind

def Check_point(Table):
	ind = []
#	print(Table['point_source_probability'])
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
	ind2 = np.where(df['PS prob'].isna())[0]
	df['PS prob'].iloc[ind2] = 'None'
	df['Redshift'] = good['transient_z']
	ind3 = np.where(df['Redshift'].isna())[0]
	df['Redshift'].iloc[ind3] = 'None'
	df['MW E(B-V)'] = good['mw_ebv']
	
	return df

def delete_rows(sheet, dlist):
	dlist.sort()
	for i in range(len(dlist)):
		row = dlist[i]+2-i
		sheet.delete_row(row)
		
def Update_sheet():
	"""
	Updates the sheet.

	"""
#	filename = './candidates.csv'
#	web = pd.read_csv(filename)

	sheet = get_gglsheet()
	st_cont = sheet.get_all_values()
	headers = st_cont.pop(0)
	web = pd.DataFrame(st_cont, columns=headers)
	df = YSE_list()

#	print(web.keys())
#	to_delete = []
#	for i in range(len(web)):
#		texas_file = lc_cfg['home_dir']+web['Name'][i]+'/'+web['Name'][i]+'_texas'
#		if not os.path.exists(texas_file+'.txt'):
#			home_dir = lc_cfg['home_dir']+web['Name'][i]+'/'
#			Save_space(home_dir)
#			print('start TEXAS on '+web['Name'][i])
#			if not texas.main([web['RA'][i], web['Dec'][i], '3', texas_file]):
#				to_delete.append(i)

#	delete_rows(sheet, to_delete)
	print(df['Name'])
	for i in range(len(df['Name'])):
		name = df['Name'][i]
		row = [df[col][i] for col in df.columns]
		home_dir = lc_cfg['home_dir']+name+'/'
		Save_space(home_dir)
		
		texas_exist = False
		if os.path.exists(home_dir+name+'_texas.txt'):
			host_exist = True
		else:
			print('start TEXAS on '+name)
			while True:
				try:
					host_exist = texas.main([df['RA'][i], df['Dec'][i], '3', home_dir + name + '_texas'])
					print(host_exist)
				except:
					time.sleep(60)					
					print('next round texas')
					continue
				break

			
		if (web['Name'] == name).any():
			ind = int(np.where(web['Name'] == name)[0][0])+2
		elif True:#host_exist:
			try:
				while True:#host_exist:
					try:
						sheet.insert_row(row, 2)
					except:
						time.sleep(100)
						print('next round insert sheet')
						continue
					print('Added ', name)
					break
			except:
				print('unable to access PanSTARRS')
	print('Updated')
	return 

if __name__ == '__main__':
	Update_sheet()

