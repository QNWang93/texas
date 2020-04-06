import matplotlib
matplotlib.use('Agg')

import coreapi
import numpy as np
import astropy
import astropy.table as at
import sys
import matplotlib.pyplot as plt
import os
import csv
from jumpssh import SSHSession
import sqlite3
from datetime import datetime as dt
import requests, re
import texas
from astropy.io import ascii
from config import lc_cfg
from astropy.time import Time

def tess_obs(ra, dec):

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
        
    time = []
    if len(sectors)>0:
        for sector in sectors:
            time.append([tess_date[int(sector)-1], tess_date[int(sector)]])
                
    return(time)


def to_array(somelist, column, start = 1):
    array_raw = np.asarray(somelist[start:])[:,column]
    array = [float(i) for i in array_raw]
    return array

def search_atlas(ra, dec):
    atlas_info = lc_cfg['atlas_info']
    gateway_session1 = SSHSession(atlas_info[0]['address'],atlas_info[0]['username'], password=atlas_info[0]['password']).open()
    if len(atlas_info)==3:
        gateway_session2 = gateway_session1.get_remote_session(atlas_info[1]['address'],atlas_info[1]['username'], password=atlas_info[2]['password'])
        remote_session = gateway_session2.get_remote_session(atlas_info[2]['address'], password=atlas_info[2]['password'])
    elif len(atlas_info)==2:
        remote_session = gateway_session1.get_remote_session(atlas_info[1]['address'], password=atlas_info[1]['password'])
    else: 
        print('wrong format for atlas_params: too many jump connection')

    today = dt.today()
    con = sqlite3.connect(":memory:")
    tdate = ' '+str(list(con.execute("select julianday('"+today.strftime("%Y-%m-%d")+"')"))[0][0]-lc_cfg['lookback_days']-2400000)

    result=remote_session.get_cmd_output('./mod_force.sh '+str(ra)+' '+ str(dec)+tdate)

    atlas_lc = at.Table(names=('jd','mag', 'mag_err','flux', 'fluxerr', 'filter','maj', 'min', 'apfit', 'sky', 'zp', 'zpsys'), dtype=('f8', 'f8', 'f8', 'f8', 'f8', 'S1' , 'f8', 'f8', 'f8', 'f8','f8', 'U8'))

    split = result.split('\r\n')

    for i in split[1:]:
        k = i.split()
#        if int(k[3])>int(k[4]):
        atlas_lc.add_row([float(k[0]), float(k[1]), float(k[2]), float(k[3]), float(k[4]), k[5], float(k[12]), float(k[13]), float(k[15]), float(k[16]),  float(k[17]), 'ab'])

    return atlas_lc


def get_ztf(ra, dec):
    ztfurl = 'https://mars.lco.global/?format=json&sort_value=jd&sort_order=desc&cone=%.7f%%2C%.7f%%2C0.0014'%(ra, dec)
    client = coreapi.Client()
    schema = client.get(ztfurl)
    return schema
    
def ztf2lc(ztf_obj):

    lc = at.Table(names=('jd','mag', 'mag_err', 'filter'), dtype=('f8', 'f8', 'f8', 'S1'))
    for i in range(len(ztf_obj['results'])):
        phot = ztf_obj['results'][i]['candidate']
        if phot['isdiffpos'] == 'f':
            continue
        lc.add_row([phot['jd'], phot['magap'], phot['sigmagap'], phot['filter']])
    
    return lc
    
def plot(ax,lc, survey):
    if survey =='atlas':
        limit =lc['flux']<lc['fluxerr']
    else:
        limit = lc['mag']>0
    
    for i in set(lc['filter']):
        index = lc['filter']==i
        jd = lc['jd'][index]
        if jd[0] > 2400000:
            jd = jd -2400000
        mag = abs(lc['mag'][index])
        err = lc['mag_err'][index]
        
        if survey =='atlas':
            flux = lc['flux'][index]
            fluxerr = lc['fluxerr'][index]

            limit =flux<lc_cfg['atlat_sig_limit']*fluxerr
        else:
            limit = mag<0

        real = [not i for i in limit]

        ax.errorbar(jd[real], mag[real], err[real], label = survey +' ' + i, fmt='o', color = lc_cfg['color'][i])
#        print(real, limit, mag[real])
        ax.scatter(jd[limit], mag[limit],c = lc_cfg['color'][i], marker = 7, s =50)

#    ax = plt.gca()

    #ax.grid(True)
    return ax

def plot_atflux(ax, lc, survey):
    
    for i in set(lc['filter']):
        jd = lc['jd'][lc['filter']==i]
        if jd[0] > 2400000:
            jd = jd -2400000
        flux = lc['flux'][lc['filter']==i]
        err = lc['fluxerr'][lc['filter']==i]
        ax.errorbar(jd, flux, err, label = survey +' '+ i, fmt='o',color = lc_cfg['color'][i])

#    ax = plt.gca()

    #ax.grid(True)
    return ax


        
def lc(ra, dec, out_fig, atlas_data_file, disc_t):
    ztf_obj = get_ztf(ra, dec)

    lc = ztf2lc(ztf_obj)
    
    fig, (ax1, ax2)=plt.subplots(2, figsize = (6,10), sharex = True)

    try:
        atlas_lc = search_atlas(ra, dec)
        mask = (abs(atlas_lc['mag'])>10)
        atlas_lc = atlas_lc[mask]
        if len(atlas_lc)>0:
            ascii.write(atlas_lc, atlas_data_file, overwrite=True)
            ax1 = plot(ax1,atlas_lc, 'atlas')

            ax2 = plot_atflux(ax2, atlas_lc, 'atlas')

            ax2.set_ylim(-0.1*max(atlas_lc['flux']), 1.2*max(atlas_lc['flux']))
            if len(lc)>0: 
                ax1 = plot(ax1,lc, 'ztf')
                ymin = min(max([max(abs(lc['mag'])), max(abs(atlas_lc['mag']))])+0.5, 21)
                ax1.set_ylim(ymin, min([min(abs(lc['mag'])), min(abs(atlas_lc['mag']))]) -0.5)
            else:
                ymin = min(max(abs(atlas_lc['mag']))+0.5,21)
                ax1.set_ylim(ymin, min(abs(atlas_lc['mag']))-0.5)
            ax1.plot(disc_t, ymin, marker = "^")
    except:
        print('unable to get atlas lc')
        if len(lc)>0:
            ax1 = plot(ax1,lc, 'ztf')
            ax1.set_ylim(min(max(abs(lc['mag']))+0.5, 21), min(abs(lc['mag']))-0.5)
    
    tess_ob = tess_obs(ra, dec)
    tess_cover = False
    
    for [t1, t2] in tess_ob:
        x1= np.arange(t1-2400000, (t1+t2)/2-2400000-1., 0.1)
        x2= np.arange((t1+t2)/2-2400000+1, t2-2400000, 0.1)
        ax1.fill_between(x1, 10, 22, facecolor='grey', alpha=0.5, label = 'TESS')
        ax1.fill_between(x2, 10, 22, facecolor='grey', alpha=0.5)
        ax2.fill_between(x1, 0, 10000, facecolor='grey', alpha=0.5, label = 'TESS')
        ax2.fill_between(x2, 0, 10000, facecolor='grey', alpha=0.5)
        if disc_t>t1 and disc_t<t2:
            tess_cover = True

    today = dt.today()
    
    con = sqlite3.connect(":memory:")
    tdate = list(con.execute("select julianday('"+today.strftime("%Y-%m-%d")+"')"))[0][0]-2400000
    ax1.axvline(x=tdate, color = 'k', label = today.strftime("%m/%d/%Y, %H:%M"))
    ax1.set_ylabel('Mag')
    ax1.axhline(y = 18, color = 'grey', linestyle = ':')
    ax1.text(tdate-20,18.2,'18 mag',fontsize = 20)
    ax2.axvline(x=tdate, color = 'k', label = 'today')
    ax2.set_xlabel('MJD')
    ax2.axhline(y = 0, color = 'k', linestyle = '-.')
    ax1.set_xlim(tdate+lc_cfg['xlim'][0], tdate+lc_cfg['xlim'][1])
    ax1.legend()
    ax2.legend()

    plt.tight_layout()

    fig.savefig(out_fig)
    return tess_cover
    
def atlas2yse(name, ra, dec, atlas_data_file):
    t = ascii.read(atlas_data_file)
    outname = atlas_data_file[:-9]+'yse.csv'
    filter_fict = {'o':'orange-ATLAS', 'c':'cyan-ATLAS'}
    
    with open(outname, 'w+') as f:
        f.write('SNID: '+name+' \nRA: '+str(ra)+'     \nDECL: '+str(dec)+' \n \nVARLIST:  MJD  FLT  FLUXCAL   FLUXCALERR MAG     MAGERR DQ \n')
        for k in t:
            flt = filter_fict[k['filter']]
            if k['mag']>0:
                flux = 10**(-0.4*(k['mag']-27.5))
                fluxerr= k['mag_err']*10**(-0.4*(k['mag']-27.5))
            else:
                flux = -10**(-0.4*(-k['mag']-27.5))
                fluxerr= k['mag_err']*10**(-0.4*(-k['mag']-27.5))
            mag = k['mag']
            magerr = k['mag_err']
            f.write('OBS: ' + str(k['jd']) +' '+ flt+' '+ str(flux)+ ' '+str(fluxerr)+' '+ str(mag)+' '+ str(magerr)+' 0 \n')

    os.system('python ./yse/uploadTransientData.py -e -s ./yse/settings.ini -i '+outname+' --instrument ACAM1 --fluxzpt 27.5')



#if __name__ == "__main__":
def main(argv):
#    print(argv)
#    ra, dec, name, date= argv[0:]
    ra = float(argv[0])
    dec = float(argv[1])
    date = argv[3]
    name = argv[2]
    disc_t = 0.

    
    home_dir = lc_cfg['home_dir']+name+'/'


    out_fig = home_dir + name + date + '_lc.'+lc_cfg['img_suffix'] 
    atlas_data_file = home_dir+name+date+'_atlas.csv'

#    print([ra, dec, '3', name + date+'_texas'])
#    if os.path.exists(home_dir+name+'_texas.txt'):
#        galcan = ascii.read(home_dir+name+'_texas.txt')
#    else:
#        try:
#            galcan = texas.main([argv[0], argv[1], '3', home_dir + name + '_texas'])
#            ascii.write(galcan, home_dir+name+'_texas.txt', overwrite=True)  
#        except:
#            print('unable to access PanSTARRS')
        
    ra = float(ra)
    dec = float(dec)  
#    tess_cover = False  
    if len(argv)>4: 
        disc_t = argv[4]
    else:
        disc_t = Time.now().jd
        
    if not os.path.exists(out_fig):
        tess_cover = lc(ra, dec, out_fig, atlas_data_file, disc_t)
        atlas2yse(name, ra, dec, atlas_data_file)
    else:    
        tess_ob = tess_obs(ra, dec)
        tess_cover = False
        for [t1, t2] in tess_ob:
            if disc_t>t1 and disc_t<t2:
                tess_cover = True

#    with open(home_dir+name+date+'_texas.txt', 'w+') as f:
#        for item in galcan:
#            f.write("%s\n" % item)
 #   if len(galcan)>0:
#        return tess_cover
 #   else:
    return tess_cover
 

if __name__ == "__main__":
    main(sys.argv[1:])

