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

def tess_obs(ra, dec):

    tess_date = [2458324.5,2458352.5,2458381.5,2458409.5,2458437.5,2458463.5,
        2458490.5,2458516.5,2458542.5,2458568.5,2458595.5,2458624.5,
        2458653.5,2458682.5,2458710.5,2458737.5,2458763.5,2458789.5,
        2458814.5,2458841.5,2458869.5,2458897.5,2458926.5,2458955.5,
        2458982.5,2459008.5,2459034.5]

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
    gateway_session1 = SSHSession('10.162.61.62','minix', password='123456').open()
    gateway_session2 = gateway_session1.get_remote_session('atlas-base-adm01.ifa.hawaii.edu', username = 'qinan' ,password='PRBNAhK93bMc88Qn')
    remote_session = gateway_session2.get_remote_session('atlas-base-db07.ifa.hawaii.edu',password='PRBNAhK93bMc88Qn')

    today = dt.today()
    con = sqlite3.connect(":memory:")
    tdate = ' '+str(list(con.execute("select julianday('"+today.strftime("%Y-%m-%d")+"')"))[0][0]-40-2400000)

    result=remote_session.get_cmd_output('./mod_force.sh '+str(ra)+' '+ str(dec)+tdate)

    atlas_lc = at.Table(names=('jd','mag', 'mag_err', 'filter'), dtype=('f8', 'f8', 'f8', 'S1'))

    split = result.split('\r\n')

    for i in split[1:]:
        k = i.split()
        if int(k[3])>int(k[4]):
            atlas_lc.add_row([float(k[0]), float(k[1]), float(k[2]), k[5]])

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
    
def plot(lc, survey):
    
    for i in set(lc['filter']):
        jd = lc['jd'][lc['filter']==i]
        if jd[0] > 2400000:
            jd = jd -2400000
        mag = lc['mag'][lc['filter']==i]
        err = lc['mag_err'][lc['filter']==i]
        plt.errorbar(jd, mag, err, label = survey + i, fmt='o')

    ax = plt.gca()
    ax.set_ylim(max(lc['mag'])+0.5, min(lc['mag'])-0.5)

    #ax.grid(True)
    return ax

def Save_space(Save):
    """
    Creates a pathm if it doesn't already exist.
    """
    try:
        if not os.path.exists(Save):
            os.makedirs(Save)
    except FileExistsError:
        pass

#if __name__ == "__main__":
def main(argv):
#    print(argv)
#    ra, dec, name, date= argv[0:]
    ra = float(argv[0])
    dec = float(argv[1])
    date = argv[3]
    name = argv[2]
    
    home_dir = "./plots/"+name+'/'
    Save_space(home_dir)

    out_fig = home_dir + name + date + '_lc'

#    print([ra, dec, '3', name + date+'_texas'])
    galcan = texas.main([argv[0], argv[1], '3', home_dir + name + date+'_texas'])

    ra = float(ra)
    dec = float(dec)
    
    ztf_obj = get_ztf(ra, dec)

    lc = ztf2lc(ztf_obj)
    
    fig, ax = plt.subplots()

    if len(lc)>0:
        ax = plot(lc, 'ztf ')

    try:
        atlas_lc = search_atlas(ra, dec)
        if len(atlas_lc)>0:
            ax = plot(atlas_lc, 'atlas ')
    except:
        print('unable to get atlas lc')

    
    tess_ob = tess_obs(ra, dec)
    i = 0
    for [t1, t2] in tess_ob:
        x= np.arange(t1-2400000, t2-2400000, 0.1)
        if i == 0:
            ax.fill_between(x, 10, 22, facecolor='grey', alpha=0.5, label = 'TESS')
            i += 1
        else:
            ax.fill_between(x, 10, 22, facecolor='grey', alpha=0.5)
    #print(tess_ob)
    ax.legend()

    today = dt.today()
    con = sqlite3.connect(":memory:")
    tdate = list(con.execute("select julianday('"+today.strftime("%Y-%m-%d")+"')"))[0][0]-2400000
    ax.axvline(x=tdate, color = 'k', label = 'today')
    ax.set_xlim(tdate-40, tdate+10)
    ax.legend()

    fig.savefig(out_fig)
    with open(home_dir+name+date+'_texas.txt', 'w+') as f:
        for item in galcan:
            f.write("%s\n" % item)
    

 

if __name__ == "__main__":
    main(sys.argv[1:])

