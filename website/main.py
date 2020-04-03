#!/usr/bin/env python
import lc_ex
import candidate_generator
from astropy.io import ascii
import datetime
import weblesniff
import sys
import gspread
from astropy.time import Time
from astropy.table import Table, Column
from oauth2client.service_account import ServiceAccountCredentials
from numpy import dtype
from config import web_cfg, lc_cfg
#import request, re
import candidate_generator


def sheet2aptable():

    sheet=candidate_generator.get_gglsheet()
    s = [sheet.col_values(1)[1:]]
    header = [sheet.col_values(1)[0]]
    t = Table(s, names = header, dtype = [web_cfg['dtype'][header[0]]])

   

    for i in range(len(sheet.row_values(1))-1):
        if sheet.row_values(1)[i+1] in web_cfg['header']:
            header = sheet.row_values(1)[i+1]
            t.add_column(Column(sheet.col_values(i+2)[1:]), name = header)

    return t
    


if __name__ == "__main__":
 #   start = sys.argv[1]
  #  end  = sys.argv[2]
    candidate_generator.Update_sheet()
#    cans = ascii.read('candidates.csv')#[start, end]
    cans = sheet2aptable()
    today = datetime.datetime.today()
    if today.month<10:
        date = '0'+str(today.month) + str(today.day)
    elif today.day<10: 
        date = str(today.month) +'0'+ str(today.day)
    else:
        date = str(today.month) + str(today.day)
    
    for h in web_cfg['header'][:-3]: 
        if not (h in cans.colnames): 
            cans.add_column(Column([web_cfg['default'][h]]*len(cans)), name = h)
    l = len(cans)
    k = 1
    for ta in cans:
        try:
            disc_date = Time(ta["Disc date"], format = 'iso', scale = 'utc')
            print('start on '+ta['Name'] + ', '+str(k)+'/'+str(l))
            k+=1
            name = ta['Name']
            ta['TESS_coverage'] = lc_ex.main([str(ta['RA']), str(ta['Dec']), ta['Name'], date, disc_date.jd])
            hosts = ascii.read(lc_cfg['home_dir']+name+'/'+ name + '_texas.txt')
            z = hosts['z'][0]
            source = hosts['source'][0]
            if z!= None and z!='null' and ta['Redshift'] =='None':
                ta['Redshift'] = str(z)
#                ta['z_source'] = source
        except:
            continue
    
#    print(cans)
    cans.sort(['Redshift','Type'])
    ascii.write(cans, 'candidates_sort.csv',names=cans.colnames, overwrite=True)

    weblesniff.main(['./web', 'candidates', 'imagelist_template.html', web_cfg['img_suffix'], cans])

    
    
