#!/usr/bin/env python
import lc_ex
import candidate_generator
from astropy.io import ascii
import datetime
import weblesniff
import sys
import gspread
from astropy.table import Table
from oauth2client.service_account import ServiceAccountCredentials



def sheet2aptable():
    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']
    creds = ServiceAccountCredentials.from_json_keyfile_name('client_secret.json', scope)
    client = gspread.authorize(creds)
    sheet = client.open("texas").sheet1
    s = [sheet.col_values(1)[1:]]

    for i in range(len(sheet.row_values(1))-1):
        s = s+[sheet.col_values(i+2)[1:]]

    t = Table(s, names = sheet.row_values(1))
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
    for ta in cans:
        try:
            z = lc_ex.main([str(ta['RA']), str(ta['Dec']), ta['Name'], date])
            if z!= None and z!='null':
                ta['Redshift'] = z
        except:
            continue

    cans.sort(['Redshift','Type'])
    ascii.write(cans, 'candidates_sort.csv',names=cans.colnames, overwrite=True)

    weblesniff.main(['./', date, 'imagelist_template.html', 'png', cans])

    
    
