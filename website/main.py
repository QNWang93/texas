import lc_ex
import candidate_generator
from astropy.io import ascii
import datetime
import weblesniff
import sys

if __name__ == "__main__":
 #   start = sys.argv[1]
  #  end  = sys.argv[2]
    candidate_generator.Update_sheet()
    cans = ascii.read('candidates.csv')#[start, end]
    today = datetime.datetime.today()
    if today.month<10:
        date = '0'+str(today.month) + str(today.day)
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

    
    
