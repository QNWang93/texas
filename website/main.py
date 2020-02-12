import lc_ex
import candidate_generator
from astropy.io import ascii
import datetime
import weblesniff

if __name__ == "__main__":
    candidate_generator.Update_sheet()
    cans = ascii.read('candidates.csv')
    today = datetime.datetime.today()
    if today.month<10:
        date = '0'+str(today.month) + str(today.day)
    else: 
        date = str(today.month) + str(today.day)
    
    for ta in cans:
        z = lc_ex.main([str(ta['RA']), str(ta['Dec']), ta['Name'], date])
        if z != None:
            ta['Redshift'] = z

    cans.sort(['Redshift','Type'])
    ascii.write(cans, 'candidates_sort.csv',names=cans.colnames, overwrite=True)

    weblesniff.main(['./', date, 'imagelist_template.html', 'png', cans])

    
    
