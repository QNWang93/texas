import lc_ex
import texas
import sys
from astropy.io import ascii

def main(filename):
#    filename = arg[0]
    f = open(filename, 'r')
    for line in f:
        para = line.split()
        print(para)
        para[3] = para[3][:4]
        lc_ex.main(para)
        
def main_texas(filename):
    t = ascii.read(filename)
    for i in t:
        try:
            galcan = texas.main([i['RA'], i['Dec'], '3', './plots/'+name + '_texas'])
            ascii.write(galcan, home_dir+name+'_texas.txt', overwrite=True)  
        else:
            print('unable to access PanSTARRS for '+i['Name'])
            
if __name__ == "__main__":
    main(sys.argv[1])
