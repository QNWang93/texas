from __future__ import print_function
from astropy.table import Tabl
import requests
from PIL import Image
from io import BytesIO
import mastcasjobs
from astropy.io import ascii
from astropy.table import Table, vstack
import sys
import os
import re
import pylab
import json
import astropy
import csv
from astropy.modeling.models import Sersic2D
import matplotlib.pyplot as lt
import numpy as np
from astropy import units as u
from astropy.modeling.models import Ellipse2D
from astropy.coordinates import Angle
import matplotlib.patches as mpatches
from astropy.io import fits
from matplotlib.colors import LogNorm
import random
from astropy.visualization import PercentileInterval, AsinhStretch
import getopt
from matplotlib.backends.backend_pdf import PdfPages

try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve

try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib

# get the WSID and password if not already defined
import getpass





def getimages(ra,dec,size=240,filters="grizy"):
    
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table


def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
    
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).from __future__ import print_function
from astropy.table import Table
import requests
from PIL import Image
from io import BytesIO
import mastcasjobs
from astropy.io import ascii
from astropy.table import Table, vstack
import sys
import os
import re
import pylab
import json
import astropy
import csv
from astropy.modeling.models import Sersic2D
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.modeling.models import Ellipse2D
from astropy.coordinates import Angle
import matplotlib.patches as mpatches
from astropy.io import fits
from matplotlib.colors import LogNorm
import random
from astropy.visualization import PercentileInterval, AsinhStretch
import getopt
from matplotlib.backends.backend_pdf import PdfPages

try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve

try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib

# get the WSID and password if not already defined
import getpass





def getimages(ra,dec,size=240,filters="grizy"):
    
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table


def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
    
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra,dec,size=size,filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size="+str(size)+"&format={format}").format(**locals())
    print(url)
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url


def getcolorim(ra, dec, size=240, output_size=None, filters="grizy", format="jpg"):
    
    """Get color image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png")
    Returns the image
    """
    
    if format not in ("jpg","png"):
        raise ValueError("format must be jpg or png")
    url = geturl(ra,dec,size=size,filters=filters,output_size=output_size,format=format,color=True)
    r = requests.get(url)
    im = Image.open(BytesIO(r.content))
    return im


def getgrayim(ra, dec, size=240, output_size=None, filter="g", format="jpg"):
    
    """Get grayscale image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filter = string with filter to extract (one of grizy)
    format = data format (options are "jpg", "png")
    Returns the image
    """
    
    if format not in ("jpg","png"):
        raise ValueError("format must be jpg or png")
    if filter not in list("grizy"):
        raise ValueError("filter must be one of grizy")
    url = geturl(ra,dec,size=size,filters=filter,output_size=output_size,format=format)
    r = requests.get(url[0])
    im = Image.open(BytesIO(r.content))
    return im

def sourcesearch_texas(ra, dec, radius):#radius in unit of arcmin
    up = int(np.ceil(dec+0.5))
    low= int(np.floor(dec-0.5))
    
    s1=ascii.read('./sky-section/dec_'+str(up-1)+'_'+str(up)+'.txt')
    s2=ascii.read('./sky-section/dec_'+str(low)+'_'+str(low+1)+'.txt')
    
    if low+2!=up:
        print("problematic")
    
    source_table=vstack([s1,s2])

    nearby_source=Table(names=('objID','ra','dec', 'raMean','decMean','raMeanErr','decMeanErr','qualityFlag','nDetections','primaryDetection','bestDetection','gSerRadius','gSerMag','gSerAb','gSerNu','gSerPhi','gSerRa','gSerDec','z'))
    for source in source_table:
        cosdec=np.cos(low*np.pi/180)
        if abs((source['raMean']-ra)*cosdec)*60<radius/1.7 and abs((source['decMean']-dec))*60<radius/1.7:
            nearby_source.add_row(source)
    
    
    return nearby_source

def sourcesearch_glade(ra, dec, radius):#radius in unit of arcmin
    up = int(np.ceil(dec+0.5))
    low= int(np.floor(dec-0.5))
    
    s1=ascii.read('../blocks/dec_'+str(up-1)+'_'+str(up)+'.txt')
    s2=ascii.read('../blocks/dec_'+str(low)+'_'+str(low+1)+'.txt')
    
    if low+2!=up:
        print("problematic")
    
    source_table=vstack([s1,s2])

    nearby_source=Table(s1[1:2])
    nearby_source.remove_row(0)
    
    for source in source_table:
        cosdec=np.cos(low*np.pi/180)
        if abs((source[6]-ra)*cosdec)*60<radius/1.7 and abs((source[7]-dec))*60<radius/1.7:
            nearby_source.add_row(source)
    
    
    return nearby_source


def nor_sep(ra, dec, raMean, decMean, a1, b1, theta1):

    d = np.sqrt((dec-decMean)*(dec-decMean) + (ra-raMean)*(ra-raMean)*np.cos(dec/180*np.pi)*np.cos(dec/180*np.pi))
    theta2=np.arctan((dec-decMean)/(ra-raMean))*180./np.pi
    if raMean>ra:
        theta=theta1+theta2
    else: 
        theta=180-theta1-theta2
    
    if theta==90. or theta==-90.:
        theta = 89.9999*abs(theta)/theta

    r=a1*b1*np.sqrt((1+np.tan(theta/180*np.pi)*np.tan(theta/180*np.pi))/(b1*b1+a1*a1*np.tan(theta/180*np.pi)*np.tan(theta/180*np.pi)))
    
    nor_d = d*3600/r
    
    return nor_d

def rearrange(s_list, column_name):
    for i in np.arange(len(s_list)):
        for j in np.arange(i):
            if s_list[j][column_name]>s_list[i][column_name]:   
                k = Table(s_list[j])
                l = Table(s_list[i])
                s_list[j]=l[0]
                s_list[i]=k[0]
    return s_list

def ser_rearrange(s_list, ra, dec):
    i=0
    if len(s_list) == 0:
        return None

    s_list['norm_dist'] = 999999.
    for s in s_list :

        i=i+1
        
        #make fitted image
        n_radius=2
        theta1 = s['gSerPhi']#rot angle
        a1= s['gSerRadius'] 
        b1= a1*s['gSerAb']

        nor_d = nor_sep(ra, dec, s['raMean'], s['decMean'], a1, b1, theta1)
        s['norm_dist'] = nor_d

        
    s_list = rearrange(s_list, 'norm_dist')
    
    return s_list

def plot_ellipse(s_list, ra, dec, size, color, ax):
    i=0
    if len(s_list) == 0:
        return None

    for s in s_list :
        x0, y0 = ((ra-s['raMean'])*4*3600*np.cos(s['decMean']/180*np.pi)+(size/2)), (s['decMean']-dec)*4*3600+(size/2)
        i=i+1
        
        y, x = np.mgrid[0:size, 0:size]# 4 pixel for 1 arcsec for PS1, here image size is set to be 20"*20", depend on your cutout image size
        #make fitted image
        n_radius=2
        theta1 = s['gSerPhi']#rot angle
        a1= s['gSerRadius'] 
        b1= a1*s['gSerAb']
        e1 = mpatches.Ellipse((x0, y0), 2*4*n_radius*a1, 2*4*n_radius*b1, theta1, edgecolor=color,
                              facecolor='none',  label='source 1')    # 4pix/arcsec * n_radius*a, 4pix/arcsec * n_radius*a*b/a 
        if 'z' in s_list.colnames:
            if x0<0 or x0>size or y0<0 or y0>size:
                print(x0, y0)
                print('missing galaxy'+str(i))
                continue
#            ax.annotate('z='+str(int(s['z']*100000)/100000.0), xy=(x0+10, y0+35), fontsize=15, ha="center", color='purple')
#            print('redshift of texas host'+str(i)+':   '+str(int(s['z']*1000)/1000.0))

        else:
            if x0<0 or x0>size or y0<0 or y0>size or s['gSerRadius']<0 or s['gSerChisq']>100:
                continue
#        ax.plot([size/2, x0], [size/2, y0], 'k', ls=':')
            
        ax.add_patch(e1)

    for j in range(len(s_list)):
        s=s_list[j]
        x0, y0 = ((ra-s['raMean'])*4*3600*np.cos(s['decMean']/180*np.pi)+(size/2)), (s['decMean']-dec)*4*3600+(size/2)
        ax.annotate(str(j+1), xy=(x0+20, y0+20), fontsize=15, ha="center", color=color)
    


def mastQuery(request):
    """Perform a MAST query.

    Parameters
    ----------
    request (dictionary): The MAST request json object

    Returns head,content where head is the response HTTP headers, and content is the returned data"""
    
    server='mast.stsci.edu'

    # Grab Python Version 
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent":"python-requests/"+version}

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)
    
    # opening the https connection
    conn = httplib.HTTPSConnection(server)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')

    # Close the https connection
    conn.close()

    return head,content


def resolve(name):
    """Get the RA and Dec for an object using the MAST name resolver
    
    Parameters
    ----------
    name (str): Name of object

    Returns RA, Dec tuple with position"""

    resolverRequest = {'service':'Mast.Name.Lookup',
                       'params':{'input':name,
                                 'format':'json'
                                },
                      }
    headers,resolvedObjectString = mastQuery(resolverRequest)
    resolvedObject = json.loads(resolvedObjectString)
    # The resolver returns a variety of information about the resolved object, 
    # however for our purposes all we need are the RA and Dec
    try:
        objRa = resolvedObject['resolvedCoordinate'][0]['ra']
        objDec = resolvedObject['resolvedCoordinate'][0]['decl']
    except IndexError as e:
        raise ValueError("Unknown object '{}'".format(name))
    return (objRa, objDec)


def fixcolnames(tab):
    """Fix column names returned by the casjobs query
    
    Parameters
    ----------
    tab (astropy.table.Table): Input table

    Returns reference to original table with column names modified"""

    pat = re.compile(r'\[(?P<name>[^[]+)\]')
    for c in tab.colnames:
        m = pat.match(c)
        if not m:
            raise ValueError("Unable to parse column name '{}'".format(c))
        newname = m.group('name')
        tab.rename_column(c,newname)
    return tab


def search_ser(ra, dec, search_size):
    query = """select o.raMean, o.decMean, s.*
    from fGetNearbyObjEq("""+str(ra)+','+str(dec)+","+str(search_size/2)+""") nb
    JOIN MeanObjectView o on o.ObjID=nb.ObjID
    JOIN StackObjectAttributes AS soa ON soa.ObjID = nb.ObjID
    JOIN StackModelFitSer s ON (s.gstackDetectID=soa.gstackDetectID AND s.ObjID=nb.ObjId)
    WHERE o.nDetections>5
    AND soa.primaryDetection>0
    AND (o.rmeanpsfmag - o.rmeankronmag > 0.05)
    """

    jobs = mastcasjobs.MastCasJobs(context="PanSTARRS_DR2")
    results = jobs.quick(query, task_name="python cone search")
    return(results)

def search_s(ra, dec, search_size):
    query = """select o.raMean, o.decMean
    from fGetNearbyObjEq("""+str(ra)+','+str(dec)+","+str(search_size/2)+""") nb
    JOIN MeanObjectView o on o.ObjID=nb.ObjID
    WHERE o.nDetections>5
    """

    jobs = mastcasjobs.MastCasJobs(context="PanSTARRS_DR2")
    results = jobs.quick(query, task_name="python cone search")
    return(results)


def plot(ra, dec, bfim, gal_list, s_list, ser_list, catalogue, search_size, filename):

    size = 240*search_size
    fig, ax = plt.subplots(1,1, figsize=(search_size*3,search_size*3))
    plt.imshow(bfim, cmap='summer')#, norm=LogNorm())

    if catalogue == 'glade':
        for j in gal_list:
            x, y = ((ra-j[6])*4*3600*np.cos(j[7]/180*np.pi)+(size/2)), (j[7]-dec)*4*3600+(size/2)
            ax.plot(x,y, 'bx', label = 'glade source')
            ax.annotate('with z '+str(j[10]), xy=(x+20, y-5), fontsize=15, ha="center", color='b')
    elif catalogue == 'texas':
        plot_ellipse(gal_list, ra, dec, size, 'k', ax)#, label = 'texas source')

    for s in s_list:
        x, y = ((ra-s['raMean'])*4*3600*np.cos(s['decMean']/180*np.pi)+(size/2)), (s['decMean']-dec)*4*3600+(size/2)
        ax.plot(x,y, 'rx', label = 'PS point source')
	    

    plot_ellipse(ser_list, ra, dec, size, 'm', ax)

    ax.plot([size/2+10, size/2+20], [size/2, size/2], 'k')
    ax.plot([size/2, size/2], [size/2+10, size/2+20], 'k')
    ax.plot([size/2-10, size/2-20], [size/2, size/2], 'k')
    ax.plot([size/2, size/2], [size/2-10, size/2-20], 'k')

    plt.savefig(filename)

def main(argv):
    try:                          
        opts, args = getopt.getopt(argv, "hi:c:n:p:d", ["help", "catalogue", "PS_WSID", "PS_password"])      
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    catalogue='glade'
    do_im = True
    os.environ['CASJOBS_WSID'] = '1827319063'
    os.environ['CASJOBS_PW'] = '123456789'

    for opt, arg in opts:                
        if opt in ("-h", "--help"):       
            usage()                     
            sys.exit()                  
        elif opt == '-d':                
            global _debug               
            _debug = 1                  
        elif opt in ("-c", "--catalogue="): 
            catalogue = arg               
        elif opt in ('-n', 'PS_WSID'):                
            os.environ['CASJOBS_WSID'] = arg
        elif opt in ('-p', 'PS_password'):                
            os.environ['CASJOBS_PW'] = arg
        elif opt in ('-i', 'do_im'):                
            do_im = arg
    ra = float(argv[0])
    dec = float(argv[1])
    search_size = int(argv[2])
    if do_im:
        filename = argv[3]
#    print(argv)
    size = 240*search_size  #PS cutout image size, 240*sidelength in arcmin
    galac_search_size = 1

    #get PS image
	 
    #downloading image part
    fitsurl = geturl(ra, dec, size=240*search_size, filters="i", format="fits")
    fh = fits.open(fitsurl[0])
    fim = fh[0].data
    fim[np.isnan(fim)] = 0.
    transform = AsinhStretch() + PercentileInterval(99)
    bfim = transform(fim)
    fig, ax = plt.subplots(1,1, figsize=(search_size*3,search_size*3))
    plt.imshow(bfim, cmap='summer')#, norm=LogNorm())
	#ax.plot(size/2,size/2, 'ko')
	#plot image

    if catalogue == 'glade':
        gal_list=sourcesearch_glade(ra,dec, size/240)#search in radius of cutout size
    elif catalogue == 'texas':
        gal_list=sourcesearch_texas(ra,dec, size/240)
    else:
        print('no this catalogue')
	#check corrsponding block of sources
	#plot sources onto image


    ser_list = search_ser(ra, dec, search_size)
    #print(ser_list)
    i=0
    if len(ser_list)>0:
        ser_list = ser_rearrange(ser_list, ra, dec)
        ser_list['z']=-999.


    gal_can = []
    #print(gal_list)
    if catalogue == 'glade':
        for j in gal_list:
            for s in ser_list:
                if abs((j[6]-s['raMean'])*np.cos(s['decMean']))<0.001 and abs(j[7]-s['decMean'])<0.001 and j[10] != 'null':
                    s['z'] = j[10]
                    gal_can.append(s)
    elif catalogue == 'texas':
        plot_ellipse(gal_list, ra, dec, size, 'k', ax)#, label = 'texas source')
    
    txt = ""
    for i in np.arange(len(ser_list)):
        k = ser_list[i]
        if k['z']<0:
            txt=txt+('"/n" normalized distance to host'+str(i+1)+':   '+str(int(k['norm_dist']*1000)/1000.0))
            print('normalized distance to host'+str(i+1)+':   '+str(int(k['norm_dist']*1000)/1000.0))
        else:
            txt=txt+('"/n" normalized distance to host'+str(i+1)+':   '+str(int(k['norm_dist']*1000)/1000.0)+' , z='+str(k['z']))
            print('normalized distance to host'+str(i+1)+':   '+str(int(k['norm_dist']*1000)/1000.0)+' , z='+str(k['z']))


    s_list = search_s(ra, dec, galac_search_size)
    
    if do_im:
        plot(ra, dec, bfim, gal_list, s_list, ser_list, catalogue, search_size, filename)


	#ax.legend()
	    
	#plt.show()
#    fig.text(.5, .05, txt, ha='center',wrap=True)
#    fig.text(-0.5, .05, '432423424214121342', ha='center',wrap=True)

    return(gal_can)
	#append necessary information onto plot, modified distance, redshift etc9


if __name__ == "__main__":
    main(sys.argv[1:])
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra,dec,size=size,filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size="+str(size)+"&format={format}").format(**locals())
    print(url)
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url


def getcolorim(ra, dec, size=240, output_size=None, filters="grizy", format="jpg"):
    
    """Get color image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png")
    Returns the image
    """
    
    if format not in ("jpg","png"):
        raise ValueError("format must be jpg or png")
    url = geturl(ra,dec,size=size,filters=filters,output_size=output_size,format=format,color=True)
    r = requests.get(url)
    im = Image.open(BytesIO(r.content))
    return im


def getgrayim(ra, dec, size=240, output_size=None, filter="g", format="jpg"):
    
    """Get grayscale image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filter = string with filter to extract (one of grizy)
    format = data format (options are "jpg", "png")
    Returns the image
    """
    
    if format not in ("jpg","png"):
        raise ValueError("format must be jpg or png")
    if filter not in list("grizy"):
        raise ValueError("filter must be one of grizy")
    url = geturl(ra,dec,size=size,filters=filter,output_size=output_size,format=format)
    r = requests.get(url[0])
    im = Image.open(BytesIO(r.content))
    return im

def sourcesearch_texas(ra, dec, radius):#radius in unit of arcmin
    up = int(np.ceil(dec+0.5))
    low= int(np.floor(dec-0.5))
    
    s1=ascii.read('./sky-section/dec_'+str(up-1)+'_'+str(up)+'.txt')
    s2=ascii.read('./sky-section/dec_'+str(low)+'_'+str(low+1)+'.txt')
    
    if low+2!=up:
        print("problematic")
    
    source_table=vstack([s1,s2])

    nearby_source=Table(names=('objID','ra','dec', 'raMean','decMean','raMeanErr','decMeanErr','qualityFlag','nDetections','primaryDetection','bestDetection','gSerRadius','gSerMag','gSerAb','gSerNu','gSerPhi','gSerRa','gSerDec','z'))
    for source in source_table:
        cosdec=np.cos(low*np.pi/180)
        if abs((source['raMean']-ra)*cosdec)*60<radius/1.7 and abs((source['decMean']-dec))*60<radius/1.7:
            nearby_source.add_row(source)
    
    
    return nearby_source

def sourcesearch_glade(ra, dec, radius):#radius in unit of arcmin
    up = int(np.ceil(dec+0.5))
    low= int(np.floor(dec-0.5))
    
    s1=ascii.read('../blocks/dec_'+str(up-1)+'_'+str(up)+'.txt')
    s2=ascii.read('../blocks/dec_'+str(low)+'_'+str(low+1)+'.txt')
    
    if low+2!=up:
        print("problematic")
    
    source_table=vstack([s1,s2])

    nearby_source=Table(s1[1:2])
    nearby_source.remove_row(0)
    
    for source in source_table:
        cosdec=np.cos(low*np.pi/180)
        if abs((source[6]-ra)*cosdec)*60<radius/1.7 and abs((source[7]-dec))*60<radius/1.7:
            nearby_source.add_row(source)
    
    
    return nearby_source


def nor_sep(ra, dec, raMean, decMean, a1, b1, theta1):

    d = np.sqrt((dec-decMean)*(dec-decMean) + (ra-raMean)*(ra-raMean)*np.cos(dec/180*np.pi)*np.cos(dec/180*np.pi))
    theta2=np.arctan((dec-decMean)/(ra-raMean))*180./np.pi
    if raMean>ra:
        theta=theta1+theta2
    else: 
        theta=180-theta1-theta2
    
    if theta==90. or theta==-90.:
        theta = 89.9999*abs(theta)/theta

    r=a1*b1*np.sqrt((1+np.tan(theta/180*np.pi)*np.tan(theta/180*np.pi))/(b1*b1+a1*a1*np.tan(theta/180*np.pi)*np.tan(theta/180*np.pi)))
    
    nor_d = d*3600/r
    
    return nor_d

def rearrange(s_list, column_name):
    for i in np.arange(len(s_list)):
        for j in np.arange(i):
            if s_list[j][column_name]>s_list[i][column_name]:   
                k = Table(s_list[j])
                l = Table(s_list[i])
                s_list[j]=l[0]
                s_list[i]=k[0]
    return s_list


def plot_ellipse(s_list, ra, dec, size, color, ax):
    i=0
    if len(s_list) == 0:
        return None

    s_list['norm_dist'] = 999999.
    for s in s_list :
        x0, y0 = ((ra-s['raMean'])*4*3600*np.cos(s['decMean']/180*np.pi)+(size/2)), (s['decMean']-dec)*4*3600+(size/2)
        i=i+1
        
        y, x = np.mgrid[0:size, 0:size]# 4 pixel for 1 arcsec for PS1, here image size is set to be 20"*20", depend on your cutout image size
        #make fitted image
        n_radius=2
        theta1 = s['gSerPhi']#rot angle
        a1= s['gSerRadius'] 
        b1= a1*s['gSerAb']
        e1 = mpatches.Ellipse((x0, y0), 2*4*n_radius*a1, 2*4*n_radius*b1, theta1, edgecolor=color,
                              facecolor='none',  label='source 1')    # 4pix/arcsec * n_radius*a, 4pix/arcsec * n_radius*a*b/a 
        if 'z' in s_list.colnames:
            if x0<0 or x0>size or y0<0 or y0>size:
                print(x0, y0)
                print('missing galaxy'+str(i))
                continue
#            ax.annotate('z='+str(int(s['z']*100000)/100000.0), xy=(x0+10, y0+35), fontsize=15, ha="center", color='purple')
            print('redshift of texas host'+str(i)+':   '+str(int(s['z']*1000)/1000.0))

        else:
            if x0<0 or x0>size or y0<0 or y0>size or s['gSerRadius']<0 or s['gSerChisq']>100:
                continue
#        ax.plot([size/2, x0], [size/2, y0], 'k', ls=':')
        nor_d = nor_sep(ra, dec, s['raMean'], s['decMean'], a1, b1, theta1)
        s['norm_dist'] = nor_d
        angle = -np.arctan((size/2-y0)/(size/2-x0))*180/np.pi
#        th1 = plt.text(size/4+x0/2,size/4+y0/2, 'norm-dist:'+str(int(nor_d*1000)/1000.0), fontsize=12,
#                    rotation=angle, rotation_mode='anchor', color = 'purple', horizontalalignment='center')
            
        ax.add_patch(e1)
        
        
    s_list = rearrange(s_list, 'norm_dist')

    for j in range(len(s_list)):
        s=s_list[j]
        x0, y0 = ((ra-s['raMean'])*4*3600*np.cos(s['decMean']/180*np.pi)+(size/2)), (s['decMean']-dec)*4*3600+(size/2)
        ax.annotate(str(j+1), xy=(x0+20, y0+20), fontsize=15, ha="center", color=color)
    
    return s_list

def mastQuery(request):
    """Perform a MAST query.

    Parameters
    ----------
    request (dictionary): The MAST request json object

    Returns head,content where head is the response HTTP headers, and content is the returned data"""
    
    server='mast.stsci.edu'

    # Grab Python Version 
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent":"python-requests/"+version}

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)
    
    # opening the https connection
    conn = httplib.HTTPSConnection(server)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')

    # Close the https connection
    conn.close()

    return head,content


def resolve(name):
    """Get the RA and Dec for an object using the MAST name resolver
    
    Parameters
    ----------
    name (str): Name of object

    Returns RA, Dec tuple with position"""

    resolverRequest = {'service':'Mast.Name.Lookup',
                       'params':{'input':name,
                                 'format':'json'
                                },
                      }
    headers,resolvedObjectString = mastQuery(resolverRequest)
    resolvedObject = json.loads(resolvedObjectString)
    # The resolver returns a variety of information about the resolved object, 
    # however for our purposes all we need are the RA and Dec
    try:
        objRa = resolvedObject['resolvedCoordinate'][0]['ra']
        objDec = resolvedObject['resolvedCoordinate'][0]['decl']
    except IndexError as e:
        raise ValueError("Unknown object '{}'".format(name))
    return (objRa, objDec)


def fixcolnames(tab):
    """Fix column names returned by the casjobs query
    
    Parameters
    ----------
    tab (astropy.table.Table): Input table

    Returns reference to original table with column names modified"""

    pat = re.compile(r'\[(?P<name>[^[]+)\]')
    for c in tab.colnames:
        m = pat.match(c)
        if not m:
            raise ValueError("Unable to parse column name '{}'".format(c))
        newname = m.group('name')
        tab.rename_column(c,newname)
    return tab


def search_ser(ra, dec, search_size):
    query = """select o.raMean, o.decMean, s.*
    from fGetNearbyObjEq("""+str(ra)+','+str(dec)+","+str(search_size/2)+""") nb
    JOIN MeanObjectView o on o.ObjID=nb.ObjID
    JOIN StackObjectAttributes AS soa ON soa.ObjID = nb.ObjID
    JOIN StackModelFitSer s ON (s.gstackDetectID=soa.gstackDetectID AND s.ObjID=nb.ObjId)
    WHERE o.nDetections>5
    AND soa.primaryDetection>0
    AND (o.rmeanpsfmag - o.rmeankronmag > 0.05)
    """

    jobs = mastcasjobs.MastCasJobs(context="PanSTARRS_DR2")
    results = jobs.quick(query, task_name="python cone search")
    return(results)

def search_s(ra, dec, search_size):
    query = """select o.raMean, o.decMean
    from fGetNearbyObjEq("""+str(ra)+','+str(dec)+","+str(search_size/2)+""") nb
    JOIN MeanObjectView o on o.ObjID=nb.ObjID
    WHERE o.nDetections>5
    """

    jobs = mastcasjobs.MastCasJobs(context="PanSTARRS_DR2")
    results = jobs.quick(query, task_name="python cone search")
    return(results)





def main(argv):
    try:                          
        opts, args = getopt.getopt(argv, "hc:n:p:d", ["help", "catalogue", "PS_WSID", "PS_password"])      
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    catalogue='glade'
    os.environ['CASJOBS_WSID'] = '1827319063'
    os.environ['CASJOBS_PW'] = '123456789'

    for opt, arg in opts:                
        if opt in ("-h", "--help"):       
            usage()                     
            sys.exit()                  
        elif opt == '-d':                
            global _debug               
            _debug = 1                  
        elif opt in ("-c", "--catalogue="): 
            catalogue = arg               
        elif opt in ('-n', 'PS_WSID'):                
            os.environ['CASJOBS_WSID'] = arg
        elif opt in ('-p', 'PS_password'):                
            os.environ['CASJOBS_PW'] = arg
    ra = float(argv[0])
    dec = float(argv[1])
    search_size = int(argv[2])
    filename = argv[3]
#    print(argv)
    size = 240*search_size  #PS cutout image size, 240*sidelength in arcmin
    galac_search_size = 1

    #get PS image
	 
    #downloading image part
    fitsurl = geturl(ra, dec, size=240*search_size, filters="i", format="fits")
    fh = fits.open(fitsurl[0])
    fim = fh[0].data
    fim[np.isnan(fim)] = 0.
    transform = AsinhStretch() + PercentileInterval(99)
    bfim = transform(fim)
    fig, ax = plt.subplots(1,1, figsize=(search_size*3,search_size*3))
    plt.imshow(bfim, cmap='summer')#, norm=LogNorm())
	#ax.plot(size/2,size/2, 'ko')
	#plot image

    if catalogue == 'glade':
        gal_list=sourcesearch_glade(ra,dec, size/240)#search in radius of cutout size
    elif catalogue == 'texas':
        gal_list=sourcesearch_texas(ra,dec, size/240)
    else:
        print('no this catalogue')
	#check corrsponding block of sources
	#plot sources onto image


    ser_list = search_ser(ra, dec, search_size)
    #print(ser_list)
    if len(ser_list)>0:
        ser_list=plot_ellipse(ser_list, ra, dec, size, 'purple', ax)#, label = 'PS sersic fitting')
        ser_list['z']=-999.

    gal_can = []
    #print(gal_list)
    if catalogue == 'glade':
        for j in gal_list:
            x, y = ((ra-j[6])*4*3600*np.cos(j[7]/180*np.pi)+(size/2)), (j[7]-dec)*4*3600+(size/2)
            for s in ser_list:
                if abs((j[6]-s['raMean'])*np.cos(s['decMean']))<0.001 and abs(j[7]-s['decMean'])<0.001 and j[10] != 'null':
                    s['z'] = j[10]
                    gal_can.append(s)
            ax.plot(x,y, 'bx', label = 'glade source')
            ax.annotate('with z '+str(j[10]), xy=(x+20, y-5), fontsize=15, ha="center", color='b')
    elif catalogue == 'texas':
        plot_ellipse(gal_list, ra, dec, size, 'k', ax)#, label = 'texas source')
    
    txt = ""
    for i in np.arange(len(ser_list)):
        k = ser_list[i]
        if k['z']<0:
            txt=txt+('"/n" normalized distance to host'+str(i+1)+':   '+str(int(k['norm_dist']*1000)/1000.0))
            print('normalized distance to host'+str(i+1)+':   '+str(int(k['norm_dist']*1000)/1000.0))
        else:
            txt=txt+('"/n" normalized distance to host'+str(i+1)+':   '+str(int(k['norm_dist']*1000)/1000.0)+' , z='+str(k['z']))
            print('normalized distance to host'+str(i+1)+':   '+str(int(k['norm_dist']*1000)/1000.0)+' , z='+str(k['z']))


    s_list = search_s(ra, dec, galac_search_size)


    for s in s_list:
        x, y = ((ra-s['raMean'])*4*3600*np.cos(s['decMean']/180*np.pi)+(size/2)), (s['decMean']-dec)*4*3600+(size/2)
        ax.plot(x,y, 'rx', label = 'PS point source')
	    

	    
    ax.plot([size/2+10, size/2+20], [size/2, size/2], 'k')
    ax.plot([size/2, size/2], [size/2+10, size/2+20], 'k')
    ax.plot([size/2-10, size/2-20], [size/2, size/2], 'k')
    ax.plot([size/2, size/2], [size/2-10, size/2-20], 'k')

	#ax.legend()
	    
	#plt.show()
#    fig.text(.5, .05, txt, ha='center',wrap=True)
#    fig.text(-0.5, .05, '432423424214121342', ha='center',wrap=True)
    plt.savefig(filename)
    return(gal_can)
	#append necessary information onto plot, modified distance, redshift etc9


if __name__ == "__main__":
    main(sys.argv[1:])
