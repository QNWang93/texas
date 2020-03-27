#!/usr/bin/env python
import sys, os, re, math, optparse, types, copy, shutil, glob, time, shutil
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
import scipy

# add the necessary python dirs to PATH
if os.environ.has_key('JWSTTOOLS_ROOTDIR'):
    sys.path.append(os.path.join(os.environ['JWSTTOOLS_ROOTDIR'],'pythonmodules'))
else:
    print 'ERROR: environment variable JWSTTOOLS_ROOTDIR is not set!'
    sys.exit(0)

# modules by A. Rest
from texttable import txttableclass
from tools import rmfile,makepath,executecommand,unique,ConfigParser_env,AnotB
from webtools import imagestring4web,addlink2string,addtag2string,htmltable,webpageclass


class qlwebclass(txttableclass):
    def __init__(self):
        txttableclass.__init__(self)
        self.stringleftalign = True

        self.cfg = ConfigParser_env()
        #self.detections = txttableclass()
        #self.catbasename = None
        #self.detectionkeys = None
        
        self.figsuffix = 'jpg'
        
        self.font = 'sans-serif'
        self.pagetablefontsize = 2.7

        self.imagetableheaderfontscale = "+1"
        self.imagetablefontscale = None

        self.bgcolor1imtable = '#%02x%02x%02x' % (255,255,220)
        self.bgcolor4tableheader = '#E0FFFF'
        self.NotAvailable = '-'

        self.jobsperline = 1000

        self.colnamemapping = {}
        self.newfileinfokeys = None
        self.allinfokeys = None
        
        self.thumbnails = False
        self.imagepathprefix = ''

        self.sortcol = 'jobN'

    def loadconfigfile(self,options):
        if not os.path.isfile(options.defaultconfigfile):
            print 'ERROR: config file %s does not exist, exiting!' % options.defaultconfigfile
            sys.exit(0)
        self.cfg.readfp(open(options.defaultconfigfile))
        if options.additionalconfigfile != None:
            self.cfg.read(options.additionalconfigfile)
        
        if options.forcenewinfotable:
            self.cfg.set("webpage","forcenewinfotable","True")

        self.cfg.setvals_nosection(options.params,allflag=True)
        self.cfg.setvals(options.params4sections)
        
        if options.add2search != None:
            s = self.cfg.get("webpage","searchpatterns4pages")
            if s =='':
                s = options.add2search
            else:
                s += ','+options.add2search
            self.cfg.set("webpage","searchpatterns4pages",s)
                

    def add_options(self, parser=None, usage=None):

        # This is directory where the default config files are located
        configdir = os.environ['QUICKLOOK_CONFIGDIR']
        jwsttest = os.environ['JWSTTEST']        
        
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")
        parser.add_option('-v', '--verbose', action="count", dest="verbose",default=0)
        parser.add_option('--defaultconfigfile', default='%s/ql.%s.cfg' % (configdir,jwsttest), type='string',
                          help='config file. (default=%default)')
#        parser.add_option('--defaultfilepattern', default='*.png,*/*.png,*/*/*.png', type="string",
#                          help='default file pattern used if no arguments are given. (default=%default)')
        parser.add_option('-c','--additionalconfigfile', default=None, type="string",
                          help='additional config file, overwriting values in defaultconfigfile. (default=%default)')
        parser.add_option('-p', '--params', action='append', default=None, nargs=2,
                          help='"param val": change parameters in config file (section independent) (default=%default)')
        parser.add_option('-s', '--params4sections', action='append', default=None, nargs=3,
                          help='"section param val". change parameters in section of config file (default=%default)')
        parser.add_option('-a', '--add2search',  default=None, type="string",
                          help='adds searchpatterns to searchpatterns4pages in the config file (default=%default)')
        parser.add_option('-f','--forcenewinfotable',default=False, action='store_true',
                          help='redo the infotable. Takes time, but makes sure that the table is up to date with all fitskeys etc')
        return(parser)

    def findfigs(self, filepatterns):
        figs = []
        os.chdir(self.cfg.getstring("webpage","webrootdir"))
        #print 'FFF',self.cfg.getstring("webpage","webrootdir"),filepatterns
        for pattern in filepatterns:
            if self.verbose: print 'Looking for files matching %s in %s' % (pattern,self.cfg.getstring("webpage","webrootdir"))
            figs.extend(glob.glob(pattern))
        figs = unique(figs)
        figs.sort()

        # strip all the suffixes to get the basename
        suffixpattern = '\.TN\.\w+\.%s|\.\w+\.%s' % (self.cfg.getstring("webpage","figsuffix"),self.cfg.getstring("webpage","figsuffix"))
        basenames = []
        for fig in figs: 
            basename =  re.sub(suffixpattern,'',fig)
            if self.verbose>2: print '%s -> %s' % (fig,basename)
            basenames.append(basename)
        basenames = unique(basenames)
        basenames.sort()
        print '%d figs found, resulting into %d basenames' % (len(figs),len(basenames))
        del figs
        return(basenames)
              
    def searchbasenames(self,basenames,searchpattern):
        goodbasenames = []
        try:
            patternobject = re.compile(searchpattern,flags=re.IGNORECASE)
        except:
            # old python
            patternobject = re.compile(searchpattern)
        for basename in basenames:
            if patternobject.search(basename):
                goodbasenames.append(basename)
        return(goodbasenames)
    
    def searchpatternname(self,searchpattern):
        name=''
        for i in xrange(len(searchpattern)):
            if not(re.search('[A-Za-z0-9]',searchpattern[i])):
                name +='_'
            else:
                name +=searchpattern[i]
        return(name)
    
    def webaddresscompliantname(self,s):
        name=''
        for i in xrange(len(s)):
            if not(re.search('[A-Z|a-z|0-9|\-]',s[i])):
                name +='_'
            else:
                name +=s[i]
        return(name)
    
    def mkBackToMainlinks(self):
        return(addlink2string('Back to Main','index.html'))

    def Noutsubdir(self):
        N = 0
        if self.cfg.get('webpage','outsubdir_depth')!=None and self.cfg.getint('webpage','outsubdir_depth')!=0:
            N = self.cfg.getint('webpage','outsubdir_depth')
            if self.cfg.getboolean('webpage','skiplastsubdir'):
                N-=1
        if self.cfg.getboolean('webpage','taskID2subdir'):
            N += 1
        return(N)
        
    def getfirstlastdate(self,keys,sortedbytime=False,datecol='DATE-OBS',timecol='TIME-OBS'):
        if len(keys)<1:
            return(['?','?','?','?'])
        if not sortedbytime:
            if self.verbose>0: print 'Sorting keys by date/time!!!'            
            keys = self.sortkeysbycols(keys,[datecol,timecol],asstring=1)
        return([self.getentry(keys[0],datecol),self.getentry(keys[0],timecol),self.getentry(keys[-1],datecol),self.getentry(keys[-1],timecol)])
        
    def getcolnamesmapping(self):
        self.colnamemapping = {}
        filename = self.cfg.getstring("webpage","colnamemappingfilename")
        if os.path.isfile(filename):
            if self.verbose>2: print 'Reading column name mapping file ',filename
            lines = open(filename,'r').readlines()
            for line in lines:
                if re.match('#',line):continue
                (col,colstring)=line.split('|',2)
                col=col.strip()
                colstring=colstring.strip()
                self.colnamemapping[col]=colstring

    def initfileinfocols(self,fitskeys=None):
        
        self.configcols(['basename'],'s','%s',visible=1)

        if fitskeys != None:
            self.configcols(fitskeys,'s','%s',visible=1)

        if self.cfg.getboolean("webpage","median2infotable"):
            self.configcols(['median0','mean0','stdev0','medianX','meanX','stdevX','delta_median'],'f','%.2f',visible=1)
            
        
        return(0)

    def getinfotablefilename(self,infofilename=None):
        if infofilename == None: 
            infofilename = self.cfg.getstring("webpage","infotablefilename")
            if infofilename.lower()=='auto':
                infofilename = '%s/fileinfo.txt' % (self.cfg.getstring("webpage","webrootdir"))
        return(infofilename)

    def saveinfofile(self,infofilename=None,keys=None):
        infofilename = self.getinfotablefilename(infofilename)
        self.save2file(infofilename,keys=keys,verbose=1)
            
    def loadinfofile(self,basenames,infofilename=None):
        # read the keywords from previous executions
        infofilename = self.getinfotablefilename(infofilename)
        if os.path.isfile(infofilename) and (not self.cfg.getboolean("webpage","forcenewinfotable")):
            if self.verbose>=1:
                print 'Loading info file %s' % infofilename
            self.loadfile(infofilename)

            newbasenames = AnotB(basenames,self.col_as_list('basename'))
            oldbasenames = AnotB(basenames,newbasenames)
        else:
            newbasenames = basenames
            oldbasenames = AnotB(basenames,newbasenames)
            
        return(newbasenames,oldbasenames)

    def getfileinfo(self, basenames):

        fitskeys = self.cfg.getstring("webpage","fitskeys4infotable").split(',')

        # load file and return the basenames that have no info entry in the infofile yet
        (newbasenames,oldbasenames) = self.loadinfofile(basenames)

        self.initfileinfocols(fitskeys=fitskeys)
        
        if self.verbose>=1:
            print 'Getting info for files from fits headers for %d of %d files...' % (len(newbasenames),len(basenames))

        unknownpattern = re.compile('UNKNOWN')            

        # fill up the new entries only!
        counter = 1
        self.newfileinfokeys=[]
        for basename in newbasenames:
            if self.verbose>=1:
                if (counter % 10)==0:
                    print '%d out of %d' % (counter,len(newbasenames))
                    
                
            # get the fits keywords
            fitsheaderfile = '%s/%s.%s.header.txt' % (self.cfg.getstring("webpage","webrootdir"),basename,self.cfg.getstring("webpage","suffix4fitsheader"))

            filekey = self.newrow({'basename':basename})
            self.newfileinfokeys.append(filekey)

            if not os.path.isfile(fitsheaderfile):
                print 'WARNING!!! file %s does not exist, cannot get file info!!!' % fitsheaderfile
                continue

            header = pyfits.Header.fromstring(''.join(open(fitsheaderfile,'r').readlines()),sep='\n')

            for fitskey in fitskeys:
                val = header.get(fitskey,None)
                
                if (type(val) is types.StringType) and unknownpattern.match(val):
                    val = None
                if (type(val) is types.BooleanType) and val == True:
                    self.data[filekey][fitskey]='True'
                elif (type(val) is types.BooleanType) and  val == False:
                    self.data[filekey][fitskey]='False'
                else:
                    self.data[filekey][fitskey]=val

            # get the median/mean info
            if self.cfg.getboolean("webpage","median2infotable"):
                medianfile = '%s/%s.medianramp.txt' % (self.cfg.getstring("webpage","webrootdir"),basename)
                if os.path.isfile(medianfile):
                    t = txttableclass()
                    t.loadfile(medianfile)
                    if len(t.allrowkeys)>0:
                        self.add2row(filekey,{'median0':t.getentry(t.allrowkeys[0],'median'),
                                              'mean0':t.getentry(t.allrowkeys[0],'mean'),
                                              'stdev0':t.getentry(t.allrowkeys[0],'stdev'),
                                              'medianX':t.getentry(t.allrowkeys[-1],'median'),
                                              'meanX':t.getentry(t.allrowkeys[-1],'mean'),
                                              'stdevX':t.getentry(t.allrowkeys[-1],'stdev')})
                        if t.getentry(t.allrowkeys[-1],'median')!=None and t.getentry(t.allrowkeys[0],'median')!=None:
                            self.setentry(filekey,'delta_median',float(t.getentry(t.allrowkeys[-1],'median'))-float(t.getentry(t.allrowkeys[0],'median')))
            
                    del t
            counter += 1
        
        if self.verbose>=1:
            print '... info from fits headers done!'

        if len(self.newfileinfokeys)>0 and self.verbose>2:
            print 'info for new files:'
            self.printtxttable(keys=self.newfileinfokeys)
        
        
        # now also get the previous keys
        if len(newbasenames)+len(oldbasenames) == len(self.allrowkeys):
            #fully up-to-date, no need to go through everything
            self.allinfokeys = self.allrowkeys[:]
        else:
            temphash = {}
            self.allinfokeys = []
            for key in self.allrowkeys:
                temphash[self.getentry(key,'basename')]=key
            for basename in oldbasenames:
                self.allinfokeys.append(temphash[basename])
            self.allinfokeys.extend(self.newfileinfokeys)

            if self.verbose>=1: print 'Using %d out of %d of all files in fileinfo table' % (len(self.allinfokeys),len(self.allrowkeys))

        
        if len(self.allinfokeys)>0 and self.verbose>2:
            print 'info for all files:'
            self.printtxttable(keys=self.allinfokeys)

        return(0)
            
    def fillup_fileinfo(self):
        # dummy function that can be overwritten by specific instruments to fill up the fileinfo table
        return(0)

    def basenames2fileinfo(self, basenames):
        self.getfileinfo(basenames)
        self.fillup_fileinfo()
        print 'test point 2.5'
        #print self.sortcol
        all_keys_sorted = self.sortkeysbycols(self.rowkeys(),self.sortcol,asstring=1)
        self.saveinfofile(keys= all_keys_sorted)
        print 'test point 3'
        
    def definepages_Nim(self,keys,sortcol,NperPage=50):
        pages = txttableclass()
        pages.configcols(['pageID','N','imin','imax'],'d','%d',visible=1)
        pages.configcols(['minval','maxval'],'s','%s',visible=1)
        pages.configcols(['keys'],'o',visible=0)

        Npages = int(len(keys)/NperPage)+1
        
        for pageID in xrange(1,Npages+1):
            indexmin = (pageID-1)*NperPage
            if indexmin>=len(keys): break            
            indexmax = (pageID)*NperPage-1
            if indexmax>=len(keys): indexmax=len(keys)-1
            valmin = self.colvalue2string(keys[indexmin],sortcol)
            valmax = self.colvalue2string(keys[indexmax],sortcol)
            N=(indexmax-indexmin)+1
            pages.newrow({'keys':keys[indexmin:indexmax+1],'pageID':pageID,'N':N,'imin':indexmin,'imax':indexmax,'minval':valmin,'maxval':valmax})
        return(pages)

    def pagename(self,basename,pageID=None):
        if pageID!=None and pageID>1:
            pagename = '%s.p%d.html' % (basename,pageID)
        else:
            pagename = '%s.html' % (basename)
       
        return(pagename)

    def colname(self,col,extracolmapping=None):
        if extracolmapping!=None and (col in extracolmapping):
            return(extracolmapping[col])
        if col in self.colnamemapping:
            return(self.colnamemapping[col])
        
        return(col)

    def sub_UTrange_placeholder(self,webpage,UTrange,suffix=''):
        if UTrange == None: UTrange=('N/A','N/A','N/A','N/A')
        for i in xrange(4):
            if UTrange[i]==None: UTrange[i]='N/A'
        webpage.substituteplaceholder('PLACEHOLDER_UTFIRSTDATE%s_PLACEHOLDER' % (suffix),UTrange[0])
        webpage.substituteplaceholder('PLACEHOLDER_UTFIRSTTIME%s_PLACEHOLDER' % (suffix),UTrange[1])
        webpage.substituteplaceholder('PLACEHOLDER_UTLASTDATE%s_PLACEHOLDER' % (suffix),UTrange[2])
        webpage.substituteplaceholder('PLACEHOLDER_UTLASTTIME%s_PLACEHOLDER' % (suffix),UTrange[3])
        

    def mkpagetable(self,pages,pagekey0,sortcol,webbasename,NperTablerow=10,NperPage=100):

        infotable = htmltable(NperTablerow,border=1,cellspacing=0,cellpadding=2,width='100px',font=self.font,fontsize=self.pagetablefontsize)
        #if len(basenames)<1: return(infotable)
        infotable.startrow()
        infotable.addcol('Page links, sorted and grouped by %s, each with %d images' % (self.colname(sortcol),NperPage),colspan=NperTablerow,bold=1, color = 'white', bgcolor = 'blue',font=self.font,fontscale="+2")
        infotable.endrow()
        
        for i in xrange(len(pages.allrowkeys)):
            pagekey=pages.allrowkeys[i]
            if (i % NperTablerow)==0:
                infotable.startrow()
            
            valrange = '%s-%s' % (pages.getentry(pagekey,'minval'),pages.getentry(pagekey,'maxval'))
            if pagekey != pagekey0:
                pagestring = '<span style="color: rgb(255, 0, 0);">%d:</span>' % (pagekey+1)
                pagename = os.path.basename(self.pagename(webbasename,pageID=pages.getentry(pagekey,'pageID')))
                #s = pagestring+addlink2string(valrange,pagename)
                s = addlink2string(valrange,pagename)
                infotable.addcol(s)
            else:
                pagestring = '<span style="color: rgb(255, 255, 255);">%d:</span>' % (pagekey+1)
                #s = pagestring+valrange
                s = valrange
                infotable.addcol(s,bgcolor='red',color='white',bold=True,fontsize=self.pagetablefontsize*1.0)

            if ((i+1) % NperTablerow)==0:
                infotable.endrow()
            
        return(infotable)

    def colval2imagetable(self,table,key,col):
        table.addcol(self.colvalue2string(key,col))            
        return(0)
        
    def mkimagetable(self,keys,p='50',cols2show1=[],cols2show2=[]):
        figs2show = self.cfg.getstring("webpage","figs2show").split(',')

        # how many subdirs?
        #self.outsubdir_depth = self.cfg.getint("webpage","outsubdir_depth")
        # if the last subdir was dropped,
        #if self.cfg.getboolean("webpage","skiplastsubdir"):
        #    self.outsubdir_depth -= 1

        #  +! for filename, the other +1 for fits header
        Ncols = len(figs2show)+1+len(cols2show1)+len(cols2show2)+1
        rampflag=False
        if 'medianramp' in figs2show or 'xyramp' in figs2show:
            rampflag = True
            Ncols +=1
            
        infotable = htmltable(Ncols,border=1,cellspacing=0,cellpadding=2,width='100px',textalign='center',verticalalign='center',fontscale=self.imagetablefontscale,font=self.font,bgcolor=self.bgcolor1imtable)

        # make the table header
        infotable.startrow()

        # first set of columns
        for col in cols2show1:
            infotable.addcol(self.colname(col),textalign='center',width=None,bold=1, color = 'white', bgcolor = self.bgcolor4tableheader,font=self.font,fontscale=self.imagetableheaderfontscale)

        for fig2show in figs2show:
            infotable.addcol(self.colname(fig2show),textalign='center',width=None,bold=1, color = 'white', bgcolor = self.bgcolor4tableheader,font=self.font,fontscale=self.imagetableheaderfontscale)
        
        # ramp ASCII
        if rampflag:
            infotable.addcol('ASCII ramp',textalign='center',width=None,bold=1, color = 'white', bgcolor = self.bgcolor4tableheader,font=self.font,fontscale=self.imagetableheaderfontscale)

        # fits header
        infotable.addcol('fits hdr',textalign='center',width=None,bold=1, color = 'white', bgcolor = self.bgcolor4tableheader,font=self.font,fontscale=self.imagetableheaderfontscale)

        # second set of columns
        for col in cols2show2:
            infotable.addcol(self.colname(col),textalign='center',width=None,bold=1, color = 'white', bgcolor = self.bgcolor4tableheader,font=self.font,fontscale=self.imagetableheaderfontscale)

        # filename
        infotable.addcol('file',colspan=self.Noutsubdir()+1,textalign='center',width=None,bold=1, color = 'white', bgcolor = self.bgcolor4tableheader,font=self.font,fontscale=self.imagetableheaderfontscale)
        infotable.endrow()
         
        if len(keys)<1: return(infotable)

        removeJNpattern = re.compile('JN\d+\.')
        # fill the table
        for key in keys:

            basename = self.getentry(key,'basename')

            infotable.startrow()
            
            for col in cols2show1:
                self.colval2imagetable(infotable,key,col)
                #infotable.addcol(self.colvalue2string(key,col),textalign='center',width=None,bold=0,verticalalign='center',fontscale=self.imagetablefontscale,font=self.font,bgcolor=self.bgcolor1imtable)            

            string4fitsheaderstring = ''
            string4ASCIIrampstring = ''
            # add the figures
            for fig2show in figs2show:
                figfile = '%s%s.%s.%s' % (self.imagepathprefix,basename,fig2show,self.cfg.getstring("webpage","figsuffix"))
                figfileTN = '%s%s.TN.%s.%s' % (self.imagepathprefix,basename,fig2show,self.cfg.getstring("webpage","figsuffix"))
                if os.path.isfile(os.path.abspath('%s/%s' % (self.cfg.getstring('webpage','webrootdir'),figfile))):
                    if self.thumbnails:
                        s = addlink2string(imagestring4web(figfileTN,width=None,height=p),figfile)
                    else:
                        s = addlink2string(self.colname(fig2show),figfile)
                else:
                    s = self.NotAvailable
                #infotable.addcol(s,textalign='center',width=None,bold=0,verticalalign='center',font=self.font,bgcolor=self.bgcolor2imtable)            
                infotable.addcol(s)            

                if fig2show in ['raw','red','slp','cds','wfs','hdr']:                
                    fitsheaderfilename = '%s%s.%s.header.txt' % (self.imagepathprefix,basename,fig2show)
                    if string4fitsheaderstring!= '': string4fitsheaderstring += ' '
                    if os.path.isfile('%s/%s' % (self.cfg.getstring('webpage','webrootdir'),fitsheaderfilename)):
                        string4fitsheaderstring +=  addlink2string(fig2show,fitsheaderfilename)
                    else:
                        string4fitsheaderstring +=  fig2show
                        

                if fig2show in ['medianramp','xyramp']:
                    rampfilename = '%s%s.%s.txt' % (self.imagepathprefix,basename,fig2show)
                    s = re.sub('ramp','',fig2show)
                    if string4ASCIIrampstring!= '': string4ASCIIrampstring += ' '
                    if os.path.isfile('%s/%s' % (self.cfg.getstring('webpage','webrootdir'),rampfilename)):
                        string4ASCIIrampstring +=  addlink2string(s,rampfilename)
                    else:
                        string4ASCIIrampstring += s 
                    

            # ASCCI ramp column
            if rampflag:
                #infotable.addcol(string4ASCIIrampstring,textalign='center',width=None,bold=0,verticalalign='center',font=self.font,bgcolor=self.bgcolor1imtable)            
                infotable.addcol(string4ASCIIrampstring)            

            # fitsheader column
            #infotable.addcol(string4fitsheaderstring,textalign='center',width=None,bold=0,verticalalign='center',fontscale=self.imagetablefontscale,font=self.font,bgcolor=self.bgcolor1imtable)            
            infotable.addcol(string4fitsheaderstring)            
                
            for col in cols2show2:
                self.colval2imagetable(infotable,key,col)
                #infotable.addcol(self.colvalue2string(key,col),textalign='center',width=None,bold=0,verticalalign='center',fontscale=self.imagetablefontscale,font=self.font,bgcolor=self.bgcolor1imtable)            

            # add the subdirs ...
            #(subdirs_withanchors,linkflag) = self.imagegroup_subdiranchors(basename,forceflag=(basename==webbasenames[0]),mkanchorflag=imagegroupanchorflag)
            #for l in subdirs_withanchors:
            #    infotable.addcol(l,verticalalign='center',textalign='left',font=self.font,bgcolor=self.bgcolor2imtable)

            s = addlink2string(removeJNpattern.sub('',os.path.basename(basename)),'%s%s.txt' % (self.imagepathprefix,basename))
            #infotable.addcol(s,verticalalign='center',textalign='left',width=None,bold=0,fontscale=self.imagetablefontscale,font=self.font,bgcolor=self.bgcolor2imtable)            
            infotable.addcol(s,textalign='left')            

            infotable.endrow()
        return(infotable)

    def mkimagewebpage(self,webpage,pages,pagekey,sortcol,webbasename,sortedbytime=False,NperPage=5,NperTablerow=5,p='50',title='',UTrange0=None):

        if len(pages.allrowkeys)>1:
            infotable = self.mkpagetable(pages,pagekey,sortcol,webbasename,NperTablerow=NperTablerow,NperPage=NperPage)
            webpage.substituteplaceholder('PLACEHOLDER_PAGETABLE_PLACEHOLDER ',infotable.gettable())
            del infotable            
        else:
            webpage.substituteplaceholder('PLACEHOLDER_PAGETABLE_PLACEHOLDER ','')
            
        # UT range for all images
        self.sub_UTrange_placeholder(webpage,UTrange0,suffix='0')
        
        
        # UT range for images in this page
        UTrange1 = self.getfirstlastdate(pages.getentry(pagekey,'keys'),sortedbytime=sortedbytime)
        self.sub_UTrange_placeholder(webpage,UTrange1,suffix='1')

        infotable = self.mkimagetable(pages.getentry(pagekey,'keys'),p=p,
                                      cols2show1=self.cfg.getstring("webpage","cols2show1").split(','),
                                      cols2show2=self.cfg.getstring("webpage","cols2show2").split(','))
        webpage.substituteplaceholder('PLACEHOLDER_IMAGETABLE_PLACEHOLDER',infotable.gettable(sortable=True))
        del infotable
        
        webpage.substituteplaceholder('PLACEHOLDER_PAGENUMBER_PLACEHOLDER','%d/%d' % (pagekey+1,len(pages.allrowkeys)))

        webpage.substituteplaceholder('PLACEHOLDER_TITLE_PLACEHOLDER',title)
        
        if self.cfg.getboolean('webpage','webpage_thumbnails'):
            # if thumbnails webpage is created as well, add the links to toggle back and forth...
            webfilename = self.pagename(webbasename,pageID=pages.getentry(pagekey,'pageID'))
            if self.thumbnails:
                thumbnailstring = addlink2string('No Thumbnails','../%s' % os.path.basename(webfilename))
            else:
                thumbnailstring = addlink2string('Thumbnails','tn/%s' % os.path.basename(webfilename))
        else:
            thumbnailstring = ''            
        webpage.substituteplaceholder('PLACEHOLDER_BACKTOMAINLINK_PLACEHOLDER',thumbnailstring+'&nbsp;&nbsp;&nbsp;&nbsp;'+self.mkBackToMainlinks())
        
        webpage.substituteplaceholder('PLACEHOLDER_LASTUPDATE_PLACEHOLDER','%s' % time.asctime())

        return(webpage)


    def mkimagelistwebpage(self,keys,sortcol,defaultwebpagename,webbasename,alreadysorted=False,sortedbytime=False,NperPage=5,NperTablerow=5,p='50',title=''):

        # sort if necessary
        if not alreadysorted:
            if self.verbose>0: print 'Sorting keys'
            if self.colinfo[sortcol]['type']=='s':
                keys = self.sortkeysbycols(keys,[sortcol],asstring=1)
            else:
                keys = self.sortkeysbycols(keys,sortcol,asstring=0)

        pages = self.definepages_Nim(keys,sortcol,NperPage=NperPage)
        pages.printtxttable()

        UTrange0 = self.getfirstlastdate(keys,sortedbytime=sortedbytime)

        for pagekey in pages.allrowkeys:
            
            webpage = webpageclass()
            webpage.loaddefaultpage(defaultwebpagename)

            webpage = self.mkimagewebpage(webpage,pages,pagekey,sortcol,webbasename,sortedbytime=sortedbytime,NperPage=NperPage,NperTablerow=NperTablerow,p=p,title=title,UTrange0=UTrange0)

            webfilename = self.pagename(webbasename,pageID=pages.getentry(pagekey,'pageID'))
            print '### Saving ',webfilename
            webpage.savepage(webfilename)
            del webpage


        del pages

    def overviewtable4imagegroups(self,infotable,primarycol,cols='allvisible',colnamemapping=None,title=None):
        if cols=='allvisible':
            cols = infotable.__cols2use__()
        webpagetable = htmltable(len(cols),border=1,cellspacing=0,cellpadding=2,textalign='center',font=self.font)

        # title
        webpagetable.settabletitle(title,fontsize_pt=22)
        #if title!=None:
        #    webpagetable.startrow()
        #    webpagetable.addcol(title,colspan=len(cols),bold=1, color = 'white', bgcolor = self.bgcolor4tableheader,fontscale="+2")
        #    webpagetable.endrow()

        # column heads
        webpagetable.startrow()
        for col in cols:
            s = self.colname(col,extracolmapping=colnamemapping)
            webpagetable.addcol(s, color = 'black', bgcolor = 'lightblue')
        webpagetable.endrow()

        for key in infotable.allrowkeys:
            webpagetable.startrow()
            for col in cols:
                if col == primarycol and infotable.getentry(key,'webpage')!=None:
                    s = addlink2string(infotable.colvalue2string(key,primarycol),os.path.basename(infotable.getentry(key,'webpage')))
                else:
                    s = infotable.colvalue2string(key,col)
                if col=='task':
                    talign='left'
                else:
                    talign=None
                webpagetable.addcol(s,textalign=talign)
            webpagetable.endrow()
        
        return(webpagetable)
        
    def overviewtable4allimages(self,name2htmlhash):
        names = name2htmlhash.keys()
        N=len(names)
        if N==0:N=1
        t = htmltable(N,border=1,cellspacing=0,cellpadding=2)
        t.settabletitle('Links to All Images, sorted lists',fontsize_pt=22)
        #t.startrow()
        #t.addcol('Links to All Images, sorted lists',colspan=N,bold=1, color = 'white', bgcolor = self.bgcolor4tableheader,font=self.font,fontscale="+2")
        #t.endrow()

        t.startrow()
        counter = 1
        for name in names:
            t.addcol(addlink2string(name,name2htmlhash[name]))
            if (counter % 10) == 0:
                t.endrow()
                t.startrow()                
            counter+=1
        t.endrow()
        return(t)

    def mksearchpatternwebpages(self,keys,sortcol,searchpatterns,alreadysorted=False,sortedbytime=False,NperPage=5,NperTablerow=5,p='50',skipUT=False):

        t = txttableclass()
        t.configcols(['searchpattern'],'s','%s',visible=1)
        t.configcols(['Nim'],'d','%d',visible=1)
        t.configcols(['minval','maxval'],'s','%s',visible=1)
        if not skipUT:
            t.configcols(['minUT','maxUT'],'s','%s',visible=1)
        t.configcols(['keys'],'o',visible=0)
        t.configcols(['webpage'],'s','%s',visible=0)

        if len(searchpatterns)<1:
            return(t)

        # sort if necessary
        if not alreadysorted:
            if self.verbose>0: print 'Sorting keys for searchpatterns'
            if self.colinfo[sortcol]['type']=='s':
                keys = self.sortkeysbycols(keys,[sortcol],asstring=1)
            else:
                keys = self.sortkeysbycols(keys,sortcol,asstring=0)

        keyhash = {}
        patternobject={}

        for searchpattern in searchpatterns:
            keyhash[searchpattern]=[]
            try:
                patternobject[searchpattern]=re.compile(searchpattern,flags=re.IGNORECASE)
            except:
                # old python
                patternobject[searchpattern]=re.compile(searchpattern)

        for key in keys:
            basename = self.getentry(key,'basename')
            for searchpattern in searchpatterns:
                if patternobject[searchpattern].search(basename):
                    keyhash[searchpattern].append(key)
        
        for searchpattern in searchpatterns:
            keys4pattern = keyhash[searchpattern][:]
            print 'search for %s: %d images found' % (searchpattern,len(keys4pattern))

            searchwebpagebasename = '%s/%s.search' % (self.cfg.getstring("webpage","webrootdir"),self.searchpatternname(searchpattern))
            # make the webpage for all keys        
            self.mkimagelistwebpage(keys4pattern,
                                    sortcol,
                                    self.cfg.getstring("webpage","defaultpage_imagelist"),
                                    searchwebpagebasename,
                                    alreadysorted=alreadysorted,
                                    sortedbytime=sortedbytime,
                                    p=p,
                                    NperPage=NperPage,
                                    NperTablerow=NperTablerow,
                                    title = '%s: "%s" search' % (self.cfg.getstring("webpage","title"),searchpattern)
                                    )
            

            if len(keys4pattern)<1:
                t.newrow({'searchpattern':searchpattern,'keys':keys4pattern,'Nim':len(keys4pattern)})
            else:
                tkey = t.newrow({'searchpattern':searchpattern,'keys':keys4pattern,'Nim':len(keys4pattern),
                                 'minval':self.colvalue2string(keys4pattern[0],sortcol),
                                 'maxval':self.colvalue2string(keys4pattern[-1],sortcol),
                                 'webpage':'%s.html' % (searchwebpagebasename)
                                 })
                if not skipUT:
                    (firstdate1,firsttime1,lastdate1,lasttime1) = self.getfirstlastdate(keys4pattern,sortedbytime=sortedbytime)
                    t.add2row(tkey,{'minUT':'%sT%s' % (firstdate1,firsttime1),
                                    'maxUT':'%sT%s' % (lastdate1,lasttime1)
                                    })
        if self.verbose>1:
            cols = t.__cols2use__()
            cols.append('webpage')
            t.printtxttable(cols=cols)
        return(t)

    def mkfitskeysearchpatternwebpages(self,keys,sortcol,fitskeysearch_filename,alreadysorted=False,sortedbytime=False,NperPage=5,NperTablerow=5,p='50',skipUT=False):

        t = txttableclass()
        t.loadfile(fitskeysearch_filename)
        t.configcols(['Nim'],'d','%d',visible=1)
        t.configcols(['minval','maxval'],'s','%s',visible=1)
        if not skipUT:
            t.configcols(['minUT','maxUT'],'s','%s',visible=1)
        t.configcols(['keys'],'o',visible=0)
        t.configcols(['compiledsearchpattern'],'o',visible=0)
        t.configcols(['webpage'],'s','%s',visible=0)

        if len(t.allrowkeys)<1:
            return(t)

        # sort if necessary
        if not alreadysorted:
            if self.verbose>0: print 'Sorting keys for searchpatterns'
            if self.colinfo[sortcol]['type']=='s':
                keys = self.sortkeysbycols(keys,[sortcol],asstring=1)
            else:
                keys = self.sortkeysbycols(keys,sortcol,asstring=0)


        keyhash = {}
        patternobject={}
        for tkey in t.allrowkeys:
            keyhash[tkey]=[]
            patternobject[tkey]=re.compile(t.getentry(tkey,'searchpattern'))

        for key in keys:
            for tkey in t.allrowkeys:
                s = self.getentry(key,t.getentry(tkey,'fitskey'))
                if s!=None and patternobject[tkey].search(s):
                    keyhash[tkey].append(key)
                    
            
        
        for tkey in t.allrowkeys:
            keys4pattern = keyhash[tkey][:]
            print 'search for %s=%s: %d images found' % (t.getentry(tkey,'fitskey'),t.getentry(tkey,'searchpattern'),len(keys4pattern))

            searchwebpagebasename = '%s/%s_%s.search' % (self.cfg.getstring("webpage","webrootdir"),t.getentry(tkey,'fitskey'),self.searchpatternname(t.getentry(tkey,'searchpattern')))
            # make the webpage for all keys        
            self.mkimagelistwebpage(keys4pattern,
                                    sortcol,
                                    self.cfg.getstring("webpage","defaultpage_imagelist"),
                                    searchwebpagebasename,
                                    alreadysorted=alreadysorted,
                                    sortedbytime=sortedbytime,
                                    p=p,
                                    NperPage=NperPage,
                                    NperTablerow=NperTablerow,
                                    title = '%s: fitskey search %s=%s' % (self.cfg.getstring("webpage","title"),t.getentry(tkey,'fitskey'),t.getentry(tkey,'searchpattern'))
                                    )
            

            if len(keys4pattern)<1:
                t.add2row(tkey,{'keys':keys4pattern,'Nim':len(keys4pattern)})
            else:
                t.add2row(tkey,{'keys':keys4pattern,'Nim':len(keys4pattern),
                                'minval':self.colvalue2string(keys4pattern[0],sortcol),
                                'maxval':self.colvalue2string(keys4pattern[-1],sortcol),
                                'webpage':'%s.html' % (searchwebpagebasename)
                                })
                if not skipUT:
                    (firstdate1,firsttime1,lastdate1,lasttime1) = self.getfirstlastdate(keys4pattern,sortedbytime=sortedbytime)
                    t.add2row(tkey,{'minUT':'%sT%s' % (firstdate1,firsttime1),
                                    'maxUT':'%sT%s' % (lastdate1,lasttime1)
                                    })
        if self.verbose>1:
            cols = t.__cols2use__()
            cols.append('webpage')
            t.printtxttable(cols=cols)
        return(t)

    def addcolvals2table(self,t,tkey,keys,colvals2table):
        for addcol in colvals2table:
            colvals = unique(self.col_asstring_list(addcol,keys=keys))
            colvals.sort()
            for i in xrange(len(colvals)-1,-1,-1):
                if colvals[i]==None:
                    del(colvals[i])
            t.setentry(tkey,addcol,', '.join(colvals))

    def mkwebpages4columngroup(self,keys,sortcol,groupcol,alreadysorted=False,sortedbytime=False,NperPage=5,NperTablerow=5,p='50',skipUT=False,colvals2table=[]):

        t = txttableclass()
        t.configcols([groupcol],'s','%s',visible=1)
        t.configcols(['Nim'],'d','%d',visible=1)
        t.configcols(['minval','maxval'],'s','%s',visible=1)
        if not skipUT:
            t.configcols(['minUT','maxUT'],'s','%s',visible=1)
        if colvals2table!=None and len(colvals2table)>0:
            for addcol in colvals2table:
                t.configcols([addcol],'s','%s',visible=1)
            
        t.configcols(['keys'],'o',visible=0)
        t.configcols(['webpage'],'s','%s',visible=0)

        if len(keys)<1:
            return(t)

        # sort if necessary
        if not alreadysorted:
            if self.verbose>0: print 'Sorting keys for searchpatterns'
            if self.colinfo[sortcol]['type']=='s':
                keys = self.sortkeysbycols(keys,[sortcol],asstring=1)
            else:
                keys = self.sortkeysbycols(keys,sortcol,asstring=0)

        keyhash = {}

        groupvals = unique(self.col_as_list(groupcol,keys=keys))
        groupvals.sort()
        
        keyhash = {}
        for groupval in groupvals:
            keyhash[groupval]=[]

        for key in keys:
            keyhash[self.getentry(key,groupcol)].append(key)
        
        for groupval in groupvals:
            if groupval==None:
                continue
            keys4group = keyhash[groupval][:]
            print '%s=%s: %d images found' % (groupcol,str(groupval),len(keys4group))

            groupwebpagebasename = '%s/%s.%s' % (self.cfg.getstring("webpage","webrootdir"),self.webaddresscompliantname(groupval),self.webaddresscompliantname(groupcol))
            # make the webpage for all keys        
            self.mkimagelistwebpage(keys4group,
                                    sortcol,
                                    self.cfg.getstring("webpage","defaultpage_imagelist"),
                                    groupwebpagebasename,
                                    alreadysorted=alreadysorted,
                                    sortedbytime=sortedbytime,
                                    p=p,
                                    NperPage=NperPage,
                                    NperTablerow=NperTablerow,
                                    title = '%s=%s' % (groupcol,str(groupval))
                                    )
            

            if len(keys4group)<1:
                t.newrow({groupcol:groupval,'keys':keys4group,'Nim':len(keys4group)})
            else:
                tkey = t.newrow({groupcol:groupval,'keys':keys4group,'Nim':len(keys4group),
                                 'minval':self.colvalue2string(keys4group[0],sortcol),
                                 'maxval':self.colvalue2string(keys4group[-1],sortcol),
                                 'webpage':'%s.html' % (groupwebpagebasename)
                                 })

                # add UT range to table
                if not skipUT:
                    (firstdate1,firsttime1,lastdate1,lasttime1) = self.getfirstlastdate(keys4group,sortedbytime=sortedbytime)
                    t.add2row(tkey,{'minUT':'%sT%s' % (firstdate1,firsttime1),
                                    'maxUT':'%sT%s' % (lastdate1,lasttime1)
                                    })
                    
                # add extra info to table
                if colvals2table!=None and len(colvals2table)>0:
                    self.addcolvals2table(t,tkey,keys4group,colvals2table)
        if self.verbose>1:
            cols = t.__cols2use__()
            cols.append('webpage')
            t.printtxttable(cols=cols)
        return(t)

    def mkwebpages(self):       

        if self.cfg.getboolean("webpage","tablesortable"):
            sourcesortablescriptname = '%s/sortable.js' % os.environ['QUICKLOOK_CONFIGDIR']
            targetsortablescriptname = '%s/sortable.js' % self.cfg.getstring("webpage","webrootdir")
            if os.path.isfile(targetsortablescriptname):
                print '%s already exists, skipping'
            else:
                if not os.path.isdir(os.path.dirname(targetsortablescriptname)):
                    makepath(os.path.dirname(targetsortablescriptname))
                shutil.copyfile(sourcesortablescriptname,targetsortablescriptname)
                if not os.path.isfile(targetsortablescriptname):
                    raise RuntimeError,"%s could not be copied to %s" % (sourcesortablescriptname,targetsortablescriptname)


        self.getcolnamesmapping()

        #self.getfileinfo(basenames)
        #self.fillup_fileinfo()
        
        #all_keys_sorted_jobN = self.sortkeysbycols(self.rowkeys(),'jobN',asstring=0)
        #self.saveinfofile(keys= all_keys_sorted_jobN)

        overviewwebpage = webpageclass()
        overviewwebpage.loaddefaultpage(self.cfg.getstring("webpage","defaultpage_overview"))

        # we sort the keys only once! 
        keys_sorted_maincol = self.sortkeysbycols(self.allinfokeys,self.sortcol,asstring=1)

        # add the UT range to webpage
        UTrange = self.getfirstlastdate(keys_sorted_maincol,sortedbytime=True)
        self.sub_UTrange_placeholder(overviewwebpage,UTrange,suffix='')

        # make the webpage for all keys        
        self.mkimagelistwebpage(keys_sorted_maincol,
                                self.sortcol,
                                self.cfg.getstring("webpage","defaultpage_imagelist"),
                                '%s/allimages' % self.cfg.getstring("webpage","webrootdir"),
                                alreadysorted=True,
                                sortedbytime=True,
                                p=self.cfg.getstring("webpage","thumbnailsize"),
                                NperPage=self.cfg.getint("webpage","NperPage"),
                                NperTablerow=self.cfg.getint("webpage","NperTablerow"),
                                title = self.cfg.getstring("webpage","title")+' All Images (sorted by job #)'
                                )
        #infotable = self.overviewtable4allimages({'job #':'jobN.all.html','task':'task.all.html'})
        webpagetable = self.overviewtable4allimages({self.sortcol:'allimages.html'})
        overviewwebpage.substituteplaceholder('PLACEHOLDER_LINK2ALLIMAGES_PLACEHOLDER ',webpagetable.gettable())
        del webpagetable
        
        # make the webpages for the tasks
        groupinfotable = self.mkwebpages4columngroup(keys_sorted_maincol,
                                                     self.sortcol,
                                                     'task',
                                                     alreadysorted=True,
                                                     sortedbytime=True,
                                                     p=self.cfg.getstring("webpage","thumbnailsize"),
                                                     NperPage=self.cfg.getint("webpage","NperPage"),
                                                     NperTablerow=self.cfg.getint("webpage","NperTablerow"),
                                                     skipUT=False,
                                                     colvals2table=['DATE-OBS']
                                                     )
        webpagetable = self.overviewtable4imagegroups(groupinfotable,'task',
                                                      colnamemapping={'minval':'first job #','maxval':'last job #'},
                                                      title='Links to Tasks')
        overviewwebpage.substituteplaceholder('PLACEHOLDER_LINKS2TASKGROUPS_PLACEHOLDER',webpagetable.gettable(sortable=True))
        del groupinfotable,webpagetable

        # make the webpages for the dates
        groupinfotable = self.mkwebpages4columngroup(keys_sorted_maincol,
                                                     self.sortcol,
                                                     'DATE-OBS',
                                                     alreadysorted=True,
                                                     sortedbytime=True,
                                                     p=self.cfg.getstring("webpage","thumbnailsize"),
                                                     NperPage=self.cfg.getint("webpage","NperPage"),
                                                     NperTablerow=self.cfg.getint("webpage","NperTablerow"),
                                                     skipUT=True,
                                                     colvals2table=['task']
                                                     )
        webpagetable = self.overviewtable4imagegroups(groupinfotable,'DATE-OBS',
                                                      colnamemapping={'minval':'first job #','maxval':'last job #'},
                                                      title='Links to each Date')
        overviewwebpage.substituteplaceholder('PLACEHOLDER_LINKS2DATEGROUPS_PLACEHOLDER',webpagetable.gettable(sortable=True))
        del groupinfotable,webpagetable
                
        # make the webpages for the detectors
        groupinfotable = self.mkwebpages4columngroup(keys_sorted_maincol,
                                                     self.sortcol,
                                                     'DETECTOR',
                                                     alreadysorted=True,
                                                     sortedbytime=True,
                                                     p=self.cfg.getstring("webpage","thumbnailsize"),
                                                     NperPage=self.cfg.getint("webpage","NperPage"),
                                                     NperTablerow=self.cfg.getint("webpage","NperTablerow"),
                                                     skipUT=True
                                                     )
        webpagetable = self.overviewtable4imagegroups(groupinfotable,'DETECTOR',
                                                      colnamemapping={'minval':'first job #','maxval':'last job #'},
                                                      title='Links to each Detector')
        overviewwebpage.substituteplaceholder('PLACEHOLDER_LINKS2DETECTORGROUPS_PLACEHOLDER',webpagetable.gettable(sortable=True))
        del groupinfotable,webpagetable

        # make the webpages for the fitskey searches
        fitskeysearch_filename = self.cfg.getstring("webpage","fitskeysearch_filename")
        if fitskeysearch_filename != '':
            groupinfotable = self.mkfitskeysearchpatternwebpages(keys_sorted_maincol,
                                                                 self.sortcol,
                                                                 fitskeysearch_filename,
                                                                 alreadysorted=True,
                                                                 sortedbytime=True,
                                                                 p=self.cfg.getstring("webpage","thumbnailsize"),
                                                                 NperPage=self.cfg.getint("webpage","NperPage"),
                                                                 NperTablerow=self.cfg.getint("webpage","NperTablerow"),
                                                                 skipUT=False
                                                                 )
            if len(groupinfotable.allrowkeys)>1:
                webpagetable = self.overviewtable4imagegroups(groupinfotable,'searchpattern',
                                                              colnamemapping={'minval':'first job #','maxval':'last job #'},
                                                              title='Links to Fitskey Searches')
                overviewwebpage.substituteplaceholder('PLACEHOLDER_LINKS2FITSKEYSEARCHES_PLACEHOLDER',webpagetable.gettable(sortable=True))
                del groupinfotable,webpagetable
            else:
                overviewwebpage.substituteplaceholder('PLACEHOLDER_LINKS2FITSKEYSEARCHES_PLACEHOLDER','')
                del groupinfotable
                
        else:
            overviewwebpage.substituteplaceholder('PLACEHOLDER_LINKS2FITSKEYSEARCHES_PLACEHOLDER','')

        # make the webpages for the filename searches
        searchpatterns = self.cfg.getstring("webpage","searchpatterns4pages").split(',')        
        if searchpatterns != ['']:
            groupinfotable = self.mksearchpatternwebpages(keys_sorted_maincol,
                                                          self.sortcol,
                                                          searchpatterns,
                                                          alreadysorted=True,
                                                          sortedbytime=True,
                                                          p=self.cfg.getstring("webpage","thumbnailsize"),
                                                          NperPage=self.cfg.getint("webpage","NperPage"),
                                                          NperTablerow=self.cfg.getint("webpage","NperTablerow"),
                                                          skipUT=False
                                                          )
            webpagetable = self.overviewtable4imagegroups(groupinfotable,'searchpattern',
                                                          colnamemapping={'minval':'first job #','maxval':'last job #'},
                                                          title='Links to Filename Searches')
            overviewwebpage.substituteplaceholder('PLACEHOLDER_LINKS2SEARCHES_PLACEHOLDER',webpagetable.gettable(sortable=True))
            del groupinfotable,webpagetable
        else:
            overviewwebpage.substituteplaceholder('PLACEHOLDER_LINKS2SEARCHES_PLACEHOLDER','')
        

        overviewwebpage.substituteplaceholder('PLACEHOLDER_TITLE_PLACEHOLDER',self.cfg.getstring("webpage","title")+' Overview')
        overviewwebpage.substituteplaceholder('PLACEHOLDER_LASTUPDATE_PLACEHOLDER','%s' % time.asctime())

        infofilename = os.path.basename(self.getinfotablefilename())
        if self.thumbnails:
            overviewwebpage.substituteplaceholder('PLACEHOLDER_THUMBNAILLINK_PLACEHOLDER',addlink2string('No Thumbnails','../index.html'))
            overviewwebpage.substituteplaceholder('PLACEHOLDER_LINKASCIITABLE_PLACEHOLDER',addlink2string(infofilename,'../'+infofilename))
        else:
            overviewwebpage.substituteplaceholder('PLACEHOLDER_THUMBNAILLINK_PLACEHOLDER',addlink2string('Thumbnails','tn/index.html'))
            overviewwebpage.substituteplaceholder('PLACEHOLDER_LINKASCIITABLE_PLACEHOLDER',addlink2string(infofilename,infofilename))
            
            

        webpagename = '%s/index.html' % self.cfg.getstring("webpage","webrootdir")
        print '### Saving ',webpagename
        overviewwebpage.savepage(webpagename)
        

if __name__ == "__main__":
    qlweb = qlwebclass()

    usagestring='USAGE: qlweb.py [filepattern1 filepattern2 ...]'
    parser=qlweb.add_options(usage=usagestring)
    options, args = parser.parse_args()

    qlweb.verbose=options.verbose
    qlweb.loadconfigfile(options)

    if len(args)<1:
        pattern = '*.png'
        for i in xrange(qlweb.Noutsubdir()): pattern = '*/%s' % pattern
        args = [pattern]
        #args = options.defaultfilepattern.split(',')

    basenames = qlweb.findfigs(args)
    qlweb.basenames2fileinfo(basenames)
    qlweb.mkwebpages()
    if qlweb.cfg.getboolean('webpage','webpage_thumbnails'):
        qlweb.thumbnails = True
        qlweb.imagepathprefix = '../'
        oldwebroot = qlweb.cfg.getstring('webpage','webrootdir')
        qlweb.cfg.set('webpage','webrootdir','%s/tn' % oldwebroot)
        qlweb.mkwebpages()
