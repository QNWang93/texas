#!/usr/bin/env python
'''
create jpg files for images
A. Rest
'''
import argparse
import glob
import sys, os, re, types, shutil,copy,random,time, math
import astropy
import astropy.io.fits as fits
from astropy.io import ascii
import numpy as np
import pylab
from   matplotlib.ticker import FormatStrFormatter,MultipleLocator
#import lc_ex
from config import web_cfg
from tools import makepath4file,rmfile

def imagestring4web(imagename,width=None,height=None):
    #imstring = '<img src="%s"' % os.path.basename(imagename)
    imstring = '<img src="%s"' % imagename
    if height != None:
        if type(height) is int: height = str(height)
        imstring += '; height=%s' % height
    if width != None:
        if type(width) is int: width = str(width)
        imstring += '; width=%s' % width
    imstring +='>'
    return(imstring)

def addlink2string(s,link,target=None):
    line = '<a '
    if target != None:
        line += 'target="%s"' % target
    line += 'href="%s">%s</a>' % (link,s)
    return(line)

def addtag2string(s,tag,target=None):
    line = '<a name="%s"></a>%s' % (tag,s)
    return(line)

def save_digit(number, num_digit):
    d = np.power(10, num_digit)
    return float(int(number*d))/d

def tab2htmltab(table, header, othercontent):
    line = '<td><table cols = %s BORDER = 1>\n <tr>' %(len(header)+1)
    for i in header:
        line += '<td>'+i + '</td> '
    line +='</tr>\n'
    for j in np.arange(len(table)):
        line += '<tr>'
        for i in header:
            if isinstance(table[j][i], float):
                if math.isnan(table[j][i]):
                    line += '<td>NaN</td> '
                else:
                    line += '<td>'+str(save_digit(table[j][i], num_digit=5)) + '</td> ' 
            else:
                line += '<td>'+str(table[j][i]) + '</td> '
         
    line += othercontent+'</tr>\n'
    line += ' </table> </td>\n'
    return line

class htmltable:
    def __init__(self,Ncols,font=None,fontscale=None,fontsize=None,color=None,bgcolor=None,cellpadding=2,cellspacing=2,border=1,
                 width='100%',height=None,textalign='center',verticalalign='top',optionalarguments=''):
        self.Ncols = Ncols
        self.font  = font
        self.fontscale  = fontscale
        self.fontsize  = fontsize
        self.color      = color
        self.bgcolor    = bgcolor
        self.cellpadding = cellpadding
        self.cellspacing = cellspacing
        self.border      = border
        self.width       = width
        self.height      = height
        self.textalign   = textalign
        self.verticalalign=verticalalign
        self.optionalarguments  = optionalarguments
        self.tabletitle = None
        self.body = []

    def startrow(self,style = ''):
        self.body.append('<tr %s>' % style)
    def endrow(self):
        self.body.append('</tr>\n')

    def addcol(self,colval,link=None, verticalalign=None, textalign=None,
               colspan=None,rowspan=None,
               bold=None, italic=None, underline = None, 
               width=None, height=None, 
               color = None, bgcolor = None, font=None, fontscale=None, fontsize=None, typ = 'td'):
        if colval is None:
            colval = '-'  # placeholder!
            
        if link != None:
            colval = '<a href="%s">%s</a>' % (link,colval)

        pre   = ''
        after = ''
        #if textalign != None:
        #    pre   += '<%s>'  % (textalign)
        #    after  = '</%s>' % (textalign) + after
        if font != None:
            pre   += '<span style="font-family: %s;">' % (font)
            after  = '</span>' + after
        if fontsize != None:
            if type(fontsize) is str: fontsize = int(fontsize)
            pre   += '<font size=%d>' % (fontsize)
            after  = '</font>' + after
        if fontscale != None:
            pre   += '<font size="%s">' % (fontscale)
            after  = '</font>' + after
        if bold != None and bold != 0:
            pre   += '<b>'
            after  = '</b>'
        if  underline != None and underline != 0:
            pre   += '<u>'
            after  = '</u>'
        if italic != None and italic != 0:
            pre   += '<b>'
            after  = '</b>'
        if color != None:
            if type(color) is int: color = str(color)
            pre   += '<font color=%s>' % (color)
            after  = '</font>' + after
            
        line = '<'+typ
        if textalign != None:
            line += ' ALIGN="%s"' % textalign            
        if width != None:
            if type(width) is int: width = str(width)
            line += ' WIDTH="%s"' % width            
        if height != None:
            if type(height) is int: height = str(height)
            line += ' HEIGHT="%s"' % height            
        if verticalalign != None:
            line += ' VALIGN="%s"' % verticalalign            
        if bgcolor != None:
            if type(bgcolor) is int: bgcolor = str(bgcolor)
            line += ' BGCOLOR="%s"' % (bgcolor)
        if colspan != None:
            line += ' colspan="%d"' % (colspan)
        if rowspan != None:
            line += ' rowspan="%d"' % (rowspan)
        #if line != '<td':
        #    line += ' NOSAVE'
        line += '>'
        line += pre + colval + after + '</'+typ+'>'
            
        self.body.append(line)

    def addtable(self, table, header, othercontent):
        self.body.append(tab2htmltab(table, header,othercontent))

    def add_sorttablescript_before_header(self):
        return('<script type="text/javascript" src="sortable.js"></script>')

    def settabletitle(self,tabletitle,align='center',fontsize_pt=None,color='white',bgcolor='blue'):
        if tabletitle==None:
            self.tabletitle = None
        
        s = '<div '
        if align!=None: s+='align="%s"; ' % align
        s+= 'style="'
        if fontsize_pt!=None: s+='font-size: %dpt;' % fontsize_pt
        if color!=None: s+='color:%s;' % color
        if bgcolor!=None: s+='background-color:%s;' % bgcolor
        s+= '">'

        self.tabletitle = [s]
        self.tabletitle.append(tabletitle)
        self.tabletitle.append('</div>')
        return(0)

    def gettable(self, sortable=False):

        tableinitstring = '<table id="t01"'
        # sortable columns?
        if sortable:
            tableinitstring += 'class="sortable" id="anyid" '
        # note: you also have to put the call 

        tableinitstring += 'style="'
        if self.textalign!=None:tableinitstring+='text-align: %s;' % (self.textalign)
        if self.width!=None:tableinitstring+='width: %s;' % (self.width)
        if self.font!=None:tableinitstring+='font-family: %s;' % (self.font)
        if self.fontscale!=None:tableinitstring+='font size: %s;>' % (self.fontscale)
        if self.fontsize!=None:
            tableinitstring+='font size: %d;>' % (int(self.fontsize))
        tableinitstring += '"'
        if self.color!=None:tableinitstring+='color="%s" ' % (self.color)
        if self.bgcolor!=None:tableinitstring+='bgcolor="%s" ' % (self.bgcolor)
        
        tableinitstring += ' COLS=%d BORDER=%d CELLSPACING=%d CELLPADDING=%d %s ' % (self.Ncols,self.border,self.cellspacing,self.cellpadding,self.optionalarguments)
        tableinitstring += '>'


        t=[]
        # Is there a title?
        if self.tabletitle != None: t.extend(self.tabletitle)

        # initialize the table
        t.append(tableinitstring)

        if sortable:
            t.append('<script type="text/javascript" src="sortable.js"></script>')

        # add the body
        t.extend(self.body)

        # close the table
        t.append('</table>')
        t.append('')

        return(t)

        
class webpageclass:
    def __init__(self):
        self.lines=[]
    def substituteplaceholder(self, pattern2find, newlines,count=0):
        import types
        patternobject = re.compile(pattern2find)

        if type(newlines) is str:
            s = newlines
        elif type(newlines) is list:
            s = '\n'.join(newlines)
        else:
            raise RuntimeError('Error: unknown type, dont know how to deal with ')
        for i in range(len(self.lines)):
            self.lines[i] = patternobject.sub(s,self.lines[i])
        
    def loaddefaultpage(self,filename):
        if not os.path.isfile(filename):
            raise RuntimeError('ERROR: could not find file '+filename)
        self.lines = open(filename).readlines()

    def savepage(self,filename):
        # Makesure the directory exists
        dir = os.path.dirname(filename)
        if not os.path.exists(dir):  os.makedirs(dir)  # This is a recursive mkdir
        rmfile(filename)
        f = open(filename,'w')
        f.writelines(self.lines)
        f.close()



class weblesniffclass:
    def __init__(self):
        self.verbose=0
        self.debug=0
        self.webdir = None
        self.figsuffix = 'jpg'
        self.usebinnedflag = False

        self.webdir = None
        self.imagelist_htmltemplate = None
        self.imagelist_htmltemplate = None

        self.imagetableheaderfontscale = "+1"
        self.target_list = astropy.table.Table()
        self.imagetablefontscale = None


        self.font = 'sans-serif'
        self.bgcolor1imtable = '#%02x%02x%02x' % (255,255,220)
        self.bgcolor4tableheader = '#E0FFFF'
        self.NotAvailable = '-'

    def define_options(self, parser=None, usage=None, conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
        parser.add_argument("webdir", help=("web directory"))
        parser.add_argument("date", help=("date of observation"))
        parser.add_argument("imagelist_htmltemplate", help=("template for the html page"))
        parser.add_argument('--figsuffix', default='jpg',
                            help=('Define the type of figure you want to save (default=%(default)s)'))
        parser.add_argument('--rootwebaddress', default=None,
                            help=('Define the root of the webaddress (default=%(default)s)'))
        parser.add_argument('--usebinned', help="Use the gis from the binned diffims",
                            action="store_true", default=False)
        parser.add_argument('--verbose', '-v', action='count')
        parser.add_argument('-d', '--debug', help="debug", action='count')

        return(parser)

    def getfiglist(self):

        figpattern = '*{name}_texas.{suffix}'.format(name=self.date, suffix=self.figsuffix)
#        figlc = '20[0-9][0-9]*[0-9][0-9][0-9][0-9]_lc.%s' % self.figsuffix
        
        self.figlist = {}
        self.figlist['target']=[]
        self.figlist['images']=[]
        self.figlist['date']=[]
        imgdir = './plots/*/'
        if len(self.target_list) ==0:
            flist = glob.glob('%s/%s' % (imgdir,figpattern))
            print(flist)
            if len(flist)>0:
                i=0
                for fname in flist:
                    fnameshort = os.path.relpath(fname,start=self.webdir)

                    objname = re.search('\d+\D+\d+_', fnameshort).group()[:-5]
                    objdate = re.search('\d+\D+\d+_', fnameshort).group()[-5:-1]
#                    m = re.search('\w+\d+\.',fnameshort)#'\w+\d+\_+\d+\.',fnameshort)

                    print(objname, objdate)
                    imname = objname+'_texas.'+self.figsuffix
                    self.figlist['target'].append(objname)
                    self.figlist['images'].append([imname, objname+objdate+'_lc.'+self.figsuffix])
                    self.figlist['date'].append(objdate)
                    i=i+1
        else: 
            for target in self.target_list:
                self.figlist['target'].append(target['Name'])
                self.figlist['date'].append(self.date)
                imname = target['Name'] + '_texas.'+self.figsuffix
                self.figlist['images'].append([imname, target['Name']+self.date+'_lc.'+self.figsuffix])

    def getwebaddress(self,field,tag=None):
        if self.rootwebaddress == None: return None
        webaddress = '%s/%s/index.html#%s' % (self.rootwebaddress,field,tag)
        return(webaddress)
    


    def makewebpage(self, cans,p='1000'):
        colors4tmpl = [self.bgcolor1imtable,'lightcyan']
        self.getfiglist()

        self.webfilename = '%s/%s.html' % (self.webdir, self.date)
        webpage = webpageclass()
        webpage.loaddefaultpage(self.imagelist_htmltemplate)

        webpage.substituteplaceholder('PLACEHOLDER_TITLE_PLACEHOLDER', os.path.basename('Candidates'))
        ggl_link = 'https://docs.google.com/spreadsheets/d/1x6DS_CmJWfhnbpwEn25DAgH194bzM7GtHkxrBhD05B4/edit?usp=sharing'
        webpage.substituteplaceholder('PLACEHOLDER_GOOGLESHEET_PLACEHOLDER', addlink2string('List',ggl_link))
        webpage.substituteplaceholder('PLACEHOLDER_BACKTOMAINLINK_PLACEHOLDER', addlink2string('BACK','..'))
        
        infotable = htmltable(3,border=1,cellspacing=0,cellpadding=2,width='1000px',textalign='center',verticalalign='center',fontscale=self.imagetablefontscale,font=self.font,bgcolor=self.bgcolor1imtable)
        field = os.path.basename(self.webdir)
        imcounter = 0

        img_dir = './plots/'
        header = web_cfg['header']
        infotable.startrow()
        for h in header:
            infotable.addcol(h, typ = 'th')
        infotable.endrow()
        if len(self.figlist['images'])>0:
            for target in self.figlist['target']:
                infotable.startrow()
                tag = target

                webaddress = self.getwebaddress(target,tag=tag)
                if webaddress != None: s+='<br><font size=5>'+webaddress+'</font>'
                img = self.figlist['images'][imcounter]
                t_info = cans[cans['Name']==target]

                for h in header[:-3]:
                    try: 
                        if h =='Name':
                            if web_cfg['link'] == 'YSE':
                                s = addlink2string(target,'https://ziggy.ucolick.org/yse/transient_detail/'+ target)
                                s = addtag2string(s,tag)
                            elif web_cfg['link'] == 'TNS':
                                s = addlink2string(target, 'https://wis-tns.weizmann.ac.il/object/'+target)
                                s = addtag2string(s,tag)
                            else: 
                                continue
                            infotable.addcol(s)
                        else:    
                            infotable.addcol(str(t_info[h][0]))
                    except:
                        infotable.addcol('')

                tmplcounter=0

                s = addlink2string(imagestring4web(img_dir+target+ '/'+img[1],width=200,height=None),img_dir+target+ '/'+ img[1])
                s = addtag2string(s,tag)
                infotable.addcol(s)  
                
                texas_info_file = './web/plots/' + target +'/'+img[0][0:-4] + '.txt'

                try:
                    texas_info_table = ascii.read(texas_info_file)
                except: 
                    texas_info_table = ascii.read('./texas_table_sample.txt')
                    
                s = addlink2string(imagestring4web(img_dir+target+ '/'+img[0],width=300,height=None),img_dir+target+ '/'+ img[0])
                s = addtag2string(s,tag)
                infotable.addcol(s)  
                texas_header = ['ra', 'dec', 'z', 'z_flag', 'norm_d', 'd']#, 'source']
                infotable.addtable(texas_info_table, texas_header, '')
                
                infotable.endrow()
                tmplcounter+=1
                imcounter+=1
                    
        webpage.substituteplaceholder('PLACEHOLDER_IMAGETABLE_PLACEHOLDER', infotable.gettable(sortable = True))
        webpage.substituteplaceholder('PLACEHOLDER_LASTUPDATE_PLACEHOLDER', '%s' % time.asctime())
       
        print('### Saving ',self.webfilename)
        webpage.savepage('./web/candidates.html') #webpage.savepage('./%s/%s' % ('web',self.webfilename))

        del webpage


def main(args):    
    weblesniff = weblesniffclass()

    weblesniff.webdir = args[0]#.webdir
    weblesniff.date = str(args[1])#.date)
    weblesniff.imagelist_htmltemplate = args[2]#.imagelist_htmltemplate

    if args[3] != None:
        weblesniff.figsuffix = args[3]#.figsuffix
    if len(args)>4:
        weblesniff.target_list = args[4]#.figsuffix

    weblesniff.rootwebaddress = None#args.rootwebaddress

    weblesniff.makewebpage(weblesniff.target_list)
    
    print("weblesniff SUCCESS")
