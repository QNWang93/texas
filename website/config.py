#!/usr/bin/env python
#import preprocessing
from numpy import dtype
texas_cfg = {'n_radius': 2,
         'casjob_id': '',
         'casjob_pw': '',
         'catalogue': 'glade',
         'do_im': True,
         'point_search_rad' : 1,
         'filters':'i',
         'img_suffix':'jpg',
         'd_limit':60.,#unit: arcsec
         'norm_d_limit':30.,
         'z_limit':0.08}
lc_cfg = {'tess_date': [2458324.5,2458352.5,2458381.5,2458409.5,2458437.5,2458463.5,
        2458490.5,2458516.5,2458542.5,2458568.5,2458595.5,2458624.5,
        2458653.5,2458682.5,2458710.5,2458737.5,2458763.5,2458789.5,
        2458814.5,2458841.5,2458869.5,2458897.5,2458926.5,2458955.5,
        2458982.5,2459008.5,2459034.5, 2459060.5, 2459087.5, 2459114.5,      
        2459143.5,2459172.5,2459200.5,2459227.5,2459254.5,2459280.5,2459306.5,
        2459332.5,2459360.5,2459389.5],
        'home_dir': './plots/',
        'xlim': [-20, 5],
        'atlas_info': [#{'address':'10.162.61.62', 'username':'minix', 'password':'123456'},
            {'address':'', 'username':'', 'password':''},
            {'address':'', 'username':'', 'password':''}],
        'lookback_days':30,
        'color' : {'o':'orange', 'c': 'cyan', 'g':'blue', 'r': 'green', 'i':'red'},
        'atlas_sig_limit':3,
        'img_suffix':texas_cfg['img_suffix']
        }
web_cfg = {'header': ['Name','RA','Dec','Disc date','Type','PS prob','Redshift', 'MW E(B-V)', 'TESS_coverage', 'lc','texas img','texas_table'],
'dtype':{'Name':dtype('U30'), 'RA':float, 'Dec':float, 'Disc date':dtype('U32'), 'Type':dtype('U10'), 'PS prob':dtype('U10'), 'Redshift':dtype('U15'),  'MW E(B-V)':float},
'default':{'Name':'', 'RA':None, 'Dec':None, 'Disc date':'','Type':'','PS prob':'','Redshift':'','z_source':'','MW E(B-V)':None,  'TESS_coverage':False},
'link': 'YSE',
'img_suffix':texas_cfg['img_suffix']
}
google_API = {'key': 'AIzaSyAg21vS6FIR8QlvrolBSGj50ubiPDqVydA'}

#preprocessing_queue = [preprocessing.scale_and_center,
#                       preprocessing.dot_reduction,
#                       preprocessing.connect_lines]
#use_anonymous = True

