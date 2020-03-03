# texas

For daily update, run:

$ ./main.py

For individuals:

to make record of a certain AT, run:

$ python3 lc_ex.py *ra* *dec* *transient_name* *date(MMDD)* 

to make a website of all AT on some certain date, run:

$ python3 weblesniff.py ./ *date(MMDD)* imagelist_template.html --figsuffix=png

a website named index.html will occur undet /plots/MMDD/


format of config.py:



#!/usr/bin/env python

#import preprocessing

texas_cfg = {'n_radius': 2,
         'casjob_id': '',
         'casjob_pw': ,
         'catalogue': 'glade',
         'do_im': True,
         'point_search_rad' : 1,
         'filters':'i'}
lc_cfg = {'tess_date': [2458324.5,2458352.5,2458381.5,2458409.5,2458437.5,2458463.5,
        2458490.5,2458516.5,2458542.5,2458568.5,2458595.5,2458624.5,
        2458653.5,2458682.5,2458710.5,2458737.5,2458763.5,2458789.5,
        2458814.5,2458841.5,2458869.5,2458897.5,2458926.5,2458955.5,
        2458982.5,2459008.5,2459034.5, 2459060.5, 2459087.5, 2459114.5,      
        2459143.5,2459172.5,2459200.5,2459227.5,2459254.5,2459280.5,2459306.5,
        2459332.5,2459360.5,2459389.5],
        
        'home_dir': './plots/', #where you put targets and plots
        
        'xlim': [-40, 10], #x-axis limit for plots
        
        'atlas_info': [{'address':'', 'username':'', 'password':''},   
            {'address':'', 'username':'', 'password':''},
            {'address':'', 'username':'', 'password':''}],#first one is optional
            
        'lookback_days':40
        }



