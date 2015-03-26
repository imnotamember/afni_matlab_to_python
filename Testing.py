__author__ = 'Joshua'

import numpy as np

someDictionary = {}
someList1 = [1,2,3]
someList2 = ['1','2','3']
someList3 = ['a','s','d','f']
someDictionary['listA'] = someList1
someDictionary['listB'] = someList2
someDictionary['listC'] = someList3
print someDictionary
########################################################################################################################
channel = np.zeros(1, dtype = [('PRN',int)])
channel['PRN'] = someList1
print channel['PRN']
'''
########################################################################################################################
SN = np.array()
SN = someDictionary
########################################################################################################################
Opt = SN
SN = ''
Opt.err = 1
Opt.zerophaseoffset = 0
if ((!Opt['Respfile']) or Opt['Respfile'] is none):
    Opt['Respfile'] = ''
    Opt['Resp_out'] = 0
    Opt['RVT_out'] = 0
if ((!Opt['Cardfile']) or Opt['Cardfile'] is none):
    Opt['Cardfile'] = ''
    Opt['Card_out'] = 0
if ((!Opt['Respfile']) or Opt['Respfile'] is none) and ((!Opt['Cardfile']) or Opt['Cardfile'] is none):
    print 'No Respfile or Cardfile\n'
########################################################################################################################
'''