__author__ = 'Joshua Zosky'

import numpy as np

someDictionary = {}
someRespfile = [1,2,3]
someList2 = ['1','2','3']
someList3 = ['a','s','d']#,'f']
someDictionary['listA'] = someRespfile
someDictionary['listB'] = someList2
someDictionary['listC'] = someList3
print someDictionary
########################################################################################################################
channel = np.zeros(3, dtype = [('list1',int), ('list2',str), ('list3',str)])
#print channel['Respfile']
#channel['Respfile'] = someRespfile
#print channel['Respfile']
channel['list2'] = someList2
print channel['list2']
channel['list3'] = someList3
print channel['list3']
print type(channel)
########################################################################################################################
SN = channel
########################################################################################################################
Opt = SN
SN = ''
#Opt.err = 1
#Opt.zerophaseoffset = 0
if Opt.any('Respfile', 0):
    Opt['Respfile'] = ''
    print 'in if'
print 'out if'
'''
    Opt['Resp_out'] = 0
    Opt['RVT_out'] = 0
if ((!Opt['Cardfile']) or Opt['Cardfile'] is none):
    Opt['Cardfile'] = ''
    Opt['Card_out'] = 0
if ((!Opt['Respfile']) or Opt['Respfile'] is none) and ((!Opt['Cardfile']) or Opt['Cardfile'] is none):
    print 'No Respfile or Cardfile\n'
########################################################################################################################
'''
