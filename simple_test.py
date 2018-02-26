# -*- coding: utf-8 -*-
'''
Simple test function test
author Owen Chapman auberginekenobi
2/19/2018
'''

def test(testname,expected, observed):
    '''
    Simple testing function. Basically checks if expected == observed.
    '''
    try:
        assert(expected==observed)
        print('passed test {}'.format(testname))
    except:
        print('failed test {}. expected {} got {}'.format(testname,expected,observed))
    
