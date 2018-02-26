# -*- coding: utf-8 -*-
'''
Test suite for
function overlap_score
author Owen Chapman auberginekenobi
2/19/2018
'''

from GlycanParser import *
from simple_test import test

def test1():
    '''
    Two small glycans
    '''
    #https://glytoucan.org/Structures/Glycans/G00055MO
    #https://glytoucan.org/Structures/Glycans/G00065MO
    lactosamine="beta-D-Galp-(1->4)-beta-D-GlcpNAc(1->"
    a=Glycan(lactosamine)
    VIM="alpha-D-NeupAc-(2->3)-beta-D-Galp-(1->4)-beta-D-GlcpNAc(1->"
    b=Glycan(VIM)
    s,o=overlap_score(a,b)
    test(1,2,s)

def test2():
    '''
    Small glycans with extra crap
    '''
    lactosamine="beta-D-Galp-(1->4)-beta-D-GlcpNAc(1->"
    VIM="alpha-D-NeupAc-(2->3)-beta-D-Galp-(1->4)-beta-D-GlcpNAc(1->"
    test21="beta-bargledook-(2->5)-"+lactosamine
    test22="beta-bargledook-(2->5)-"+VIM
    a=Glycan(test21)
    b=Glycan(test22)
    s,o=overlap_score(a,b)
    test(2,2,s)

def test3():
    '''
    Identity test, size 3
    '''
    VIM="alpha-D-NeupAc-(2->3)-beta-D-Galp-(1->4)-beta-D-GlcpNAc(1->"
    b=Glycan(VIM)
    s,o=overlap_score(b,b)
    test(3,3,s)
    
def test4():
    '''
    Small glycans with extra crap
    '''
    lactosamine="beta-D-Galp-(1->4)-beta-D-GlcpNAc(1->"
    VIM="alpha-D-NeupAc-(2->3)-beta-D-Galp-(1->4)-beta-D-GlcpNAc(1->"
    test21="alpha-D-NeupAc-(2->3)-"+lactosamine
    test22="beta-bargledook-(2->5)-"+VIM
    a=Glycan(test21)
    b=Glycan(test22)
    s,o=overlap_score(a,b)
    test(4,3,s)

def test5():
    '''
    Null test
    '''
    a=Glycan()
    s,o=overlap_score(a,a)
    test(5,0,s)

def test6():
    '''
    Null test
    '''
    a=Glycan()
    VIM="alpha-D-NeupAc-(2->3)-beta-D-Galp-(1->4)-beta-D-GlcpNAc(1->"
    b=Glycan(VIM)
    s,o=overlap_score(a,b)
    test(6,0,s)
    s,o=overlap_score(b,a)
    test(6.1,0,s)

def test7():
    '''
    Singleton test
    '''
    a=Glycan('beta-D-GlcpNAc(1->')
    VIM="alpha-D-NeupAc-(2->3)-beta-D-Galp-(1->4)-beta-D-GlcpNAc(1->"
    b=Glycan(VIM)
    s,o=overlap_score(a,b)
    test(7,1,s)
    s,o=overlap_score(b,a)
    test(7.1,1,s)

def test8():
    '''
    Long test glycan
    '''
    # https://glytoucan.org/Structures/Glycans/G27293OK
    iupac = '''alpha-D-Manp-(1->3)beta-D-GlcpNAc-(1->4)[beta-D-Galp-(1->4)[alpha-L-Fucp-(1->3)]-beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->6)]-beta-D-Manp-(1->4)-beta-D-GlcpNAc-(1->4)[alpha-L-Fucp-(1->6)]-beta-D-GlcpNAc(1->'''
    a=Glycan(iupac)
    s,o=overlap_score(a,a)
    test(8,10,s)

def test9():
    '''
    Long against short subtree
    '''
    iupac = '''alpha-D-Manp-(1->3)beta-D-GlcpNAc-(1->4)[beta-D-Galp-(1->4)[alpha-L-Fucp-(1->3)]-beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->6)]-beta-D-Manp-(1->4)-beta-D-GlcpNAc-(1->4)[alpha-L-Fucp-(1->6)]-beta-D-GlcpNAc(1->'''
    a=Glycan(iupac)
    lactosamine="beta-D-Galp-(1->4)-beta-D-GlcpNAc(1->"
    b=Glycan(lactosamine)
    s,o=overlap_score(a,b)
    test(9,2,s)
    s,o=overlap_score(b,a)
    test(9.1,2,s)
    
def test10():
    '''
    Long against short not subtree
    '''
    iupac = '''alpha-D-Manp-(1->3)beta-D-GlcpNAc-(1->4)[beta-D-Galp-(1->4)[alpha-L-Fucp-(1->3)]-beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->6)]-beta-D-Manp-(1->4)-beta-D-GlcpNAc-(1->4)[alpha-L-Fucp-(1->6)]-beta-D-GlcpNAc(1->'''
    a=Glycan(iupac)
    VIM="alpha-D-NeupAc-(2->3)-beta-D-Galp-(1->4)-beta-D-GlcpNAc(1->"
    b=Glycan(VIM)
    s,o=overlap_score(a,b)
    test(10,2,s)
    s,o=overlap_score(b,a)
    test(10.1,2,s)
    
test1()
test2()
test3()
test4()
test5()
test6()
test7()
test8()
test9()
test10()
