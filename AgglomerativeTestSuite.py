# jtsorren

from GlycanParser import *
from agglomerative_assemble import *

def Test(name, expecteds, observeds):
    '''
    Test suite for the Agglomerative algorithm. 
    This function loops through possible accepted 
    solutions and if an observed solution is viable - pass.
    '''
    failed = 0
    true = 0
    if len(expecteds) == 0:
        if expecteds == observeds:
            print('Passed test {}, yay!'.format(name))
            print()
        else:
            print('Failed test {}. expected {}'.format(name, expecteds))
            print()
            print('got {}, boo!'.format(observed))
            print()
    elif len(expecteds) > 1:
        for expected in expecteds:
            for expect in expected:
                for observed in observeds:
                    try:
                        assert(Glycan.equals(expect, observed))
                        true += 1
                    except:
                        pass
        if true == sum(len(i) for i in expecteds):
            print('Passed test {}, yay!'.format(name))
            print()
        else:
            print('Failed test {}. expected {}'.format(name, expecteds))
            print()
            print('got {}, boo!'.format(observed))
            print()
    else:
        for expected in expecteds:
            for observed in observeds:
                try:
                    assert(Glycan.equals(expected, observed))
                except:
                    failed += 1
        if failed == len(expecteds):
            print('Failed test {}. expected {}'.format(name, expecteds))
            print()
            print('got {}, boo!'.format(observed))
            print()
        else:
            print('Passed test {}, yay!'.format(name))
            print()
        
def TestOne():
    '''
    Empty Test
    '''
    fragments = []
    expecteds = []
    Test(1, expecteds, agglomerative_assemble(fragments))
    return

def TestTwo():
    '''
    Non-Overlapping Fragments
    '''
    fragments = [Glycan(['L-Fucp', 'D-Galp'], [{1: ('alpha', 1, 2)}, {0: ('alpha', 2, 1)}]), 
                 Glycan(['L-Fucp', 'D-GlcpNAc'], [{1: ('alpha', 1, 3)}, {0: ('alpha', 3, 1)}])]
    expecteds = [[Glycan(['L-Fucp', 'D-Galp'], [{1: ('alpha', 1, 2)}, {0: ('alpha', 2, 1)}]), 
                  Glycan(['L-Fucp', 'D-GlcpNAc'], [{1: ('alpha', 1, 3)}, {0: ('alpha', 3, 1)}])],
                 [Glycan(['L-Fucp', 'D-GlcpNAc'], [{1: ('alpha', 1, 3)}, {0: ('alpha', 3, 1)}]), 
                  Glycan(['L-Fucp', 'D-Galp'], [{1: ('alpha', 1, 2)}, {0: ('alpha', 2, 1)}])]]
    Test(2, expecteds, agglomerative_assemble(fragments))
    return

def TestThree():
    '''
    Identical Fragments
    '''
    fragments = [Glycan(['L-Fucp', 'D-Galp'], [{1: ('alpha', 1, 2)}, {0: ('alpha', 2, 1)}]), 
                 Glycan(['L-Fucp', 'D-Galp'], [{1: ('alpha', 1, 2)}, {0: ('alpha', 2, 1)}])]
    expecteds = [Glycan(['L-Fucp', 'D-Galp'], [{1: ('alpha', 1, 2)}, {0: ('alpha', 2, 1)}])]
    Test(3, expecteds, agglomerative_assemble(fragments))
    return

def TestFour():
    '''
    Linear Glycan
    '''
    fragments = [Glycan(['L-Fucp', 'D-GlcpNAc'], [{1: ('alpha', 1, 4)}, {0: ('alpha', 4, 1)}]), 
                 Glycan(['D-GlcpNAc', 'D-GlcpNAc'], [{1: ('beta', 1, 4)}, {0: ('beta', 4, 1)}]), 
                 Glycan(['D-GlcpNAc', 'D-GlcpNAc', 'D-GlcpNAc'], [{1: ('beta', 1, 4)}, {0: ('beta', 4, 1), 2: ('beta', 1, 4)}, {1: ('beta', 4, 1)}]), 
                 Glycan(['D-Manp', 'D-GlcpNAc'], [{1: ('alpha', 4, 1)}, {0: ('alpha', 1, 4)}])]
    expecteds = [Glycan(['L-Fucp', 'D-GlcpNAc', 'D-GlcpNAc', 'D-GlcpNAc', 'D-Manp'], [{1: ('alpha', 1, 4)}, {0: ('alpha', 4, 1), 2: ('beta', 1, 4)}, {1: ('beta', 4, 1), 3: ('beta', 1, 4)}, {2: ('beta', 4, 1), 4: ('alpha', 1, 4)}, {3: ('alpha', 4, 1)}])]
    Test(4, expecteds, agglomerative_assemble(fragments))
    return

def TestFive():
    '''
    Branched Glycan
    '''
    fragments = [Glycan(['L-Fucp', 'D-GlcpNAc', 'L-Fucp'], [{1: ('alpha', 1, 3)}, {0: ('alpha', 3, 1), 2: ('alpha', 1, 3)}, {1: ('alpha', 3, 1)}]), 
                 Glycan(['L-Fucp', 'D-GlcpNAc', 'D-Manp'], [{1: ('alpha', 3, 1)}, {0: ('alpha', 1, 3), 2: ('beta', 4, 1)}, {1: ('beta', 1, 4)}]), 
                 Glycan(['L-Fucp', 'D-Manp'], [{1: ('beta', 1, 1)}, {0: ('beta', 1, 1)}])]
    expecteds = [Glycan(['L-Fucp', 'D-GlcpNAc', 'L-Fucp', 'D-Manp', 'D-Manp'], [{1: ('alpha', 1, 3)}, {0: ('alpha', 3, 1), 2: ('alpha', 1, 3), 4: ('beta', 4, 1)}, {1: ('alpha', 3, 1), 3: ('beta', 1, 1)}, {2: ('beta', 1, 1)}, {1: ('beta', 1, 4)}])]
    Test(5, expecteds, agglomerative_assemble(fragments))
    return

TestOne()
TestTwo()
TestThree()
TestFour()
TestFive()