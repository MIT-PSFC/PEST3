#!/usr/bin/env python

"""
run test cases and compare results
$Id: $
"""

import os

pest3x = '../build/pest3/pest3x'

tests = {
    'case1_NEW': '-i2 -feqdsk.cdf -k"70 100 140 200" -l10', 
    'case2_NEW': '-i3 -fgeqdsk -k"70 100 140 200" -l10',
    'case3_NEW': '-i1 -finp1.cdf -k"70 100 140 200" -l10',
}

for case in tests:
    os.system(pest3x + ' ' + tests[case] + (' >& pest3_%s.txt' % case))
    os.system('cmp pest3.cdf pest3_%s.cdf' % case)


#cp ../build/pest3/pest3x pest3x
#./pest3x -i2 -feqdsk.cdf -k"70 100 140 200" -l10
#./pest3x -i2 -feqdsk.cdf -k"70 100 140 200" -l10 pest3_case1_NEWs.txt 

#./pest3x -i3 -fgeqdsk -k"70 100 140 200" -l10 
#./pest3x -i3 -fgeqdsk -k"70 100 140 200" -l10 pest3_case2_NEWs.txt 

#./pest3x -i1 -finp1.cdf -k"70 100 140 200" -l10 
#./pest3x -i1 -finp1.cdf -k"70 100 140 200" -l10 pest3_case3_NEWs.txt 
