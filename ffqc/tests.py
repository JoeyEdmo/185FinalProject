import quality
import pandas as pd
import numpy as np
import os


def testMain():
    '''
    test testing that our df is created correctly.
    '''
    r0 = [34	,34,	34,	37,	37]
    r1 = [26	,26	,26	,26	,26]
    r2 = [32	,33,	34,	35,	36]
    data = np.array([r0,r1,r2])
    df = pd.DataFrame(data, columns=[0,1,2,3,4])
    assert(quality.dfScores('test_files/smallTest.fq').equals(df))
    print('test 1 passes!')
    
    
testMain()