import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
import matplotlib.pyplot as plt
import os







def mainPlot(df, args):
    
    '''
    this takes in a df of the fastq data organized and makes the main plot plotting position vs quality score boxplot.

    ===========params==========
    df: DataFrame with columns as read position and rows as quality of read at that position
    
    ===========result==========
    
    outputs graph to specified -d directory or working directory if none 
    '''
    
    modval = 5 #TODO: make variable based on num bases
    redstart = 0
    redend = 20
    yellowend = 32
    greenend = 41.1
    hsize = 101
    vsize = 50
    
    directory = handleWrite(args)
    
    
    sns.set(font_scale = 4)
    plt.figure(figsize=(hsize,vsize))
    plot = sns.boxplot(data=df, showfliers = False)
    plt.title('FlimsyQC', fontsize=100)
    plt.xlabel('Base Position', fontsize=80)
    plt.ylabel('Quality Score', fontsize=80)
    plot.set_ylim(ymin=0)
    for ind, label in enumerate(plot.get_xticklabels()):
        if ind % modval == 0:  # every 10th label is kept
            label.set_visible(True)
        else:
            label.set_visible(False)
            
    plt.axhspan(redstart,redend, facecolor='red', alpha=.1, zorder = -1)
    plt.axhspan(redend,yellowend, facecolor='yellow', alpha=.1, zorder = -1)
    plt.axhspan(yellowend,greenend, facecolor='green', alpha=.1, zorder = -1)
    plt.savefig(directory + 'mainfig.png', facecolor=plot.get_facecolor(), edgecolor='none')



def dfScores(IN_FILE):
    '''
    ===========params==========
    IN_FILE: fq file we read from
    
    ===========return==========
        
    df:  our data frame where each column is a position and the rows are just given reads
    '''
    qscores = {}
    for record in SeqIO.parse(IN_FILE, "fastq"):
        qualities = record.letter_annotations["phred_quality"]
        qscores = storeScores(qscores, record.seq, qualities)

    df = pd.DataFrame(qscores)
    return df







#=========================== Helpers for this file only ==============================

def handleWrite(args):
    
    '''
    makes sure your specified directory exists and if so makes our output go there
    '''
    
    dir = args.directory
    specified = False
    if(not dir):
        dir = ''
    elif(not os.path.exists(dir)):
        print('output directory does not exist. check and see if ' + dir + ' is the correct directory.')
        exit()
    else:
        dir += '/'
    return dir
        


def storeScores(dictionary, sequence, qualities):
    '''
    helps create the df storing our scores/positions by building a dictionary
    '''
    
    for n,c in enumerate(sequence):
        try:
            dictionary[n].append(qualities[n])
        except:
            dictionary[n] = []
            dictionary[n].append(qualities[n])
    return dictionary