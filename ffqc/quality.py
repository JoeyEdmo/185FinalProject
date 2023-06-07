import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
import matplotlib.pyplot as plt
from PIL import Image
import os

def fullReadQualityDistribution(df, args):
    
    '''
    this takes in a df of the fastq data organized and makes the main plot plotting the avg read quality for a read to read
    ===========params==========
    df: DataFrame with columns as read position and rows as quality of read at that position
    args: the arguments given to the program
    ===========result==========
    
    outputs graph to specified -d directory or working directory if none 
    '''
    directory = handleWrite(args)
    writefile = directory + 'readScoreDistribution.png'
    plt.figure(figsize=(13,7))
    plt.rcdefaults()
    plt.clf()
    valuecount = {0:0}
    for reads in df.values:
        key =int(sum(reads)/len(reads))
        if(valuecount.get(key) == None):
            valuecount[key] = 1
        else:
            valuecount[key] = valuecount[key] + 1
    valuedf = pd.DataFrame([valuecount])
    valuedf = valuedf.reindex(sorted(valuedf.columns), axis=1)
    valuedf
    plt.plot(valuedf.columns, valuedf.iloc[0])
    plt.title('Number of reads with Average Quality Score', fontsize=20)
    plt.xlabel('Average Quality Per Read', fontsize=15)
    plt.ylabel('number of Reads', fontsize=15)
    plt.savefig(writefile)


def qualityDistribution(df, args):
    
    '''
    this takes in a df of the fastq data organized and makes the quality score and number of reads
    assigned to it on a graph.
    ===========params==========
    df: DataFrame with columns as read position and rows as quality of read at that position
    args: the arguments given to the program
    ===========result==========
    
    outputs graph to specified -d directory or working directory if none 
    '''
    directory = handleWrite(args)
    writefile = directory + 'baseScoreDistribution.png'
    plt.figure(figsize=(13,7))
    plt.rcdefaults()
    plt.clf()
    valuecount = {0:0}
    for reads in df.values:
        for baseval in reads:
            if(valuecount.get(baseval) == None):
                valuecount[baseval] = 1
            else:
                valuecount[baseval] = valuecount[baseval] + 1
    valuedf = pd.DataFrame([valuecount])
    valuedf = valuedf.reindex(sorted(valuedf.columns), axis=1)
    valuedf
    plt.plot(valuedf.columns, valuedf.iloc[0])
    plt.title('Number of bases with Quality Score', fontsize=20)
    plt.xlabel('Quality Score', fontsize=15)
    plt.ylabel('Number of Bases', fontsize=15)
    plt.savefig(writefile)



def mainPlot(df, args):
    
    '''
    this takes in a df of the fastq data organized and makes the main plot plotting position vs quality score boxplot.

    ===========params==========
    df: DataFrame with columns as read position and rows as quality of read at that position
    args: the arguments given to the program
    ===========result==========
    
    outputs graph to specified -d directory or working directory if none 
    '''
    
    modval = 5 
    redstart = 0
    redend = 20
    yellowend = 32
    greenend = 41.1
    hsize = 101
    vsize = 50
    
    directory = handleWrite(args)
    writefile = directory + 'mainfig.png'
    
    
    sns.set(font_scale = 4)
    plt.figure(figsize=(hsize,vsize))
    plot = sns.boxplot(data=df, showfliers = False)
    plt.title('ffqc', fontsize=100)
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
    plt.savefig(writefile, facecolor=plot.get_facecolor(), edgecolor='none')
    im = Image.open(writefile)
    im.thumbnail((1600,1600), Image.ANTIALIAS)
    im.save(writefile)


def dfScores(IN_FILE, args):
    '''
    ===========params==========
    IN_FILE: fq file we read from
    
    ===========return==========
        
    df:  our data frame where each column is a position and the rows are just given reads
    
    ===========effect============
    outputs a txt file that contains basic stats on our data.
    
    '''
    records = list(SeqIO.parse(IN_FILE, "fastq"))
    totalseqs = len(records)
    filename = IN_FILE
    sequenceLength = len(records[0].seq)
    numGC=0
    numtotal=0
    numBadRead = 0
    qscores = {}
    for record in records:
        qualities = record.letter_annotations["phred_quality"]
        qscores = storeScores(qscores, record.seq, qualities)
        s = str(record.seq)
        numGC += s.count('G') + s.count('C') 
        numtotal += len(s)
        if(sum(qualities)/len(qualities) < 5):
            numBadRead +=1
    ratioGC = int(numGC/numtotal * 100)
    
    dir = handleWrite(args)
    writeloc = dir + 'summary.txt'
    
    f = open(writeloc, 'w')
    broken = filename.split('/')
    f.write('Filename : '+ broken[len(broken)-1] + "\n")
    f.write('totalSequences : '+  str(totalseqs) +  "\n")
    f.write('Sequences Flagged as Poor Quality : '+  str(numBadRead) +  "\n")
    f.write('SequenceLength : '+  str(sequenceLength) +  "\n")
    f.write('GC Ratio : '+  str(ratioGC) +  "\n")
    f.close()
        
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