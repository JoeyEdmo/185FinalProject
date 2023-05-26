import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
import matplotlib.pyplot as plt







def mainPlot(file):
    IN_FILE = file
    qscores = {}
    for record in SeqIO.parse(IN_FILE, "fastq"):
        qualities = record.letter_annotations["phred_quality"]
        qscores = storeScores(qscores, record.seq, qualities)

    df = pd.DataFrame(qscores)
    sns.set(font_scale = 4)
    plt.figure(figsize=(101,50))
    plot = sns.boxplot(data=df, showfliers = False)
    plt.title('FlimsyQC', fontsize=100)
    plt.xlabel('Base Position', fontsize=80)
    plt.ylabel('Quality Score', fontsize=80)
    plot.set_ylim(ymin=0)
    for ind, label in enumerate(plot.get_xticklabels()):
        if ind % 5 == 0:  # every 10th label is kept
            label.set_visible(True)
        else:
            label.set_visible(False)
            
    plt.axhspan(20,20.1, facecolor='red', alpha=.1)
    plt.axhspan(28,28.1, facecolor='yellow', alpha=.1)
    plt.axhspan(41,41.1, facecolor='green', alpha=.1)
    plt.show(block=True)











#=========================== Helpers for this file only ==============================

def storeScores(dictionary, sequence, qualities):
    for n,c in enumerate(sequence):
        try:
            dictionary[n].append(qualities[n])
        except:
            dictionary[n] = []
            dictionary[n].append(qualities[n])
    return dictionary