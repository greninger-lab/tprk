
import pandas as pd
import sys
import os
import argparse
from Bio.Seq import Seq

def translate_nucs(read_seq):
	coding_dna = Seq(read_seq)
	translation = str(coding_dna.translate())
	return translation

parser = argparse.ArgumentParser()

parser.add_argument('allreads_filtered')

parser.add_argument('metadata')

args = parser.parse_args()

#df = pd.read_csv('/Users/administrator/Downloads/allreads_filtered.csv')

#meta = pd.read_csv('/Users/administrator/Desktop/Illumina_PacBio/metadata_trimmed_pacbio_Ill_test_copy.csv')

df = pd.read_csv(args.allreads_filtered)

#Removes Zeros
df = df.loc[(df.sum(axis=1) != 0)]

meta = pd.read_csv(args.metadata)

df['Read AA'] = ""

for q in range(len(df)):
    if (pd.isnull(df.iloc[q,1]) != True):
        df.iloc[q,len(df.columns)-1] = translate_nucs(df.iloc[q,1])

#We wont need this for pipeline:
df = df[df.columns.drop(list(df.filter(regex='_RelativeFreq')))]

df2 = df[['Region','Read','Read AA' ]]

#df2.to_csv('test3.csv', index = False, header = True)

sampleNum = int(len(df.columns))

for i in range(sampleNum - 2):

    for j in range(len(meta.index)):

        #check for sample pairs
        if(df.columns[i+2] == ("Ill_" + meta.iloc[j,0] + "_Count") and (meta.iloc[j,4] == "A")):
            sampleA = "Ill_" + meta.iloc[j,0] + "_Count"
            #print(sampleA)

            for k in range(len(meta.index)):
                if(meta.iloc[j,3] == meta.iloc[k,3] and meta.iloc[k,4] == "B"):
                    sampleB = "Ill_" + meta.iloc[k,0] + "_Count"
                    #print(sampleB)

            #set null to zero
            for l in range(len(df.index)):
                if(pd.isnull(df.loc[l,sampleA]) or pd.isnull(df.loc[l,sampleB])):
                    df.loc[l, sampleA] = 0
                    df.loc[l, sampleB] = 0

            #Combine and average columns
            sum = pd.DataFrame((df[sampleA] + df[sampleB])/2)
            sum.columns=[meta.iloc[j,3]]

            df2 = pd.concat([df2, sum], axis=1)

#sum variable regions

df2.to_csv('tprk_count_technical_rep.csv', index = False, header = True)

sampleCol = int(len(df2.columns))

test = len(df2.index)

d = 0

newCount = len(df2.index)

V1 = [0] * (sampleCol-2)
V2 = [0] * (sampleCol-2)
V3 = [0] * (sampleCol-2)
V4 = [0] * (sampleCol-2)
V5 = [0] * (sampleCol-2)
V6 = [0] * (sampleCol-2)
V7 = [0] * (sampleCol-2)

for m in range(len(df2.index)):
    for n in range(sampleCol - 3):
        if(df2.loc[m,"Region"] == 'V1'):
            V1[n] = V1[n] + df2.iloc[m,n+3]

        if(df2.loc[m, "Region"] == 'V2'):
            V2[n] = V2[n] + df2.iloc[m, n+3]

        if(df2.loc[m, "Region"] == 'V3'):
            V3[n] = V3[n] + df2.iloc[m, n+3]

        if(df2.loc[m, "Region"] == 'V4'):
            V4[n] = V4[n] + df2.iloc[m, n+3]

        if(df2.loc[m, "Region"] == 'V5'):
            V5[n] = V5[n] + df2.iloc[m, n+3]

        if(df2.loc[m, "Region"] == 'V6'):
            V6[n] = V6[n] + df2.iloc[m, n+3]

        if (df2.loc[m, "Region"] == 'V7'):
            V7[n] = V7[n] + df2.iloc[m, n+3]

for o in range(len(df.index)):
    for p in range(sampleCol - 3):
        if(df2.iloc[o,p+3] != 0):
            if(df2.loc[o,"Region"] == 'V1'):
                df2.iloc[o,p+3] = (df2.iloc[o,p+3]/(V1[p])*100)

            if(df2.loc[o, "Region"] == 'V2'):
                df2.iloc[o,p+3] = (df2.iloc[o,p+3]/V2[p])*100

            if(df2.loc[o, "Region"] == 'V3'):
                df2.iloc[o,p+3] = (df2.iloc[o,p+3]/V3[p])*100

            if(df2.loc[o, "Region"] == 'V4'):
                df2.iloc[o,p+3] = (df2.iloc[o,p+3]/V4[p])*100

            if(df2.loc[o, "Region"] == 'V5'):
                df2.iloc[o,p+3] = (df2.iloc[o,p+3]/V5[p])*100

            if(df2.loc[o, "Region"] == 'V6'):
                df2.iloc[o,p+3] = (df2.iloc[o,p+3]/V6[p])*100

            if (df2.loc[o, "Region"] == 'V7'):
                df2.iloc[o,p+3] = (df2.iloc[o,p+3]/V7[p])*100

df2.to_csv('tprk_percent_technical_rep.csv', index = False, header = True)
