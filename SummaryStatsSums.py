import pandas as pd
import sys
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('all_summary_stats')

parser.add_argument('metadata')

args = parser.parse_args()

meta = pd.read_csv(args.metadata)

df = pd.read_csv(args.all_summary_stats)

#df = df.sort_values(by=['Sample'])
#df = pd.read_csv('./all_summary_stats.csv')

df['V_region_summary'] = ''
df['% input'] = ''

sampleNum = int(len(df.index)-1)

regionSum = [0]*sampleNum
percentInput = [0]*sampleNum

for i in range(len(df.index)-1):
    for j in range(7):
        regionSum[i] = regionSum[i] + df.iloc[i,2+j]
    percentInput[i] =  ( regionSum[i] / df.iloc[i, 1] ) * 100

for q in range(len(df.index)-1):
    df.iloc[q, 9] = regionSum[q]
    df.iloc[q, 10] = percentInput[q]

df.to_csv('More_summary_stats.csv', index = False, header = True)

df2 = pd.read_csv('./all_summary_stats.csv')

sampleNum = int(len(df2.index)-1)

regionSum = [0]*sampleNum
percentInput = [0]*sampleNum

for i in range(len(df.index)-1):
    for j in range(7):
        if(df2.iloc[i,2+j] < 5000):
            print(df2.iloc[i,2+j])
        else: df2.iloc[i,2+j] = ""

df2.to_csv('less_than_5000_all_summary_stats.csv', index = False, header = True)

#meta = pd.read_csv(args[2])

df['Rep Ratio'] = ''

sampleNum = int(len(meta.index))

for i in range(sampleNum):
    if(meta.iloc[i,4] == "A"):

        id = "Ill_" + (meta.iloc[i, 0])

        test = df[df['Sample'].str.contains(id)].index

        for j in range(sampleNum):
            if(meta.iloc[j,3] == meta.iloc[i,3] and meta.iloc[j,4] == 'B'):
                id2 = "Ill_" + (meta.iloc[j, 0])

                test2 = df[df['Sample'].str.contains(id2)].index

        if(int(df.iloc[test2, 9]) != 0):

            ratio = int(df.iloc[test, 9])/int(df.iloc[test2, 9])

            df.iloc[test, 11] = ratio

            df.iloc[test2, 11] = ratio

df.to_csv('technicalRep_all_summary_stats.csv', index = False, header = True)

for k in range(sampleNum):
    if(df.iloc[k,11] != ""):
         if (df.iloc[k, 9] < 70):
             for m in range(12):
                df.iloc[k, m] = ""
    if (df.iloc[k, 11] != ""):
         if(df.iloc[k,11] < 0.8 or df.iloc[k,11] > 1.25):
             for l in range(12):
                 df.iloc[k, l] = ""
    # if(df.iloc[k,11] != ""):
    #     if(df.iloc[k,11] < 0.8 or df.iloc[k,11] > 1.25):
    #         for l in range(12):
    #             df.iloc[k, l] = ""

df = df.drop_duplicates()

df.to_csv('technicalRepCheck_all_summary_stats.csv', index = False, header = True)

numSamp = (sum(df['Rep Ratio'] != "?")) - 1

df2 = pd.DataFrame(index=range(numSamp),columns=range(10))

df2.columns = ["Sample","Total Input Reads","V1_Reads","V2_Reads","V3_Reads","V4_Reads","V5_Reads","V6_Reads","V7_Reads","V_region_summary"]

for n in range(numSamp):
    if (df.iloc[n, 11] != ""):
        #print(df.iloc[n,0])
        #for o in range(sampleNum):

        test

        test2

        for o in range(sampleNum):
            if (("Ill_" + meta.iloc[o,0]) == df.iloc[n,0]):
                for p in range(sampleNum):
                    if(meta.iloc[o,4] == "A"):
                        if (meta.iloc[p, 3] == meta.iloc[o, 3] and meta.iloc[p, 4] == 'B'):
                            test = "Ill_" + meta.iloc[p,0]
                            test2 = meta.iloc[p, 3]
                    if(meta.iloc[o, 4] == "B"):
                        if (meta.iloc[p, 3] == meta.iloc[o, 3] and meta.iloc[p, 4] == 'A'):
                            test = "Ill_" + meta.iloc[p,0]
                            test2 = meta.iloc[p, 3]
        check = 0

        for r in range(numSamp):
            if(df.iloc[r, 0] == test):
                df2.iloc[n, 0] = test2
                for q in range(9):
                    print(df.iloc[n,0])
                    print(df.iloc[r,0])
                    df2.iloc[n,q+1] = (df.iloc[n,q+1] + df.iloc[r,q+1])/2
                    check = 1
                    print("here")


        if(check == 0):
            df2.iloc[n, 0] = test2
            for s in range(9):
                df2.iloc[n,s+1] = df.iloc[n,s+1]

#df2 = df2.drop_duplicates()

df2.to_csv('average_technicalRepCheck_all_summary_stats.csv', index = False, header = True)
