#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#baseline correction package
from pybaselines import Baseline, utils
from scipy.integrate import trapezoid
from functools import reduce

#open hypertable
tri=pd.read_csv("21_Dic_iter0_TCA_qc_hypertable.csv")
#file from BO iteration, just a label to save files accordingly
itera="0"

#check unique precursor/fragments masses"
sci=pd.unique(tri["scan_type"])

#compounds names are in this table, it was exported from the MS method
values=pd.read_csv("SRM Table_TCA.csv")

#sort by polarity/precursor
val2=values.sort_values(["Polarity","Precursor (m/z)"])

#separate between negative and positive polarity
valpol=val2.loc[val2["Polarity"]=="Positive"]
valneg=val2.loc[val2["Polarity"]=="Negative"]

#list for unique compounds in positive mode
poslist=list(pd.unique(valpol["Compound"]))

#list for unique compounds in negative mode
neglist=list(pd.unique(valneg["Compound"]))

#add lists
lista=poslist+neglist

#check if number of unique precursor/fragments in the hypertable match with
#the number of unique compounds in the SRM exported table
print(len(sci))
print(len(lista))

#make dictionary to transform the names into compounds
dictio=dict(zip(sci,lista))

#do the name conversion
tri["scan_type"]=tri["scan_type"].apply(lambda x: dictio[x] if x in list(dictio.keys()) else "file")

#open file with list of compounds whose peaks are horrible by visual inspection
bad=pd.read_csv("bad_quality.txt", header=None)

#list of these compounds
badlist=bad[0].to_list()

#compounds can be analysed further
goodtri=tri.loc[~tri["scan_type"].isin(badlist)]

#small copy of goodtri to separate code
df=goodtri.copy()

#open a list of compounds that don't need correction
uncorr=pd.read_csv("uncorrected.txt",header=None)
uncolist=uncorr[0].to_list()

#open a list of compounds that need correction
tocorr=pd.read_csv("to_correct.txt",header=None)
toclist=tocorr[0].to_list()

#create separated datframe for compounds that don't need correction and those who need it
df_good=df.loc[df["scan_type"].isin(uncolist)].copy()
df_to_correct=df.loc[df["scan_type"].isin(toclist)].copy()


savefilt=[]
#uncorrected
for i in pd.unique(df_good["scan_type"]):   
    scan_type_to_plot = i

    # filter the dataframe to include only the desired scan_type
    df_filtered = df_good[df_good['scan_type'] == scan_type_to_plot].copy()
    #no correction is needed
    #integration using trapezoid rule
    integral=df_filtered.groupby("file").apply(lambda x: trapezoid(x["tic"],x=x["rt"])).reset_index()
    integral.columns=["file","integral"]
    integral["molecule"]=i
    #get minimum and maximum area value
    minimo=integral["integral"].min()
    maximo=integral["integral"].max()
    #append dataframe result to savefilt 
    savefilt.append(integral)
    print(minimo,maximo,i)
     
    #plot the tic vs rt as a line plot, colored by file
    plt.figure(figsize=(16, 12))
    sns.lineplot(data=df_filtered, x='rt', y='tic', hue='file')
    
    # add plot title and axis labels
    plt.title(f"TIC vs RT for {scan_type_to_plot}")
    plt.xlabel("Retention Time (RT)")
    plt.ylabel("Total Ion Current (TIC)")
    plt.savefig(scan_type_to_plot+"_file_plot_iter_"+itera+"_uncorrected.png",dpi=300,bbox_inches="tight")
    # show the plot
    #plt.show()
    plt.close()
    
#corrected
for i in pd.unique(df_to_correct["scan_type"]):
    
    scan_type_to_plot = i

    # filter the dataframe to include only the desired scan_type
    df_filtered = df_to_correct[df_to_correct['scan_type'] == scan_type_to_plot].copy()
    #split the dataframe by file, creating a new dataframe for each
    #list to add dataframes per file
    juy=[]
    for filo in pd.unique(df_filtered["file"]):
        #select file on dataframe
        intor=df_filtered.loc[df_filtered["file"]==filo].copy()
        #declare the baseline correction filter
        baseline_fitter = Baseline(x_data=intor["rt"])
        #fit the asls algorithm to get the baseline
        #bkg_1 = baseline_fitter.modpoly(df_filtered["tic"], poly_order=3)[0]
        bkg_2 = baseline_fitter.asls(intor["tic"], lam=1e9, p=0.01)[0]
        #remove the baseline from the original intensity
        intor["tic"]=intor["tic"]-bkg_2
        #append file dataframe to empty list
        juy.append(intor)
    #reconstruct the dataframe, now eith the corrected values
    df_new=reduce(lambda x,y:pd.concat([x,y]),juy)
    #perform the peak integration (area under the curve step), same as for uncorrected
    integral=df_new.groupby("file").apply(lambda x: trapezoid(x["tic"],x=x["rt"])).reset_index()
    integral.columns=["file","integral"]
    integral["molecule"]=i
    #get minimum and maximum area value
    minimo=integral["integral"].min()
    maximo=integral["integral"].max()
    #append dataframe result to savefilt 
    savefilt.append(integral)
    print(minimo,maximo,i)
    # add plot title and axis labels
    plt.figure(figsize=(16, 12))
    sns.lineplot(data=df_filtered, x='rt', y='tic', hue='file')
    plt.title(f"TIC vs RT for {scan_type_to_plot}")
    plt.xlabel("Retention Time (RT)")
    plt.ylabel("Total Ion Current (TIC)")
    plt.savefig(scan_type_to_plot+"_file_plot_iter_"+itera+"_byfile.png",dpi=300,bbox_inches="tight")
    #plt.show()
    plt.close()
    sns.lineplot(data=df_new, x='rt', y='tic', hue='file')
    plt.title(f"TIC vs RT for {scan_type_to_plot}")
    plt.xlabel("Retention Time (RT)")
    plt.ylabel("Total Ion Current (TIC)")
    plt.savefig(scan_type_to_plot+"_file_plot_iter_"+itera+"_corrected_byfile.png",dpi=300,bbox_inches="tight")
    #plt.show()
    plt.close()

#savefilt contains peak area (area under the curve) data for each compound and each file
info=reduce(lambda x,y:pd.concat([x,y]),savefilt)
info=info.sort_values("file")
lic=pd.unique(info["file"])