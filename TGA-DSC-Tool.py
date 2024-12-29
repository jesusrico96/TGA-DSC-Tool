import sys
from tkinter import *
import pandas as pd
import numpy as np
import re
from tkinter import filedialog
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid, simpson
import math
import os

def browseFiles1():
    global filename1, folderpath, experiments, df2, cells, df1, experimentname
    filename1 = filedialog.askopenfilename(initialdir=sys.path[0],
                                           title="Choose Sample File",
                                           filetypes=(("Excel Files",
                                                       "*.xls*"),
                                                      ("Text Files",
                                                       "*.txt*"),
                                                      ("All Files",
                                                       "*.*")))

    # Used file label
    label_file_explorer.configure(text=filename1)
    filelabel = os.path.basename(filename1)
    experimentname = str(filelabel).split("_")
    print("Filename is: "+str(filelabel))

    print("Folderpath is: "+sys.path[0])
    folderpath = sys.path[0]

    print("Reading now: "+str(filename1))
    df1 = pd.read_excel(filename1)

    # Locate columns names by looking for a cell containing
    # "Time (s)" and split in two dataframes: Header and rest
    d = dict(zip(df1.columns, range(len(df1.columns))))
    s = df1.rename(columns=d).stack()
    a = (s == "Time (s)").idxmax()
    dfh = df1.iloc[:a[0],:]
    df2 = df1.iloc[a[0]:,:]

    # Update the column names to be correct
    df2.columns = df2.iloc[0]
    df2 = df2[1:]

    # Count the number of experiments in file
    ex = (dfh == "# Column mapping:").idxmax()
    experiments = int(ex.iloc[0]/9)
    print("Number of experiments: "+str(experiments))
    
    # Replace missing data with NaNs
    df2 = df2.replace("-",np.nan).infer_objects(copy=False)

    # Determination of size of the dataframe to cover
    cells = range(1, df2.shape[1], 1)


def hfcorr():
    #Plot heat flow graphs vs temperature all at once, QRAD-CORRECTED

    if "20mg" in filename1:
        mass = 20.0
    else:
        mass = 5.0
    print("Filename mass is " + str(mass) + " mg")

    if "10HR" in filename1:
        hr = 10
    elif "2,5HR" in filename1:
        hr = 2.5
    else:
        hr = 5
    print("Detected heating rate is "+str(hr)+" K/min")


    df3 = pd.DataFrame(df2[:].values)
    df4 = df2.columns
    dfheat = pd.DataFrame()
    dfheattemp = pd.DataFrame()
    dfderiv = pd.DataFrame()
    dfmasstemp = pd.DataFrame()
    dfmass = pd.DataFrame()
    for i in cells:
        if "(mW)" in df4[i]:
            dfheat[str(df4[i])] = df3[i]
            dfheattemp[str(df4[i-1])] = df3[i-1]
        if "(mg/s)" in df4[i]:
            dfderiv[str(df4[i])] = df3[i]
        if "(mg)" in df4[i]:
            dfmass[str(df4[i])] = df3[i]
            dfmasstemp[str(df4[i-1])] = df3[i-1]
        else:
            continue

    dfheat["heataverage"] = dfheat.mean(axis=1)

    dfheattemp["tempaverage"] = dfheattemp.mean(axis=1)
    if dfheattemp.at[0, "tempaverage"]>250:
        dfheattemp = dfheattemp-273.15
    dfderiv["derivaverage"] = dfderiv.mean(axis=1)
    dfmasstemp["masstempaverage"] = dfmasstemp.mean(axis=1)
    if dfmasstemp.at[0, "masstempaverage"]>250:
        dfmasstemp = dfmasstemp-273.15
    dfmass["massaverage"] = dfmass.mean(axis=1)
    
    dfheattemp["tempaverage"] = dfheattemp.mean(axis=1)
    dfderiv["derivaverage"] = dfderiv.mean(axis=1)
    dfmasstemp["masstempaverage"] = dfmasstemp.mean(axis=1)
    dfmass["massaverage"] = dfmass.mean(axis=1)

    hf = dfheat["heataverage"]
    heattemp = dfheattemp["tempaverage"]
    deriv = dfderiv["derivaverage"]
    masstemp = dfmasstemp["masstempaverage"]
    initialmass = dfmass["massaverage"].drop(range(25, dfmass.shape[0], 1))
    initialmassavg = initialmass.mean()
    finalmass = dfmass["massaverage"].drop(range(0, dfmass.shape[0]-25, 1))
    finalmassavg = finalmass.mean()
    print("Detected initial average mass is: "+str(initialmassavg)+" mg")
    print("Detected final average mass is: "+str(finalmassavg)+" mg")

    dfrad = pd.DataFrame()
    listalpha = []
    listq = []

    # Determination of size of the dataframe to cover
    dimension = range(1, dfmass.shape[0], 1)

    # Correction based on radiation
    for i in dimension:
        alpha_x = ((dfmass.at[0, "massaverage"]-dfmass.at[i, "massaverage"])/(dfmass.at[0, "massaverage"]-dfmass.at[dfmass.shape[0]-1, "massaverage"]))
        listalpha.append(alpha_x)

        # Surface area of the particle, diminishes with time due to the effect of pyrolysis
        As = alpha_x*math.pi*pow((0.003),2)+(1-alpha_x)*math.pi*pow((0.003),2)

        # Emissivity of the particle, higher as it blackens due to pyrolysis
        eps = alpha_x*0.95+(1-alpha_x)*0.6

        # Stefan-Boltzmann constant
        sigma = 5.67e-8

        # Wall temperature
        # tw = dfmasstemp.at[i, "masstempaverage"] + 220
        tw = alpha_x*(dfmasstemp.at[i, "masstempaverage"] + 10) + (1-alpha_x)*(dfmasstemp.at[i, "masstempaverage"] + 20)

        Qrad = 1000*eps*As*sigma*(pow((tw+273.15),4)-pow((dfmasstemp.at[i, "masstempaverage"]+273.15),4))

        listq.append(Qrad)

    dfrad["alpha"] = listalpha
    dfrad["Qrad"] = listq
    Q = dfrad["Qrad"]
    zeros = [0]*len(hf)
    heattempQ = heattemp.drop(heattemp.tail(1).index)
    hfcorr = (hf-Q)/initialmassavg
    dev = (np.std(dfheat, axis=1))/initialmassavg


    # ----------------CORRECTION OF HEAT DUE TO WOOD AND CHAR HEATING---------------------------

    dfcpwood = pd.DataFrame()
    listcpwood = []
    dfcpchar = pd.DataFrame()
    listcpchar = []
    lenheat = range(0, heattemp.shape[0], 1)
    for i in lenheat:
        cpwood = 1113.68+4.8567*(dfheattemp.at[i, "tempaverage"]) # Formula from Rath 2003 (https://doi.org/10.1016/S0016-2361(02)00138-2)
        listcpwood.append(cpwood)
    
    dfcpwood["cpwood"] = listcpwood
    cp = dfcpwood["cpwood"]
    listqwood = [] 

    for i in dimension:
        alpha_x = ((dfmass.at[0, "massaverage"]-dfmass.at[i, "massaverage"])/(dfmass.at[0, "massaverage"]-dfmass.at[dfmass.shape[0]-1, "massaverage"]))
        cpwood = 1113.68+4.8567*(dfheattemp.at[i, "tempaverage"]) # Formula from Rath 2003 (https://doi.org/10.1016/S0016-2361(02)00138-2)
        Qwood = 1000*(0.95-alpha_x)*(initialmassavg/1000000)*cpwood*(hr/60) # This is prepared to give mW, Formula from Rath 2003 (https://doi.org/10.1016/S0016-2361(02)00138-2)
        listqwood.append(Qwood)

    for i in lenheat:
        cpchar = (8314/5.75)*(pow(math.e, 380/(dfheattemp.at[i, "tempaverage"]+273.15))*pow(((pow(math.e, 380/(dfheattemp.at[i, "tempaverage"]+273.15))-1)/(380/(dfheattemp.at[i, "tempaverage"]+273.15))), -2) + 2*pow(math.e, 1800/(dfheattemp.at[i, "tempaverage"]+273.15))*pow(((pow(math.e, 1800/(dfheattemp.at[i, "tempaverage"]+273.15))-1)/(1800/(dfheattemp.at[i, "tempaverage"]+273.15))), -2)) # Formula from Rath 2003 (https://doi.org/10.1016/S0016-2361(02)00138-2)
        listcpchar.append(cpchar)
    
    dfcpchar["cpchar"] = listcpchar
    cpchar = dfcpchar["cpchar"]
    listqchar = [] 

    for i in dimension:
        alpha_x = ((dfmass.at[0, "massaverage"]-dfmass.at[i, "massaverage"])/(dfmass.at[0, "massaverage"]-dfmass.at[dfmass.shape[0]-1, "massaverage"]))
        cpchar = (8314/5.75)*(pow(math.e, 380/(dfheattemp.at[i, "tempaverage"]+273.15))*pow(((pow(math.e, 380/(dfheattemp.at[i, "tempaverage"]+273.15))-1)/(380/(dfheattemp.at[i, "tempaverage"]+273.15))), -2) + 2*pow(math.e, 1800/(dfheattemp.at[i, "tempaverage"]+273.15))*pow(((pow(math.e, 1800/(dfheattemp.at[i, "tempaverage"]+273.15))-1)/(1800/(dfheattemp.at[i, "tempaverage"]+273.15))), -2)) # Formula from Rath 2003 (https://doi.org/10.1016/S0016-2361(02)00138-2)
        Qchar = 1000*(alpha_x)*(initialmassavg/1000000)*cpchar*(hr/60) # This is prepared to give mW, Formula from Rath 2003 (https://doi.org/10.1016/S0016-2361(02)00138-2)
        listqchar.append(Qchar)
    
    dfqwood = pd.DataFrame()
    dfqwood["qwood"] = listqwood

    # Qwoodcp, in mW
    Qwoodcp = dfqwood["qwood"]

    dfqchar = pd.DataFrame()
    dfqchar["qchar"] = listqchar

    # Qwcharcp, in mW
    Qcharcp = dfqchar["qchar"]

    # ----------------END OF CORRECTION OF HEAT DUE TO WOOD AND CHAR HEATING---------------------------

    zeros = [0]*len(hf)
    hfcorr = (hf-Q-Qwoodcp-Qcharcp)/initialmassavg
    hfcorr[np.isnan(hfcorr)] = 0
    dev = (np.std(dfheat, axis=1))/initialmassavg

    diff = len(heattemp)-len(hfcorr)
    if diff < 0:
        hfcorr = hfcorr.drop(hfcorr.tail(abs(diff)).index)
    else:
        heattemp = heattemp.drop(heattemp.tail(abs(diff)).index)

    
    # ----------------CONVERSION OF THE TEMP. AXIS INTO A TIME AXIS---------------------------
    if hr==10:
        tailcut = 2100
    elif hr==2.5:
        tailcut = 2900
    else:
        tailcut = 2100


    hfcorr2 = hfcorr.drop(hfcorr.tail(tailcut).index)
    heattemp2 = heattemp.drop(heattemp.tail(tailcut).index) # Check visually that range of temperatures: 150-550 C

    
    # Converting temperature axis to time axis. Example, if HR=10 K/min, then 1/6 K/s
    hrsec = hr/60
    heattemp2 = pd.DataFrame([x for x in heattemp2 if x > 150]) # Only integrating for temperatures greater than 150 C
    heattemp2 = heattemp2[0]

    hfcorr2 = hfcorr2.drop(hfcorr.head(hfcorr2.shape[0]-heattemp2.shape[0]).index) # Equalizing the size of x and y axes

    reftemp = heattemp2[0]
    heattemp2 = (heattemp2-reftemp)/hrsec # Now heattemp2 is a time axis
    heattempt = (heattemp2*hrsec)+reftemp # Temperature axis is calculated as well for representation purposes

    # ----------------END OF THE CONVERSION OF THE TEMP. AXIS INTO A TIME AXIS---------------------------

    zeros2 = [0]*len(hfcorr2)

    # Integration of the area with the time axis to calculate heat of reaction
    area = trapezoid(y = -hfcorr2-zeros2, x = heattemp2)
    print("Area with trapezoids is "+str(area)+" kJ/kg")
    area2 = simpson(y = -hfcorr2-zeros2, x = heattemp2)
    print("Area with Simpson is "+str(area2)+" kJ/kg")

    # --------------------------------PLOTTING EXPLANATORY FIGURES------------------------------------------
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Temperature (°C)')
    ax1.set_ylabel('Heat flow (mW/g)', color="xkcd:blue")
    ax1.spines['left'].set_color("xkcd:blue")
    ax1.tick_params(axis='y', labelcolor="xkcd:blue", color="xkcd:blue")
    ax1.plot(heattemp, -hfcorr, label=experimentname[0]+" Heat Flow", linewidth=0.75)
    ax1.plot(heattemp, -hfcorr+dev, color="black", linewidth=0.5)
    ax1.plot(heattemp, -hfcorr-dev, color="black", linewidth=0.5)
    ax1.plot(heattemp, hf/initialmassavg, color="red", linewidth=1, label = "Uncorrected HF")
    ax1.plot(heattempQ, Q/initialmassavg, color="green", linewidth=0.5, label = "Correction line")
    ax1.plot(heattemp.drop(heattemp.tail(1).index), Qwoodcp/initialmassavg, color="xkcd:military green", linewidth=1.5, label = "Correction line: Qwood")
    ax1.plot(heattemp.drop(heattemp.tail(1).index), Qcharcp/initialmassavg, color="xkcd:bright violet", linewidth=1.5, label = "Correction line: Qchar")
    ax1.plot(heattemp, zeros, color="black", linewidth=0.5)
    
    ax1.fill_between(dfheattemp["tempaverage"], -hfcorr+dev, -hfcorr-dev, color="xkcd:powder blue")
    ax1.fill_between(heattempt, -hfcorr2, zeros2, color="pink") # While the integration has been made with the time x-axis, representation is much clearer with the temperature x-axis
    ax1.set_ylim(-9,9)
    ax1.set_xlim(100,550)

    ax2 = ax1.twinx()
    ax2.set_ylabel('dTG (mg/s)', color="xkcd:burnt orange")
    ax2.set_xlabel('Temperature (°C)')
    ax2.spines['left'].set_color("xkcd:blue")
    ax2.spines['right'].set_color("xkcd:burnt orange")
    ax2.tick_params(axis='y', labelcolor="xkcd:burnt orange", color="xkcd:burnt orange")
    ax2.plot(masstemp, -deriv, label=experimentname[0]+" Mass derivative", color="xkcd:burnt orange", linewidth=0.75)
    plt.xlabel('Temperature (°C)')
    # ax1.set_ylim(-0.025,0.005)
    ax2.set_ylim(bottom=0)

    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=2)

    fig.set_size_inches(10, 5)
    plt.title(str(os.path.basename(filename1)) + " plot corrected with modeled radiation EXPLANATORY FIGURE")

    plt.show()


    # --------------------------------PLOTTING CLEAN FIGURES------------------------------------------
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Temperature (°C)')
    ax1.set_ylabel('Heat flow (mW/g)', color="xkcd:blue")
    ax1.spines['left'].set_color("xkcd:blue")
    ax1.tick_params(axis='y', labelcolor="xkcd:blue", color="xkcd:blue")
    ax1.plot(heattemp, -hfcorr, label=experimentname[0]+" Heat Flow", linewidth=0.75)
    ax1.plot(heattemp, -hfcorr+dev, color="black", linewidth=0.5)
    ax1.plot(heattemp, -hfcorr-dev, color="black", linewidth=0.5)
    # ax1.plot(heattemp, hf/initialmassavg, color="red", linewidth=1, label = "Uncorrected HF")
    # ax1.plot(heattempQ, Q/initialmassavg, color="green", linewidth=0.5, label = "Correction line")
    ax1.plot(heattemp, zeros, color="black", linewidth=0.5)
    
    ax1.fill_between(dfheattemp["tempaverage"], -hfcorr+dev, -hfcorr-dev, color="xkcd:powder blue")
    # ax1.fill_between(heattempt, -hfcorr2, zeros2, color="pink") # While the integration has been made with the time x-axis, representation is much clearer with the temperature x-axis
    ax1.set_ylim(-6,6)
    ax1.set_xlim(100,550)

    ax2 = ax1.twinx()
    ax2.set_ylabel('dTG (mg/s)', color="xkcd:burnt orange")
    ax2.set_xlabel('Temperature (°C)')
    ax2.spines['left'].set_color("xkcd:blue")
    ax2.spines['right'].set_color("xkcd:burnt orange")
    ax2.tick_params(axis='y', labelcolor="xkcd:burnt orange", color="xkcd:burnt orange")
    ax2.plot(masstemp, -deriv, label=experimentname[0]+" Mass derivative", color="xkcd:burnt orange", linewidth=0.75)
    plt.xlabel('Temperature (°C)')
    # ax1.set_ylim(-0.025,0.005)
    ax2.set_ylim(bottom=0)

    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=2)

    fig.set_size_inches(10, 5)
    plt.title(str(os.path.basename(filename1)) + " plot corrected with modeled radiation")

    plt.show()


def hfcorr_char():
    #Plot heat flow graphs vs temperature, Char+Rad-CORRECTED, calculates heat of reaction in the range 150-550 ºC
    # Initial mass of the sample


    if "20mg" in filename1:
        mass = 20.0
    else:
        mass = 5.0
    print("Filename mass is " + str(mass) + " mg")

    if "10HR" in filename1:
        hr = 10
    elif "2,5HR" in filename1:
        hr = 2.5
    else:
        hr = 5

    print("Detected heating rate is "+str(hr)+" K/min")

    df3 = pd.DataFrame(df2[:].values)
    df4 = df2.columns
    dfheat = pd.DataFrame()
    dfheattemp = pd.DataFrame()
    dfderiv = pd.DataFrame()
    dfmasstemp = pd.DataFrame()
    dfmass = pd.DataFrame()
    for i in cells:
        if "(mW)" in df4[i]:
            dfheat[str(df4[i])] = df3[i]
            dfheattemp[str(df4[i-1])] = df3[i-1]
        if "(mg/s)" in df4[i]:
            dfderiv[str(df4[i])] = df3[i]
        if "(mg)" in df4[i]:
            dfmass[str(df4[i])] = df3[i]
            dfmasstemp[str(df4[i-1])] = df3[i-1]
        else:
            continue

    dfheat["heataverage"] = dfheat.mean(axis=1)
    
    dfheattemp["tempaverage"] = dfheattemp.mean(axis=1)
    if dfheattemp.at[0, "tempaverage"]>250:
        dfheattemp = dfheattemp-273.15
    dfderiv["derivaverage"] = dfderiv.mean(axis=1)
    dfmasstemp["masstempaverage"] = dfmasstemp.mean(axis=1)
    if dfmasstemp.at[0, "masstempaverage"]>250:
        dfmasstemp = dfmasstemp-273.15
    dfmass["massaverage"] = dfmass.mean(axis=1)

    hf = dfheat["heataverage"]
    heattemp = dfheattemp["tempaverage"]
    deriv = dfderiv["derivaverage"]
    masstemp = dfmasstemp["masstempaverage"]
    initialmass = dfmass["massaverage"].drop(range(25, dfmass.shape[0], 1))
    initialmassavg = initialmass.mean()
    finalmass = dfmass["massaverage"].drop(range(0, dfmass.shape[0]-25, 1))
    finalmassavg = finalmass.mean()
    print("Detected initial average mass is: "+str(initialmassavg)+" mg")
    print("Detected final average mass is: "+str(finalmassavg)+" mg")

    # ----------------START OF THE CHAR BASELINE PART---------------------------
    # Now, it will ask to provide a char file to serve as baseline

    filename_char = filedialog.askopenfilename(initialdir=sys.path[0],
                                           title="Choose Char File",
                                           filetypes=(("Excel Files",
                                                       "*.xls*"),
                                                      ("Text Files",
                                                       "*.txt*"),
                                                      ("All Files",
                                                       "*.*")))

    # Used file label
    filelabel_char = os.path.basename(filename_char)
    experimentname_char = str(filelabel_char).split("_")
    print("Char filename is: "+str(filelabel_char))

    print("Char folderpath is: "+sys.path[0])
    folderpath_char = sys.path[0]

    print("Reading now char: "+str(filename_char))
    df1_char = pd.read_excel(filename_char)

    # Locate columns names by looking for a cell containing
    # "Time (s)" and split in two dataframes: Header and rest
    d_char = dict(zip(df1_char.columns, range(len(df1_char.columns))))
    s_char = df1_char.rename(columns=d_char).stack()
    a_char = (s_char == "Time (s)").idxmax()
    dfh_char = df1_char.iloc[:a_char[0],:]
    df2_char = df1_char.iloc[a_char[0]:,:]

    # Update the column names to be correct
    df2_char.columns = df2_char.iloc[0]
    df2_char = df2_char[1:]

    # Count the number of experiments in file
    ex_char = (dfh_char == "# Column mapping:").idxmax()
    experiments_char = int(ex_char.iloc[0]/9)
    print("Number of experiments in char file: "+str(experiments_char))
    
    # Replace missing data with NaNs
    df2_char = df2_char.replace("-",np.nan).infer_objects(copy=False)

    # Determination of size of the dataframe to cover
    cells_char = range(1, df2_char.shape[1], 1)

    df3_char = pd.DataFrame(df2_char[:].values)
    df4_char = df2_char.columns
    dfheat_char = pd.DataFrame()
    dfheattemp_char = pd.DataFrame()
    dfmass_char = pd.DataFrame()

    for i in cells_char:
        if "(mW)" in df4_char[i]:
            dfheat_char[str(df4_char[i])] = df3_char[i]+(0.5*initialmassavg) # <--------------------------------Adjust compensation parameter so "uncorrected HF" and "correction line: char" coincide at max.
            dfheattemp_char[str(df4_char[i-1])] = df3_char[i-1]
        if "(mg)" in df4[i]:
            dfmass_char[str(df4[i])] = df3[i]
        else:
            continue


    dfheat_char["heataverage"] = dfheat_char.mean(axis=1)
    
    dfheattemp_char["tempaverage"] = dfheattemp_char.mean(axis=1)
    if dfheattemp_char.at[0, "tempaverage"]>250:
        dfheattemp_char = dfheattemp_char-273.15
    dfmass_char["charmassaverage"] = dfmass_char.mean(axis=1)

    hf_char = dfheat_char["heataverage"]
    heattemp_char = dfheattemp_char["tempaverage"]

    dfhf = pd.DataFrame()
    listhf = []

    # Determination of size of the dataframe to cover
    dimension = range(0, min(dfheat_char.shape[0], dfmass.shape[0]), 1)

    # Correction based on char
    for i in dimension:
        alpha_x = ((dfmass.at[0, "massaverage"]-dfmass.at[i, "massaverage"])/(dfmass.at[0, "massaverage"]-dfmass.at[dfmass.shape[0]-1, "massaverage"]))
        hfalpha = alpha_x*dfheat_char.at[i, "heataverage"]
        listhf.append(hfalpha)
    
    dfhf["hfalpha"] = listhf



    # ----------------END OF THE CHAR BASELINE PART---------------------------

    # ----------------START OF THE RADIATION BASELINE PART---------------------------
    dfrad = pd.DataFrame()
    listalpha = []
    listq = []

    # Determination of size of the dataframe to cover
    dimension = range(0, dfmass.shape[0], 1)

    # Correction based on radiation
    for i in dimension:
        alpha_x = ((dfmass.at[0, "massaverage"]-dfmass.at[i, "massaverage"])/(dfmass.at[0, "massaverage"]-dfmass.at[dfmass.shape[0]-1, "massaverage"]))
        listalpha.append(alpha_x)

        # Surface area of the particle, decreases with time due to the effect of pyrolysis (m2)
        As = alpha_x*math.pi*pow((0.003),2)+(1-alpha_x)*math.pi*pow((0.003),2)

        # Emissivity of the particle, higher as it blackens due to pyrolysis
        eps = alpha_x*0.95+(1-alpha_x)*0.6

        # Stefan-Boltzmann constant (W/m2K4)
        sigma = 5.67e-8

        # Wall temperature
        # tw = dfmasstemp.at[i, "masstempaverage"] + 220
        tw = ((1-alpha_x)*(dfmasstemp.at[i, "masstempaverage"] + 35) + (alpha_x*(dfmasstemp.at[i, "masstempaverage"] + 0)))
        
        # Qrad, in mW. First line is out because I don't understand what the "100" is doing there
        # Qrad = (0.95-alpha_x)*100*eps*As*sigma*(pow((tw+273.15),4)-pow((dfmasstemp.at[i, "masstempaverage"]+273.15),4))
        Qrad = 1000*(0.95-alpha_x)*eps*As*sigma*(pow((tw+273.15),4)-pow((dfmasstemp.at[i, "masstempaverage"]+273.15),4)) # <-----------------------------------
        
        listq.append(Qrad)
    
    dfrad["alpha"] = listalpha
    dfrad["Qrad"] = listq
    Q = dfrad["Qrad"]

    # ----------------END OF THE RADIATION BASELINE PART---------------------------

    # ----------------CORRECTION OF HEAT DUE TO WOOD HEATING---------------------------

    dfcpwood = pd.DataFrame()
    listcpwood = []
    lenheat = range(0, heattemp.shape[0], 1)
    for i in lenheat:
        cpwood = 1113.68+4.8567*(dfheattemp.at[i, "tempaverage"]) # Formula from Rath 2003 (https://doi.org/10.1016/S0016-2361(02)00138-2)
        listcpwood.append(cpwood)
    
    dfcpwood["cpwood"] = listcpwood
    cp = dfcpwood["cpwood"]
    listqwood = [] 

    for i in dimension:
        alpha_x = ((dfmass.at[0, "massaverage"]-dfmass.at[i, "massaverage"])/(dfmass.at[0, "massaverage"]-dfmass.at[dfmass.shape[0]-1, "massaverage"]))
        cpwood = 1113.68+4.8567*(dfheattemp.at[i, "tempaverage"]) # Formula from Rath 2003 (https://doi.org/10.1016/S0016-2361(02)00138-2)
        Qwood = 1000*(0.95-alpha_x)*(initialmassavg/1000000)*cpwood*(hr/60) # This is prepared to give mW, Formula from Rath 2003 (https://doi.org/10.1016/S0016-2361(02)00138-2)
        listqwood.append(Qwood)
    
    dfqwood = pd.DataFrame()
    dfqwood["qwood"] = listqwood

    # Qwoodcp, in mW
    Qwoodcp = dfqwood["qwood"]

    # ----------------END OF CORRECTION OF HEAT DUE TO WOOD HEATING---------------------------

    zeros = [0]*len(hf)
    hfcorr = (hf-dfhf["hfalpha"]-Q-Qwoodcp)/initialmassavg
    hfcorr[np.isnan(hfcorr)] = 0
    dev = (np.std(dfheat, axis=1))/initialmassavg

    diff = len(heattemp)-len(hfcorr)
    if diff < 0:
        hfcorr = hfcorr.drop(hfcorr.tail(abs(diff)).index)
    else:
        heattemp = heattemp.drop(heattemp.tail(abs(diff)).index)

    diff_char = len(heattemp_char)-len(dfhf)
    if diff_char < 0:
        dfhf = dfhf.drop(dfhf.tail(abs(diff_char)).index)
    else:
        heattemp_char = heattemp_char.drop(heattemp_char.tail(abs(diff_char)).index)

    # ----------------CONVERSION OF THE TEMP. AXIS INTO A TIME AXIS---------------------------
    if hr==10:
        tailcut = 2100
    elif hr==2.5:
        tailcut = 2900
    else:
        tailcut = 2100


    hfcorr2 = hfcorr.drop(hfcorr.tail(tailcut).index)
    heattemp2 = heattemp.drop(heattemp.tail(tailcut).index) # Check visually that range of temperatures: 150-550 C

    
    # Converting temperature axis to time axis. Example, if HR=10 K/min, then 1/6 K/s
    hrsec = hr/60
    heattemp2 = pd.DataFrame([x for x in heattemp2 if x > 150]) # Only integrating for temperatures greater than 150 C
    heattemp2 = heattemp2[0]

    hfcorr2 = hfcorr2.drop(hfcorr.head(hfcorr2.shape[0]-heattemp2.shape[0]).index) # Equalizing the size of x and y axes

    reftemp = heattemp2[0]
    heattemp2 = (heattemp2-reftemp)/hrsec # Now heattemp2 is a time axis
    heattempt = (heattemp2*hrsec)+reftemp # Temperature axis is calculated as well for representation purposes

    # ----------------END OF THE CONVERSION OF THE TEMP. AXIS INTO A TIME AXIS---------------------------

    zeros2 = [0]*len(hfcorr2)


    # ----------------Plotting explanatory Figure----------------

    
    SMALL_SIZE = 14
    MEDIUM_SIZE = 16
    BIGGER_SIZE = 12

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=SMALL_SIZE)  # fontsize of the figure title
    
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Temperature (°C)')
    ax1.set_ylabel('Heat flow (W/g)', color="xkcd:blue") # The math behind the problem results in mW of heat flow divideed by the initial mass in mg, i.e. mW/mg == W/g
    ax1.spines['left'].set_color("xkcd:blue")
    ax1.tick_params(axis='y', labelcolor="xkcd:blue", color="xkcd:blue")
    ax1.plot(heattemp, -hfcorr, label=experimentname[0]+" Heat Flow", linewidth=0.75)

    # Integration of the area with the time axis to calculate heat of reaction
    area = trapezoid(y = -hfcorr2-zeros2, x = heattemp2)
    print("Area with trapezoids is "+str(area)+" kJ/kg")
    area2 = simpson(y = -hfcorr2-zeros2, x = heattemp2)
    print("Area with Simpson is "+str(area2)+" kJ/kg")

    # # Used for the raw heat flow drawings
    # zeros_raw = [0]*len(hf)
    # heattemp = dfheattemp["tempaverage"]
    # reftemp_raw = heattemp[0]
    # heattemp_raw = (heattemp-reftemp_raw)/hrsec
    # raw_area = trapezoid(y = -hf-zeros_raw, x = heattemp_raw)
    # print("RAW Area with trapezoids is "+str(raw_area)+" kJ/kg")



    ax1.plot(heattemp, -hfcorr+dev, color="black", linewidth=0.5)
    ax1.plot(heattemp, -hfcorr-dev, color="black", linewidth=0.5)
    ax1.plot(heattemp, hf/initialmassavg, color="red", linewidth=1, label = "Uncorrected HF")
    ax1.plot(heattemp_char, dfhf["hfalpha"]/initialmassavg, color="green", linewidth=0.5, label = "Correction line: Char")
    ax1.plot(heattemp, Q/initialmassavg, color="xkcd:vivid blue", linewidth=0.5, label = "Correction line: Rad")
    ax1.plot(heattemp, Qwoodcp/initialmassavg, color="xkcd:military green", linewidth=1.5, label = "Correction line: Qwood")
    ax1.plot(heattemp, zeros, color="black", linewidth=0.5)

    # ax1.plot(heattemp, zeros_raw, color="black", linewidth=0.5) # Used for the raw heat flow drawings
    
    ax1.fill_between(dfheattemp["tempaverage"], -hfcorr+dev, -hfcorr-dev, color="xkcd:powder blue")
    ax1.fill_between(heattempt, -hfcorr2, zeros2, color="pink") # While the integration has been made with the time x-axis, representation is much clearer with the temperature x-axis
    # ax1.plot(heattemp_raw, hf, color="red", linewidth=0.5) # Used for the raw heat flow drawings
    # ax1.fill_between(heattemp_raw, hf, zeros_raw, color="pink") # Used for the raw heat flow drawings

    if hr == 10:
        ax1.set_ylim(-2,2)
    elif hr == 2.5:
        ax1.set_ylim(-2,2)
    else:
        ax1.set_ylim(-2,2)

    ax1.set_xlim(100,550)

    ax2 = ax1.twinx()
    ax2.set_ylabel('dTG (mg/s)', color="xkcd:burnt orange")
    ax2.set_xlabel('Temperature (°C)')
    ax2.spines['left'].set_color("xkcd:blue")
    ax2.spines['right'].set_color("xkcd:burnt orange")
    ax2.tick_params(axis='y', labelcolor="xkcd:burnt orange", color="xkcd:burnt orange")
    ax2.plot(masstemp, -deriv, label=experimentname[0]+" Mass derivative", color="xkcd:burnt orange", linewidth=0.75)
    plt.xlabel('Temperature (°C)')
    ax2.set_ylim(bottom=0)

    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=2)
    fig.set_size_inches(10, 5)
    plt.title(str(os.path.basename(filename1)) + " plot corrected with char measurements EXPLANATORY GRAPH")

    plt.show()

    # ----------------Plotting Clean Figure----------------

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Temperature (°C)')
    ax1.set_ylabel('Heat flow (W/g)', color="xkcd:blue")
    ax1.spines['left'].set_color("xkcd:blue")
    ax1.tick_params(axis='y', labelcolor="xkcd:blue", color="xkcd:blue")
    # ax1.plot(heattemp, -hfcorr, label=experimentname[0]+" Heat Flow", linewidth=0.75, zorder = 10)
    ax1.plot(heattemp, -hfcorr, label="Heat Flow", linewidth=0.75, zorder = 10)
    ax1.plot(heattemp, -hfcorr+dev, color="black", linewidth=0.5, zorder = 9)
    ax1.plot(heattemp, -hfcorr-dev, color="black", linewidth=0.5, zorder = 9)
    ax1.plot(heattemp, zeros, color="black", linewidth=0.5, zorder = 11)
    
    ax1.fill_between(dfheattemp["tempaverage"], -hfcorr+dev, -hfcorr-dev, color="xkcd:powder blue", zorder = 5)
    # ax1.fill_between(heattempt, -hfcorr2, zeros2, color="pink", zorder = 0, label=experimentname[0]+" Reaction area") # While the integration has been made with the time x-axis, representation is much clearer with the temperature x-axis

    print(hr)
    if hr == 10:
        ax1.set_ylim(-3.3,3.3)
    elif hr == 2.5:
        ax1.set_ylim(-1,1)
    else:
        ax1.set_ylim(-2,2)

    ax1.set_xlim(100,550)

    ax2 = ax1.twinx()
    ax2.set_ylabel('dTG (mg/s)', color="xkcd:burnt orange")
    ax2.set_ylim(0, 0.0081)
    ax2.set_xlabel('Temperature (°C)')
    ax2.spines['left'].set_color("xkcd:blue")
    ax2.spines['right'].set_color("xkcd:burnt orange")
    ax2.tick_params(axis='y', labelcolor="xkcd:burnt orange", color="xkcd:burnt orange")
    # ax2.plot(masstemp, -deriv, label=experimentname[0]+" Mass derivative", color="xkcd:burnt orange", linewidth=0.75)
    ax2.plot(masstemp, -deriv, label="Mass derivative", color="xkcd:burnt orange", linewidth=0.75)

    plt.xlabel('Temperature (°C)')
    ax2.set_ylim(bottom=0)
    
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=4)
    fig.set_size_inches(6, 3)
    fig.tight_layout()
    plt.title(str(os.path.basename(filename1)) + " plot corrected with char measurements")

    plt.show()


#The rest of the code is the setup and positioning of the GUI

# Main window
window = Tk()

# Main window title
window.title('TGA-DSC-Tool')

# Main window size
window.geometry("200x330")

# Main window color
window.config(background="white")

# Background frames
top_frame = Frame(window, bg="black", width=200, height=50, pady=3)
top_frame2 = Frame(window, bg="alice blue", width=400, height=100, pady=20)
graphing_frame = Frame(window, bg="antiquewhite", width=400, height=50, pady=6)
low_frame = Frame(window, bg="gray84", width=400, height=50, pady=3)

# Background frame positioning
window.grid_rowconfigure(1, weight=1)
window.grid_columnconfigure(0, weight=1)

top_frame.grid(row=0, sticky="ew")
top_frame2.grid(row=1, sticky="ew")
graphing_frame.grid(row=2, sticky="ew")
low_frame.grid(row=3, sticky="ew")

# File Explorer labels
label_file_explorer = Label(top_frame,
                            text="TGA/DSC analysis tool",
                            width=30, height=4,
                            fg="blue")

# File Explorer button
button_explore1 = Button(top_frame2,
                         text="Load Excel file",
                         command=browseFiles1)

button_exit = Button(low_frame,
                     text="Exit",
                     fg="red",
                     command=exit)

# Graphing frame label
grf_label = Label(graphing_frame, 
                  width=30, height=2, 
                  text='Available graphs:')


# Heat flow corrected graph, mass detected from filename
button_hfcorr20 = Button(graphing_frame,
                       text="HF Qrad corr.",
                       command=hfcorr)


# Heat flow corrected with char file, mass detected from filename
button_hfuncorr20 = Button(graphing_frame,
                       text="HF char corr.",
                       command=hfcorr_char)

# Grid method is chosen for placing the widgets at respective positions
# in a table like structure by specifying rows and columns. Place method is
# also used because I´m a lazy person sometimes and it´s easier to place by hand.

label_file_explorer.grid(column=1, row=1)

button_explore1.grid(column=1, row=1, padx=70)

grf_label.grid(column=2, row=1, pady=10, padx=5)

button_hfcorr20.grid(column=2, row=4, pady=10)

button_hfuncorr20.grid(column=2, row=5, pady=10)

button_exit.grid(column=3, row=9, padx=5)

# Let the window wait for any events
window.mainloop()
