import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

from src.vLab.RamanAnalytics.Preprocessing import msc, detrend, sav_gol, snv
from src.vLab.RamanAnalytics.Analysis import RR, LA, EN, PCR, PLS

if __name__ == '__main__':

    ###############################################
    ######### Data Import from .csv file###########
    ###############################################
    df = pd.read_csv('data/RamanRaw.csv')
       
    X = df[df.columns[301:1702]]
    y = df[df.columns[-4:]]
    # pd.to_numeric(X.columns)
       
    ########################################
    ########## Preprocessing ###############
    ########################################
    # MSC
    MSC = msc()
    MSC.fit(X)
    X_msc = MSC.transform(X)
    
    plt.figure()
    X.transpose().plot(legend=False, title = "Raw spectra")
    plt.ylim(bottom = 0, top = 2.5 * 10**6)
    plt.xlabel('RamanShift(cm-1)')
    plt.ylabel('Counts')
    plt.show()
    
        
    plt.figure()
    X_msc.transpose().plot(legend=False, title = "MSC spectra")
    plt.ylim(bottom = 0, top = 2.5 * 10**6)
    plt.xlabel('RamanShift(cm-1)')
    plt.ylabel('Counts')
    plt.show()

    # SNV and Detrend
    SNV = snv()
    X_snv = SNV.fit_transform(X)

    Detr = detrend()
    X_detrend = Detr.fit_transform(spc=X_snv, wave=np.array(X_snv.columns.astype('float64')), deg = 1)
    
    
    plt.figure()
    X.transpose().plot(legend=False, title = "Raw spectra")
    plt.ylim(bottom = 0, top = 2.5 * 10**6)
    plt.xlabel('RamanShift(cm-1)')
    plt.ylabel('Counts')
    plt.show()
    
        
    plt.figure()
    X_snv.transpose().plot(legend=False, title = "SNV spectra")
    #plt.ylim(bottom = 0, top = 2.5 * 10**6)
    plt.xlabel('RamanShift(cm-1)')
    plt.ylabel('Counts')
    plt.show()

    plt.figure()
    X_detrend.transpose().plot(legend=False, title = "SNV+ Detrend spectra")
    #plt.ylim(bottom = 0, top = 2.5 * 10**6)
    plt.xlabel('RamanShift(cm-1)')
    plt.ylabel('Counts')
    plt.show()


    # Savitzky - Golay
    SAV_gol = sav_gol()
    X_SAV = SAV_gol.transform(X, deriv=2)

    plt.figure()
    X.transpose().plot(legend=False, title = "Raw spectra")
    plt.ylim(bottom = 0, top = 2.5 * 10**6)
    plt.xlabel('RamanShift(cm-1)')
    plt.ylabel('Counts')
    plt.show()
    
        
    plt.figure()
    X_snv.transpose().plot(legend=False, title = "Savitzky golay smmothing spectra")
    #plt.ylim(bottom = 0, top = 2.5 * 10**6)
    plt.xlabel('RamanShift(cm-1)')
    plt.ylabel('Counts')
    plt.show()

    ##################################
    #########Split The Dataset########
    ##################################
    
    # split the dataset into training (70%) and testing (30%) sets
    X_train, X_test, y_all_train, y_all_test = train_test_split(X_detrend, y, test_size=0.3, random_state=0)

    Metabolites = {0:'Glucose', 1: 'Lactate', 2: 'Glutamine', 3: 'NH4'}
    #0:GLuc; 1: Lac; 2:Gln; 3: NH4;
    index = 0
    y_train, y_test =  y_all_train[y.columns[index]].to_frame(), y_all_test[y.columns[index]].to_frame()

    ##################################
    ####Raman Spectra Data Analysis###
    ##################################
    
    
    ################Ridge Regression ################## 
    rr = RR()
    rr.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
    
    
    ################Lasso ##################
    la = LA()
    la.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
    
    
    ################Elastic Net ##################
    en = EN()
    en.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
    
        
    ######## Principal Components Regression ###########
    pcr = PCR()
    pcr.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
    
    
    ################PLS ##################
    pls = PLS()
    pls.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
