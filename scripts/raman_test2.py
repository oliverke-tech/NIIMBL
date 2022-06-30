import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from sklearn import model_selection
from sklearn.model_selection import RepeatedKFold
from sklearn.model_selection import train_test_split
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import mean_squared_error

from src.vLab.RamanAnalytics.ReadSpc import read_spc
from src.vLab.RamanAnalytics.ReadDx import ReadDx
from src.vLab.RamanAnalytics.Preprocessing import msc, detrend, sav_gol, snv
from sklearn.linear_model import Ridge, RidgeCV, Lasso, LassoCV, ElasticNet, ElasticNetCV, LinearRegression
from src.vLab.RamanAnalytics.Analysis import RR, LA, EN, PCR, PLS

if __name__ == '__main__':

    ######### Data Import from .spc file###########
    spc = read_spc('data/VIAVI/JDSU_Phar_Rotate_S06_1_20171009_1540.spc')
    spc.plot()
    plt.xlabel("nm")
    plt.ylabel("Abs")
    plt.grid(True)
    print(spc.head())

    # Instantiate an object
    Foss_single = ReadDx()
    # Run  read method
    df = Foss_single.read(file='data/FOSS/FOSS.dx')
    df.transpose().plot(legend=False)

    Y = Foss_single.Samples['29179']['X']

    ########## Preprocessing ################
    # MSC
    MSC = msc()
    MSC.fit(df)
    df_msc = MSC.transform(df)

    f, ax = plt.subplots(2, 1, figsize=(14, 8))
    ax[0].plot(df.transpose())
    ax[0].set_title("Raw spectra")

    ax[1].plot(df_msc.transpose())
    ax[1].set_title("MSC spectra")
    plt.show()

    # SNV and Detrend
    SNV = snv()
    df_snv = SNV.fit_transform(df)

    Detr = detrend()
    df_detrend = Detr.fit_transform(spc=df_snv, wave=np.array(df_snv.columns))

    f, ax = plt.subplots(3, 1, figsize=(18, 8))
    ax[0].plot(df.transpose())
    ax[0].set_title("Raw spectra")

    ax[1].plot(df_snv.transpose())
    ax[1].set_title("SNV spectra")

    ax[2].plot(df_detrend.transpose())
    ax[2].set_title("SNV+ Detrend spectra")

    plt.tight_layout()
    plt.show()

    # Savitzky - Golay
    SAV_gol = sav_gol()
    df_SAV = SAV_gol.transform(df, deriv=2)

    f, ax = plt.subplots(2, 1, figsize=(14, 8))
    ax[0].plot(df.transpose())
    ax[0].set_title("Raw spectra")

    ax[1].plot(df_SAV.transpose())
    ax[1].set_title("Savitzky golay smmothing spectra")
    plt.show()

    ######## PCA unsupervised transformation ###########
    pca = PCA()
    pca.fit(df_msc)
    plt.figure(figsize=(18, 8))
    plt.plot(range(1, len(pca.explained_variance_) + 1),
             100 * pca.explained_variance_.cumsum() / pca.explained_variance_.sum())
    plt.grid(True)
    plt.xlabel("Number of components")
    plt.ylabel(" cumulative % of explained variance")

    df_pca = pd.DataFrame(pca.transform(df_msc))
    plt.figure(figsize=(18, 8))
    plt.plot(df_pca.loc[:, 0:25].transpose())

    plt.title("Transformed spectra PCA")
    plt.ylabel("Response feature")
    plt.xlabel("Principal component")
    plt.grid(True)
    plt.show()
    
    
    

    ######### Data Import from .csv file###########
    df = pd.read_csv('data/RamanRaw.csv')
       
    X = df[df.columns[301:1702]]
    y = df[df.columns[-4:]]
    pd.to_numeric(X.columns)
       
    ########## Preprocessing ################
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


    # split the dataset into training (70%) and testing (30%) sets
    #X_train, X_test, y_train, y_test = train_test_split(X, y[y.columns[0]].to_frame(), test_size=0.3, random_state=0)
    X_train, X_test, y_all_train, y_all_test = train_test_split(X_detrend, y, test_size=0.3, random_state=0)

    Metabolites = {0:'GLucose', 1: 'Lactate', 2: 'Glutamine', 3: 'NH4'}
    #0:GLuc; 1: Lac; 2:Gln; 3: NH4;
    index = 3
    y_train, y_test =  y_all_train[y.columns[index]].to_frame(), y_all_test[y.columns[index]].to_frame()

    rr = RR()
    rr.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
    
    la = LA()
    la.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
    
    en = EN()
    en.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
    
    pcr = PCR()
    pcr.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
    
    pls = PLS()
    pls.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
    ################Ridge Regression ##################   
    #Cross-validation
    ridgecv = RidgeCV(scoring = 'neg_mean_squared_error')
    ridgecv.fit(X_train, y_train)
    #ridgecv.alpha_

    
    ridge = Ridge(alpha = ridgecv.alpha_)
    ridge.fit(X_train, y_train)
    mean_squared_error(y_test, ridge.predict(X_test))

    ridge_coef = ridge.coef_
    
    plt.figure()
    g = plt.scatter(pd.to_numeric(X_detrend.columns), ridge_coef)
    g.axes.set_title('Regression Coefficients (RR)')
    g.axes.set_xlabel('Raman Shift (cm-1)')
    g.axes.set_ylabel('Coefficient')
    plt.show()
    

    predicted_RR_train = ridge.predict(X_train)

    plt.figure()
    g = plt.scatter(y_train, predicted_RR_train)
    g.axes.set_title(str(Metabolites[index]) + ' Concentration Train (RR)')
    g.axes.set_xlabel('True Values (g/L)')
    g.axes.set_ylabel('Predictions (g/L)')
    g.axes.axis('equal')
    g.axes.axis('square')
    g.axes.axline([0, 0], [1, 1], color='r')
    #plt.plot(predicted_value_train,predicted_value_train, 'r-' )
    plt.show()
    
    predicted_RR_test = ridge.predict(X_test)
    
    plt.figure()
    g = plt.scatter(y_test, predicted_RR_test)
    g.axes.set_title(str(Metabolites[index]) + ' Concentration Test (RR)')
    g.axes.set_xlabel('True Values (g/L)')
    g.axes.set_ylabel('Predictions (g/L)')
    g.axes.axis('equal')
    g.axes.axis('square')
    g.axes.axline([0, 0], [1, 1], color='r')
    #plt.plot([0,1],[0,np.max(g.axes.get_xlim(),g.axes.get_ylim()),predicted_value_test)], 'r-' )
    plt.show()
    
    
    ################Lasso ##################
    #Cross-validation
    lassoCV = LassoCV(max_iter = 10000)
    lassoCV.fit(X_train, y_train)

    lasso = Lasso(alpha = lassoCV.alpha_ , normalize = True)
    lasso.fit(X_train, y_train)
    mean_squared_error(y_test, lasso.predict(X_test))
    
    lasso_coef = lasso.coef_
    
    plt.figure()
    g = plt.scatter(pd.to_numeric(X_detrend.columns), lasso_coef)
    g.axes.set_title('Regression Coefficients (LA)')
    g.axes.set_xlabel('Raman Shift (cm-1)')
    g.axes.set_ylabel('Coefficient')
    plt.show()
    
        
    predicted_LA_train = lasso.predict(X_train)

    plt.figure()
    g = plt.scatter(y_train, predicted_LA_train)
    g.axes.set_title(str(Metabolites[index]) + ' Concentration Train (LA)')
    g.axes.set_xlabel('True Values (g/L)')
    g.axes.set_ylabel('Predictions (g/L)')
    g.axes.axis('equal')
    g.axes.axis('square')
    g.axes.axline([0, 0], [1, 1], color='r')
    plt.show()
    
    
    predicted_LA_test = lasso.predict(X_test)
    
    
    plt.figure()
    g = plt.scatter(y_test, predicted_LA_test)
    g.axes.set_title(str(Metabolites[index]) + ' Concentration Test (LA)')
    g.axes.set_xlabel('True Values (g/L)')
    g.axes.set_ylabel('Predictions (g/L)')
    g.axes.axis('equal')
    g.axes.axis('square')
    g.axes.axline([0, 0], [1, 1], color='r')
    plt.show()
    
    
    ################Elastic Net ##################
    #Cross-validation
    elasticCV = ElasticNetCV(max_iter = 10000)
    elasticCV.fit(X_train, y_train)

    elastic = ElasticNet(alpha = elasticCV.alpha_, normalize = True)
    elastic.fit(X_train, y_train)
    mean_squared_error(y_test, elastic.predict(X_test))

    elastic_coef = elastic.coef_
    
    plt.figure()
    g = plt.scatter(pd.to_numeric(X_detrend.columns), elastic_coef)
    g.axes.set_title('Regression Coefficients (EN)')
    g.axes.set_xlabel('Raman Shift (cm-1)')
    g.axes.set_ylabel('Coefficient')
    plt.show()

    predicted_EN_train = elastic.predict(X_train)

    plt.figure()
    g = plt.scatter(y_train, predicted_EN_train)
    g.axes.set_title(str(Metabolites[index]) + ' Concentration Train (EN)')
    g.axes.set_xlabel('True Values (g/L)')
    g.axes.set_ylabel('Predictions (g/L)')
    g.axes.axis('equal')
    g.axes.axis('square')
    g.axes.axline([0, 0], [1, 1], color='r')
    plt.show()
        
    
    predicted_EN_test = elastic.predict(X_test)
    
    plt.figure()
    g = plt.scatter(y_test, predicted_EN_test)
    g.axes.set_title(str(Metabolites[index]) + ' Concentration Test (EN)')
    g.axes.set_xlabel('True Values (g/L)')
    g.axes.set_ylabel('Predictions (g/L)')
    g.axes.axis('equal')
    g.axes.axis('square')
    g.axes.axline([0, 0], [1, 1], color='r')
    plt.show()
    

    ######## Principal Components Regression ###########
    # pca = PCA()
    # pca.fit(X_train)
    # plt.figure(figsize=(18, 8))
    # plt.plot(range(1, len(pca.explained_variance_) + 1),
    #          100 * pca.explained_variance_.cumsum() / pca.explained_variance_.sum())
    # plt.grid(True)
    # plt.xlabel("Number of components")
    # plt.ylabel(" cumulative % of explained variance")

    # df_pca = pd.DataFrame(pca.transform(df_msc))
    # plt.figure(figsize=(18, 8))
    # plt.plot(df_pca.loc[:, 0:25].transpose())

    # plt.title("Transformed spectra PCA")
    # plt.ylabel("Response feature")
    # plt.xlabel("Principal component")
    # plt.grid(True)
    # plt.show()    
        
    pca = PCA()

    # Scale the data
    X_reduced_train = pca.fit_transform(scale(X_train))
    n = len(X_reduced_train)

    # 10-fold CV, with shuffle
    kf_10 = model_selection.KFold(n_splits=10, shuffle=True, random_state=1)

    mse = []
    regr = LinearRegression()
    # Calculate MSE with only the intercept (no principal components in regression)
    score = -1*model_selection.cross_val_score(regr, np.ones((n,1)), y_train, cv=kf_10, scoring='neg_mean_squared_error').mean()    
    mse.append(score)

    # Calculate MSE using CV for the 20 principle components, adding one component at the time.
    for i in np.arange(1, 21):
        score = -1*model_selection.cross_val_score(regr, X_reduced_train[:,:i], y_train, cv=kf_10, scoring='neg_mean_squared_error').mean()
        mse.append(score)
        
    plt.plot(np.array(mse), '-v')
    plt.xlabel('Number of principal components in regression')
    plt.ylabel('MSE')
    plt.title(str(Metabolites[index]))
    plt.xlim(xmin=-1);

    pca_component = mse.index(np.min(mse))+1

    # Train regression model on training data 
    regr.fit(X_reduced_train[:,:(pca_component)], y_train)
    
    
    pcr_coef = np.matmul(pca.components_.T[:,0:pca_component], regr.coef_.T)
    
    plt.figure()
    g = plt.scatter(pd.to_numeric(X_detrend.columns), pcr_coef)
    g.axes.set_title('Regression Coefficients (PCR)')
    g.axes.set_xlabel('Raman Shift (cm-1)')
    g.axes.set_ylabel('Coefficient')
    plt.show()

    # test = np.matmul(scale(X_train), np.matmul(pca.components_.T[:,0:PCA_component], regr.coef_.T)) + np.mean(y_train)[0]
    # np.matmul(scale(X_train), pca.components_.T[:,0:PCA_component])
    

    predicted_PCR_train = regr.predict(X_reduced_train[:,:(pca_component)])
    
    plt.figure()
    g = plt.scatter(y_train, predicted_PCR_train)
    g.axes.set_title(str(Metabolites[index]) + ' Concentration Train (PCR)')
    g.axes.set_xlabel('True Values (g/L)')
    g.axes.set_ylabel('Predictions (g/L)')
    g.axes.axis('equal')
    g.axes.axis('square')
    g.axes.axline([0, 0], [1, 1], color='r')
    plt.show()
     

    # Prediction with test data
    X_reduced_test = pca.transform(scale(X_test))[:,:(pca_component)]
    predicted_PCR_test = regr.predict(X_reduced_test)
    mean_squared_error(y_test, predicted_PCR_test)
    
    
    plt.figure()
    g = plt.scatter(y_test, predicted_PCR_test)
    g.axes.set_title(str(Metabolites[index]) + ' Concentration Test (PCR)')
    g.axes.set_xlabel('True Values (g/L)')
    g.axes.set_ylabel('Predictions (g/L)')
    g.axes.axis('equal')
    g.axes.axis('square')
    g.axes.axline([0, 0], [1, 1], color='r')
    plt.show()
    
    
    
    ################PLS ##################
    n = len(X_train)

    # 10-fold CV, with shuffle
    kf_10 = model_selection.KFold(n_splits=10, shuffle=True, random_state=1)

    mse = []

    for i in np.arange(1, 20):
        pls = PLSRegression(n_components=i)
        score = -1*model_selection.cross_val_score(pls, scale(X_train), y_train, cv=kf_10, scoring='neg_mean_squared_error').mean()
        mse.append(score)

    # Plot results
    plt.plot(np.arange(1, 20), np.array(mse), '-v')
    plt.xlabel('Number of principal components in regression')
    plt.ylabel('MSE')
    plt.title(str(Metabolites[index]))
    plt.xlim(xmin=-1)


    num_component = mse.index(np.min(mse))+1
    pls = PLSRegression(n_components=num_component)
    pls.fit(scale(X_train), y_train)
    
    
    pls_coef = pls.coef_
    
    plt.figure()
    g = plt.scatter(pd.to_numeric(X_detrend.columns), pls_coef)
    g.axes.set_title('Regression Coefficients (PLS)')
    g.axes.set_xlabel('Raman Shift (cm-1)')
    g.axes.set_ylabel('Coefficient')
    plt.show()
    

    predicted_PLS_train = pls.predict(scale(X_train))
    
    plt.figure()
    g = plt.scatter(y_train, predicted_PLS_train)
    g.axes.set_title(str(Metabolites[index]) + ' Concentration Train (PLS)')
    g.axes.set_xlabel('True Values (g/L)')
    g.axes.set_ylabel('Predictions (g/L)')
    g.axes.axis('equal')
    g.axes.axis('square')
    g.axes.axline([0, 0], [1, 1], color='r')
    plt.show()
     
   
    
    # Prediction with test data
    predicted_PLS_test = pls.predict(scale(X_test))
    mean_squared_error(y_test, pls.predict(scale(X_test)))
    
    plt.figure()
    g = plt.scatter(y_test, predicted_PLS_test)
    g.axes.set_title(str(Metabolites[index]) + ' Concentration Test (PLS)')
    g.axes.set_xlabel('True Values (g/L)')
    g.axes.set_ylabel('Predictions (g/L)')
    g.axes.axis('equal')
    g.axes.axis('square')
    g.axes.axline([0, 0], [1, 1], color='r')
    plt.show()
    