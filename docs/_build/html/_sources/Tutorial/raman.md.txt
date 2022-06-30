# Raman
Welcome to Raman. <br>
This package is created for NIIMBL4.1 Case2. <br>

Currently supported input files are:
* .spc
* .dx
* .csv


Currently supported a set of routines to execute spectral pre-processing like:<br>
* MSC
* SNV
* Detrend
* Savitzky - Golay
* Derivatives
* ..


```python
#Import basic libraries
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
```

# Read .spc file
## Read a single file


```python
spc = read_spc('VIAVI/JDSU_Phar_Rotate_S06_1_20171009_1540.spc')
spc.plot()
plt.xlabel("nm")
plt.ylabel("Abs")
plt.grid(True)
print(spc.head())
```

    gx-y(1)
    908.100000    0.123968
    914.294355    0.118613
    920.488710    0.113342
    926.683065    0.108641
    932.877419    0.098678
    dtype: float64
    

![image](https://user-images.githubusercontent.com/56982084/163854218-c5a8407e-2d36-4ab2-8a2c-13aef19c5341.png)


# Read .dx spectral files
Raman also allows to read the most common .dx file formats

## Read a single .dx file

```python
Foss_single = read_dx()
# Run read method
df = Foss_single.read(file='../../data/FOSS/FOSS.dx')
df.transpose().plot(legend=False)
```



![image](https://user-images.githubusercontent.com/56982084/163855034-65caff70-aec4-4fef-b022-9acf08699bbe.png)


```python
for c in Foss_single.Samples['29179'].keys():
    print(c)
```
    y
    Conc
    TITLE
    JCAMP_DX
    DATA TYPE
    CLASS
    DATE
    DATA PROCESSING
    XUNITS
    YUNITS
    XFACTOR
    YFACTOR
    FIRSTX
    LASTX
    MINY
    MAXY
    NPOINTS
    FIRSTY
    CONCENTRATIONS
    XYDATA
    X
    Y


# Read .csv files
Raman also allows to read .csv file

```python
df = pd.read_csv('data/RamanRaw.csv')
       
    X = df[df.columns[301:1702]]
    y = df[df.columns[-4:]]
```



# Spectra preprocessing
A set of built in classes to perform spectra pre-processing like: <br>
* MSC: Multiplicative scattering correction
* SNV: Standard normal variate
* Detrend
* n order derivative
* Savitzky golay smmothing

## MSC
```python
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
```
![image](https://user-images.githubusercontent.com/56982084/171699029-aa8b4746-bfcb-4247-b46c-9fc097735107.png)

![image](https://user-images.githubusercontent.com/56982084/171699048-d6d2243c-df1e-41e0-a98c-288108596dbb.png)

## SNV and Detrend

```python
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
```
![image](https://user-images.githubusercontent.com/56982084/171699288-78f47aa5-fc83-4f61-8f45-9963e89a6d2c.png)
![image](https://user-images.githubusercontent.com/56982084/171699300-667279eb-b849-4e97-a84b-30908a5237e5.png)
![image](https://user-images.githubusercontent.com/56982084/171699311-cdd63bb0-43b4-429a-8fd8-57a413ca28db.png)


## Savitzky - Golay
```python
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
```
![image](https://user-images.githubusercontent.com/56982084/171699461-bb22831b-eae5-4ee0-8455-58e1e65f19e2.png)
![image](https://user-images.githubusercontent.com/56982084/171699468-c040a4d3-1036-40b4-8639-2a9ddfadfdea.png)


# Split The Dataset
## split the dataset into training (70%) and testing (30%) sets

We now split the samples into a training set and a test set in order to estimate the test error from different methods.
In this part, we use the pre-processed raman data via SNV and Detrend as example. The concentration of glucose is target variable.

```python
X_train, X_test, y_all_train, y_all_test = train_test_split(X_detrend, y, test_size=0.3, random_state=0)

Metabolites = {0:'Glucose', 1: 'Lactate', 2: 'Glutamine', 3: 'NH4'}
#0:GLuc; 1: Lac; 2:Gln; 3: NH4;
index = 0
y_train, y_test =  y_all_train[y.columns[index]].to_frame(), y_all_test[y.columns[index]].to_frame()
```

# Raman Spectra Data Analysis


## Ridge Regression


We fit a ridge regression model on the training set. Without given any value of Regularization strength. The function will implement cross validation to find optimal value of it with the smallest validation error.

The output includings: 
1. Visualization of regression coefficients 
2. Regression results on the training data
3. Regression results on the test data


```python
rr = RR()
rr.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
```

![image](https://user-images.githubusercontent.com/56982084/171699875-c9d280a8-b8de-4b46-b68d-f9075ab17ea5.png)
![image](https://user-images.githubusercontent.com/56982084/171699890-4a911f36-ea0b-446a-983e-a359679c6046.png)
![image](https://user-images.githubusercontent.com/56982084/171699900-b358d97d-b353-455b-bf2d-cfbee1c67da8.png)



## Lasso Regression

We fit a Lasso regression model on the training set. Without given any value of Regularization strength. The function will implement cross validation to find optimal value of it with the smallest validation error.

The output includings: 
1. Visualization of regression coefficients 
2. Regression results on the training data
3. Regression results on the test data


```python
la = LA()
la.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
```

![image](https://user-images.githubusercontent.com/56982084/171700170-4778cca4-6c95-44b9-9049-df7f8193b5fe.png)
![image](https://user-images.githubusercontent.com/56982084/171700186-7e2483e8-bb68-4822-8d20-950ed3688cdf.png)
![image](https://user-images.githubusercontent.com/56982084/171700203-dd8e2a40-2502-4a7a-a0a2-fa4564bc3cd9.png)

We fit a Elastic Net model on the training set. Without given any value of Regularization strength. The function will implement cross validation to find optimal value of it with the smallest validation error.

The output includings: 
1. Visualization of regression coefficients 
2. Regression results on the training data
3. Regression results on the test data


## Elastic Net
```python
en = EN()
en.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
```
![image](https://user-images.githubusercontent.com/56982084/171700458-1a1085eb-feed-4d5f-b2d1-5480fe51ebe2.png)
![image](https://user-images.githubusercontent.com/56982084/171700466-75d99417-87ce-48d4-b1e0-c8f611c85d32.png)
![image](https://user-images.githubusercontent.com/56982084/171700478-999adb98-abaf-468c-893f-9fe2d8bc235c.png)





## Principal Components Regression

We fit Principal Components Regression on the training set. Without given any number of principal component to keep. The function will implement cross validation to find optimal value of it with the smallest validation error.

The output includings: 
1. Visualization of regression coefficients 
2. Regression results on the training data
3. Regression results on the test data


```python
pcr = PCR()
pcr.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
```
![image](https://user-images.githubusercontent.com/56982084/171700543-24cef100-b61c-48ef-b715-651426fb1b4f.png)
![image](https://user-images.githubusercontent.com/56982084/171700551-f124ca12-f5b9-4e1b-8ab6-a8fd2521b1b9.png)
![image](https://user-images.githubusercontent.com/56982084/171700560-5a2304df-a26b-4377-8657-22709af790c6.png)





## Partial least-squares (PLS) regression

We fit Partial least-squares (PLS) regression on the training set. Without given any number of principal component to keep. The function will implement cross validation to find optimal value of it with the smallest validation error.

The output includings: 
1. Visualization of regression coefficients 
2. Regression results on the training data
3. Regression results on the test data


```python
pls = PLS()
pls.analysis(X_train, X_test, y_train, y_test, Metabolites[index])
```
![image](https://user-images.githubusercontent.com/56982084/171700723-4786f9f5-85f4-49d6-9280-8845aff8fc38.png)
![image](https://user-images.githubusercontent.com/56982084/171700734-a312b931-dc3f-4029-a67b-5082ff67bad6.png)
![image](https://user-images.githubusercontent.com/56982084/171700746-16a423cc-a0e6-4d33-8ee6-57634aee433d.png)
