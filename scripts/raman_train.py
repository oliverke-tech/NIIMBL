import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

from src.vLab.RamanAnalytics.Preprocessing import msc, detrend, sav_gol, snv
from sklearn.linear_model import LinearRegression
import pickle


df_all = pd.read_csv('data/Ambr250.csv')

# df = df_all.iloc[0:75]
df = df_all
concen = df[df.columns[-5:-1]]
raman = df[df.columns[312:-430]]


# SNV and Detrend
SNV = snv()
raman_snv = SNV.fit_transform(raman)

Detr = detrend()
raman_detrend = Detr.fit_transform(spc=raman_snv, wave=np.array(raman_snv.columns.astype('float64')), deg = 1)


plt.figure()
raman.transpose().plot(legend=False, title = "Raw spectra")
plt.xlabel('RamanShift(cm-1)')
plt.ylabel('Counts')
plt.show()

    
plt.figure()
raman_snv.transpose().plot(legend=False, title = "SNV spectra")
plt.xlabel('RamanShift(cm-1)')
plt.ylabel('Counts')
plt.show()

plt.figure()
raman_detrend.transpose().plot(legend=False, title = "SNV+ Detrend spectra")
plt.xlabel('RamanShift(cm-1)')
plt.ylabel('Counts')
plt.show()


y = raman_detrend-raman_detrend.iloc[0]
y = y.iloc[1:]

X = concen-concen.iloc[0]
X = X.iloc[1:]

# X = pd.concat([X, -X.sum(axis = 1)],axis=1)
# X.columns = ['Glucose', 'Lactate', 'Glutamine', 'Ammonia', 'Glutamate', 'Others']
# X.columns = ['Glucose', 'Lactate', 'Glutamine', 'Ammonia', 'Others']


# split the dataset into training (70%) and testing (30%) sets
X_train_diff, X_test_diff, y_train_diff, y_test_diff = train_test_split(X, y, test_size=0.3, random_state=0)


regr = LinearRegression()
regr.fit(X_train_diff, y_train_diff)
# coef = regr.coef_

# regr.score(X_train_diff, y_train_diff)

predicted_test_diff = pd.DataFrame(regr.predict(X_test_diff))
predicted_test_diff.columns = raman.columns

predicted_train_diff = pd.DataFrame(regr.predict(X_train_diff))
predicted_train_diff.columns = raman.columns


y_train = y_train_diff + raman_detrend.iloc[0]
y_test = y_test_diff + raman_detrend.iloc[0]

predicted_train = predicted_train_diff + raman_detrend.iloc[0]
predicted_test = predicted_test_diff + raman_detrend.iloc[0]


#True spectra vs simulation spectra in every 10 records at train set 
# for i in range(10, 200,10):
#     fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1)
#     fig.tight_layout(h_pad=3)
#     np.transpose(y_train).iloc[:,(i-10):i].plot(legend=False, title = "True_Train", ax=ax1)
#     ax1.set_ylim(bottom = -0.3, top = 0.6)
# # plt.xlim(left = -100, right = 2800)
#     ax1.set_xlabel('RamanShift(cm-1)')
#     ax1.set_ylabel('Counts')
#     np.transpose(predicted_train).iloc[:,(i-10):i].plot(legend=False, title = "Pred_Train", ax=ax2)
# # plt.plot(,np.array(predicted_test))
#     ax2.set_ylim(bottom = -0.3, top = 0.6)
# # plt.xlim(left = -100, right = 2800)
#     ax2.set_xlabel('RamanShift(cm-1)')
#     ax2.set_ylabel('Counts')
#     plt.show()


#True spectra vs simulation spectra in every 10 records at test set 
# for i in range(10, 200,10):
#     fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1)
#     fig.tight_layout(h_pad=3)
#     np.transpose(y_test).iloc[:,(i-10):i].plot(legend=False, title = "True_test", ax=ax1)
#     ax1.set_ylim(bottom = -0.3, top = 0.6)
# # plt.xlim(left = -100, right = 2800)
#     ax1.set_xlabel('RamanShift(cm-1)')
#     ax1.set_ylabel('Counts')
#     np.transpose(predicted_test).iloc[:,(i-10):i].plot(legend=False, title = "Pred_test", ax=ax2)
# # plt.plot(,np.array(predicted_test))
#     ax2.set_ylim(bottom = -0.3, top = 0.6)
# # plt.xlim(left = -100, right = 2800)
#     ax2.set_xlabel('RamanShift(cm-1)')
#     ax2.set_ylabel('Counts')
#     plt.show()
    
    
#True spectra vs simulation spectra at train set 
fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1)
fig.tight_layout(h_pad=3)
np.transpose(y_train).plot(legend=False, title = "True_Train", ax=ax1)
# ax1.set_ylim(bottom = -0.3, top = 1)
# plt.xlim(left = -100, right = 2800)
ax1.set_xlabel('RamanShift(cm-1)')
ax1.set_ylabel('Counts')
np.transpose(predicted_train).plot(legend=False, title = "Pred_Train", ax=ax2)
# plt.plot(,np.array(predicted_test))
# ax2.set_ylim(bottom = -0.3, top = 1)
# plt.xlim(left = -100, right = 2800)
ax2.set_xlabel('RamanShift(cm-1)')
ax2.set_ylabel('Counts')
plt.show()


#True spectra vs simulation spectra at train set 
fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1)
fig.tight_layout(h_pad=3)
np.transpose(y_test).plot(legend=False, title = "True_test", ax=ax1)
# ax1.set_ylim(bottom = -0.3, top = 0.6)
# plt.xlim(left = -100, right = 2800)
ax1.set_xlabel('RamanShift(cm-1)')
ax1.set_ylabel('Counts')
np.transpose(predicted_test).plot(legend=False, title = "Pred_test", ax=ax2)
# plt.plot(,np.array(predicted_test))
# ax2.set_ylim(bottom = -0.3, top = 0.6)
# plt.xlim(left = -100, right = 2800)
ax2.set_xlabel('RamanShift(cm-1)')
ax2.set_ylabel('Counts')
plt.show()

# save the model
filename = 'src/vLab/RamanSimulator/simulator_wo_Glu.sav'
pickle.dump(regr, open(filename, 'wb'))


# save the model
# filename = 'src/vLab/RamanSimulator/simulator_wi_Glu.sav'
# pickle.dump(regr, open(filename, 'wb'))

