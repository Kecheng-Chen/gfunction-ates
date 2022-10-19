import joblib
import numpy as np
import pandas as pd

def myfunc2(row):
	rhoc_f=1000*4000
	alpha = row['Q']*rhoc_f/(4*np.pi*row['H']*row['hea_cond'])
	D=row['hea_cond']/(row['hea_qua']*1e6)
	t_star = row['Year']*360*24*60*60/(row['H']**2/(4*D))
	dl_star = row['H']/((D*row['Year']*360*24*60*60)**0.5)
	hcap_star = row['hea_qua']/2.05
	hcond_star= row['hea_cond']/3.77
	return pd.Series([alpha,t_star,dl_star,hcap_star,hcond_star])

def result(X,name):
	knn_from_joblib = joblib.load(name)
	df = pd.DataFrame(columns=['bound','mag_a','Q','hea_qua','hea_cond','H','Year'])
	df.loc[0] = [X[0],X[1],X[2],X[3],X[4],X[5],X[6]]
	df[['alpha','t_star','dl_star','hcap_star','hcond_star']] = df.apply(myfunc2,axis=1)
	gg=df[['bound','mag_a','alpha','t_star','dl_star','hcap_star','hcond_star']].values[0]
	return (knn_from_joblib(gg[0],gg[1],gg[2],gg[3],gg[4],gg[5],gg[6]))*(1/df['H'].values[0]*(2*df['hea_cond'].values[0]))

def classify(X):
	knn_from_joblib = joblib.load('model/classify.pkl')  
	return knn_from_joblib.predict([[X[0],X[1],X[2],X[3],X[4],X[5]]])[0]

def extract1(X):
	return result(X,'model/pred31.pkl')

def extract2(X):
	return result(X,'model/pred32.pkl')

def inject10(X):
	return result(X,'model/pred11.pkl')

def inject11(X):
	return result(X,'model/pred21.pkl')

def inject20(X):
	return result(X,'model/pred12.pkl')

def inject21(X):
	return result(X,'model/pred22.pkl')

e1=extract1(X)
e2=extract2(X)
i10=inject10(X)
i11=inject11(X)
i20=inject20(X)
i21=inject21(X)
c=classify(X)