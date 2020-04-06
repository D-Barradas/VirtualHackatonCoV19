#!/usr/bin/env python
# coding: utf-8

# In[1]:


import requests
import json
import pandas as pd
url = "https://glyconnect.expasy.org/api/glycosylations"


# In[2]:


## send the correct params to query the api
params = {'taxonomy':'Severe acute respiratory syndrome coronavirus 2 (2019-nCoV)', 'protein': 'Recombinant Spike glycoprotein (HEK293) - DRAFT DATA'}

# Severe acute respiratory syndrome coronavirus2  (2019-nCoV)&protein=Recombinant Spike glycoprotein (HEK293)
response = requests.get(url ,params=params)


# In[3]:


my_response = response.json()
df_dump = pd.DataFrame()
for r in range(len(my_response['results'])):
    df_results_uniprots = pd.DataFrame(my_response['results'][r]['protein']['uniprots'],index=[r])
    df_results_site = pd.DataFrame(my_response['results'][r]['site'],index=[r])
    df_results_composition = pd.DataFrame(my_response['results'][r]["composition"],index=[r])
    df_temp = pd.concat([df_results_composition,df_results_site, df_results_uniprots], sort=False, axis=1)
    df_dump = pd.concat([df_dump,df_temp],sort=True,axis=0)
df_dump.drop_duplicates(inplace=True)


# In[4]:


df_dump.to_csv("../data/HEK293_glycosilations_COVID19.csv",index=False)


# In[5]:


## params for the second protein 
params = {'taxonomy':'Severe acute respiratory syndrome coronavirus 2 (2019-nCoV)', 'protein': "Recombinant Spike glycoprotein (BTI-Tn-5B1-4) - DRAFT DATA"}

response = requests.get(url ,params=params)
my_response = response.json()


# In[6]:


df_dump = pd.DataFrame()
for r in range(len(my_response['results'])):
    df_results_uniprots = pd.DataFrame(my_response['results'][r]['protein']['uniprots'],index=[r])
    df_results_site = pd.DataFrame(my_response['results'][r]['site'],index=[r])
    df_results_composition = pd.DataFrame(my_response['results'][r]["composition"],index=[r])
    df_temp = pd.concat([df_results_composition,df_results_site, df_results_uniprots], sort=False, axis=1)
    df_dump = pd.concat([df_dump,df_temp],sort=True,axis=0)
df_dump.drop_duplicates(inplace=True)


# In[7]:


df_dump.to_csv("../data/BTI-Tn-5B1-4_glycosilations_COVID19.csv",index=False)


# In[ ]:




