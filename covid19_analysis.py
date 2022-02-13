#!/usr/bin/env python3
# coding: utf-8

# In[47]:


import pandas as pd
import sqlite3
import matplotlib.pyplot as plt
import datetime as dt


# In[73]:


with sqlite3.connect('covid.db') as conn:
    euro = (
        "SELECT Date_reported, SUM(New_cases) as new_EURO "
        "FROM data WHERE WHO_region='EURO' "
        "GROUP BY Date_reported"
    )
    
    emro = (
        "SELECT Date_reported, SUM(New_cases) as new_EMRO "
        "FROM data WHERE WHO_region='EMRO' "
        "GROUP BY Date_reported"
    )
    
    afro = (
        "SELECT Date_reported, SUM(New_cases) as new_AFRO "
        "FROM data WHERE WHO_region='AFRO' "
        "GROUP BY Date_reported"
    )
    
    amro = (
        "SELECT Date_reported, SUM(New_cases) as new_AMRO "
        "FROM data WHERE WHO_region='AMRO' "
        "GROUP BY Date_reported"
    )
    
    searo = (
        "SELECT Date_reported, SUM(New_cases) as new_SEARO "
        "FROM data WHERE WHO_region='SEARO' "
        "GROUP BY Date_reported"
    )
    
    wpro = (
        "SELECT Date_reported, SUM(New_cases) as new_WPRO "
        "FROM data WHERE WHO_region='WPRO' "
        "GROUP BY Date_reported"
    )
    
dfeuro = pd.read_sql_query(euro, conn, parse_dates=["Date_reported"])
dfemro = pd.read_sql_query(emro, conn, parse_dates=["Date_reported"])
dfafro = pd.read_sql_query(afro, conn, parse_dates=["Date_reported"])
dfamro = pd.read_sql_query(amro, conn, parse_dates=["Date_reported"])
dfsearo = pd.read_sql_query(searo, conn, parse_dates=["Date_reported"])
dfwpro = pd.read_sql_query(wpro, conn, parse_dates=["Date_reported"])
data = pd.concat([dfeuro,dfemro['new_EMRO'],dfafro['new_AFRO'],dfamro['new_AMRO'],dfsearo['new_SEARO'],
                  dfwpro['new_WPRO']],axis=1)
data.index = data['Date_reported']
data.index.names=['date']
data = data.drop(['Date_reported'],axis=1)
data


# In[78]:


datemin = dt.datetime(2021,1,21)
datemax = dt.datetime(2022,1,31);
date_range = [datemin, datemax]
fig, ax = plt.subplots(2,3)
ax[0,0].bar(data.index.values,dfeuro['new_EURO'])
ax[0,0].plot(data.index.values,data['new_EURO'].rolling(window=7).mean(),color='red',label='7-day MA')
ax[0,0].set_xlim(date_range)
#ax[0,0].set_xlabel('Date Reported')
ax[0,0].set_ylabel('New case count')
ax[0,0].set_ylim(0,1.75e6)
ax[0,0].set_title('EURO')
ax[0,0].legend(loc='upper left')

ax[0,1].bar(data.index.values,dfemro['new_EMRO'])
ax[0,1].plot(data.index.values,data['new_EMRO'].rolling(window=7).mean(),color='red',label='7-day MA')
ax[0,1].set_xlim(date_range)
ax[0,1].set_title('EMRO')
ax[0,1].legend(loc='upper left')
#ax[0,1].xaxis_date()
#ax[0,1].set_xlabel('Date Reported')
#ax[0,1].set_ylabel('New case count')

ax[0,2].bar(dfafro['Date_reported'],dfafro['new_AFRO'])
ax[0,2].plot(data.index.values,data['new_AFRO'].rolling(window=7).mean(),color='red',label='7-day MA')
ax[0,2].set_xlim(date_range)
ax[0,2].set_title('AFRO')
ax[0,2].legend(loc='upper left')
#ax[0,2].set_xlabel('Date Reported')
#ax[0,2].set_ylabel('New case count')

ax[1,0].bar(dfamro['Date_reported'],dfamro['new_AMRO'])
ax[1,0].plot(data.index.values,data['new_AMRO'].rolling(window=7).mean(),color='red',label='7-day MA')
ax[1,0].set_xlim(date_range)
ax[1,0].set_xlabel('Date Reported')
ax[1,0].set_ylabel('New case count')
ax[1,0].set_title('AMRO')
ax[1,0].legend(loc='upper left')

ax[1,1].bar(dfsearo['Date_reported'],dfsearo['new_SEARO'])
ax[1,1].plot(data.index.values,data['new_SEARO'].rolling(window=7).mean(),color='red',label='7-day MA')
ax[1,1].set_xlim(date_range)
ax[1,1].set_xlabel('Date Reported')
ax[1,1].set_title('SEARO')
ax[1,1].legend(loc='upper left')
#ax[1,1].set_ylabel('New case count')

ax[1,2].bar(dfwpro['Date_reported'],dfwpro['new_WPRO'])
ax[1,2].plot(data.index.values,data['new_WPRO'].rolling(window=7).mean(),color='red',label='7-day MA')
ax[1,2].set_xlim(date_range)
ax[1,2].set_xlabel('Date Reported')
ax[1,2].set_title('WPRO')
ax[1,2].legend(loc='upper left')
#ax[1,2].set_ylabel('New case count')
plt.gcf().set_size_inches(30,15)
plt.savefig('worldcases.png',dpi=600)


# In[58]:


with sqlite3.connect('covid.db') as conn:
    euro = (
        "SELECT Date_reported, SUM(New_deaths) as new_EURO "
        "FROM data WHERE WHO_region='EURO' "
        "GROUP BY Date_reported"
    )
    
    emro = (
        "SELECT Date_reported, SUM(New_deaths) as new_EMRO "
        "FROM data WHERE WHO_region='EMRO' "
        "GROUP BY Date_reported"
    )
    
    afro = (
        "SELECT Date_reported, SUM(New_deaths) as new_AFRO "
        "FROM data WHERE WHO_region='AFRO' "
        "GROUP BY Date_reported"
    )
    
    amro = (
        "SELECT Date_reported, SUM(New_deaths) as new_AMRO "
        "FROM data WHERE WHO_region='AMRO' "
        "GROUP BY Date_reported"
    )
    
    searo = (
        "SELECT Date_reported, SUM(New_deaths) as new_SEARO "
        "FROM data WHERE WHO_region='SEARO' "
        "GROUP BY Date_reported"
    )
    
    wpro = (
        "SELECT Date_reported, SUM(New_deaths) as new_WPRO "
        "FROM data WHERE WHO_region='WPRO' "
        "GROUP BY Date_reported"
    )
    
dfeuro = pd.read_sql_query(euro, conn, parse_dates=["Date_reported"])
dfemro = pd.read_sql_query(emro, conn, parse_dates=["Date_reported"])
dfafro = pd.read_sql_query(afro, conn, parse_dates=["Date_reported"])
dfamro = pd.read_sql_query(amro, conn, parse_dates=["Date_reported"])
dfsearo = pd.read_sql_query(searo, conn, parse_dates=["Date_reported"])
dfwpro = pd.read_sql_query(wpro, conn, parse_dates=["Date_reported"])
data = pd.concat([dfeuro,dfemro['new_EMRO'],dfafro['new_AFRO'],dfamro['new_AMRO'],dfsearo['new_SEARO'],
                  dfwpro['new_WPRO']],axis=1)
data.index = data['Date_reported']
data.index.names=['date']
data = data.drop(['Date_reported'],axis=1)


# In[59]:


fig, ax = plt.subplots(2,3)
ax[0,0].bar(dfeuro['Date_reported'],dfeuro['new_EURO'])
ax[0,0].plot(data.index.values,data['new_EURO'].rolling(window=7).mean(),color='red',label='7-day MA')
ax[0,0].set_xlim(date_range)
#ax[0,0].set_xlabel('Date Reported')
ax[0,0].set_ylabel('New death count')
ax[0,0].legend(loc='upper left')
ax[0,0].set_title('EURO')

ax[0,1].bar(dfemro['Date_reported'],dfemro['new_EMRO'])
ax[0,1].plot(data.index.values,data['new_EMRO'].rolling(window=7).mean(),color='red',label='7-day MA')
ax[0,1].set_xlim(date_range)
ax[0,1].legend(loc='upper left')
ax[0,1].set_title('EMRO')
#ax[0,1].xaxis_date()
#ax[0,1].set_xlabel('Date Reported')
#ax[0,1].set_ylabel('New case count')

ax[0,2].bar(dfafro['Date_reported'],dfafro['new_AFRO'])
ax[0,2].plot(data.index.values,data['new_AFRO'].rolling(window=7).mean(),color='red',label='7-day MA')
ax[0,2].set_xlim(date_range)
ax[0,2].legend(loc='upper left')
ax[0,2].set_title('AFRO')
#ax[0,2].set_xlabel('Date Reported')
#ax[0,2].set_ylabel('New case count')

ax[1,0].bar(dfamro['Date_reported'],dfamro['new_AMRO'])
ax[1,0].plot(data.index.values,data['new_AMRO'].rolling(window=7).mean(),color='red',label='7-day MA')
ax[1,0].set_xlim(date_range)
ax[1,0].set_xlabel('Date Reported')
ax[1,0].set_ylabel('New death count')
ax[1,0].legend(loc='upper left')
ax[1,0].set_title('AMRO')

ax[1,1].bar(dfsearo['Date_reported'],dfsearo['new_SEARO'])
ax[1,1].plot(data.index.values,data['new_SEARO'].rolling(window=7).mean(),color='red',label='7-day MA')
ax[1,1].set_xlim(date_range)
ax[1,1].set_xlabel('Date Reported')
ax[1,1].legend(loc='upper left')
ax[1,1].set_title('SEARO')
#ax[1,1].set_ylabel('New case count')

ax[1,2].bar(dfwpro['Date_reported'],dfwpro['new_WPRO'])
ax[1,2].plot(data.index.values,data['new_WPRO'].rolling(window=7).mean(),color='red',label='7-day MA')
ax[1,2].set_xlim(date_range)
ax[1,2].set_xlabel('Date Reported')
ax[1,2].legend(loc='upper left')
ax[1,2].set_title('WPRO')
#ax[1,2].set_ylabel('New case count')
plt.gcf().set_size_inches(30,15)
plt.savefig('worlddeaths.png',dpi=600)


# In[60]:


with sqlite3.connect('covid.db') as conn:
    phnew = (
        "SELECT Date_reported, SUM(New_cases) as new_cases "
        "FROM data WHERE Country='Philippines' "
        "GROUP BY Date_reported"
    )
    
    phdeath = (
        "SELECT Date_reported, SUM(New_deaths) as new_deaths "
        "FROM data WHERE Country='Philippines' "
        "GROUP BY Date_reported"
    )
dfphnew = pd.read_sql_query(phnew, conn, parse_dates=["Date_reported"])
dfphdeath = pd.read_sql_query(phdeath, conn, parse_dates=["Date_reported"])
data = pd.concat([dfphnew,dfphdeath['new_deaths']],axis=1)
data.index = data['Date_reported']
data.index.names=['date']
data = data.drop(['Date_reported'],axis=1)
data


# In[64]:


fig, (ax,ax1) = plt.subplots(1,2)
ax.bar(data.index.values,data['new_cases'])
ax.plot(data.index.values,data['new_cases'].rolling(window=7).mean(),color='red',label='7-day MA')
ax.set_xlim(date_range)
ax.set_ylim(0,50000)
#ax[0,0].set_xlabel('Date Reported')
ax.set_ylabel('New case count')
ax.set_xlabel('Date Reported')
ax.set_title('PH New Cases')
ax.legend(loc='upper left')

ax1.bar(data.index.values,data['new_deaths'])
ax1.plot(data.index.values,data['new_deaths'].rolling(window=7).mean(),color='red',label='7-day MA')
ax1.set_xlim(date_range)
ax1.set_ylim(0,500)
ax1.set_ylabel('New death count')
ax1.set_xlabel('Date Reported')
ax1.set_title('PH New Deaths')
ax1.legend(loc='upper left')
plt.gcf().set_size_inches(30,15)
plt.savefig('PHdata.png',dpi=600)
#ax[0,1].xaxis_date()
#ax[0,1].set_xlabel('Date Reported')
#ax[0,1].set_ylabel('New case count')


# In[ ]:




