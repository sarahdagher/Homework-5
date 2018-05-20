#dependencies
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

#just in case!
from scipy.stats import sem

#read in data
clinical_trial_df = pd.read_csv('clinicaltrial_data.csv')
mouse_drug_df = pd.read_csv('mouse_drug_data.csv')


#clinical_trial_df.head()
#mouse_drug_df.head()

#find duplicate mice
drop_dup_mouse_id = clinical_trial_df.loc[clinical_trial_df.duplicated(subset=['Mouse ID', 'Timepoint',]),'Mouse ID'].unique()

clean_clinical_trial_df = clinical_trial_df[clinical_trial_df['Mouse ID'].isin(drop_dup_mouse_id)==False]
clean_mouse_df = mouse_drug_df[mouse_drug_df['Mouse ID'].isin(drop_dup_mouse_id)==False]

#combine the dataframes using Mouse ID
combined_data = pd.merge(clean_clinical_trial_df, clean_mouse_df, on = "Mouse ID")
#combined_data.head()
#sort by Timepoint 
sort_by_time = combined_data.sort_values("Timepoint", ascending= True) 

#to reset index
all_sort_by_time = sort_by_time.reset_index()
del all_sort_by_time['index']
all_sort_by_time.head()

#Tumor Response to Treatment

#delete Metastatic Sites column
tumor_volume = all_sort_by_time.drop('Metastatic Sites', axis=1)

#make 'Drug' index
drug_groups = tumor_volume.pivot_table(tumor_volume, index=['Drug','Timepoint'], aggfunc='mean')

#find standard error of mean
sem_tumor_volume = tumor_volume.pivot_table(tumor_volume, index=['Drug','Timepoint'], aggfunc='sem')

#Make a new df for tumor volume to delta
change_tumor_volume = drug_groups
change_tumor_volume['Change in Tumor Vol. (mm3)'] = change_tumor_volume.groupby('Drug')['Tumor Volume (mm3)'].diff()

#remove tumor volume column
new_change_tumor_vol = change_tumor_volume.drop('Tumor Volume (mm3)', axis= 1)

#correct NaNs
new_change_tumor_vol.fillna(0, inplace=True)
new_change_tumor_vol.head(100)

#pivot
sem_table = sem_tumor_volume.pivot_table('Tumor Volume (mm3)', ['Timepoint'], 'Drug')
tumor_vol_table = drug_groups.pivot_table('Tumor Volume (mm3)', ['Timepoint'], 'Drug')
tumor_vol_table.head()   

# Set the x_axis, colors, makers, and xlim for the line graphs with drugs in the legend
x_axis= np.arange(0,50,5)
drugs = tumor_volume["Drug"].unique()
count = np.arange(0,len('drugs'))
plt.xlim(0,50) 

colors = ['green','orange','blue','red','black','yellow','pink','magenta','brown','purple' ]
markers = ['d','v','o','^','*','2','3','4','8','.']

from scipy import stats
plt.style.use('seaborn-whitegrid')
plt.title("Tumor Response to Treatment")
plt.xlabel("Time (Days)")
plt.ylabel("Tumor Volume (mm3)")
plt.grid(alpha = 0.5)
for i in count:
    standard_errors = stats.sem(sem_table[drugs[i]])
    plt.errorbar(x_axis, tumor_vol_table[drugs[i]], yerr = standard_errors, marker= markers[i], color= colors[i], alpha = 0.5, label = drugs[i])
plt.legend(bbox_to_anchor=(1.05,1),loc= 2, borderaxespad = 0.)
plt.xticks(np.arange(min(x_axis)-5, max(x_axis)+5, 5.0))
plt.show()

#Metastatic Response to Treatment
metastatic_response = all_sort_by_time.drop('Tumor Volume (mm3)', axis = 1)
#metastatic_response.head()
#make drug the index and get sem
metastatic_mean = metastatic_response.pivot_table(metastatic_response, index = ['Drug','Timepoint',], aggfunc='mean')
metastatic_sem = metastatic_response.pivot_table(metastatic_response, index = ['Drug','Timepoint',], aggfunc='sem')

#tables
metastatic_mean_table = metastatic_mean.pivot_table('Metastatic Sites', ['Timepoint'],'Drug')
metastatic_sem_table = metastatic_sem.pivot_table('Metastatic Sites', ['Timepoint'],'Drug')
metastatic_mean_table.head()

#graph
x_axis = np.arange(0,50,5)
drugs = all_sort_by_time["Drug"].unique()
count = np.arange(0,len(drugs))
colors = ['green','orange','blue','red','black','yellow','pink','magenta','brown','purple' ]
markers = ['d','v','o','^','*','2','3','4','8','.']

plt.style.use('seaborn-whitegrid')
plt.title("Metastatic Spread During Treatment")
plt.xlabel("Treatment Duration(Days)")
plt.ylabel("Met. Sites")
plt.grid(alpha = 0.5)


for i in count:
    graph_data = stats.sem(metastatic_sem_table[drugs[i]])
    plt.errorbar(x_axis, metastatic_mean_table[drugs[i]], yerr = graph_data, marker= markers[i], color= colors[i], alpha = 0.5, label = drugs[i])
plt.legend(bbox_to_anchor=(1.05,1),loc= 2, borderaxespad = 0.)
plt.xticks(np.arange(min(x_axis)-5, max(x_axis)+5, 5.0))
plt.show()

#Survival Rate
#Survival Rate
mouse_survive = metastatic_response.drop('Metastatic Sites', axis = 1)

#make Drug the index and show means
mouse_survival = mouse_survive.pivot_table(mouse_survive, index=['Drug','Timepoint'], aggfunc='count')
#rename Mouse ID --> Mouse Count
mouse_survival_rename = mouse_survival.rename(columns={"Mouse ID":"Mouse Count"})

#make Timepoint index
mouse_survival_tbl = mouse_survival_rename.pivot_table('Mouse Count',['Timepoint'],'Drug')

#rate
percent_surviving = (1-(mouse_survival_tbl.iloc[0]- mouse_survival_tbl)/mouse_survival_tbl.iloc[0])*100
percent_surviving

#graph
x_axis = np.arange(0,50,5)
drugs = all_sort_by_time["Drug"].unique()
count = np.arange(0,len(drugs))
colors = ['green','orange','blue','red','black','yellow','pink','magenta','brown','purple' ]
markers = ['d','v','o','^','*','2','3','4','8','.']

plt.style.use('seaborn-whitegrid')
plt.title("Survival During Treatment")
plt.xlabel("Time(Days)")
plt.ylabel("Survival Rate(%)")
plt.grid(alpha = 0.5)

#points
for i in count:
    plt.errorbar(x_axis, percent_surviving[drugs[i]], yerr = graph_data, marker= markers[i], color= colors[i], alpha = 0.5, label = drugs[i])
plt.legend(bbox_to_anchor=(1.05,1),loc= 2, borderaxespad = 0.)
plt.xticks(np.arange(min(x_axis)-5, max(x_axis)+5, 5.0))
plt.show()
                                      
#summary bar graph
#first and last row of tumor vol. table
summ_tumor_vol = tumor_vol_table.iloc[[0,-1]]
#tumor volume change
percent_change_tumor_vol= (((summ_tumor_vol -tumor_vol_table.iloc[0])/tumor_vol_table.iloc[0]))*100

#last row
percent_changes = percent_change_tumor_vol.iloc[1:]
percent_changes.sum()                                     
