# Homework-5
#dependencies
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

#just in case!
from scipy.stats import sem

#read in data
clinical_trial_data = pd.read_csv('clinicaltrial_data.csv')
mouse_drug_data = pd.read_csv('mouse_drug_data.csv')

#clinical_trial_data.head()
#mouse_drug_data.head()
