#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np


# In[4]:


cohort_name = 'UKBB'
gene_name = 'SLC25A46'


# In[5]:


with open(f'{cohort_name}.VEPannotated') as f:
    for i, line in enumerate(f):
        if line.startswith("#Uploaded_variation"):
            cols = line.lstrip("#").rstrip("\n").split("\t")
            header_row = i
            break

# Read only the variant table (everything below that header)
VEP_df = pd.read_csv(
    f'{cohort_name}.VEPannotated',
    sep="\t",
    skiprows=header_row + 1,  # drop meta-header lines
    names=cols,               # assign correct column names
    comment="#",              # ignore any stray comment lines
    low_memory=False,
)

extra_df = (
    VEP_df["Extra"]
      .str.split(';')
      .apply(lambda items: dict(s.split('=', 1) for s in items))
      .apply(pd.Series)          # turn each row-dict into columns
)

# 2) merge with the main table and drop the original Extra column
VEP_df = pd.concat([VEP_df.drop(columns="Extra"), extra_df], axis=1)

VEP_df


# In[6]:


############################################
# 1) Loss-of-function (LOF) variants
############################################
lof_terms = [
    'stop_gained',
    'frameshift_variant',
    'start_lost',
    'stop_lost',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'transcript_ablation',
    'exon_loss_variant'
]

lof_mask = VEP_df['Consequence'].str.contains('|'.join(lof_terms), regex=True)
lof_df   = VEP_df[lof_mask]

############################################
# 2) Non-synonymous, non-LOF variants
############################################
nonsyn_terms = [
    'missense_variant',
    'inframe_insertion',
    'inframe_deletion',
    'stop_retained_variant',
    'start_retained_variant',
    'protein_altering_variant',
    'coding_sequence_variant',      # generic “changes CDS” term
    'splice_region_variant'         # outside the AG/GT ±1-2 sites
]

nonsyn_mask = (
      VEP_df['Consequence'].str.contains('|'.join(nonsyn_terms), regex=True)
  & ~lof_mask                                               # not LOF
  & ~VEP_df['Consequence'].str.contains('synonymous_variant')  # not silent
)

nonsynonymous_df = VEP_df[nonsyn_mask]

############################################
# 3) CADD ≥ 20 variants (any consequence)
############################################
VEP_df['CADD_PHRED'] = pd.to_numeric(
    VEP_df['CADD_PHRED'], errors='coerce'
)

cadd_mask = VEP_df['CADD_PHRED'] >= 20
cadd_df   = VEP_df[cadd_mask]

############################################
# 4) Alpha Missense Likely Pathogenic variants
############################################
alphamissense_df = VEP_df[VEP_df['am_class'] == 'likely_pathogenic']


############################################
# 5) All variants (for completeness)
############################################
all_df = VEP_df.copy()



# In[7]:


nonsynonymous_df


# In[ ]:


list_of_dfs = {'ALL': all_df,'CADD': cadd_df, 'LOF': lof_df, 'nonsyn': nonsynonymous_df, 'AM': alphamissense_df}


# In[ ]:


setid_filename = f'{cohort_name}.SETID'

with open(setid_filename, 'w') as output_file:
    for df_key, df_value in list_of_dfs.items():
      for variantID in df_value['Uploaded_variation']:
          output_file.write(
              f'{gene_name}_{df_key}\t{variantID}\t{gene_name}\n')

output_file.close()

