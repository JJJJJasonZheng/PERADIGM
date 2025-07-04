# train_icd10_embeddings.py
"""
Train word embeddings using medical descriptions to generate ICD-10 code embeddings.
"""

import pandas as pd
import numpy as np
import string
from gensim.models import Word2Vec
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
import nltk
from itertools import chain

# Download NLTK resources
nltk.download('punkt')
nltk.download('stopwords')

DATA_DIR = '../data'

# 1. Load Data
hesin_descrip = pd.read_table(f'{DATA_DIR}/hesin_200k_descrip_icd10.txt')
icd_mapping = pd.read_csv(f'{DATA_DIR}/ICD10_mapping.csv')
icd_200k_codes = pd.read_csv(f'{DATA_DIR}/icd_200k_codes.csv')  # should have a column 'code'

# 2. Aggregate Descriptions Per Patient (eid)
df_aggregated = hesin_descrip.groupby('eid')['meaning'].apply(' '.join).reset_index()

# 3. Text Preprocessing
def clean_tokenize(text):
    text = text.lower()
    text = text.translate(str.maketrans('', '', string.punctuation))
    tokens = word_tokenize(text)
    tokens = [word for word in tokens if word not in stopwords.words('english')]
    return tokens

df_aggregated['tokens'] = df_aggregated['meaning'].apply(clean_tokenize)
sentences = df_aggregated['tokens'].tolist()

# 4. Word2Vec Training
print(f"Training Word2Vec on {len(sentences)} sentences...")
model = Word2Vec(sentences, vector_size=100, window=5, min_count=1, workers=4)

# 5. Map ICD-10 Codes to Descriptions and Tokenize
icd_mapping['tokens'] = icd_mapping['meaning'].apply(clean_tokenize)

# 6. Compute Average Embeddings for Each ICD-10 Code
def average_embedding(tokens, model):
    vectors = [model.wv[token] for token in tokens if token in model.wv]
    if vectors:
        return np.mean(vectors, axis=0)
    else:
        return np.zeros(model.vector_size)

icd_mapping['average_embedding'] = icd_mapping['tokens'].apply(lambda tokens: average_embedding(tokens, model))

# 7. Filter to 200k Codes of Interest
filtered_icd = icd_mapping[icd_mapping['coding'].isin(icd_200k_codes['code'])]

# 8. Convert Embeddings to DataFrame
embeddings_df = pd.DataFrame(filtered_icd['average_embedding'].tolist(), index=filtered_icd['coding'])
embeddings_df.columns = [f'Dimension_{i+1}' for i in range(embeddings_df.shape[1])]

# 9. Save Embeddings
embeddings_df.to_csv(f'{DATA_DIR}/icd10_embeddings_200k.csv')
print(f"Embeddings for {len(embeddings_df)} ICD-10 codes saved to icd10_embeddings_200k.csv")
