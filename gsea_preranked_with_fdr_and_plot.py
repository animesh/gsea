import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re

# 1. Load ranked gene list
ranked_file = "/home/ash022/gsea/leukemia_ranked_gene_list.txt"
rnk = pd.read_csv(ranked_file, sep='\t', header=None)
rnk.columns = ['GENE_SYMBOL', 'SCORE']
rnk = rnk.dropna()
rnk = rnk[rnk['GENE_SYMBOL'].str.lower() != 'null']
rnk = rnk[rnk['SCORE'].apply(lambda x: str(x).lower() != 'null')]
genes = rnk['GENE_SYMBOL'].tolist()
scores = rnk['SCORE'].astype(float).values

# 2. Load gene sets from GMT
def load_gmt(gmt_path):
    gene_sets = {}
    with open(gmt_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            set_name = parts[0]
            gene_list = set(parts[2:])  # skip description
            gene_sets[set_name] = gene_list
    return gene_sets

gmt_file = "/mnt/z/Download/c2.symbols.gmt"
gene_sets = load_gmt(gmt_file)

# 2b. Filter gene sets to those mentioned in the article
keywords = [
    'IL6', 'MYC', 'P53', 'NFkB', 'Ras', 'E2F', 'RB', 'TGF', 'WNT', 'APC', 'KRAS', 'EGFR', 'AKT',
    'MAPK', 'ERK', 'JAK', 'STAT', 'BCR', 'ABL', 'oncogenic', 'leukemia', 'ALL', 'AML', 'pathway', 'signature'
]
pattern = re.compile(r"(" + "|".join(keywords) + ")", re.IGNORECASE)
filtered_gene_sets = {k: v for k, v in gene_sets.items() if pattern.search(k)}

print(f"Found {len(filtered_gene_sets)} gene sets matching article keywords.")

# 3. GSEA Enrichment Score calculation (from the paper)
def enrichment_score(ranked_genes, gene_set, scores, return_running_sum=False):
    N = len(ranked_genes)
    Nh = len(gene_set)
    tag_indicator = np.array([g in gene_set for g in ranked_genes], dtype=bool)
    correl_vector = np.abs(scores)
    NR = np.sum(correl_vector[tag_indicator])
    if NR == 0:
        if return_running_sum:
            return 0.0, [0.0]*N
        return 0.0
    running_sum = []
    sum_val = 0.0
    for i in range(N):
        if tag_indicator[i]:
            sum_val += correl_vector[i] / NR
        else:
            sum_val -= 1.0 / (N - Nh)
        running_sum.append(sum_val)
    ES = max(running_sum, key=abs)
    if return_running_sum:
        return ES, running_sum
    return ES

# 4. Permutation test for significance
def permutation_test(ranked_genes, scores, gene_set, n_perm=100, seed=42):
    observed_es = enrichment_score(ranked_genes, gene_set, scores)
    np.random.seed(seed)
    null_dist = []
    for _ in range(n_perm):
        permuted = np.random.permutation(ranked_genes)
        es = enrichment_score(permuted, gene_set, scores)
        null_dist.append(es)
    null_dist = np.array(null_dist)
    if observed_es >= 0:
        pval = np.sum(null_dist >= observed_es) / n_perm
    else:
        pval = np.sum(null_dist <= observed_es) / n_perm
    return observed_es, pval, null_dist

# 5. Run GSEA for filtered gene sets with progress output
results = []
all_nulls = []
gene_set_items = list(filtered_gene_sets.items())
for i, (set_name, gene_set) in enumerate(gene_set_items):
    print(f"Processing gene set {i+1}/{len(gene_set_items)}: {set_name}")
    es, pval, null_dist = permutation_test(genes, scores, gene_set, n_perm=100)
    results.append({'GeneSet': set_name, 'ES': es, 'p-value': pval})
    all_nulls.append(null_dist)

results_df = pd.DataFrame(results)
results_df = results_df.sort_values('p-value')

# 6. FDR/q-value calculation (Benjamini-Hochberg)
def benjamini_hochberg(pvalues):
    pvalues = np.array(pvalues)
    n = len(pvalues)
    sorted_indices = np.argsort(pvalues)
    sorted_pvalues = pvalues[sorted_indices]
    qvalues = np.zeros(n)
    prev_q = 1.0
    for i in reversed(range(n)):
        rank = i + 1
        q = sorted_pvalues[i] * n / rank
        q = min(q, prev_q)
        qvalues[sorted_indices[i]] = q
        prev_q = q
    return qvalues

results_df['q-value'] = benjamini_hochberg(results_df['p-value'].values)
results_df.to_csv('gsea_preranked_results_with_fdr.csv', index=False)
print("GSEA complete. Results with FDR saved to gsea_preranked_results_with_fdr.csv")

# 7. Plotting top N gene sets
def plot_enrichment(ranked_genes, scores, gene_set, set_name, outdir='gsea_plots'):
    os.makedirs(outdir, exist_ok=True)
    _, running_sum = enrichment_score(ranked_genes, gene_set, scores, return_running_sum=True)
    plt.figure(figsize=(8,4))
    plt.plot(running_sum, color='blue')
    plt.title(f'Enrichment Score: {set_name}')
    plt.xlabel('Rank in Ordered Dataset')
    plt.ylabel('Running Enrichment Score')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{outdir}/{set_name}_enrichment.png")
    plt.close()

# Plot all processed gene sets
top_sets = results_df['GeneSet']
for set_name in top_sets:
    plot_enrichment(genes, scores, filtered_gene_sets[set_name], set_name)

print("Enrichment plots saved in 'gsea_plots' directory.")
