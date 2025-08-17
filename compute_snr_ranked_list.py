import pandas as pd
import numpy as np

# Load expression data
gct_file = "/mnt/z/Download/Leukemia_hgu95av2.gct"
with open(gct_file) as f:
    lines = f.readlines()
# Skip first two lines, use third as header
data = pd.read_csv(gct_file, sep='\t', skiprows=2)

# Load class labels
cls_file = "/mnt/z/Download/Leukemia.cls"
with open(cls_file) as f:
    cls_lines = f.readlines()
labels = cls_lines[2].strip().split()

# Identify ALL and AML columns
all_cols = [i for i, l in enumerate(labels) if l == 'ALL']
aml_cols = [i for i, l in enumerate(labels) if l == 'AML']
# The data columns start from the third column (after NAME, Description)
all_colnames = data.columns[2:][all_cols]
aml_colnames = data.columns[2:][aml_cols]

# Compute SNR for each gene
def signal_to_noise(a, b):
    a, b = np.array(a), np.array(b)
    mean_a, mean_b = np.mean(a), np.mean(b)
    std_a, std_b = np.std(a, ddof=1), np.std(b, ddof=1)
    # Avoid division by zero
    denom = std_a + std_b if (std_a + std_b) != 0 else 1e-6
    return (mean_a - mean_b) / denom

snr_list = []
for idx, row in data.iterrows():
    gene = row['NAME']
    all_expr = row[all_colnames].values
    aml_expr = row[aml_colnames].values
    snr = signal_to_noise(all_expr, aml_expr)
    snr_list.append((gene, snr))

# Sort by SNR (descending)
snr_list_sorted = sorted(snr_list, key=lambda x: x[1], reverse=True)

# Save as tab-delimited file (no header), as expected by GSEA script
with open("/home/ash022/gsea/leukemia_ranked_gene_list.txt", "w") as f:
    for gene, snr in snr_list_sorted:
        f.write(f"{gene}\t{snr}\n")

print("Ranked gene list saved to /home/ash022/gsea/leukemia_ranked_gene_list.txt")
