import re

iqtree_file = "Artiodactyla.aln.iqtree"
msa_file = "Artiodactyla.aln"
tree_file = "Artiodactyla.aln.treefile"
out_root = "neutral"

with open(iqtree_file) as f:
    content = f.read()

# Extract best-fit model
best_model = re.search(r'Best-fit model according to BIC: (\S+)', content).group(1)

# Extract frequencies
pi_match = re.search(r'pi\(A\) = ([0-9.]+)\s+pi\(C\) = ([0-9.]+)\s+pi\(G\) = ([0-9.]+)\s+pi\(T\) = ([0-9.]+)', content)
A, C, G, T = pi_match.groups()

# Extract rates
rates_match = re.search(r'Rate parameter R:\s+A-C: ([0-9.]+)\s+A-G: ([0-9.]+)\s+A-T: ([0-9.]+)\s+C-G: ([0-9.]+)\s+C-T: ([0-9.]+)\s+G-T: ([0-9.]+)', content)
rates = " ".join(rates_match.groups())

# Extract gamma alpha
alpha_match = re.search(r'Gamma shape alpha: ([0-9.]+)', content)
alpha = alpha_match.group(1)

# Construct phyloFit command
cmd = f"""phyloFit \
--msa {msa_file} \
--tree {tree_file} \
--subst-mod REV \
--rates {rates} \
--pi {A} {C} {G} {T} \
--gamma --alpha {alpha} --ncat 4 \
--out-root {out_root} \
--msa-format FASTA"""

print(cmd)
