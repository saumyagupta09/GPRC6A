alt_file=LRT_test_values.txt
null_file=LRT_null_value.txt
omega_file="branchmodel_results.txt"
species="$1"                              # Species name passed as argument

# === Get alternative lnL for the species ===
alt_lnL=$(awk -v sp="$species" '
  $0 == sp { getline; print $5 }
' "$alt_file")

# === Get null lnL ===
null_lnL=$(awk '{ print $5 }' "$null_file")

# === Get omega values ===
read bg_omega fg_omega <<< $(awk -v sp="$species" '
  $0 == sp { getline
    for (i=1; i<=NF; i++) {
      if ($i ~ /^[0-9.]+$/) omega[++n] = $i
    }
    if (n >= 2) print omega[1], omega[2]
  }
' "$omega_file")


# === Check if values were extracted properly ===
if [[ -z "$alt_lnL" || -z "$null_lnL" ]]; then
  echo "Error: lnL not found for species or null model."
  exit 1
fi

# === Compute LRT and p-value ===
python3 - <<EOF
import math
from scipy.stats import chi2

lnL_alt = float("$alt_lnL")
lnL_null = float("$null_lnL")
lrt = 2 * (lnL_alt - lnL_null)
p = chi2.sf(lrt, df=1)

print("Species\tlnL_null\tlnL_alt\tLRT\tp-value\tbg_omega\tfg_omega")
print(f"{'$species'}\t{lnL_null}\t{lnL_alt}\t{lrt:.4f}\t{p:.4g}\t{'$bg_omega' if '$bg_omega' else 'NA'}\t{'$fg_omega' if '$fg_omega' else 'NA'}")
EOF
