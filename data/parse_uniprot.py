import json
import pandas as pd

# Load JSON data
# r"string" : raw string, making Python treat slash as a normal character
json_path = r"C:\Users\User\Desktop\Lab\druggable_dynomics\parseUniprot.json"
with open(json_path, "r") as f:
    data = json.load(f)

# Extract UniProt ID and determine presence of cofactor
output_data = []
for entry in data:
    uniprot_id = entry.get("id", "")
    has_cofactor = "Yes" if entry.get("cofactor") else "No"
    output_data.append({"UniProt_id": uniprot_id, "cofactor": has_cofactor})

# Convert to DataFrame and save to CSV
df = pd.DataFrame(output_data)
csv_path = r"C:\Users\User\Desktop\Lab\druggable_dynomics\if_cofactors.csv"
df.to_csv(csv_path, index=False)
