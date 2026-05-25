import os

script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(script_dir)
pdb_list_path = os.path.join(base_dir, "data", "existing_pdbs.txt")

with open(pdb_list_path, 'r') as f:
    existing_pdbs = set(line.strip().lower() for line in f if line.strip())

test_pdb = "1ajx"
print(f"Is '{test_pdb}' in set? {test_pdb in existing_pdbs}")
print(f"Total PDBs loaded: {len(existing_pdbs)}")
print(f"Sample from set: {list(existing_pdbs)[:5]}")

bindingdb_sample = "1W5Y,1W5X,1W5W,1W5V,2FDE,7UPJ,6UWC,6UWB,6D0E,6D0D,5TYS,5TYR,4I8Z,4I8W,4HLA,4FE6,3T11,3PSU,3PHV,3GGX,3GGV,3GGA,3CKT,3BHE,3BGC,3BGB,2ZGA,2WKZ,2UY0,2UXZ,2UPJ,2QNQ,2QNP,2QNN,2PWR,2PWC,2PQZ,2CEN,2CEM,2CEJ,2BQV,2BPZ,2BPY,2BPX,2BPW,2BPV,2BBB,2BB9,2AQU,2A4F,1ZSR,1ZSF,1YT9,1XL5,1XL2,1WBM,1WBK,1UPJ,1U8G,1T7K,1SP5,1OHR,1NPW,1NPV,1NPA,1NH0,1MUI,1M0B,1IIQ,1HVL,1HVK,1HVJ,1HVI,1HTG,1HTF,1HTE,1HSG,1HPX,1HPV,1HPS,1HOS,1HIH,1HHP,1HBV,1GNO,1G35,1G2K,1FQX,1EC3,1EC2,1EC1,1EC0,1EBZ,1EBW,1DIF,1D4J,1D4I,1D4H,1C70,1AJX,1AJV"
parts = bindingdb_sample.lower().split(',')
matches = [p for p in parts if p in existing_pdbs]
print(f"Matches found in sample: {matches}")
