# Description of `parseUniprot_881.json`

## 1. First level

| Key            | Description |
|----------------|-------------|
| **id**         | UniProt accession ID (e.g., `P12345`) |
| **name**       | The recommended full name of the protein|
| **protein**    | List of all names of the protein, from commonly used to obsolete, to allow unambiguous identification of a protein. |
| **gene**       | Name(s) of the gene(s) that code for the protein sequence(s) |
| **organism**   | Name(s) of the organism that is the source of the protein sequence. |
| **sequence**   | UniProt protein sequence |
| **cell_location** | Look for membrane protein |
| **cofactor**   | Provided in Cofactor section |
| **binding_site** | Feature subsection under Function section |
| **active_site** | Feature subsection under Function section |
| **dna_binding** | Feature subsection under Function section |
| **zinc_finger** | Feature subsection under Function section |
| **pdb**        | Information of from PDB (Structure section) and so on |
| **alphafold**  | AlphaFold ID |

## 2. Second level 

| Key                     | Description |
|-------------------------|-------------|
| **PDB_ID**              | 4-letter PDB entry ID (e.g., `2JLE`). |
| **method**              | Information of from PDB (Structure section) and so on |
| **resolution**          | Information of from PDB (Structure section) and so on |
| **chains**              | Information of from PDB (Structure section) and so on |
| **resrange**            | Information of from PDB (Structure section) and so on |
| **ligand**              | List of ligands with `comp_id`, `name`, and `pdbx_description`. |
| **initial_release_date**| First release date (ISO format). |
| **ls_R_factor_R_free**  | R-Value Free (Depositor)  |
| **ls_R_factor_R_work**  | R-Value Work (Depositor) |
| **ls_R_factor_obs**     | R-Value Observed (Depositor) |
| **seq_annot**           | Per-chain sequence info. |
| **unobs_res**           | Missing residues. |
| **unobs_atom**          | Missing atoms. |
| **coverage**            | Fraction missing (residues/atoms). |
| **modeled_residue_count** | Count of modeled residues. Experimentally resolved residues |
| **pdbx_mutation**       | Mutation details (if any). |
| **rcsb_mutation_count** | Mutation count. |
