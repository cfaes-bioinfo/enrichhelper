# run_ora() requires exactly one of 'term_map'/'OrgDb'

    Code
      run_ora(focal_genes = c("g1", "g2"))
    Condition
      Error in `run_ora()`:
      ! Supply either a 'term_map', an 'OrgDb' (for GO), or an 'organism' (for KEGG) -- none was given

---

    Code
      run_ora(focal_genes = c("g1", "g2"), term_map = make_term_map(), OrgDb = "org.Hs.eg.db")
    Condition
      Error in `run_ora()`:
      ! Supply either a 'term_map' or an 'OrgDb'/'organism', not both

# run_ora() requires 'organism' for KEGG and 'OrgDb' for GO

    Code
      run_ora(focal_genes = c("g1", "g2"), ontology_type = "KEGG", OrgDb = "org.Hs.eg.db")
    Condition
      Error in `run_ora()`:
      ! ontology_type = 'KEGG' requires an 'organism' (KEGG species code, e.g. 'hsa')

---

    Code
      run_ora(focal_genes = c("g1", "g2"), ontology_type = "GO", organism = "hsa")
    Condition
      Error in `run_ora()`:
      ! ontology_type = 'GO' requires an 'OrgDb'

