# run_ora() requires exactly one of 'term_map'/'OrgDb'

    Code
      run_ora(focal_genes = c("g1", "g2"))
    Condition
      Error in `run_ora()`:
      ! Supply either a 'term_map' or an 'OrgDb' (neither was given)

---

    Code
      run_ora(focal_genes = c("g1", "g2"), term_map = make_term_map(), OrgDb = "org.Hs.eg.db")
    Condition
      Error in `run_ora()`:
      ! Supply either a 'term_map' or an 'OrgDb', not both

