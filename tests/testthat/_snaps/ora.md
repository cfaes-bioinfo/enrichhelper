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

# run_ora() rejects a KEGG keyType unsupported by KEGG's own API

    Code
      run_ora(focal_genes = c("g1", "g2"), ontology_type = "KEGG", organism = "hsa",
      keyType = "ENSEMBL", verbose = FALSE)
    Condition
      Error in `resolve_kegg_keytype()`:
      ! keyType = 'ENSEMBL' is not supported by KEGG -- KEGG's API only recognizes: kegg, ncbi-geneid, ncbi-proteinid, uniprot (or 'ENTREZID', auto-translated to 'ncbi-geneid'). Convert your gene IDs to one of these first, e.g. via AnnotationDbi::mapIds(OrgDb, keys, keytype = 'ENSEMBL', column = 'ENTREZID').

