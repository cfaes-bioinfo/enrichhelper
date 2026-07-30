# run_gsea() requires exactly one of 'term_map'/'OrgDb'/'organism'

    Code
      run_gsea(df = fixture$de_df, contrast = "c1")
    Condition
      Error in `run_gsea()`:
      ! Supply either a 'term_map', an 'OrgDb' (for GO), or an 'organism' (for KEGG) -- none was given

---

    Code
      run_gsea(df = fixture$de_df, contrast = "c1", term_map = fixture$term_map,
      OrgDb = "org.Hs.eg.db")
    Condition
      Error in `run_gsea()`:
      ! Supply either a 'term_map' or an 'OrgDb'/'organism', not both

# run_gsea() requires 'organism' for KEGG and 'OrgDb' for GO

    Code
      run_gsea(df = fixture$de_df, contrast = "c1", ontology_type = "KEGG", OrgDb = "org.Hs.eg.db")
    Condition
      Error in `run_gsea()`:
      ! ontology_type = 'KEGG' requires an 'organism' (KEGG species code, e.g. 'hsa')

---

    Code
      run_gsea(df = fixture$de_df, contrast = "c1", ontology_type = "GO", organism = "hsa")
    Condition
      Error in `run_gsea()`:
      ! ontology_type = 'GO' requires an 'OrgDb'

# run_gsea() rejects a KEGG keyType unsupported by KEGG's own API

    Code
      run_gsea(df = fixture$de_df, contrast = "c1", ontology_type = "KEGG", organism = "hsa",
      keyType = "ENSEMBL")
    Output
      Contrast: c1
    Condition
      Error in `resolve_kegg_keytype()`:
      ! keyType = 'ENSEMBL' is not supported by KEGG -- KEGG's API only recognizes: kegg, ncbi-geneid, ncbi-proteinid, uniprot (or 'ENTREZID', auto-translated to 'ncbi-geneid'). Convert your gene IDs to one of these first, e.g. via AnnotationDbi::mapIds(OrgDb, keys, keytype = 'ENSEMBL', column = 'ENTREZID').

