# cdotplot() errors with a clear message when a required column is missing

    Code
      cdotplot(df)
    Condition
      Error in `cdotplot()`:
      ! df is missing required column(s): term

# cdotplot() errors when DE_dirs is given but DE_direction is missing

    Code
      cdotplot(df, DE_dirs = "up")
    Condition
      Error in `cdotplot()`:
      ! df is missing the 'DE_direction' column, needed because 'DE_dirs' was specified

# cdotplot() errors when x_var/fill_var/facet_var names a missing column

    Code
      cdotplot(make_enrich_df(), facet_var1 = "ontology")
    Condition
      Error in `cdotplot()`:
      ! df is missing column(s) requested via x_var/fill_var/label_var/facet_var1/facet_var2: ontology

# cdotplot() errors when nothing is left to plot after filtering

    Code
      cdotplot(df)
    Condition
      Error in `cdotplot()`:
      ! No rows left to plot after filtering for sig == TRUE and the given contrasts/DE_dirs

