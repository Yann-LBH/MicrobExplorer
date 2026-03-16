def filter_tsv_contigs(path_in, path_out, path_report, min_size, sample_name):
    df.to_csv('output.tsv', 
          sep='\t', 
          index=False, 
          na_rep='NA', 
          encoding='utf-8',
          quoting=3)