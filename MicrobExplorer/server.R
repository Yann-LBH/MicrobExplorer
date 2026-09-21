function(input, output, session) {
  
  # Initialize module servers
  mod_home_server("home")
  mod_benchmarks_server("benchmarks")
  mod_qc_server("qc")
  
  # Store settings returned by modules
  reads_cfg   <- mod_reads_server("reads")
  contigs_cfg <- mod_contigs_server("contigs")
  kegg_cfg    <- mod_kegg_server("kegg")
  
  # Pass settings to the export module
  mod_export_server("export_config", 
                    reads_config   = reads_cfg, 
                    contigs_config = contigs_cfg, 
                    kegg_config    = kegg_cfg)
}