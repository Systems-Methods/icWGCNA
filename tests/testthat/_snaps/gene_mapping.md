# gene_mapping success defaults

    Code
      gene_mapping(tmp_expr_data, tmp_mapping_file)
    Message
      3 exprs_data rows could not be linked using the mapping file, resulting in 7 rows. 
      These rows link to 6 distinct gene symbols. 
      Will compress duplicate rows using the mean method with no transformation.
    Output
                 ID_1      ID_2        ID_3
      Gene_1 2.666667  8.333333    5.000000
      Gene_2 1.000000 10.000000    2.000000
      Gene_3 4.333333  6.666667    3.666667
      Gene_4 2.000000  9.000000    3.000000
      Gene_5 3.000000  8.000000    6.000000
      Gene_6 6.000000  5.000000 1000.000000

# gene_mapping success log mean

    Code
      gene_mapping(tmp_expr_data, tmp_mapping_file, "highest_mean", "log", verbose = FALSE)
    Output
             ID_1 ID_2 ID_3
      Gene_1    5    6   10
      Gene_2    1   10    2
      Gene_3    4    7    7
      Gene_4    2    9    3
      Gene_5    3    8    6
      Gene_6    6    5 1000

# gene_mapping pca

    Code
      gene_mapping(tmp_expr_data, tmp_mapping_file, compress_fun = "pc1")
    Message
      3 exprs_data rows could not be linked using the mapping file, resulting in 7 rows. 
      These rows link to 6 distinct gene symbols. 
      Will compress duplicate rows using the pc1 method with no transformation.
    Output
                     ID_1        ID_2      ID_3
      Gene_1  0.001982281    1.578949 -1.580932
      Gene_2  1.000000000    2.000000  1.000000
      Gene_3 -0.321982279    1.817832 -1.495850
      Gene_4  2.000000000    3.000000  1.000000
      Gene_5  3.000000000    6.000000  1.000000
      Gene_6  6.000000000 1000.000000 10.000000

# gene_mapping no dups

    Code
      gene_mapping(tmp_expr_data, tmp_mapping_file)
    Message
      0 exprs_data rows could not be linked using the mapping file, resulting in 10 rows. 
      These rows link to 10 distinct gene symbols. 
      Will compress duplicate rows using the mean method with no transformation.
    Output
             ID_1 ID_2 ID_3
      Gene_0    1   10    2
      Gene_1    2    9    3
      Gene_2    3    8    6
      Gene_3    4    7    7
      Gene_4    5    6   10
      Gene_5    6    5 1000
      Gene_6    7    4    1
      Gene_7    8    3    1
      Gene_8    9    2    1
      Gene_9   10    1    1

