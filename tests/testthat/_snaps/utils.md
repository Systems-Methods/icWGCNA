# dropModuels logFlag and Kurts

    Code
      dropModuels(testing_data[2:10, 1:10], Kurts = 1:9, corCut = 0.01, logFlag = TRUE)
    Message
      eigegenes trimmed to 2 due to correlation > 0.01 max eigenCor = 0
    Output
            TCGA-05-4389-01 TCGA-91-6829-01 TCGA-69-7763-01 TCGA-05-4410-01
      GTSE1        2.163415        2.129160        1.866310        2.162300
      ITGA8        1.718113        2.183678        2.194833        2.183374
            TCGA-50-6595-01 TCGA-49-6743-01 TCGA-44-2657-01 TCGA-93-A4JP-01
      GTSE1        2.273465        2.338734        2.123255        2.165011
      ITGA8        1.921559        2.048131        2.166846        2.371075
            TCGA-53-7624-01 TCGA-78-7167-01
      GTSE1       2.3361746        1.601386
      ITGA8       0.6467891        2.033267

# dropModuels large corCut

    Code
      dropModuels(testing_data[2:10, 1:10], corCut = 0.01, logFlag = TRUE)
    Message
      eigegenes trimmed to 4 due to correlation > 0.01 max eigenCor = 0
    Output
              TCGA-05-4389-01 TCGA-91-6829-01 TCGA-69-7763-01 TCGA-05-4410-01
      HIF3A          2.271002        2.325295        2.375454        1.912294
      RTN4RL2        2.081577        1.932014        1.994265        2.386154
      GTSE1          2.163415        2.129160        1.866310        2.162300
      ITGA8          1.718113        2.183678        2.194833        2.183374
              TCGA-50-6595-01 TCGA-49-6743-01 TCGA-44-2657-01 TCGA-93-A4JP-01
      HIF3A          1.987929        1.612892        2.307523        2.366967
      RTN4RL2        1.826161        2.287025        2.065355        1.759374
      GTSE1          2.273465        2.338734        2.123255        2.165011
      ITGA8          1.921559        2.048131        2.166846        2.371075
              TCGA-53-7624-01 TCGA-78-7167-01
      HIF3A         2.3450985        2.414583
      RTN4RL2       2.1684135        2.222156
      GTSE1         2.3361746        1.601386
      ITGA8         0.6467891        2.033267

