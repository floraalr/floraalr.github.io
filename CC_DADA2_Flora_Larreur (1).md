R Notebook
================

# Tutoriel DADA2

Dans ce tutoriel, nous suivons les étapes principales du pipeline DADA2
afin de traiter des données de séquençage Illumina. À partir de fichiers
FASTQ appariés déjà démultiplexés et nettoyés, nous allons produire une
table d’ASV et attribuer une taxonomie aux séquences.

## I. Se préparer à l’analyse des données

``` r
#La toute première étape consiste à charger le package DADA2, à installer préalablement sur R

library(dada2)
```

    ## Loading required package: Rcpp

``` r
#On commence par indiquer le répertoire de travail où se trouvent les fichiers FASTQ issus du séquençage (après décompression). La commande list.files() permet de vérifier que les fichiers ont bien été localisés.

path <- "~/MiSeq_SOP" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

``` r
#Cette étape identifie et classe les fichiers FASTQ avant (R1) et arrière (R2) présents dans le répertoire de travail. Les fonctions list.files() et sort() permettent de repérer les fichiers correspondant aux deux types de lectures et de s’assurer qu’ils sont listés dans le même ordre pour chaque échantillon.

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
```

``` r
#On extrait ici les noms d’échantillons à partir des fichiers FASTQ avant (fnFs).La fonction basename() récupère uniquement le nom du fichier (sans le chemin), puis strsplit() le découpe en plusieurs parties selon le caractère “_”.L’expression [ , 1] sélectionne la première partie du nom, qui correspond à l’identifiant de l’échantillon.

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

## II. Inspecter les profils de qualité de lecture des fichiers FASTQ

``` r
#Nous générons ici un profil de qualité pour les deux premiers fichiers FASTQ avant (R1).Le graphique produit par la fonction plotQualityProfile() affiche la qualité moyenne des bases le long des lectures, ce qui aide à déterminer où les séquences doivent être tronquées lors de l’étape de filtrage.

plotQualityProfile(fnFs[1:2])
```

![](CC_DADA2_Flora_Larreur_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
Le fond gris illustre la distribution des scores de qualité, la ligne
verte la moyenne, les lignes orange les quartiles, et la ligne rouge la
proportion de lectures atteignant chaque position. La qualité reste
bonne jusqu’à environ 245 pb.

``` r
#Après avoir examiné les lectures avant, nous observons ici la qualité des lectures arrière (R2). Le graphique obtenu permet d’évaluer la dégradation éventuelle de la qualité vers la fin des séquences, ce qui guidera également le choix des paramètres de troncature lors du filtrage.

plotQualityProfile(fnRs[1:2])
```

![](CC_DADA2_Flora_Larreur_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
Contrairement aux lectures forward, la qualité chute fortement après 160
pb, ce qui est souvent observé avec la méthode de séquençage Illumina.
Ainsi, les lectures restent de bonne qualité jusqu’à 160 pb dans ce cas.

## III. Filtrer et rogner

``` r
#Cette commande définit le chemin et le nom des fichiers FASTQ filtrés pour les lectures avant (R1).

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
```

``` r
# Nous faisons de même pour les lectures arrières (R2)
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#On attribue ici les noms d’échantillons comme identifiants aux fichiers filtrés avant (filtFs) et arrière (filtRs). Cela permet à DADA2 de suivre plus facilement les paires de fichiers correspondant à chaque échantillon tout au long du pipeline.

names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

``` r
#Cette commande effectue le filtrage et la troncature des lectures avant (R1) et arrière (R2).

out <- filterAndTrim(fnFs,#fichiers brut "avant"
                     filtFs, #fichiers filtrés "avant" de sortie
                     fnRs, #fichiers brut "arrières" 
                     filtRs,#fichiers "arrières" filtrés de sortie 
 truncLen=c(240,160),#tronque les lectures à 240 pb pour les forward et 160 pb pour les reverse. 
              maxN=0,#élimine toute lecture contenant un N (base indéterminée)
              maxEE=c(2,2),#conserve les lectures dont le nombre d’erreurs attendues est ≤ 2.
             truncQ=2, #tronque une lecture dès qu’un score de qualité ≤ 2 est rencontré.
             rm.phix=TRUE,#supprime les séquences PhiX de contrôle.
              compress=TRUE, #compresse les fichiers de sortie.
              multithread=FALSE)#empêche d’utiliser plusieurs cœurs, ce qui pourrait accélérer le traitement. Dans ce cas, on veut éviter de faire "planter" la VM.
head(out)#affiche un tableau, résultante du filtrage 
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

Comme on peut le voir, un très grande partie des séquences est conservée
par ce filtrage !

## IV. Constat des taux d’erreur de séquençage

Cette étape de la pipeline a pour but de constater les erreurs de
séquençage (Illumina). DADA2 est un très bon outil pour les distinguer

``` r
#Apprentissage du modèle d’erreur pour les lectures forward filtrées, utilisé pour distinguer erreurs et vrais variants.
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
##Apprentissage du modèle d’erreur pour les lectures reverse filtrées, utilisé pour distinguer erreurs et vrais variants.
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
#affichage graphique qui visualise le modèle d’erreur appris pour les lectures forward (R1). On ne n'utilise pas R2 qui lui est de moins bonne qualité

plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    ## Transformation introduced infinite values in continuous y-axis

![](CC_DADA2_Flora_Larreur_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->
Le graphique présente les taux d’erreur pour chaque transition (A→C,
A→G, …) en fonction du score de qualité. Les points noirs sont les
erreurs observées, la ligne noire les erreurs estimées par l’algorithme,
et la ligne rouge les erreurs théoriques attendues selon le Q score. On
voit que les points suivent bien la ligne noire, indiquant que les taux
observés correspondent aux taux estimés, et que les erreurs diminuent
lorsque le Q score augmente, comme prévu. Le modèle d’erreur est donc
fiable.

## V. Inférence des variants de séquence (ASV)

``` r
#Cette étape applique l’algorithme DADA2 aux lectures forward filtrées (filtFs) en utilisant le modèle d’erreur appris (errF).

dadaFs <- dada(filtFs,#Fichiers forward filtrés
               err=errF,#Modèle d'erreurs appris précédemment
               multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

exemple :

Sample 1 – 7113 reads in 1979 unique sequences.

Cela signifie que l’échantillon 1 contenait 7113 lectures de bonne
qualité après filtrage.

Parmi elles, DADA2 a identifié 1979 séquences distinctes, c’est-à-dire
des variantes uniques avant la correction finale des erreurs et la
fusion des lectures.

Le nombre de lectures (“reads”) varie selon les échantillons (de ~2900 à
plus de 18 000), ce qui est normal : la profondeur de séquençage n’est
jamais parfaitement égale.

Le nombre de séquences uniques reflète la diversité brute apparente de
chaque échantillon avant la détection des artefacts.

Des valeurs plus élevées (ex. Sample 13 : 3707 séquences uniques)
peuvent indiquer une diversité microbienne plus riche,ou simplement un
plus grand nombre de lectures, ce qui augmente les chances de détecter
de nouvelles variantes.

Ces chiffres vont diminuer après les étapes suivantes (fusion et
suppression des chimères), car DADA2 va éliminer les séquences
redondantes ou erronées,ne conserver que les ASV réels, c’est-à-dire les
séquences biologiquement distinctes.

``` r
#On applique ici aussi l’algorithme DADA2 aux lectures reverse filtrées (filtRs) en utilisant le modèle d’erreur appris (errR).

dadaRs <- dada(filtRs, #fichiers reverse filtrés
               err=errR,#modèle d'erreur appris précedemment 
               multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

``` r
#ici, on ffiche le résultat du débruitage du premier échantillon pour vérifier les variantes de séquences (ASV) détectées.
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

Avant débruitage : 1979 séquences uniques Après débruitage : 128 ASV

Le débruitage réduit le bruit de séquençage et produit un jeu de données
fiable pour les analyses taxonomiques et écologiques. Les paramètres
affichés montrent que le pipeline a utilisé des seuils stricts,
garantissant peu de faux positifs.

## VI. Fusionner les lectures appariées

A présent, cette étape assemble les lectures avant et arrière en
séquences complètes grâce à leur zone de chevauchement. Seules les
paires cohérentes et de bonne qualité sont conservées, garantissant des
séquences fiables pour la construction de la table d’ASV.

``` r
#Nous combinons les lectures « forward » et « reverse » après débruitage afin de former des séquences consensus, tout en vérifiant leur cohérence et leur chevauchement.

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 6540 paired-reads (in 107 unique pairings) successfully merged out of 6891 (in 197 pairings) input.

    ## 5028 paired-reads (in 101 unique pairings) successfully merged out of 5190 (in 157 pairings) input.

    ## 4986 paired-reads (in 81 unique pairings) successfully merged out of 5267 (in 166 pairings) input.

    ## 2595 paired-reads (in 52 unique pairings) successfully merged out of 2754 (in 108 pairings) input.

    ## 2553 paired-reads (in 60 unique pairings) successfully merged out of 2785 (in 119 pairings) input.

    ## 3646 paired-reads (in 55 unique pairings) successfully merged out of 4109 (in 157 pairings) input.

    ## 6079 paired-reads (in 81 unique pairings) successfully merged out of 6514 (in 198 pairings) input.

    ## 3968 paired-reads (in 91 unique pairings) successfully merged out of 4388 (in 187 pairings) input.

    ## 14233 paired-reads (in 143 unique pairings) successfully merged out of 15355 (in 352 pairings) input.

    ## 10528 paired-reads (in 120 unique pairings) successfully merged out of 11165 (in 278 pairings) input.

    ## 11154 paired-reads (in 137 unique pairings) successfully merged out of 11797 (in 298 pairings) input.

    ## 4349 paired-reads (in 85 unique pairings) successfully merged out of 4802 (in 179 pairings) input.

    ## 17431 paired-reads (in 153 unique pairings) successfully merged out of 17812 (in 272 pairings) input.

    ## 5850 paired-reads (in 81 unique pairings) successfully merged out of 6095 (in 159 pairings) input.

    ## 3716 paired-reads (in 86 unique pairings) successfully merged out of 3894 (in 147 pairings) input.

    ## 6865 paired-reads (in 99 unique pairings) successfully merged out of 7191 (in 187 pairings) input.

    ## 4426 paired-reads (in 67 unique pairings) successfully merged out of 4603 (in 127 pairings) input.

    ## 4576 paired-reads (in 101 unique pairings) successfully merged out of 4739 (in 174 pairings) input.

    ## 6092 paired-reads (in 109 unique pairings) successfully merged out of 6315 (in 173 pairings) input.

    ## 4269 paired-reads (in 20 unique pairings) successfully merged out of 4281 (in 28 pairings) input.

``` r
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                       sequence
    ## 1 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACAGG
    ## 2 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAGG
    ## 3 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATGTCGGGGCTCAACCCCGGCCTGCCGTTGAAACTGGCGGCCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ## 4 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTTTTAAGTCAGCGGTAAAAATTCGGGGCTCAACCCCGTCCGGCCGTTGAAACTGGGGGCCTTGAGTGGGCGAGAAGAAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCCTTCCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCGAACAGG
    ## 5 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGACTCTCAAGTCAGCGGTCAAATCGCGGGGCTCAACCCCGTTCCGCCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAAGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCCTACCGGCGCGCAACTGACGCTCATGCACGAAAGCGTGGGTATCGAACAGG
    ## 6 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1       579       1       1    148         0      0      1   TRUE
    ## 2       470       2       2    148         0      0      2   TRUE
    ## 3       449       3       4    148         0      0      1   TRUE
    ## 4       430       4       3    148         0      0      2   TRUE
    ## 5       345       5       6    148         0      0      1   TRUE
    ## 6       282       6       5    148         0      0      2   TRUE

## VII. Construire une table de séquence

``` r
#Nous créons une table résumant toutes les séquences obtenues pour chaque échantillon, étape essentielle avant le filtrage et l’assignation taxonomique.
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]  20 293

``` r
#Cette commande calcule la fréquence des différentes longueurs de séquences après fusion et construction de la table d’ASV. Utile pour détecter des séquences anormalement courtes ou longues qui pourraient provenir d’erreurs de séquençage.
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  88 196   6   2

## VIII. Supprimer les chimères

A ce moment du tutoriel, nous identifions et supprimons les séquences
chimériques, c’est-à-dire celles formées artificiellement lors de la PCR
par la combinaison de deux séquences réelles. Cela permet d’obtenir une
table d’ASV purifiée, contenant uniquement des séquences fiables
représentatives des vrais organismes présents dans les échantillons.

``` r
##Cette commande filtre les ASV suspectes d’être des chimères en comparant chaque séquence aux autres plus abondantes. La table résultante, seqtab.nochim, contient uniquement les séquences considérées comme fiables pour les analyses ultérieures.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 61 bimeras out of 293 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]  20 232

``` r
##Cette commande compare l’abondance totale des séquences avant et après suppression des chimères, fournissant une mesure quantitative de la perte de données due au filtrage chimérique.

sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.9640374

## IX. Suivi des lectures tout au long du pipeline

Cette étape consiste à suivre l’évolution du nombre de lectures pour
chaque échantillon à chaque phase du pipeline DADA2 : après le filtrage,
le débruitage, la fusion des lectures et la suppression des chimères.
Cela permet de contrôler la qualité et l’efficacité du traitement,
d’identifier d’éventuelles pertes importantes de séquences et de
s’assurer que les données finales sont représentatives et fiables pour
les analyses ultérieures.

``` r
#Ce code crée un tableau de suivi qui montre le nombre de lectures pour chaque échantillon à chaque étape du pipeline DADA2 : lectures initiales, après filtrage, après débruitage avant et arrière, après fusion, et après suppression des chimères. La fonction getN calcule le nombre de lectures uniques, sapply l’applique aux différents objets DADA2, et rowSums(seqtab.nochim) donne les lectures finales non-chimériques. Ce suivi permet de visualiser les pertes de lectures et de vérifier la qualité des données finales.

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##        input filtered denoisedF denoisedR merged nonchim
    ## F3D0    7793     7113      6976      6979   6540    6528
    ## F3D1    5869     5299      5227      5239   5028    5017
    ## F3D141  5958     5463      5331      5357   4986    4863
    ## F3D142  3183     2914      2799      2830   2595    2521
    ## F3D143  3178     2941      2822      2868   2553    2519
    ## F3D144  4827     4312      4151      4228   3646    3507

Les séquences sont filtrées, corrigées, fusionnées et purifiées,

Chaque étape élimine le “bruit” et les erreurs,

Les valeurs finales (nonchim) représentent les lectures fiables qui
serviront à créer la table d’ASV.

## X. Attribution taxonomique

Une fois les séquences filtrées, débruitées, fusionnées et débarrassées
des chimères, il est possible d’attribuer une classification taxonomique
à chaque variante de séquence (ASV). Cette étape consiste à comparer les
ASV à une base de référence de séquences connues pour déterminer leur
règne, phylum, classe, ordre, famille, genre, voire espèce. DADA2
propose la fonction assignTaxonomy, qui utilise un classificateur
bayésien naïf pour fournir ces attributions avec un niveau de confiance,
ce qui permet d’interpréter biologiquement les données et de préparer
les analyses écologiques ou statistiques.

``` r
#À ce stade du pipeline, après avoir filtré, débruité, fusionné les lectures et supprimé les chimères, on dispose de nos ASV fiables dans seqtab.nochim. Cette ligne de code permet d’attribuer une taxonomie à chacune de ces séquences en les comparant à une base de référence (ici Silva). Le but est de savoir à quel règne, phylum, genre ou espèce correspond chaque ASV, ce qui est essentiel pour interpréter biologiquement les résultats et préparer les analyses écologiques ou statistiques. L’option multithread=TRUE accélère simplement le processus en utilisant plusieurs cœurs de calcul.

taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr99_v138.2_toSpecies_trainset.fa.gz?download=1", multithread=TRUE)
```

``` r
#Cette ligne complète l’attribution taxonomique en ajoutant, lorsque possible, le niveau espèce aux ASV déjà classifiées dans taxa, en comparant les séquences à une base de référence spécifique (silva_species_assignment_v132.fa.gz).

taxa <- addSpecies(taxa, "~/silva_species_assignment_v132.fa.gz")
```

``` r
##Crée une copie de la table taxonomique taxa, supprime les noms de lignes pour faciliter l’affichage, et montre les premières lignes.

taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum         Class         Order           Family          
    ## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      Genus         Species      Species
    ## [1,] NA            NA           NA     
    ## [2,] NA            NA           NA     
    ## [3,] NA            NA           NA     
    ## [4,] NA            NA           NA     
    ## [5,] "Bacteroides" "caecimuris" NA     
    ## [6,] NA            NA           NA

Toutes les séquences appartiennent au règne Bacteria, au phylum
Bacteroidota, et à la classe Bacteroidia :

Ce sont donc des bactéries anaérobies strictes, très fréquentes dans le
microbiote intestinal des mammifères.

L’ordre Bacteroidales regroupe plusieurs familles :

Muribaculaceae : bactéries communes dans le microbiote de la souris,
impliquées dans la dégradation des polysaccharides complexes.

Bacteroidaceae (avec le genre Bacteroides) : très abondantes dans le
microbiote intestinal humain et animal.

Le seul ASV identifié jusqu’à l’espèce est :

Bacteroides caecimuris C’est une espèce typique du microbiote
intestinal, souvent utilisée dans les modèles expérimentaux.

Les autres ASV restent non assignés au genre ou à l’espèce (NA) car la
similarité avec les bases de données n’était pas suffisante pour une
identification plus fine — ce qui est tout à fait normal en
métagénomique.

## XI. Evaluation de la précision

``` r
#seqtab.nochim["Mock",] récupère toutes les variantes de séquences (ASV) présentes dans l’échantillon contrôle « Mock ». sort(unqs.mock[unqs.mock>0], decreasing=TRUE) sélectionne uniquement les ASV effectivement présentes et les classe par abondance décroissante. La commande cat() affiche ensuite combien de séquences différentes DADA2 a identifiées dans ce contrôle, ce qui permet de vérifier la précision du pipeline sur un échantillon de référence.

unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) 
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

    ## DADA2 inferred 20 sample sequences present in the Mock community.

DADA2 a détecté 20 séquences différentes dans l’échantillon témoin
(“Mock community”), et que ces séquences correspondent aux 20 taxons
réellement présents dans cette communauté de référence.

``` r
#Cette étape évalue la proportion de séquences identifiées qui correspondent exactement aux séquences attendues dans l’échantillon témoin, offrant un indicateur direct de la performance du débruitage et de la suppression des erreurs/chimères.

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

    ## Of those, 20 were exact matches to the expected reference sequences.

parmi les 20 séquences détectées dans la communauté témoin (Mock), les
20 correspondent exactement aux séquences de référence attendues

Grâce au tutoriel DADA2, on a pu nettoyer et affiner mes données de
séquençage pour obtenir des ASV fiables et prêtes à être analysées. La
suite logique est maintenant d’explorer ces résultats avec phyloseq,
afin de mieux visualiser et comprendre la composition et la diversité
des communautés microbiennes.
