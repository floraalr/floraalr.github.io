# Analyses de données métabarcodes, Flora LARREUR

Les analyses de données métabarcodes permettent d’étudier la composition des communautés microbiennes à partir de séquences d’ADN extraites d’un échantillon environnemental. Le gène de l’ARN ribosomique 16S est souvent ciblé car il est présent chez toutes les bactéries et contient des régions conservées et variables, idéales pour l’identification taxonomique. Cette approche offre une vision globale de la diversité microbienne sans culture préalable. Elle est essentielle en écologie microbienne, santé ou agroalimentaire pour comparer la structure et la dynamique des microbiotes.

## Tutoriel de la pipeline DADA2

Vous trouverez ici le tutoriel DADA2 annoté et expliquéé à chaque étape via ce lien : [Mon Tutoriel](CC_DADA2_Flora_Larreur%20(1).md)

https://github.com/floraalr/floraalr.github.io/blob/main/CC_DADA2_Flora_Larreur%20(1).md

La pipeline DADA2 traite les séquences brutes issues du séquençage d’ADN microbien pour en extraire une image fidèle de la diversité réelle d’un échantillon. Elle commence par un contrôle de qualité rigoureux, éliminant les lectures trop courtes ou comportant des erreurs. Ensuite, elle modélise statistiquement les erreurs de séquençage propres à la machine afin de corriger les séquences plutôt que de simplement les regrouper. Cette étape permet d’identifier avec une très grande précision les ASV (Amplicon Sequence Variants), c’est-à-dire les séquences uniques réellement présentes dans l’échantillon. Contrairement aux anciennes méthodes basées sur les OTU, DADA2 offre une résolution quasi nucléotidique, essentielle pour distinguer finement les espèces ou les souches proches au sein des communautés microbiennes.

