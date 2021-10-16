!======================================================================!
!=== VERSION 1
!=== date     : 30/11/2016
!=== author   : Christophe Mathé
!=== new      : dynamical allocation to matrix C et C0
!===            control on the method used (method_invert)
!===            inverse by Cholesky decomposition
!=== ------------------------------------------------------------------!
!=== VERSION 2.0
!=== date     : 09/12/2016
!=== objectif : passage en implicit none
!=== ------------------------------------------------------------------!
!=== VERSION 3.0
!=== date     : 14/12/2016
!=== objectif : restructuration et allocation dynamique
!===          - conv* : mise en format select case
!===          - suppression des boucles implicites
!===          - readline2
!=== ------------------------------------------------------------------!
!=== VERSION 3.1
!=== date     : 02/02/2017
!=== objectif : optimisation des allocations dynamiques
!===            si refprofT à l'envers => à l'endroit (sol:ciel)
!=== ------------------------------------------------------------------!
!=== VERSION 3.2
!=== date     : 16/06/2017
!=== objectif : ajout de idist = 4 (lecture de fichier en input)
!===          - remettre au propre les variables des codes
!===          - gérer les formats de lecture dynamiques et d'écriture
!======================================================================!
