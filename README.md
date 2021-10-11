# long-reads
AMR detection in uncorrected long reads

Input is uncorrected(?)/corrected nanopore sequenced reads. 

Simulated reads: 

* https://github.com/bcgsc/NanoSim
* https://github.com/LomanLab/mockcommunity 

Just do Ryanwick for error correction (ML?)


## Method

locate amr genes. 

fetch long reads. 

consensus. 

## Output

output is the sytenic block of ARGs
mutuations in those genes. 

plasmid typing on reads associated with the ARGs. 

Plot options - Artemis looking output. 

## Previous & related tools: 

plasmid typing from long reads 
https://github.com/andrewjpage/tiptoft

Database of plasmid
https://ccb-microbe.cs.uni-saarland.de/plsdb/

plasmid typing
https://github.com/aldertzomer/RFPlasmid

plasmid detection + typing
mobsuite

## Output libraries 

* https://flanker.readthedocs.io/en/latest/
* https://github.com/gamcil/clinker
* https://cran.r-project.org/web/packages/gggenes/vignettes/introduction-to-gggenes.html
 
