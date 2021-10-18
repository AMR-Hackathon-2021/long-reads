[![AMR Long Reads](docs/amr-long-reads.png)](https://github.com/AMR-Hackathon-2021/long-reads#readme)

# long-reads

AMR detection in uncorrected long reads. Input is uncorrected/corrected nanopore sequenced reads. 

See [Docs](docs/README.md)

## A possible approach

* Locate AMR genes. 

* Fetch long reads. 

* Consensus. 

## Output

* Sytenic block of ARGs, mutuations in those genes. 

* Plasmid typing on reads associated with the ARGs. 

* Plot options - Artemis looking output. 

* Circularity of reads with detected AMR gene?

## Previous & related tools: 

* plasmid typing from long reads: [andrewjpage/tiptoft](https://github.com/andrewjpage/tiptoft)

* Database of plasmid: [plsdb](https://ccb-microbe.cs.uni-saarland.de/plsdb/)

* Database of plasmid: [itsmeludo/COMPASS](https://github.com/itsmeludo/COMPASS)

* Plasmid typing: [aldertzomer/RFPlasmid](https://github.com/aldertzomer/RFPlasmid)

* Plasmid classifier: [santirdnd/COPLA](https://github.com/santirdnd/COPLA)

* Plasmid detection + typing: [MobSuite](mobsuite)

## Output libraries 

* https://flanker.readthedocs.io/en/latest/
* https://github.com/gamcil/clinker
* https://cran.r-project.org/web/packages/gggenes/vignettes/introduction-to-gggenes.html
 


