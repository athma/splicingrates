# splicingrates

CODE FOR PAI et al 2018 eLIFE

1) Simulations to assess approachs to calculate rates of intron splicing from nascent 4sU-seq data (to be run over various parameters as inputs)

     get reads/value for all 3 methods: *simulation.R* <br>
     simulation of PSI decrease method (inputs of MISO summary PSI values run with output from simulation.R): *miso_modeling.R* <br>
     simulation of junction dynamics method (inputs are output from simulation.R): *junction_modeling.R* <br>


2) Code to calculate and analyze mRNA splicing rates from Drosophila 4sU-seq data
     
     *SplicingRates.Rmd*

3) Code to recreate figures in Pai et al. 2018
     
     *SplicingRates_figures.R*
