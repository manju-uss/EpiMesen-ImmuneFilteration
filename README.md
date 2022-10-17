# EpiMesen-ImmuneFilteration

<p align="center">
  <img src="https://user-images.githubusercontent.com/28807444/146238344-07163b85-5450-49ac-ab7e-7d5c89cb03c7.jpg" align="left" width="35%" height='35%'/>
</p>

This script categorize cancer tissues/samples into epithelial and mesenchymal state using normalized RNAseq dataset downloaded from http://gdac.broadinstitute.org. Categorization is based on the expression of CDH1, DSP, TJP1, FZR1 genes for epithelial and VIM, CDH2, FOXC2, SNAI1, SNAI2, TWIST1, GSC, FN1, ITBG6, MMP2, MMP3, MMP9 and SOX10 genes for mesenchymal tissues. In the futher step, it explore immune-infiltration markers, infilteration of CD8+ cell and expression of immunogenic checkpoints in epithelial and mesechymal tissues. Young Kwang Chae et al (2018) (https://www.nature.com/articles/s41598-018-21061-1). 


# Citation

Tiwari J, Negi S, Kashyap M, Nizamuddin S, Singh A, Khattri A. Pan-cancer analysis shows enrichment of macrophages, overexpression of checkpoint molecules, inhibitory cytokines and immune exhaustation signatures in EMT-high tumors. Frontiers in oncology 2021
https://www.frontiersin.org/articles/10.3389/fonc.2021.793881/full


REQUIREMENTS:

        R
        This script will automatically install other requirements 
        in its first run if computer/laptop connected with internet.
        
USAGE:

    Rscript run.R [options]
    
Options:

        -f CHARACTER, --file=CHARACTER
                input file

        --output=CHARACTER
                output folder

        --script_folder=CHARACTER
                location of folder 'EpiMesen-ImmuneFilteration'
                
        -h, --help
                Show this help message and exit
                

