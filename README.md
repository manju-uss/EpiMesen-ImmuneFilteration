# EpiMesen-ImmuneFilteration

This script categorize cancer tissues/samples into epithelial and mesenchymal state using normalized RNAseq dataset (data_RNA_Seq_v2_expression_median.txt) downloaded from https://www.cbioportal.org/. Categorization is based on the expression of CDH1, DSP, TJP1, FZR1 genes for epithelial and VIM, CDH2, FOXC2, SNAI1, SNAI2, TWIST1, GSC, FN1, ITBG6, MMP2, MMP3, MMP9 and SOX10 genes for mesenchymal tissues. In the futher step, it explore immune-infiltration markers, infilteration of CD8+ cell and expression of immunogenic checkpoints in epithelial and mesechymal tissues.


REQUIREMENTS:

        R, this script will automatically install requirements in its first run, if computer/laptop connected with internet.
        
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
                
