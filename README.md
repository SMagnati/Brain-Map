# Brain-Map 🧠🗺
_A survey on the expression of the UPS components Hect-, RBR-E3 Ub-ligase, E2 Ub-conjugating, and E1 Ub-activating enzymes during the human brain development_ (link)

### Introduction
Welcome to our repository! 
For the first time 🥇, we embarked on a journey to map the human brain gene expression 🧬 of specific proteins called 'ubiquitins' during development using data analysis and machine learning techniques. This endeavor aims to shed light on the intricate mechanisms underlying brain development and the role of ubiquitins in this process. Through comprehensive analysis and advanced methodologies, we strive to unravel the complexities of gene expression dynamics in the human brain. 

### Data Collection 
We utilized a publicly available dataset from open-access repositories: _Brain Span_ (https://www.brainspan.org/). Brain Span provides comprehensive gene expression datasets generated using mRNA sequencing technology. The dataset comprises 504 observations and features such as age group, brain regions, gender, and ethnicity of donors. The values are in log2 RPKM (reads per kilobase per million).

### Feature Engineering
The Brain Span dataset includes expression data for 89 genes, categorized into nine representative age groups: early prenatal, early mid-prenatal, late-mid prenatal, late prenatal, infancy, early childhood, late childhood, adolescence, and adulthood. These genes consist of nine E1, thirty-seven E2, fourteen E3-RBR, and twenty-nine E3-HECT ubiquitin components.
The 'ethnicity' and 'gender' variables were manually inserted following the donor ID as indicated in the Brain Span documentation found in the 'Brain Span Documentation' file. If you wish to download genes of interest from Brain Span, you can copy and paste the 'gender' and 'ethnicity' columns.

### Language and Envoironment
All analyses and implementation were conducted using Python (www.python.org). We recommend using Jupyter notebook for running the codes.

### Functions and libraries
The 'Magnati_functions.py' file contains all the functions codes and packages used in the main files described below.

### Main files
    0.Import_Data 📨
    1.EDA_Stats_and_Cluster_Analysis.ipynb 📊📚
    2.Volcano_Plot.ipynb 🌋
    3.Machine_Learning_Classification.ipynb 💻
Within each main file, you will find code cells ready for use. We have ensured that the code is easy to understand and highly flexible. It is recommended to execute one code cell at a time. If a cell does not finish its run, it could be due to excessive computational calculation.


### Citation 📝
If you find this repository useful for your work, please consider citing:

[Insert citation details here]

### Acknowledgement 🙏🏻
We would like to thank the medical student Giada Gaveglio 👩🏻‍⚕️ for her dedication and assistance in developing the pipeline figure. Thank you!
It is also necessary to thank Brain Span for this great data collection and OpenAI (https://openai.com/) and its AI language, ChatGPT (version 3.5), which facilitated the development of various parts of the code. 

Last but not least, thank you to all for your interest in our research! 😊
                                                                                                                            
__Stefano Magnati, Eugenio Alladio and Professor Enrico Bracco__
