# Installation

### Install requirements

For running the Cognitive Tomography inference algorithms, you will need a Python environment. From the repository's root folder, open a shell and run

`pip install -r cogtom_requirements.txt`

Installation should normally take only a few minutes.

### Compile stan models

This step may take from a few minutes to about 10 minutes depending on the computer.

`cd Python && python compile_stan.py`



Now you are ready to run model inference on data. This may take a few minutes depending on the settings in the ini file. The `nums` value governs the number of samples. The demo should run in 10 minutes maximum. However, the runs used in the manuscript, which include 850 trials and take 1600 samples may take 18-48 hrs.

`python cogtom.py <ct_parameters.ini> <asrt_data.pkl>`

For example with the included data sets:

`python cogtom.py demo.ini artificial_demo_ASRT_0_blocks_1_4.pkl`

To run the Markov model run (this is needed to have all analyses in later sections). The Markov model runs much faster, should finish within a few minutes.

`python cogtom.py markov demo.ini artificial_demo_ASRT_0_blocks_1_4.pkl`

The output of these are iHMM model samples put into `Python/Output`



### Run analyses

Before running the analyses, you need to first run `python organize.py` which puts the model files in the correct folders.

To run analyses, you need a template file. An example can be found at `Data/demo.template`

Then you need to select which commit's models to run using which analysis template. Example:

`python analysis.py -c NOGTHSH -t demo -p1`

**Output** 
After the analysis has run, the console should state where the saved outputs are located.



### Running on human data

First, you need to download the human dataset and unzip it into the Data folder.

Then you need to convert the raw experiment data into the input structure of the task.

Run `python data_handling.py --convert` from the Python folder.

It will ask you for which experiment to take: `elarasztas`

Participants: e.g. `102,103`

Blocks_from: e.g. `14`

Blocks_to: e.g. `16`

Then it generates a .pkl file that you can use with CT model inference.

With the above setting you will be able to run

`python cogtom.py demo.ini elarasztas_102_blocks_14_16.pkl`



## System requirements

Linux: Ubuntu 16.04+

Python 3.7+

Conda version: 4.7.12+

Required python packages are included in `cogtom_requirements.txt`

For installation see the "Installation" section.



## Experiment

The Eprime2 files for the experiment can be found in the `Experiment` folder.


## Generate figures
In the following, I will assume you are located in the root directory of this repository.
1. Download dataset from [figshare](https://figshare.com/articles/dataset/Tracking_the_contribution_of_inductive_bias_to_individualized_internal_models_--_figure_datasets/19620285) into `Data` folder.
1. Checkout this repository.
1. Create output folder for figures: `mkdir Figures`
1. Run the following docker container with a previously installed R environment and set memory limit to 3G (if docker on your machine is setup to limit 2GB of memory for each container, you need to extend that.)
    ```
    docker run --rm -it \
        -v $(pwd)/R:/asrt-beamsampling/R \
        -v $(pwd)/Data:/data \
        -v $(pwd)/Figures:/asrt-beamsampling/Figures/figlist \
        --memory="3g" \
        -e DATA_DIR=/data \
        mzperix/asrt-beamsampling-r:0.1.0
    ```