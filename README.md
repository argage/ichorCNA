# ichorCNA WDL
A revised version of [this WDL](https://github.com/GavinHaLab/wdl-repo/tree/main/processes/ichorCNA) so that it outputs plots.

## Instructions
This is just a rough explanation of what I've tested this with and how you can do it yourself, from a beginner's perspective.
You can run this on your own if you have [cromwell installed](https://github.com/broadinstitute/cromwell/releases/tag/85) and want to use the 'run' mode. To do this, enter this into your terminal:

    java -jar cromwell-XX.jar run your_workflow.wdl -i your_input.json -o your_configuration.conf

with `XX` being the version of cromwell you have. Make sure all of your files (WDL, input, options) 
