# Use the official Debian Buster base image
FROM r-base:latest

## Install additional R packages
RUN R -e "install.packages(c('BiocManager', 'readxl', 'ggplot2', 'shiny', 'cowplot'), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('HDCytoData', 'CATALYST', 'FlowSOM', 'ConsensusClusterPlus', 'ComplexHeatmap', 'scater', 'diffcyt'), force=TRUE)"

# Copy your R script into the container
COPY Cytof2.R /usr/src/app/

# Expose port 3838 for Shiny app
EXPOSE 3838

# Run your script as a Shiny app
CMD ["R", "-e", "shiny::runApp('/usr/src/app/Cytof2.R', host = '0.0.0.0', port = 3838)"]
