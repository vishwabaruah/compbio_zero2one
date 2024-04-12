# How to create Dockerfile for your (developed) R/Shiny-app

You created a R/Shiny-app and now want to deploy it. The best method is to create a Dockerfile and use prebuilt R/Shiny docker images. 
That will reduce your time and hard work so that you can focus on building your app image.

-  Create a file named `Dockerfile`. Open it and write the commands to execute. Remove the comments if needed.

```
# Use this *FROM* syntax to use a prebuild (shiny) base image. Here I am using rocker/shiny and keeping the version same as the one I am using in my Shiny-app. The image takes care of the ubuntu and r package dependencies for shiny version also. 
FROM rocker/shiny:4.3.2

# Set the Environment for R repository
ENV RENV_CONFIG_REPOS_OVERRIDE https://packagemanager.rstudio.com/cran/latest

# Install any additional required dependencies that your app demands 
RUN apt-get update && apt-get install -y \
    python3 \
    make \
    zlib1g-dev
  
# OPTIONAL step: If your R/shiny-app requires additional r packages, download and install specific r packages and any other dependencies required to run your R/shiny_app 
# It is recommended to have version control over the r packages so that issues related to upgrades do not lead to a crash. that you install the correct verCopy your renv.lock to the Docker image. Prior to this create a shiny_renv.lock file and keep it in the same folder system as the Dockerfile. The shiny_renv.lock basically contain r packages and corresponding version. It is necessary for (other) r package version control. 
RUN Rscript -e "install.packages('devtools')"
RUN R -e "install.packages('rmarkdown')"

# Copy the R/Shiny `app` directory (preferably on the working directory) into the Docker image. Preferebly, keep the `app` in same directory as in the `Dockerfile`.
COPY ./app /srv/shiny-server/app

# Expose the port to `3838`. You may also define it in `app/global.R`. This version of rocker/shiny by default expose the shiny to `3838` of `localhost`
EXPOSE 3838

# Set the working directory to default shiny server. Shiny will search for the `index.html` or `server.R` file in this location.
WORKDIR /srv/shiny-server/

# Run the Shiny app
CMD ["/usr/bin/shiny-server"]                         
```

-   Modify/edit the `Dockerfile` as per your requirement.
-   Finally, save and close the `Dockerfile`.
-   Build and run docker image using `Dockerfile`
-   Open internet browser and go to this link - `[http://localhost:3838/app]`. If everything is correct than you should be able to see your R/Shiny running. 
  
Voila! your app is ready and running
