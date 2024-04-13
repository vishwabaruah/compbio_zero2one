## What is a Dockerfile

`Dockerfile` is a human-readable flat file (without extension), which uses very specific keywords that supervise specific docker image building process.

## Why use a `Dockerfile`

You created an app. You want to run it or you want to deploy it. The best (and arguably the simplest) way is to create a Dockerfile for your app. Using this `Dockerfile` you can build (create) a Docker image of your app and then run that image for the app to come alive.

## Pre-requisites

-   A text (visual) editor (preferably) e.g. *gedit*, *vim*, *VScode*, etc.
-   Docker installed in the computer where docker image is to be build.
-   Knowledge on non-interactive command-lines software(s) installation process in unix-based systems.

## Docker keywords that are used in `Dockerfile`

Use these docker commands keywords to build your `Dockerfile`. The specific keywords you can use in a file are:

| Keywords   | Function                                                                                           | Usage                                      |
|------------------|---------------------------|---------------------------|
| ADD        | Copies the files from a remote location via a URL into the destination directory of container’s filesystem. Recognizes and extract compressed files from local directories before copying to container. | `ADD app.gzip /root/myapp`                                      |
| ARG        | Defines the name of a parameter and its default value. | `ARG WEB=https://github.com/`                                      |
| CMD        | Can be used for executing a specific command within the container.                                 | `CMD ["bash"]`                             |
| COPY       | Copies the files from a source on the host into the container’s filesystem at the set destination.                                                                                                   | `COPY app_dir /home/john/new_app_dir`      |
| ENTRYPOINT | Sets a default application to be used every time a container is created with the image.            | `ENTRYPOINT javac Hello.java && java Hello`                               |
| ENV        | Sets environment variables.                                                                        | `ENV TZ UTC`                               |
| EXPOSE     | Associates a specific port to enable networking between the container and the outside world.       | `EXPOSE 8080`                              |
| FROM       | Defines the base image used to start the build process.                                            | `FROM ubuntu:jammy`                        |
| LABEL      | Allows you to add a label to your docker image.                                                    | `LABEL website="compbio_zer2one"`                                    |
| MAINTAINER | Defines the full name and email address of the image creator.                                      | `MAINTAINER Vishwa Jyoti Baruah "baruahvishwa@fakeorg.com"` |
| RUN        | Central executing directive for Dockerfiles, also known as the run Dockerfile command.      | `RUN apt-get update`                       |
| USER       | Sets the UID (or username) which is to run the container.                                          | `USER root`                                     |
| VOLUME     | Enables access from the container to a directory on the host machine.                    | `VOLUME /test`                                   |
| WORKDIR    | Sets the path where the command, defined with CMD, is to be executed.                              | `WORKDIR /src`                      |

Table: Available Keywords for creating a `Dockerfile`.

## Creating a basic `Dockerfile` for your app with keywords

Here is a step-by-step protocol for creating a basic `Dockerfile`.

-   Create a file named **Dockerfile**. Put it in the location where your app folder/directory is present
-   Make sure that file doesn't have an extension. Docker recognizes this file(name) and will execute the commands provided in it to create a docker image.
-   Open the `Dockerfile` and start writing/editing it.

```         
# Use *FROM* to use a prebuild base image if any required.
FROM Pre-built_base_image

# Optional: Provide the maintainers informations. 
MAINTAINER user_name "email@institute.edu"

# Set the Environment for running unattended scripts and disabeling interactive xprompts so thatinstallation does not get stuck at an interactive prompt.
ENV DEBIAN_FRONTEND=noninteractive

# Install any additional dependencies that your app requires.
RUN apt-get update && apt-get install -y \
    package1 \
    package2 \
    package3-dev
  
#   Set the working directory
WORKDIR /path/to/be/used/as/working/directory/

# Copy your app or/and other files from relative path of the target directory/file to full destination directory/file of the container’s filesystem.
COPY relative/path/of/target/app_dir /full/path/to/destination/app_dir

# Run any additional app preparation commands or setting configuration files.
RUN command_to_execute

# Expose the port to empty port like `3838` or `8080` or any specific port (if your app requires an empty port or if it has a html frontend).
EXPOSE 8000

# When all the parameters are set, provide the command to run the app
CMD ["npm", "start"]
```

-   You can remove the comments from the final `Dockerfile`.
-   Use the extended docker commands section to edit the `Dockerfile`.
-   Finally, save and close the `Dockerfile`.

Voila!! Your `Dockerfile` is ready and your **App** is ready to build image and run it. You should be able to see the app running at `[http://localhost:8000/]` 

## Demo

Here is a demo of creating a [`Dockerfile` for R/Shiny app](https://github.com/vishwabaruah/compbio_zero2hero/blob/main/tutorial/docker/creating_dockerfile/dockerfile4shiny.md).
Here is a demo of creating a [`Dockerfile` for node.js application](https://github.com/vishwabaruah/compbio_zero2hero/blob/main/tutorial/docker/creating_dockerfile/dockerfile4shiny.md).


