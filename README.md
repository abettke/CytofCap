## Build

`docker build -t cyto2 .`
* `-t cyto2` will name the image `cyto2` (you can name it whatever you want)
* `.` is the path to the Dockerfile. This points to the current directory.

## Run

`docker run -p 3838:3838 cyto2`
* `-p 3838:3838` maps the port 3838 from the container to the port 3838 on your machine