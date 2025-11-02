# ENR_NRB

# Docker Image

## Build

`docker build -t asard .`
## Run

`docker run -d -v [host_path_to_work_dir]:[container_path_to_work_dir] -e FTP_USER=[user_for_the_dem_ftp] -e FTP_PASS=[password_for_the_dem_ftp] asard asard -c [path_to_config_file]`

- The [host_path_to_work_dir] path should contain all necessary inputs and output folders. 
  Paths written in the config file should resolve to the [container_path_to_work_dir].


