
## File structure:
* main.c
  * Build UDP connection to local host
  * Run FDM generated from simulink
  * Pack heartbeat and attitude messages
  * Send via UDP
  * Qgroundstation can receive the messages and connect automatically.
* mavlink - MAVLink library
* UAV_Dynamics.c   UAV_Dynamics.h   
   FDM file generated from Simulink.  The input of FDM is 'UAV_Dynamics_U' and the output is UAV_Dynamics_Y
* rtw_continuous.h    rtw_solver.h     rtwtypes.h
   Headfile genrated/copied from Simulink 
* udp_io.c udp_io.h
  Code to handle UDP communicaiton to local host

## build and run
~~~
gcc -std=c99 main.c udp_io.c UAV_Dynamics.c     -I.     -I~/testcode/hellodemo/mavlink/common    -o fdm_publisher     -lm
./fdm_publisher
~~~


## for snap
* use 'snapcraft pack' to generate snap file
* use 'sudo snap install fdm-publisher-snap_1.0_amd64.snap --dangerous --devmode' to intall locally
* use 'fdm-publisher-snap.fdm-publisher' to run the application
* More actions to improve the yaml settings and publish to snap store.
