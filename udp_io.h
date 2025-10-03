#ifndef UDP_IO_H
#define UDP_IO_H

#include <stdint.h>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

// Standard MAVLink GCS/QGroundControl listening port
#define TARGET_PORT 14550
// Standard GCS IP address (change this if your GCS is on a different machine)
#define TARGET_IP "127.0.0.1" 

// Public function prototypes
int udp_init();
int udp_send_bytes(const uint8_t *buf, uint16_t len);
void udp_close();

#endif // UDP_IO_H