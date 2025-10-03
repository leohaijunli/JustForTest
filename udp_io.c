#include "udp_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int sockfd = -1;
static struct sockaddr_in target_addr;

/**
 * @brief Initializes the UDP socket connection.
 * @return 0 on success, -1 on failure.
 */
int udp_init() {
    // 1. Create UDP socket
    if ((sockfd = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP)) < 0) {
        perror("Error: socket creation failed");
        return -1;
    }

    // 2. Setup target (GCS) address structure
    memset(&target_addr, 0, sizeof(target_addr));
    target_addr.sin_family = AF_INET;
    target_addr.sin_port = htons(TARGET_PORT);

    // Convert IP address string to binary form
    if (inet_pton(AF_INET, TARGET_IP, &target_addr.sin_addr) <= 0) {
        perror("Error: Invalid address/Address not supported");
        close(sockfd);
        sockfd = -1;
        return -1;
    }

    printf("UDP Initialized: Sending to %s:%d\n", TARGET_IP, TARGET_PORT);
    return 0;
}

/**
 * @brief Sends a buffer of bytes over the initialized UDP socket.
 * This is the 'send_bytes' function required by MAVLink.
 * @param buf The raw MAVLink packet buffer.
 * @param len The length of the buffer.
 * @return Number of bytes sent, or -1 on error.
 */
int udp_send_bytes(const uint8_t *buf, uint16_t len) {
    if (sockfd < 0) {
        fprintf(stderr, "Error: UDP socket not initialized.\n");
        return -1;
    }
    
    // Use sendto to dispatch the packet
    int bytes_sent = sendto(sockfd, buf, len, 0, 
                            (struct sockaddr *)&target_addr, 
                            sizeof(target_addr));
    
    if (bytes_sent < 0) {
        perror("Error: sendto failed");
    }
    return bytes_sent;
}

/**
 * @brief Cleans up the UDP socket.
 */
void udp_close() {
    if (sockfd >= 0) {
        close(sockfd);
        sockfd = -1;
        printf("UDP socket closed.\n");
    }
}