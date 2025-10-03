#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <stdint.h>
#include "mavlink/common/mavlink.h" 

// --- FDM MODEL HEADERS (Replace 'model_name' with your actual model name) ---
#include "UAV_Dynamics.h" 
// The FDM runs at a fixed rate, e.g., 100 Hz (10000 microseconds)
#define FDM_SAMPLE_TIME_US 10000 
// --------------------------------------------------------------------------

// --- MAVLINK HEADERS (Ensure this path is correct for your mavgen output) ---
#include "mavlink/common/mavlink.h" 
#include "udp_io.h"
// --------------------------------------------------------------------------

// Define System IDs
#define SYSTEM_ID 1    // The ID of this vehicle/system
#define COMPONENT_ID 1 // The ID of the FDM component (MAV_COMP_ID_AUTOPILOT1)

// External FDM I/O and data structure definitions (Assuming ERT target)
// These typically align with the root Inport/Outport blocks in your Simulink model


// Global state for MAVLink
static mavlink_system_t mavlink_system = {
    .sysid = SYSTEM_ID,
    .compid = COMPONENT_ID,
};

// Global sequence number for MAVLink packets
static uint8_t mavlink_seq = 0;

/**
 * @brief Publishes a MAVLink HEARTBEAT message.
 * This should run at ~1Hz to maintain the GCS connection.
 */
void mavlink_send_heartbeat() {
    mavlink_message_t msg;
    uint8_t buf[MAVLINK_MAX_PACKET_LEN];
    uint16_t len;

    // --- FDM Data to MAVLink mapping (Static/placeholder values) ---
    uint8_t type = MAV_TYPE_FIXED_WING; // Change based on your FDM vehicle
    uint8_t autopilot = MAV_AUTOPILOT_GENERIC;
    uint8_t base_mode = MAV_MODE_FLAG_SAFETY_ARMED | MAV_MODE_FLAG_CUSTOM_MODE_ENABLED;
    uint32_t custom_mode = 0;
    uint8_t system_status = MAV_STATE_ACTIVE;
    // ---------------------------------------------------------------

    mavlink_msg_heartbeat_pack(
        mavlink_system.sysid, mavlink_system.compid, &msg,
        type, autopilot, base_mode, custom_mode, system_status
    );

    len = mavlink_msg_to_send_buffer(buf, &msg);
    udp_send_bytes(buf, len);
}

/**
 * @brief Publishes a MAVLink ATTITUDE message using FDM outputs.
 * This should run at a higher rate (e.g., 50Hz).
 */
void mavlink_send_attitude() {
    mavlink_message_t msg;
    uint8_t buf[MAVLINK_MAX_PACKET_LEN];
    uint16_t len;
    uint64_t time_us = (uint64_t)clock() * 1000000L / CLOCKS_PER_SEC;

    // --- FDM Data to MAVLink mapping (Example using FDM outputs) ---
    float roll_rad = 0;    // Assuming FDM output exists
    float pitch_rad = 45;  // Assuming FDM output exists
    float yaw_rad = UAV_Dynamics_Y.Acc[2];      // Assuming FDM output exists
    // ---------------------------------------------------------------
    
    // Pack the ATTITUDE message
    mavlink_msg_attitude_pack(
        mavlink_system.sysid, mavlink_system.compid, &msg, 
        time_us, 
        roll_rad, pitch_rad, yaw_rad, 
        0.0f, 0.0f, 0.0f // Roll/Pitch/Yaw angular velocities (set to 0 for simplicity)
    );

    len = mavlink_msg_to_send_buffer(buf, &msg);
    udp_send_bytes(buf, len);
}

int main(void) {
    long long loop_counter = 0;

    printf("Starting FDM-MAVLink Application...\n");

    // 1. Initialize FDM Model
    UAV_Dynamics_initialize();
    printf("FDM Initialized.\n");

    // 2. Initialize UDP Communication
    if (udp_init() != 0) {
        fprintf(stderr, "Failed to initialize UDP. Exiting.\n");
        return EXIT_FAILURE;
    }

    // Main Simulation Loop
    while (1) {
        // --- A. Update FDM Inputs (e.g., control surfaces, throttle) ---
        UAV_Dynamics_U.PWMInputs[0] = 0.2;
        UAV_Dynamics_U.PWMInputs[1] = 0.2;
        UAV_Dynamics_U.PWMInputs[2] = 0.2;
        UAV_Dynamics_U.PWMInputs[3] = 0.2;
        // ... (Update other inputs)

        // --- B. Execute FDM Step ---
        UAV_Dynamics_step();

        // --- C. MAVLink Publishing ---
        // Heartbeat at 1 Hz (1000ms / 10ms step = 100 steps)
        if (loop_counter % 100 == 0) {
            mavlink_send_heartbeat();
            printf("HEARTBEAT Sent. Loop: %lld\n", loop_counter);
        }

        // Attitude at 50 Hz (20ms / 10ms step = 2 steps)
        if (loop_counter % 2 == 0) {
             mavlink_send_attitude();
        }
        
        loop_counter++;

        // --- D. Time Step Management ---
        usleep(FDM_SAMPLE_TIME_US); 
    }

    // Cleanup (This code is typically unreachable in an infinite loop)
    UAV_Dynamics_terminate();
    udp_close();
    return EXIT_SUCCESS;
}