#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <stdint.h>
#include "mavlink/common/mavlink.h" 
#include <netinet/in.h>
#include <arpa/inet.h>
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

//Structure to hold static vehicle data
typedef struct{
    uint32_t icao_address;
    int32_t latitude;     // degE7
    int32_t longitude;    // degE7
    int32_t altitude_msl; // mm
    char callsign[9];
} vehicle_t;

// Stationary Vehicles Data
const vehicle_t vehicles[] = {
    {0x100001, 486801413, -1233990053, 500000, "UAVSIM03"}, // Los Angeles (Example)
    {0x100002, 486102410, -1233990050, 500000, "UAVSIM01"}, // Los Angeles (Example)
    {0x100003, 486403412, -1233990052, 500000, "UAVSIM02"}, // Los Angeles (Example)

};
#define NUM_VEHICLES (sizeof(vehicles) / sizeof(vehicles[0]))
uint16_t adsb_flags_my = ADSB_FLAGS_VALID_COORDS | 
                         ADSB_FLAGS_VALID_ALTITUDE | 
                         ADSB_FLAGS_VALID_HEADING | 
                         ADSB_FLAGS_VALID_VELOCITY | 
                         ADSB_FLAGS_VALID_CALLSIGN | 
                         ADSB_FLAGS_SIMULATED;

void send_adsb_vehicle(const vehicle_t *vehicle) {
    mavlink_message_t msg;
    uint8_t buf[MAVLINK_MAX_PACKET_LEN];

    // Encode the ADSB_VEHICLE message
    mavlink_msg_adsb_vehicle_pack(
        SYSTEM_ID, 
        COMPONENT_ID, 
        &msg, 
        vehicle->icao_address, // ICAO Address
        vehicle->latitude,     // Latitude in degE7
        vehicle->longitude,    // Longitude in degE7
        ADSB_ALTITUDE_TYPE_PRESSURE_QNH, // Altitude type
        vehicle->altitude_msl, // Altitude in mm (MSL)
        5,                     // Heading (0 for stationary)
        1000000,                     // Horizontal Velocity (cm/s)
        0,                     // Vertical Velocity (cm/s)
        vehicle->callsign,     // Callsign
        MAV_TYPE_QUADROTOR, // Emitter type
        0,                     // Squawk (0)
        adsb_flags_my, // Flags
        0                      // Time since last communication (s)
    );

    uint16_t len = mavlink_msg_to_send_buffer(buf, &msg);
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
            for (int i = 0; i < 3; i++) {
            // Note: Mavlink pack functions automatically update the sequence number 'g_sequence'
                send_adsb_vehicle(&vehicles[i]);
            }
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