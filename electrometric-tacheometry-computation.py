# Istanbul Technical University
# Department of Geomatics Engineering
# FUNDAMENTAL OF PROGRAMMING - TERM PROJECT
# Part 2: Electrometric Tacheometry Computation
# StudentID_010930605_2.py

import math
import sys
import os

# Clearing Terminal Screen
if sys.platform.startswith('win'):
    os.system('cls')
else:
    os.system('clear')

# Constants
GRAD_TO_RAD_0605 = math.pi / 200
RAD_TO_GRAD_0605 = 200 / math.pi

print("This code for Electrometric Tacheometry Computation")
print("*" * 50)

# Given data from the image
station_traverse_ID_0605 = input("Enter the stationary traverse ID:")
ref_traverse_ID_0605 = input("Enter the referenced traverse ID:")

ref_traverse_Y_0605 = float(input(f"Enter the Y coordinates of {station_traverse_ID_0605} (m):"))
ref_traverse_X_0605 = float(input(f"Enter the X coordinates of {station_traverse_ID_0605} (m):"))

station_traverse_Y_0605 = float(input(f"Enter the Y coordinates of {ref_traverse_ID_0605} (m):"))
station_traverse_X_0605 = float(input(f"Enter the X coordinates of {ref_traverse_ID_0605} (m):"))

station_height_0605 = float(input(f"Enter the Height of {station_traverse_ID_0605} (m):"))

# Detail point measurements
detail_point_ID_0605 = input("Enter the point ID of detail point:")
detail_point_hor_dir_0605 = float(input(f"Enter the horizontal direction of point {detail_point_ID_0605} (grad):"))
detail_point_vert_angle_0605 = float(input(f"Enter the vertical angle of point {detail_point_ID_0605} (grad):"))
slope_dist_0605 = float(input(f"Enter the slope distance between {station_traverse_ID_0605} and {detail_point_ID_0605} (m):"))

instrument_height_0605 = float(input("Enter the height of the instrument (m):"))
reflector_height_0605 = float(input("Enter the height of the reflector (m):"))

# --- Calculations ---

# Compute delta Y and delta X
delta_Y_0605 = station_traverse_Y_0605 - ref_traverse_Y_0605
delta_X_0605 = station_traverse_X_0605 - ref_traverse_X_0605

# Compute azimuth angle
if delta_X_0605 == 0:
    azimuth_angle_0605 = 100.0 if delta_Y_0605 > 0 else 300.0 if delta_Y_0605 < 0 else 0.0
else:
    azimuth_angle_0605 = math.atan(abs(delta_Y_0605 / delta_X_0605)) * RAD_TO_GRAD_0605

    # Determine quadrant
    if delta_Y_0605 > 0 and delta_X_0605 < 0:
        azimuth_angle_0605 = 200 - azimuth_angle_0605
    elif delta_Y_0605 < 0 and delta_X_0605 < 0:
        azimuth_angle_0605 += 200
    elif delta_Y_0605 < 0 and delta_X_0605 > 0:
        azimuth_angle_0605 = 400 - azimuth_angle_0605

# Compute horizontal distance
horizontal_distance_0605 = slope_dist_0605 * math.sin(detail_point_vert_angle_0605 * GRAD_TO_RAD_0605)

# Compute height difference
# Formula from image calculation block (slope_dist * cos(vertical_angle) + (reflector_height - instrument_height))
height_difference_0605 = (slope_dist_0605 * math.cos(detail_point_vert_angle_0605 * GRAD_TO_RAD_0605)) + (reflector_height_0605 - instrument_height_0605)

# Compute elevation of detail point
elevation_detail_point_0605 = station_height_0605 + height_difference_0605

# Compute total direction (azimuth to detail point)
# This is typically the sum of the initial azimuth and the horizontal direction read by the instrument,
# adjusted to be within the 0-400 grad range.
total_direction_0605 = (azimuth_angle_0605 + detail_point_hor_dir_0605) % 400

# Compute X and Y coordinates of the detail point
# These formulas assume the horizontal distance and total direction are relative to the station.
X_detail_point_0605 = station_traverse_X_0605 + (horizontal_distance_0605 * math.cos(total_direction_0605 * GRAD_TO_RAD_0605))
Y_detail_point_0605 = station_traverse_Y_0605 + (horizontal_distance_0605 * math.sin(total_direction_0605 * GRAD_TO_RAD_0605))

# --- Formatting for the table ---
horizontal_distance_formatted_0605 = f"{horizontal_distance_0605:.2f}"
height_difference_formatted_0605 = f"{height_difference_0605:.3f}"
elevation_detail_point_formatted_0605 = f"{elevation_detail_point_0605:.3f}"
Y_detail_point_formatted_0605 = f"{Y_detail_point_0605:.3f}"
X_detail_point_formatted_0605 = f"{X_detail_point_0605:.3f}"

# --- Table formatting with proper column widths ---
print("\n") # Add a newline for separation

# Define column widths
col_widths = [12, 15, 24, 10, 13, 14, 14]
total_width = sum(col_widths) + len(col_widths) + 1  # +1 for each separator + 1 for final border

print("-" * total_width)
print(f"|{'Electronic Tacheometry Calculation Table':^{total_width-2}}|")
print("-" * total_width)

# Header row with proper spacing
header = (
    f"| {'Station ID':^{col_widths[0]-1}}"
    f"| {'Target ID':^{col_widths[1]-1}}"
    f"| {'Horizontal Distance(m)':^{col_widths[2]-1}}"
    f"| {'ΔH(m)':^{col_widths[3]-1}}"
    f"| {'Elevation(m)':^{col_widths[4]-1}}"
    f"| {'Y(m)':^{col_widths[5]-1}}"
    f"| {'X(m)':^{col_widths[6]-1}}|"
)
print(header)
print("-" * total_width)

# Data row with proper spacing
data = (
    f"| {station_traverse_ID_0605:^{col_widths[0]-1}}"
    f"| {detail_point_ID_0605:^{col_widths[1]-1}}"
    f"| {horizontal_distance_formatted_0605:^{col_widths[2]-1}}"
    f"| {height_difference_formatted_0605:^{col_widths[3]-1}}"
    f"| {elevation_detail_point_formatted_0605:^{col_widths[4]-1}}"
    f"| {Y_detail_point_formatted_0605:^{col_widths[5]-1}}"
    f"| {X_detail_point_formatted_0605:^{col_widths[6]-1}}|"
)
print(data)
print("-" * total_width)
