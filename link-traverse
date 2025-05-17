# Istanbul Technical University
# Department of Geomatics Engineering
# FUNDAMENTAL OF PROGRAMMING - TERM PROJECT
# Part 1: Link Traverse Computation
# StudentID_010930605_1.py 

import math
import sys 
import os

#Clearing Terminal Screen 
print("\033c", end="")  # ANSI escape sequence for clearing the screen

# Execute the 'cls' command in Windows Command Prompt
os.system('cls') #clear the console screen
os.system("clear")  # Works on Mac and Linux


STUDENT_ID_SUFFIX = "0605" #to use in the output file name

sample_data_1234 = {
    "id_p1": "C2001", "y_p1": 417545.886, "x_p1": 4552506.610, 
    "id_p2": "C2002", "y_p2": 417740.231, "x_p2": 4552274.937,
    "id_p3": "C2003", "y_p3": 418245.208, "x_p3": 4553186.147, 
    "id_p4": "C2004", "y_p4": 418574.862, "x_p4": 4553238.143,
    "num_unknown_points": 3,
    "unknown_point_ids": ["C3001", "C3002", "C3003"],
    "measured_angles": [103.1912, 128.4724, 302.3489, 105.0024, 295.4614],
    "measured_distances": [418.825, 333.576, 265.547, 289.230]
}


# Function to convert grads to radians
def grads_to_radians_1234(grads_1234):
    return grads_1234 * (math.pi / 200.0)

# Function to convert radians to grads
def radians_to_grads_1234(rads_1234):
    return rads_1234 * (200.0 / math.pi)

# Function to calculate azimuth in grads using atan2() function in python
'''
def calculate_azimuth_1234(y1_1234, x1_1234, y2_1234, x2_1234):

    delta_y_1234 = y2_1234 - y1_1234
    delta_x_1234 = x2_1234 - x1_1234
    # Check for zero division 
    if abs(delta_x_1234) < 1e-9 and abs(delta_y_1234) < 1e-9: return 0.0 
    azimuth_rad_1234 = math.atan2(delta_y_1234, delta_x_1234)
    azimuth_grads_1234 = radians_to_grads_1234(azimuth_rad_1234)
    if azimuth_grads_1234 < 0:
        azimuth_grads_1234 += 400.0
    return azimuth_grads_1234

'''
# Function to calculate delta Y

def cal_delta_y_1234(y2_1234, y1_1234):
    return y2_1234 - y1_1234

# Function to calculate delta X

def cal_delta_x_1234(x2_1234, x1_1234):
    return x2_1234 - x1_1234

# Function to calculate azimuth in grads using atan() function in python

# Function to calculate azimuth in grads

def calculate_azimuth_1234(y1_1234, x1_1234, y2_1234, x2_1234):
    
    delta_y_1234 = cal_delta_y_1234(y2_1234, y1_1234)
    delta_x_1234 = cal_delta_x_1234(x2_1234, x1_1234)
    

    if abs(delta_y_1234) < 1e-9 and abs(delta_x_1234) < 1e-9: return 0.0
    
    elif delta_y_1234 > 0 and delta_x_1234 > 0:
        azimuth = math.atan(abs(delta_y_1234) / abs(delta_x_1234))
        azimuth_grads_1234 = azimuth * (200 / math.pi)

    elif delta_y_1234 > 0 and delta_x_1234 < 0:
        azimuth = math.atan(abs(delta_y_1234) / abs(delta_x_1234))
        azimuth_grads_1234 = 200 - (azimuth * (200 / math.pi))

    elif delta_y_1234 < 0 and delta_x_1234 < 0:
        azimuth = math.atan(abs(delta_y_1234) / abs(delta_x_1234))
        azimuth_grads_1234 = 200 + (azimuth * (200 / math.pi))

    elif delta_y_1234 < 0 and delta_x_1234 > 0:
        azimuth = math.atan(abs(delta_y_1234) / abs(delta_x_1234))
        azimuth_grads_1234 = 400 - (azimuth * (200 / math.pi))

    elif delta_y_1234 == 0 and delta_x_1234 > 0:
        azimuth_grads_1234 = 0 

    elif delta_y_1234 == 0 and delta_x_1234 < 0:
        azimuth_grads_1234 = 200

    elif delta_y_1234 > 0 and delta_x_1234 == 0:
        azimuth_grads_1234 = 100

    elif delta_y_1234 < 0 and delta_x_1234 == 0:
        azimuth_grads_1234 = 300

    elif delta_y_1234 == 0 and delta_x_1234 == 0:
        azimuth_grads_1234 = 0

    return azimuth_grads_1234   
    

# Function to calculate distance using Pythagorean theorem

def calculate_distance_1234(y1_1234, x1_1234, y2_1234, x2_1234):
    return math.sqrt((y2_1234 - y1_1234)**2 + (x2_1234 - x1_1234)**2)


# Tee class to write to multiple files
class Tee(object):
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f_obj in self.files: 
            f_obj.write(obj)
            f_obj.flush() 
    def flush(self):
        for f_obj in self.files: 
            f_obj.flush()

def main_1234():
    
    original_stdout_1234 = sys.stdout 
    output_file_path_1234 = None
    file_out_1234 = None 

    save_to_file_1234 = input("Do you want to save the output to a file? (y/n): ").strip().lower()
    if save_to_file_1234 == 'y':
        default_folder_1234 = os.path.join(os.getcwd(), "TraverseResults") 
        default_filename_1234 = f"LinkTraverse_Output_{STUDENT_ID_SUFFIX}.txt"
        
        folder_path_1234 = input(f"Enter folder path to save output (default: {default_folder_1234}): ").strip()
        if not folder_path_1234: folder_path_1234 = default_folder_1234
        
        filename_1234 = input(f"Enter filename (default: {default_filename_1234}): ").strip()
        if not filename_1234: filename_1234 = default_filename_1234

        if not os.path.exists(folder_path_1234):
            try:
                os.makedirs(folder_path_1234)
            except OSError as e:
                original_stdout_1234.write(f"Error creating folder {folder_path_1234}: {e}. Saving to current directory instead.\n")
                folder_path_1234 = os.getcwd() 
        
        output_file_path_1234 = os.path.join(folder_path_1234, filename_1234)
        
        try:
            file_out_1234 = open(output_file_path_1234, 'w', encoding='utf-8')
            sys.stdout = Tee(original_stdout_1234, file_out_1234) 
            print(f"Output will also be saved to: {output_file_path_1234}\n")
        except Exception as e:
            sys.stdout = original_stdout_1234 
            print(f"Error opening file for writing: {e}. Output will only be to console.")
            output_file_path_1234 = None 
            if file_out_1234: file_out_1234.close() 
            file_out_1234 = None
    
    try:
        print("This code calculates the coordinates in a link traverse series.")
        print("******************************************************************")
        
        load_sample_1234 = input("Do you want to load sample data? (y/n): ").strip().lower()

        if load_sample_1234 == 'y':
            data_source = sample_data_1234
            id_p1_1234 = data_source["id_p1"]
            y_p1_1234 = data_source["y_p1"]
            x_p1_1234 = data_source["x_p1"]
            id_p2_1234 = data_source["id_p2"]
            y_p2_1234 = data_source["y_p2"]
            x_p2_1234 = data_source["x_p2"]
            id_p3_1234 = data_source["id_p3"] 
            y_p3_1234 = data_source["y_p3"]
            x_p3_1234 = data_source["x_p3"]
            id_p4_1234 = data_source["id_p4"] 
            y_p4_1234 = data_source["y_p4"]
            x_p4_1234 = data_source["x_p4"]
            num_unknown_points_1234 = data_source["num_unknown_points"]
            point_ids_chain_for_angles_1234 = [id_p2_1234] + data_source["unknown_point_ids"] + [id_p3_1234]
            point_ids_for_legs_1234 = [id_p2_1234] + data_source["unknown_point_ids"] + [id_p3_1234]
            measured_angles_1234 = data_source["measured_angles"]
            measured_distances_1234 = data_source["measured_distances"]
        else:
            print("\nProceeding with manual data entry...")
            print("\nEnter details for the STARTING known control points:")
            id_p1_1234 = input(f"Enter the point ID of first known point (e.g., C2001): ")
            y_p1_1234 = float(input(f"Enter the Y coordinate of {id_p1_1234} (m): ")) 
            x_p1_1234 = float(input(f"Enter the X coordinate of {id_p1_1234} (m): "))
            
            id_p2_1234 = input(f"Enter the point ID of second known point (start of traverse, e.g., C2002): ")
            y_p2_1234 = float(input(f"Enter the Y coordinate of {id_p2_1234} (m): "))
            x_p2_1234 = float(input(f"Enter the X coordinate of {id_p2_1234} (m): "))

            print("\nEnter details for the ENDING known control points:")
            id_p3_1234 = input(f"Enter the point ID of third known point (end of traverse, e.g., C2003): ")
            y_p3_1234 = float(input(f"Enter the Y coordinate of {id_p3_1234} (m): "))
            x_p3_1234 = float(input(f"Enter the X coordinate of {id_p3_1234} (m): "))

            id_p4_1234 = input(f"Enter the point ID of fourth known point (e.g., C2004): ")
            y_p4_1234 = float(input(f"Enter the Y coordinate of {id_p4_1234} (m): "))
            x_p4_1234 = float(input(f"Enter the X coordinate of {id_p4_1234} (m): "))

            num_unknown_points_1234 = int(input("\nEnter the number of unknown traverse points: "))
            
            unknown_point_ids_list_1234 = [] # list to store unknown point IDs
            for i in range(num_unknown_points_1234): # loop to get unknown point IDs
                uid_1234 = input(f"Enter the point ID of unknown point {i + 1}: ") # e.g., C3001 , i+1 to get exact number
                if uid_1234 in unknown_point_ids_list_1234: # check if the ID already exists
                    print(f"Error: Point ID {uid_1234} already exists. Please enter a unique ID.")
                    uid_1234 = input(f"Enter a unique point ID for unknown point {i + 1}: ")
                unknown_point_ids_list_1234.append(uid_1234) # add to the list of the unknown point IDs
            
            # construct the point IDs for angles and legs
            point_ids_chain_for_angles_1234 = [id_p2_1234] + unknown_point_ids_list_1234 + [id_p3_1234]
            point_ids_for_legs_1234 = [id_p2_1234] + unknown_point_ids_list_1234 + [id_p3_1234]

            # create empty lists to store measured angles and distances
            measured_angles_1234 = []
            measured_distances_1234 = []
                                 
            #for i_1234 in range(num_unknown_points_1234 + 2):
            for i_1234 in range(len(point_ids_chain_for_angles_1234)): # 5 points reason e.g C2002, C3001, C3002, C3003, C2003
                angle_1234 = float(input(f"Enter the traverse angle of {point_ids_chain_for_angles_1234[i_1234]} (grad): "))
                measured_angles_1234.append(angle_1234)

            #for i_1234 in range(num_unknown_points_1234 + 1):
            for i_1234 in range(len(point_ids_for_legs_1234) - 1):   # 5 points 4 legs reason for -1
                dist_1234 = float(input(f"Enter the horizontal distance between {point_ids_for_legs_1234[i_1234]} and {point_ids_for_legs_1234[i_1234+1]} (m): "))
                measured_distances_1234.append(dist_1234)
        
        # Calculate azimuth and distance for the known first and last legs
        t_P1P2_1234 = calculate_azimuth_1234(y_p1_1234, x_p1_1234, y_p2_1234, x_p2_1234) 
        s_P1P2_1234 = calculate_distance_1234(y_p1_1234, x_p1_1234, y_p2_1234, x_p2_1234)
        t_P3P4_1234 = calculate_azimuth_1234(y_p3_1234, x_p3_1234, y_p4_1234, x_p4_1234) 
        s_P3P4_1234 = calculate_distance_1234(y_p3_1234, x_p3_1234, y_p4_1234, x_p4_1234)

        num_traverse_legs_1234 = num_unknown_points_1234 + 1
        num_measured_angles_1234 = len(measured_angles_1234)
        sum_beta_1234 = sum(measured_angles_1234) 
        
        raw_leg_azimuths_1234 = []
        current_raw_azimuth_1234 = t_P1P2_1234
        
        for i in range(num_traverse_legs_1234):
            raw_leg_azimuth_1234 = (current_raw_azimuth_1234 + measured_angles_1234[i] - 200.0) % 400.0
            if raw_leg_azimuth_1234 < 0: raw_leg_azimuth_1234 += 400.0
            raw_leg_azimuths_1234.append(raw_leg_azimuth_1234)
            current_raw_azimuth_1234 = raw_leg_azimuth_1234 
        final_calc_raw_az_to_P4_1234 = (current_raw_azimuth_1234 + measured_angles_1234[num_measured_angles_1234 - 1] - 200.0) % 400.0
        if final_calc_raw_az_to_P4_1234 < 0: final_calc_raw_az_to_P4_1234 += 400.0
        angular_misclosure_g_1234 = t_P3P4_1234 - final_calc_raw_az_to_P4_1234
        if angular_misclosure_g_1234 > 200.0: angular_misclosure_g_1234 -= 400.0
        elif angular_misclosure_g_1234 < -200.0: angular_misclosure_g_1234 += 400.0
        angular_misclosure_cc_1234 = angular_misclosure_g_1234 * 10000.0

        delta_beta_cc_list_1234 = [0] * num_measured_angles_1234
        total_correction_to_apply_cc_1234 = angular_misclosure_cc_1234 
        if num_measured_angles_1234 > 0:
                base_corr_cc = math.floor(total_correction_to_apply_cc_1234 / num_measured_angles_1234)
                delta_beta_cc_list_1234 = [base_corr_cc] * num_measured_angles_1234
                residual_cc = round(total_correction_to_apply_cc_1234 - sum(delta_beta_cc_list_1234))
                if residual_cc != 0:
                    indexed_angles = sorted([(measured_angles_1234[i], i) for i in range(num_measured_angles_1234)], reverse=True)
                    op = 1 if residual_cc > 0 else -1
                    for k in range(abs(int(residual_cc))):
                        idx_to_corr = indexed_angles[k % num_measured_angles_1234][1]
                        delta_beta_cc_list_1234[idx_to_corr] += op
        adjusted_measured_angles_1234 = [measured_angles_1234[i] + (delta_beta_cc_list_1234[i] / 10000.0) for i in range(num_measured_angles_1234)]
        adjusted_leg_azimuths_1234 = []
        current_adj_azimuth_1234 = t_P1P2_1234
        for i in range(num_traverse_legs_1234):
            adj_leg_azimuth_1234 = (current_adj_azimuth_1234 + adjusted_measured_angles_1234[i] - 200.0) % 400.0
            if adj_leg_azimuth_1234 < 0: adj_leg_azimuth_1234 += 400.0
            adjusted_leg_azimuths_1234.append(adj_leg_azimuth_1234)
            current_adj_azimuth_1234 = adj_leg_azimuth_1234
        angularly_adjusted_delta_y_list_1234 = []
        angularly_adjusted_delta_x_list_1234 = []
        for i in range(num_traverse_legs_1234):
            s_1234 = measured_distances_1234[i]
            az_rad = grads_to_radians_1234(adjusted_leg_azimuths_1234[i])
            angularly_adjusted_delta_y_list_1234.append(s_1234 * math.sin(az_rad))
            angularly_adjusted_delta_x_list_1234.append(s_1234 * math.cos(az_rad))

        y_p3_calc_after_angular_adj_1234 = y_p2_1234 + sum(angularly_adjusted_delta_y_list_1234)
        x_p3_calc_after_angular_adj_1234 = x_p2_1234 + sum(angularly_adjusted_delta_x_list_1234)
        linear_misclosure_y_m_1234 = y_p3_1234 - y_p3_calc_after_angular_adj_1234 
        linear_misclosure_x_m_1234 = x_p3_1234 - x_p3_calc_after_angular_adj_1234 
        linear_misclosure_total_m_1234 = math.sqrt(linear_misclosure_x_m_1234**2 + linear_misclosure_y_m_1234**2)
        total_traverse_dist_1234 = sum(measured_distances_1234)
        delta_delta_y_cm_list_1234 = [0] * num_traverse_legs_1234
        delta_delta_x_cm_list_1234 = [0] * num_traverse_legs_1234
        sum_s_val_for_linear_adj_1234 = sum(measured_distances_1234)
        if sum_s_val_for_linear_adj_1234 > 1e-9:
            for i in range(num_traverse_legs_1234):
                s_i = measured_distances_1234[i]
                delta_delta_y_cm_list_1234[i] = round(linear_misclosure_y_m_1234 * (s_i / sum_s_val_for_linear_adj_1234) * 100)
                delta_delta_x_cm_list_1234[i] = round(linear_misclosure_x_m_1234 * (s_i / sum_s_val_for_linear_adj_1234) * 100)
            current_sum_ddy_cm = sum(delta_delta_y_cm_list_1234)
            residual_ddy_cm = round(linear_misclosure_y_m_1234 * 100 - current_sum_ddy_cm)
            idx = 0
            while residual_ddy_cm != 0 and num_traverse_legs_1234 > 0 : 
                op = 1 if residual_ddy_cm > 0 else -1
                delta_delta_y_cm_list_1234[idx % num_traverse_legs_1234] += op
                residual_ddy_cm -= op; idx += 1
            current_sum_ddx_cm = sum(delta_delta_x_cm_list_1234)
            residual_ddx_cm = round(linear_misclosure_x_m_1234 * 100 - current_sum_ddx_cm)
            idx = 0
            while residual_ddx_cm != 0 and num_traverse_legs_1234 > 0 : 
                op = 1 if residual_ddx_cm > 0 else -1
                delta_delta_x_cm_list_1234[idx % num_traverse_legs_1234] += op
                residual_ddx_cm -= op; idx += 1

        output_rows_1234 = []
        final_adjusted_coords_for_print_1234 = [(id_p1_1234, x_p1_1234, y_p1_1234), (id_p2_1234, x_p2_1234, y_p2_1234)]
        output_rows_1234.append({"ID": id_p1_1234, "Y_m": f"{y_p1_1234:.3f}", "X_m": f"{x_p1_1234:.3f}"})
        output_rows_1234.append({"t_grad": f"{t_P1P2_1234:.4f}", "S_m": f"{s_P1P2_1234:.3f}"})
        current_final_y_1234 = y_p2_1234
        current_final_x_1234 = x_p2_1234
        s_sum_for_image_footer_1234 = s_P1P2_1234 
        sum_ang_adj_dy_for_table_1234 = 0.0
        sum_ang_adj_dx_for_table_1234 = 0.0
        sum_dbeta_cc_for_table_1234 = 0 
        sum_ddy_cm_for_table_1234 = 0
        sum_ddx_cm_for_table_1234 = 0

        for i in range(num_traverse_legs_1234):
            from_station_id_1234 = point_ids_chain_for_angles_1234[i]
            to_station_id_1234 = point_ids_chain_for_angles_1234[i+1] 
            if i == 0: 
                row1_data = {
                    "ID": from_station_id_1234, "beta": f"{measured_angles_1234[i]:.4f}",
                    "d_beta_cc": str(int(delta_beta_cc_list_1234[i])),
                    "Y_m": f"{current_final_y_1234:.3f}", "X_m": f"{current_final_x_1234:.3f}"
                }
                output_rows_1234.append(row1_data)
            sum_dbeta_cc_for_table_1234 += delta_beta_cc_list_1234[i]
            s_leg = measured_distances_1234[i]
            s_sum_for_image_footer_1234 += s_leg
            ang_adj_dy = angularly_adjusted_delta_y_list_1234[i]
            ang_adj_dx = angularly_adjusted_delta_x_list_1234[i]
            sum_ang_adj_dy_for_table_1234 += ang_adj_dy
            sum_ang_adj_dx_for_table_1234 += ang_adj_dx
            sum_ddy_cm_for_table_1234 += delta_delta_y_cm_list_1234[i]
            sum_ddx_cm_for_table_1234 += delta_delta_x_cm_list_1234[i]
            row2_data = {
                "t_grad": f"{adjusted_leg_azimuths_1234[i]:.4f}", "S_m": f"{s_leg:.3f}",
                "dY_m": f"{ang_adj_dy:.3f}", "dX_m": f"{ang_adj_dx:.3f}", 
                "d_dY_cm": str(int(delta_delta_y_cm_list_1234[i])),
                "d_dX_cm": str(int(delta_delta_x_cm_list_1234[i]))
            }
            output_rows_1234.append(row2_data)
            current_final_y_1234 += (ang_adj_dy + (delta_delta_y_cm_list_1234[i] / 100.0))
            current_final_x_1234 += (ang_adj_dx + (delta_delta_x_cm_list_1234[i] / 100.0))
            final_adjusted_coords_for_print_1234.append((to_station_id_1234, current_final_x_1234, current_final_y_1234))
            row1_data_next_st = {
                "ID": to_station_id_1234, "beta": f"{measured_angles_1234[i+1]:.4f}",
                "d_beta_cc": str(int(delta_beta_cc_list_1234[i+1])),
                "Y_m": f"{current_final_y_1234:.3f}", "X_m": f"{current_final_x_1234:.3f}"
            }
            if i+1 == num_measured_angles_1234 -1 : 
                sum_dbeta_cc_for_table_1234 += delta_beta_cc_list_1234[i+1]
            if to_station_id_1234 == id_p3_1234:
                row1_data_next_st["Y_m"] = f"{y_p3_1234:.3f}" 
                row1_data_next_st["X_m"] = f"{x_p3_1234:.3f}"
                for k_idx, (p_id_k, _, _) in enumerate(final_adjusted_coords_for_print_1234):
                    if p_id_k == id_p3_1234:
                        final_adjusted_coords_for_print_1234[k_idx] = (id_p3_1234, x_p3_1234, y_p3_1234)
                        break
            output_rows_1234.append(row1_data_next_st)
        output_rows_1234.append({"t_grad": f"{t_P3P4_1234:.4f}", "S_m": f"{s_P3P4_1234:.3f}"})
        s_sum_for_image_footer_1234 += s_P3P4_1234 
        output_rows_1234.append({"ID": id_p4_1234, "Y_m": f"{y_p4_1234:.3f}", "X_m": f"{x_p4_1234:.3f}"})

        cols_print_config_1234 = [
            ("ID", 12, '<', '<'), ("beta", 11, '>', '>'), ("d_beta_cc", 8, '>', '>'),
            ("t_grad", 11, '>', '>'), ("S_m", 11, '>', '>'), ("dY_m", 11, '>', '>'), 
            ("dX_m", 11, '>', '>'), ("d_dY_cm", 9, '>', '>'), ("d_dX_cm", 9, '>', '>'),
            ("Y_m", 13, '>', '>'), ("X_m", 13, '>', '>')
        ]
        w = {c[0]: c[1] for c in cols_print_config_1234} 
        
        num_cols = len(cols_print_config_1234)
        separator_width = 3 
        calculated_total_width = sum(c[1] for c in cols_print_config_1234) + (num_cols - 1) * separator_width
        
        table_title_1234 = "Link Traverse Computation Table"
        print("\n" + "-" * calculated_total_width) 
        print(f"{table_title_1234:^{calculated_total_width}}") 
        
        header_str_1234 = (
            f"{'Station ID':^{w['ID']}} | "
            f"{'β(grad)':^{w['beta']}} | "
            f"{'δβ(cc)':^{w['d_beta_cc']}} | "
            f"{'t(grad)':^{w['t_grad']}} | "
            f"{'S(m)':^{w['S_m']}} | "
            f"{'ΔY(m)':^{w['dY_m']}} | "
            f"{'ΔX(m)':^{w['dX_m']}} | "
            f"{'δΔY(cm)':^{w['d_dY_cm']}} | "
            f"{'δΔX(cm)':^{w['d_dX_cm']}} | "
            f"{'Y(m)':^{w['Y_m']}} | "
            f"{'X(m)':^{w['X_m']}}"
        )

        print("-" * calculated_total_width)
        print(header_str_1234)
        print("-" * calculated_total_width)

        for row_data_1234 in output_rows_1234:
            parts_1234 = []
            for key_cfg, width_cfg, data_align_cfg, _ in cols_print_config_1234:
                value_str_1234 = str(row_data_1234.get(key_cfg, "")) 
                parts_1234.append(f"{value_str_1234:{data_align_cfg}{width_cfg}}")
            print(" | ".join(parts_1234))

        print("-" * calculated_total_width) 
        
        sum_values_map_1234 = {
            "ID": "", "beta": f"Σ={sum_beta_1234:.4f}", 
            "d_beta_cc": f"Σ={int(sum_dbeta_cc_for_table_1234)}",
            "S_m": f"Σ={s_sum_for_image_footer_1234:.3f}", 
            "dY_m": f"Σ={sum_ang_adj_dy_for_table_1234:.3f}", 
            "dX_m": f"Σ={sum_ang_adj_dx_for_table_1234:.3f}", 
            "d_dY_cm": f"Σ={int(sum_ddy_cm_for_table_1234)}", 
            "d_dX_cm": f"Σ={int(sum_ddx_cm_for_table_1234)}",
            "Y_m": f"ΔY={y_p3_1234 - y_p2_1234:.3f}", 
            "X_m": f"ΔX={x_p3_1234 - x_p2_1234:.3f}"  
        }
        
        footer_parts_1234 = []
        for key_cfg, width_cfg, _, sum_align_cfg in cols_print_config_1234:
            value_str_1234 = str(sum_values_map_1234.get(key_cfg, ""))
            footer_parts_1234.append(f"{value_str_1234:{sum_align_cfg}{width_cfg}}")
        print(" | ".join(footer_parts_1234))
        print("-" * calculated_total_width)

    finally:
        if output_file_path_1234 and file_out_1234 is not None: 
            sys.stdout = original_stdout_1234 
            file_out_1234.close()
            original_stdout_1234.write(f"\nOutput also saved to: {output_file_path_1234}\n") 
        elif 'file_out_1234' in locals() and file_out_1234 is not None: 
             sys.stdout = original_stdout_1234
             file_out_1234.close()


        print(t_P1P2_1234)
        print(s_P1P2_1234)
        print(t_P3P4_1234) 
        print(s_P3P4_1234)

    print(sum_beta_1234)
    print(measured_angles_1234)
    print(delta_beta_cc_list_1234)
    print(num_measured_angles_1234)
    print(num_unknown_points_1234)
    print(len(point_ids_for_legs_1234))
    print(point_ids_for_legs_1234)
    print(point_ids_chain_for_angles_1234)

if __name__ == "__main__":
    main_1234()
