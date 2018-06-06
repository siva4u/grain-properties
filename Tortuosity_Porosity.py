# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 10:41:14 2018

@author: Siva
"""
from numpy import arange, float16
import time, os, itertools, shutil, cv2, gc
import multiprocessing
import matplotlib.image as mpimg
import numpy as np
from skimage.external.tifffile import imsave as skimtiffsave
from subprocess import Popen
#------------
def func_input_section():
    # Input section
    # Directory containing positions files
    Directory = 'C:\\tortuosity'
    # File name characteristics
    left_strip = "final_"
    right_strip = ".posi"
    # Flags for different options
    flags = { 'plot_3D_figure' :  False ,\
              'calc_porosity_flag' :  False, \
              'calc_tortuosity_flag' : False , \
              'clean_up_required_flag' : True \
             }
    # End of input section
    return Directory, left_strip, right_strip, flags

def func_3D_fig(data_file_name, x_coord, y_coord, z_coord, particle_diameter):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection =  '3d')
    ax.scatter(x_coord, y_coord, z_coord, s= particle_diameter)
    plt.savefig('3D visual of ' + data_file_name + 'particle_positions.jpg', dpi= 300)
    plt.close()
    del fig, Axes3D, ax
#----------
#--- func_read_file section
# Does 2 things
# 1) get tortuosity box extremes
# 2) filter out particles not within the tortuosity box
def func_read_file(position_file_name, three_D_figure_needed):
    # --- Input Section
    # Enter number of cross sectional slices needed
    number_of_slices = 256
    # Enter width of tortuosity box
    tortuosity_box_width = 0.0
    #--- End of input section
    #-------------------
    #--- Read input box dimension from file
    with open(position_file_name) as lines:
        box_dims = np.loadtxt(list(itertools.islice(lines,5,8)))
    #--- Calculater center of simulation box
    x_bottom, x_top = box_dims[0,0], box_dims[0,1]
    y_bottom, y_top = box_dims[1,0], box_dims[1,1]
    z_bottom, z_top = box_dims[2,0], box_dims[2,1]
    
    x_center = x_bottom + (x_top - x_bottom) * 0.5
    y_center = y_bottom + (y_top - y_bottom) * 0.5
    z_center = z_bottom + (z_top - z_bottom) * 0.5
    #----------------------------------------------
    #--- Getting particle data
    data = np.genfromtxt(position_file_name, unpack = True, skip_header = 9).astype(float16)
    #---
    print("# of particles before filtering: ", np.shape(data))
    #--- Ignore columns of data not required for tortuosity estimation
    data = data[1:5]
    #--- may need to convert to microns
    particle_diameter = data[0]
    #--- Largest radius of particle
    largest_rad = 0.5 *  np.max(particle_diameter)
    #---
    #--- Getting half width of tortuosity box
    if tortuosity_box_width <= 0.0:
        x_thick = x_top - x_bottom
        y_thick = y_top - y_bottom
        z_thick = z_top - z_bottom
    
    tortuosity_box_width = min(x_thick, y_thick , z_thick) - np.max(particle_diameter)
    
    box_half_width = 0.5 * tortuosity_box_width
    #--- Make cubic tortuosity box
    box_x_bottom, box_x_top = x_center -  box_half_width, x_center + box_half_width
    box_y_bottom, box_y_top = y_center -  box_half_width, y_center + box_half_width
    box_z_bottom, box_z_top = z_center -  box_half_width, z_center + box_half_width
    
    tortuosity_box_dimensions = np.array([box_x_bottom, box_x_top, \
                                          box_y_bottom, box_y_top, \
                                          box_z_bottom, box_z_top, \
                                          ])
    
    #--- Make an array of Z-coordinates at which slices are made
    z_plane_slice_coord = box_z_bottom + (box_z_top - box_z_bottom) / \
                          number_of_slices * arange(number_of_slices)
    #--- Plot 3D figure
    if three_D_figure_needed:
        x_coord, y_coord, z_coord = data[1], data[2], data[3]
        func_3D_fig(position_file_name, x_coord, y_coord, z_coord, particle_diameter)
    # End of plotting 3D figure
    #-----
    #--- Calculating position of farthesy from edge plane away from box_center
    extreme_particle = box_half_width + largest_rad
    x_bottom_plane   = x_center - extreme_particle
    x_top_plane      = x_center + extreme_particle
    y_bottom_plane   = y_center - extreme_particle
    y_top_plane      = y_center + extreme_particle
    z_bottom_plane   = z_center - extreme_particle
    z_top_plane      = z_center + extreme_particle
    #--- Removing particle falling oufside tortuosity box
    data = data[:, np.logical_and(data[1] > x_bottom_plane, data[1] < x_top_plane)]
    data = data[:, np.logical_and(data[2] > y_bottom_plane, data[2] < y_top_plane)]
    data = data[:, np.logical_and(data[3] > z_bottom_plane, data[3] < z_top_plane)]
    #---
    print("# of particles outside tortuosity box : ", np.shape)
    return data, tortuosity_box_dimensions, z_plane_slice_coord

# Multiprocessing function
# Function for slicing different Z planes
def func_z_plane_slice(q, data, slice_counter, box_dimensions):
    from matplotlib.pyplot  import figure, gca, savefig
    from matplotlib.patches import Circle
    figure(slice_counter)
    ax1 = gca(aspect = 'equal', \
              frameon = False,  \
              xlim = [box_dimensions[0], box_dimensions[1]], \
              ylim = [box_dimensions[2], box_dimensions[3]]  \
              )
    
    ax1.xaxis._visible = False
    ax1.yaxis._visible = False
    
    for particle_projected_radius, x_coord, y_coord in zip(  \
                          data[0], data[1], data[2]) :
        ax1.add_artist(Circle((x_coord, y_coord), particle_projected_radius, \
                              color = 'white', fill = True))
    
    savefig(str(slice_counter).zfill(3) +'.jpg', dpi = 120, bbox_inches = 'tight',\
               axes = 'off', pad_inches=0.0, facecolor="black")
    q.put(('S'))
    q.close()

#---  End of multiprocessing

def func_trim_n_tiff(posi_name):
    tiff_file_dir= ".\\CS_tiff_images"
    if not os.path.exists(tiff_file_dir):
        os.mkdir(tiff_file_dir)

    jpg_file_list = list(file_name for file_name in os.listdir(os.getcwd())\
                         if file_name.endswith(".jpg"))

    for jpg_file_name in jpg_file_list:
        img = cv2.imread( jpg_file_name, 0)
        cropped_img = img[3:447, 3:447]
        resized_img = cv2.resize(cropped_img  ,(256,256))
        tiff_file_name = tiff_file_dir + "\\" +jpg_file_name.rstrip(".jpg") + ".tiff"
        skimtiffsave(tiff_file_name, resized_img, compress = 0)

    print("\n Trimmed & Tiffed the posi data", posi_name, time.ctime(), "\n")
# Porosity estimation
def func_dark_pixel_count(image_file, white_intensity):
    img = mpimg.imread(image_file)
    total_pixel_count = img.size
    dark_pixel_count = img[(img < white_intensity)].size
    del img
    return dark_pixel_count / total_pixel_count
def func_calc_porosity(posi_name):
    white_intensity = 3
    porosity        = 0.0
    file_counter    = 0
    tiff_files_in_dir = [item for item in os.listdir(".\\CS_tiff_images") \
                         if item.endswith(('.tiff'))]
    out_file = open('.\\' + posi_name + '.porosity', 'w')
    for tiff_file in tiff_files_in_dir:
        
        figure_porosity = func_dark_pixel_count(".CS_tiff_images\\"+ tiff_file,\
                                                white_intensity)
        out_file.write(str(tiff_file) + '    ' + str(figure_porosity) + '\n')
        porosity += figure_porosity
        file_counter+=1
    out_file.write('\n Mean Porosity Z- direction = ', porosity / file_counter)
    out_file.close()
# Rename results
def func_rename_results(posi_name):
    Results_dir = '.\\Results_'+ posi_name
    if not os.path.exists(Results_dir):
        os.mkdir(Results_dir)
    
    try:
        os.renames("Graph.pdf", str(Results_dir +"\\Tortuosity_"+posi_name+ ".pdf"))
        os.renames("Rwalk_x.csv", str(Results_dir +"\\Tortuosity_steps_"+posi_name+".csv"))
        os.renames("raw_data_on_gradients_of_msd.tif", \
                   str(Results_dir + "\\gradients_msd_"+ posi_name + ".tiff"))
        os.renames("C_Dellipsoid_x_Graph2D.pdf", \
                   str(Results_dir + "C_Dellipsoid_" + posi_name + "_Graph2D.pdf"))
        os.renames("C_Dellipsoid_x_log.txt", \
                   str(Results_dir + "\\C_Dellipsoid_"+ posi_name + "_log.txt"))
        os.renames(".\\graph_3D_s\\index.html",\
                   str( Results_dir + "3D_unwrapd_dir_diff_" + posi_name + ".html"))
        os.renames(".\\graph_3D_e\\index.html",\
                   str( Results_dir + "3D_wrapd_dir_diff_" + posi_name + ".html"))
        os.renames(".\\graph_3D_m\\index.html",\
                   str( Results_dir + "3D_map_walkers_" + posi_name + ".html"))
    except:
        print("Error *** while renaming results files", posi_name)
        pass
#--- Run tortuosity estimator program
def func_calc_tortuosity(posi_name):
    print("Tortuosity estimation started at time : ", time.ctime())
    if not os.path.exists(".\\perlRscript.bat"):
        perlRscript = open(".\\perlRscript.bat", "w")
        commands = [ "copy ..\\source_files\\*.*  .\\"  , \
                 "mkdir .\labeled",  \
                 "perl ClabelM6.pl",   \
                 "perl PreRwalkM6.pl", \
                 "perl C_Rwalk_x.pl",  \
                 "Rscript .\C_Dellipsoid_x.R"  \
                ]
        perlRscript.writelines('\n'.join(commands))
        perlRscript.close()
      
    try:
        p = Popen("perlRscript.bat")
        stdout, stderr = p.communicate()
        func_rename_results(posi_name)
    except:
        try:
            print(stderr)
        except:
            pass
        raise SystemExit("Error : Running Perl R scripts", posi_name)
    print("tortuosity estimation ended at time : ", time.ctime())
    return None
# Begin clean up section
def func_clean_up(posi_name, posi_file_name):
    files_to_remove = [item for item in os.listdir(os.getcwd()) if item.endswith   \
    ((".jpg", ".exe", ".pl", ".txt", ".csv", ".bin", ".R", ".bat", ".dll"))]
    
    for item in files_to_remove:
        os.remove(item)
    
    shutil.make_archive(posi_name + '-CS_tiff_Images', 'zip', "CS_tiff_images")
    
    shutil.copy("..\\" + posi_file_name, ".\\")
    
    try:
        shutil.rmtree(".\\CS_tiff_images", ignore_errors = True, onerror = None)
        shutil.rmtree(".\\labeled", ignore_errors = True, onerror = None)
        print("Clean-up Completed @ :", time.ctime())
    except:
        print("Error : Cleaning up", posi_name)
        pass
# End of Cleanup Section

def func_make_cross_sectional_images( left_strip, right_strip, flags):
    posi_files_list = [file_name for file_name in os.listdir(os.getcwd())\
                       if file_name.endswith(right_strip)   ]
    if not posi_files_list:
        raise SystemExit("No posi files")
    else:
        print("number of posi files = ", len(posi_files_list))
    processed_posi_file_count = 0
    for posi_file_name in posi_files_list:
        print("Reading file", posi_file_name, "@ time : ", time.ctime())
        data, box_dimensions, z_plane_slice_coord = \
                         func_read_file(posi_file_name, flags['3D_fig_needed'])
        particle_diameter = data[0]
        largest_particle_rad = np.max(particle_diameter) * 0.5
        z_coord = data[3]
        # Get tortuosity box dimensions
        x_bottom, x_top = box_dimensions[0], box_dimensions[1]
        y_bottom, y_top = box_dimensions[2], box_dimensions[3]
        # Get position of particle at extremities
        x_bottom_most_particle = x_bottom - largest_particle_rad
        x_top_most_particle    = x_top + largest_particle_rad
        y_bottom_most_particle = y_bottom - largest_particle_rad
        y_top_most_particle    = y_top + largest_particle_rad
        
        particle_radius_squared = particle_diameter * particle_diameter * 0.25
        slice_counter = arange(len(z_plane_slice_coord))
        print(" multiprocessing started @", time.ctime())
        
        posi_name = (posi_file_name.lstrip(left_strip)).rstrip(right_strip)
        path = ".\\" + posi_name
        
        if not os.path.exists(path): os.mkdir(path)
        
        os.chdir(path)
        q = multiprocessing.Queue()
        # --- Create worker processes empty list
        workers = None * len(z_plane_slice_coord)

        for p in range(len(z_plane_slice_coord)):
            part_proj_depth = z_plane_slice_coord[p] - z_coord
            # Eliminate particles with 'NaN' as projected depth.
            # This filters out particles not in the plane
            data_slice = np.copy(data[0:3]).astype(float16)
            data_slice[0] = np.sqrt( particle_radius_squared -   \
                           (part_proj_depth * part_proj_depth)).astype(float16)

            data_slice = data_slice[:, np.isfinite(data_slice[0])]

            # Eliminate particles outside of x_bottom and x_top
            data_slice = data_slice[:, np.logical_and(\
                         data_slice[1] >= x_bottom_most_particle, \
                         data_slice[1] <= x_top_most_particle)]

            # Eliminate particles outside of y_bottom and y_top
            data_slice = data_slice[:, np.logical_and(\
                         data_slice[2] >= y_bottom_most_particle, \
                         data_slice[2] <= y_top_most_particle)]

            # Create worker processes
            workers[p] = multiprocessing.Process(target = func_z_plane_slice, \
                   args = (q, data_slice, slice_counter[p], box_dimensions))

        #---End of workers creation
        Num_of_CPU_cores = multiprocessing.cpu_count()
        scaling_factor = 2.
        Num_of_jobs_per_batch = int(Num_of_CPU_cores * scaling_factor)
        Num_of_batches = len(workers) // Num_of_jobs_per_batch
        ceil_job = 0

        for batch_count in range(1, (Num_of_batches + 1)):
            floor_job = (batch_count - 1) * Num_of_jobs_per_batch
            ceil_job = floor_job + Num_of_jobs_per_batch
            for p in workers[floor_job : ceil_job]:
                p.start()
            for p in workers[floor_job : ceil_job]:
                p.join()
                while not q.empty(): q.get()
        for p in workers[ceil_job : ]:
            p.start()
        for p in workers[ceil_job : ]:
            p.join()
            while not q.empty(): q.get()
        # End of job submission
        for p in multiprocessing.active_children():
            p.terminate()
        print("multiprocessing ended at ", time.ctime())

        # Trim and Tiff
        func_trim_n_tiff(posi_name)
        if flags['calc_porosity_flag']:
            func_calc_porosity(posi_name)

        # Start of tortuosity estimation
        if flags['calc_tortuosity_flag']:
            func_calc_tortuosity(posi_name)
            if flags['clean_up_required_flag']:
                func_clean_up(posi_name, posi_file_name)

        # End of tortuosity estimation
        os.chdir('..\\')

        processed_posi_file_count += 1
        print("Tortuosity estimated file count = ", processed_posi_file_count)
        print("Ended at ", time.ctime())
        gc.collect()
    print("\n Total number of posi files processed", processed_posi_file_count)
#--g
def main():
    Dir_containing_posi_files, left_strip, right_strip, flags = func_input_section()
    os.chdir("Dir_containing_posi_files")
    func_make_cross_sectional_images(left_strip, right_strip, flags)

if __name__ == "__main__":
    start_time = time.time()
    multiprocessing.freeze_support()
    main()
    end_time = time.time()
    print(" Total execution time = ", end_time - start_time)