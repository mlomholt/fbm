When running the fbm_skel_*.m files first run

addpath('nestedsampling','fbm_util')

Files:

fbm_skel_artificial.m
Runs nested sampling for a single artificial trajectory. Un-remark the relevant trajectory data file in the beginning of the program.

fbm_skel_CHO.m
Runs nested sampling for several trajectories. Un-remark the relevant trajectory data file in the beginning of the program.

generate_subdiffusive_track.m
Generates a single subdiffusive track

generate_superdiffusive_track.m
Generates a single superdiffusive track

generate_random_tracks.m
Generates several random tracks

fbm_load_CHO_data.m
Loads a bunch of data into memory to use for fbm_plot_CHO_figs.m. Note: first run the fbm_skel_CHO.m file to generate the data.

fbm_plot_CHO_figs.m
Plots a number of figures.

fbm_print_table.m
Prints a table for use in an article written in Latex


