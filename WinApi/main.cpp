#include <iostream>
#include <regex>
#include <windows.h>
#include <tchar.h>
#include <gtk/gtk.h>
#include <gtkgraph.h>
#include <sstream>
#include <string>
#include <cmath>
#include <vector> 
// Application Includes
#include <buttons.h>



//int main( int argc, char **argv );
//int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine,int nCmdShow){   
//    return main (__argc, __argv);     
//}

int main( int argc, char **argv ) {

    // create the console
  //  if(AllocConsole()) {
  //      freopen("CONOUT$", "w", stdout);
   
   SetConsoleTitle("Three-phase 1D Well Simulator Console");
  //      SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_RED);  
 //   }
    gtk_init( &argc, &argv );
    // set std::cout to use my custom streambuf
  //  outbuf ob;
  //  std::streambuf *sb = std::cout.rdbuf(&ob);       
   
    GtkWidget *window;
    GtkWidget *vbox;
    GtkWidget *menubar;
    GtkWidget *filemenu;
    GtkWidget *file;
    GtkWidget *quit;
    GtkWidget *helpmenu;
    GtkWidget *help;
    GtkWidget *about;


    // Setting Main Window
    
    window = gtk_window_new( GTK_WINDOW_TOPLEVEL );     
    gtk_window_set_title(GTK_WINDOW(window), "WellDrift");
    gtk_window_set_default_size(GTK_WINDOW(window), 700, 400);
    gtk_window_set_position(GTK_WINDOW(window), GTK_WIN_POS_CENTER);
    gtk_window_set_icon(GTK_WINDOW(window), create_pixbuf(".\\images\\icon.png"));
    gtk_container_set_border_width(GTK_CONTAINER(window), 5);

    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(window), vbox);

    menubar = gtk_menu_bar_new();
    filemenu = gtk_menu_new();
    helpmenu = gtk_menu_new();

    file = gtk_menu_item_new_with_label("File");
    quit = gtk_menu_item_new_with_label("Quit");
    help = gtk_menu_item_new_with_label("Help");
    about = gtk_menu_item_new_with_label("About");

    gtk_menu_item_set_submenu   (GTK_MENU_ITEM(file), filemenu);
    gtk_menu_item_set_submenu   (GTK_MENU_ITEM(help), helpmenu);
    gtk_menu_shell_append       (GTK_MENU_SHELL(filemenu), quit);
    gtk_menu_shell_append       (GTK_MENU_SHELL(helpmenu), about);
    gtk_menu_shell_append       (GTK_MENU_SHELL(menubar), file);
    gtk_menu_shell_append       (GTK_MENU_SHELL(menubar), help);
    gtk_box_pack_start          (GTK_BOX(vbox), menubar, FALSE, FALSE, 3); 

    // Main Table
    GtkWidget *main_table;
    main_table = gtk_table_new(2, 3, FALSE); 
    gtk_box_pack_start(GTK_BOX(vbox), main_table, TRUE, TRUE, 0);

    // Reading Initial Data from File
    SharedPointer<WellInitialData> well_initial_data(new WellInitialData);
    std::string well_name_path("..\\WellData\\Setup.txt");        
    std::ifstream InFile( well_name_path.c_str() );


    InFile.ignore(50,':');	InFile >> well_initial_data->m_number_of_nodes;  
    InFile.ignore(50,':');	InFile >> well_initial_data->m_well_inclination;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_well_length;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_delta_t;		
    InFile.ignore(50,':');	InFile >> well_initial_data->m_initial_pressure;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_initial_mixture_velocity;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_initial_gas_vol_frac;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_initial_oil_vol_frac;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_well_diameter;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_toe_mixture_velocity;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_tolerance;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_max_delta_t;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_final_time;          
    // InFile.ignore(50,':');	InFile >> p_initial_data.PROFILE_PARAM_C_0;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_heel_pressure;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_ref_temperature;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_oil_density;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_oil_viscosity;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_water_density;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_water_viscosity;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_gas_ref_density;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_gas_ref_pressure;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_gas_sound_speed;
    InFile.ignore(50,':');	InFile >> well_initial_data->m_gas_viscosity;
    
    
    // ------------------------------------------------------------------------------
    // Well Dimension Layout
    GtkWidget *frame_well_dimensions;
    frame_well_dimensions = gtk_frame_new("Well Dimensions");
    gtk_table_attach(GTK_TABLE(main_table), frame_well_dimensions, 0, 1, 0, 1, GTK_SHRINK, GTK_SHRINK, 5, 5);
   // gtk_box_pack_start(GTK_BOX(vbox), frame_well_dimensions, TRUE, TRUE, 0);
    GtkWidget *table_well_dimensions;
    table_well_dimensions = gtk_table_new(4, 2, FALSE);  
    gtk_container_add(GTK_CONTAINER(frame_well_dimensions), table_well_dimensions);
    GtkWidget *label_number_of_nodes; label_number_of_nodes     = gtk_label_new("Number of points"); 
    GtkWidget *label_well_length  ; label_well_length           = gtk_label_new("Well Length [m]");
    GtkWidget *label_well_diameter; label_well_diameter         = gtk_label_new("Well Diameter [m]");
    GtkWidget *label_well_inclination; label_well_inclination   = gtk_label_new("Well Inclination [deg]");  
    GtkObject *adjustment_number_of_nodes;   
    GtkObject *adjustment_well_length;
    GtkObject *adjustment_well_diameter;
    GtkObject *adjustment_well_inclination;
    GtkWidget *spin_button_number_of_nodes;
    GtkWidget *spin_button_well_length;
    GtkWidget *spin_button_well_diameter;
    GtkWidget *spin_button_well_inclination; 
    adjustment_number_of_nodes      = gtk_adjustment_new(well_initial_data->m_number_of_nodes, 0, 1.0e10, 1, 1, 0);
    adjustment_well_length          = gtk_adjustment_new(well_initial_data->m_well_length, -1.0e50, 1.0e50, 1, 1, 0);
    adjustment_well_diameter        = gtk_adjustment_new(well_initial_data->m_well_diameter, -1.0e50, 1.0e50, 1, 1, 0);
    adjustment_well_inclination     = gtk_adjustment_new(well_initial_data->m_well_inclination, -1.0e50, 1.0e50, 1, 1, 0);

    spin_button_number_of_nodes     = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_number_of_nodes), 1,0);
    spin_button_well_length         = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_well_length), 0.1,1);
    spin_button_well_diameter       = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_well_diameter), 0.1,2);    
    spin_button_well_inclination    = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_well_inclination), 0.1,1);
    gtk_table_attach(GTK_TABLE(table_well_dimensions), label_number_of_nodes, 0, 1, 0, 1, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_well_dimensions), label_well_length, 0, 1, 1, 2, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_well_dimensions), label_well_diameter, 0, 1, 2, 3, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_well_dimensions), label_well_inclination, 0, 1, 3, 4, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_well_dimensions), spin_button_number_of_nodes, 1, 2, 0, 1, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_well_dimensions), spin_button_well_length, 1, 2, 1, 2, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_well_dimensions), spin_button_well_diameter, 1, 2, 2, 3, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_well_dimensions), spin_button_well_inclination, 1, 2, 3, 4, GTK_SHRINK, GTK_SHRINK, 5, 5);
    // ------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------
    // Boundary Conditions Layout
    GtkWidget *frame_boundary_conditions;
    frame_boundary_conditions = gtk_frame_new("Boundary Conditions");
    //gtk_box_pack_start(GTK_BOX(vbox), frame_boundary_conditions, TRUE, TRUE, 0);
    gtk_table_attach(GTK_TABLE(main_table), frame_boundary_conditions, 0, 1, 1, 2, GTK_SHRINK, GTK_SHRINK, 5, 5); 
    GtkWidget *table_boundary_conditions;
    table_boundary_conditions = gtk_table_new(4, 2, FALSE);  
    gtk_container_add(GTK_CONTAINER(frame_boundary_conditions), table_boundary_conditions);
    GtkWidget *label_toe_velocity  ; label_toe_velocity    = gtk_label_new("Toe Velocity [m/s]");
    GtkWidget *label_heel_pressure ; label_heel_pressure   = gtk_label_new("Heel Pressure [Pa]");
    GtkWidget *label_file_chooser_button_inflow ; label_file_chooser_button_inflow   = gtk_label_new("Inflow File [m3/s]");
    GtkObject *adjustment_toe_velocity;
    GtkObject *adjustment_heel_pressure;
    GtkWidget *spin_button_toe_velocity;
    GtkWidget *spin_button_heel_pressure;
    adjustment_toe_velocity   = gtk_adjustment_new(well_initial_data->m_toe_mixture_velocity, -1.0e50, 1.0e50, 1, 1, 0);
    adjustment_heel_pressure = gtk_adjustment_new (well_initial_data->m_heel_pressure, -1.0e50, 1.0e50, 1, 1, 0);
    spin_button_toe_velocity = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_toe_velocity), 0.1,1);
    spin_button_heel_pressure = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_heel_pressure), 0.1,1);
    gtk_table_attach(GTK_TABLE(table_boundary_conditions), label_toe_velocity, 0, 1, 0, 1, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_boundary_conditions), label_heel_pressure, 0, 1, 1, 2, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_boundary_conditions), label_file_chooser_button_inflow, 0, 2, 2, 3, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_boundary_conditions), spin_button_toe_velocity, 1, 2, 0, 1, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_boundary_conditions), spin_button_heel_pressure, 1, 2, 1, 2, GTK_SHRINK, GTK_SHRINK, 5, 5);
    GtkWidget *file_chooser_button_inflow;  file_chooser_button_inflow = gtk_file_chooser_button_new("Select an Inflow File", GTK_FILE_CHOOSER_ACTION_OPEN);
    gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER (file_chooser_button_inflow), "..\\WellData\\");
    gtk_file_chooser_button_set_width_chars(GTK_FILE_CHOOSER_BUTTON (file_chooser_button_inflow), 15);
    gtk_table_attach(GTK_TABLE(table_boundary_conditions), file_chooser_button_inflow, 0, 2, 3, 4, GTK_SHRINK, GTK_SHRINK, 5, 5);
    // ------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------
    // Initial Conditions Layout
    GtkWidget *frame_initial_conditions;
    frame_initial_conditions = gtk_frame_new("Initial Conditions");
    //gtk_box_pack_start(GTK_BOX(vbox), frame_initial_conditions, TRUE, TRUE, 0);
    gtk_table_attach(GTK_TABLE(main_table), frame_initial_conditions, 1, 2, 0, 1, GTK_SHRINK, GTK_SHRINK, 5, 5); 
    GtkWidget *table_initial_conditions;
    table_initial_conditions = gtk_table_new(4, 2, FALSE);  
    gtk_container_add(GTK_CONTAINER(frame_initial_conditions), table_initial_conditions);
    GtkWidget *label_initial_pressure ; label_initial_pressure   = gtk_label_new("Initial Pressure [Pa]");
    GtkWidget *label_initial_velocity ; label_initial_velocity   = gtk_label_new("Initial Velocity [m/s]");
    GtkWidget *label_initial_gas_vol_frac ; label_initial_gas_vol_frac   = gtk_label_new("Initial Gas Volume Fraction [-]");
    GtkWidget *label_initial_oil_vol_frac ; label_initial_oil_vol_frac   = gtk_label_new("Initial Oil Volume Fraction [-]");

    GtkObject *adjustment_initial_pressure; adjustment_initial_pressure         = gtk_adjustment_new(well_initial_data->m_initial_pressure, -1.0e50, 1.0e50, 1, 1, 0);
    GtkObject *adjustment_initial_velocity; adjustment_initial_velocity         = gtk_adjustment_new(well_initial_data->m_initial_mixture_velocity, -1.0e50, 1.0e50, 1, 1, 0);
    GtkObject *adjustment_initial_gas_vol_frac; adjustment_initial_gas_vol_frac = gtk_adjustment_new(well_initial_data->m_initial_gas_vol_frac, -1.0e50, 1.0e50, 1, 1, 0);
    GtkObject *adjustment_initial_oil_vol_frac; adjustment_initial_oil_vol_frac = gtk_adjustment_new(well_initial_data->m_initial_oil_vol_frac, -1.0e50, 1.0e50, 1, 1, 0);

    GtkWidget *spin_button_initial_pressure;    spin_button_initial_pressure        = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_initial_pressure), 0.1,1);
    GtkWidget *spin_button_initial_velocity;    spin_button_initial_velocity        = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_initial_velocity), 0.1,1);
    GtkWidget *spin_button_initial_gas_vol_frac;spin_button_initial_gas_vol_frac    = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_initial_gas_vol_frac), 0.1,3);
    GtkWidget *spin_button_initial_oil_vol_frac;spin_button_initial_oil_vol_frac    = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_initial_oil_vol_frac), 0.1,3);         
    
    gtk_table_attach(GTK_TABLE(table_initial_conditions), label_initial_pressure, 0, 1, 0, 1, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_initial_conditions), label_initial_velocity, 0, 1, 1, 2, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_initial_conditions), label_initial_gas_vol_frac, 0, 1, 2, 3, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_initial_conditions), label_initial_oil_vol_frac, 0, 1, 3, 4, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_initial_conditions), spin_button_initial_pressure, 1, 2, 0, 1, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_initial_conditions), spin_button_initial_velocity, 1, 2, 1, 2, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_initial_conditions), spin_button_initial_gas_vol_frac, 1, 2, 2, 3, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_initial_conditions), spin_button_initial_oil_vol_frac, 1, 2, 3, 4, GTK_SHRINK, GTK_SHRINK, 5, 5);
    // ------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------
    // Simulation Settings Layout
    GtkWidget *frame_simulation_settings;
    frame_simulation_settings = gtk_frame_new("Simulation Settings");
    //gtk_box_pack_start(GTK_BOX(vbox), frame_simulation_settings, TRUE, TRUE, 0);
    gtk_table_attach(GTK_TABLE(main_table), frame_simulation_settings, 1, 2, 1, 2, GTK_SHRINK, GTK_SHRINK, 5, 5); 
    GtkWidget *table_simulation_settings;
    table_simulation_settings = gtk_table_new(4, 2, FALSE);  
    gtk_container_add(GTK_CONTAINER(frame_simulation_settings), table_simulation_settings);
    GtkWidget *label_timestep ; label_timestep      = gtk_label_new("Initial Timestep [s]");
    GtkWidget *label_max_delta_t ; label_max_delta_t  = gtk_label_new("Maximum Timestep [s]");
    GtkWidget *label_final_time ; label_final_time  = gtk_label_new("Final Time [s]");
    GtkWidget *label_tolerance ; label_tolerance    = gtk_label_new("Tolerance [-]");

    GtkObject *adjustment_timestep; adjustment_timestep     = gtk_adjustment_new(well_initial_data->m_delta_t, 0, 1.0e50, 1, 1, 0);
    GtkObject *adjustment_max_delta_t; adjustment_max_delta_t = gtk_adjustment_new(well_initial_data->m_max_delta_t, 0, 1.0e50, 1, 1, 0);
    GtkObject *adjustment_final_time; adjustment_final_time = gtk_adjustment_new(well_initial_data->m_final_time, -1.0e50, 1.0e50, 1, 1, 0);
    GtkObject *adjustment_tolerance; adjustment_tolerance   = gtk_adjustment_new(well_initial_data->m_tolerance, -1.0e50, 1.0e50, 1, 1, 0);

    GtkWidget *spin_button_timestep;    spin_button_timestep   = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_timestep), 0.1,3);
    GtkWidget *spin_button_max_delta_t;  spin_button_max_delta_t = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_max_delta_t), 0.1,1);
    GtkWidget *spin_button_final_time;  spin_button_final_time = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_final_time), 0.1,1);
    GtkWidget *spin_button_tolerance;   spin_button_tolerance  = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_tolerance), 0.1,10);  

    gtk_table_attach(GTK_TABLE(table_simulation_settings), label_timestep, 0, 1, 0, 1, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_simulation_settings), label_max_delta_t, 0, 1, 1, 2, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_simulation_settings), label_final_time, 0, 1, 2, 3, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_simulation_settings), label_tolerance, 0, 1, 3, 4, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_simulation_settings), spin_button_timestep, 1, 2, 0, 1, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_simulation_settings), spin_button_max_delta_t, 1, 2, 1, 2, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_simulation_settings), spin_button_final_time, 1, 2, 2, 3, GTK_SHRINK, GTK_SHRINK, 5, 5);
    gtk_table_attach(GTK_TABLE(table_simulation_settings), spin_button_tolerance, 1, 2, 3, 4, GTK_SHRINK, GTK_SHRINK, 5, 5);
    // ------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------
    // Drift Flux Parameters Layout
    //GtkWidget *frame_drift_flux_parameters;
    //frame_drift_flux_parameters = gtk_frame_new("Drift Flux Parameters");
    ////gtk_box_pack_start(GTK_BOX(vbox), frame_simulation_settings, TRUE, TRUE, 0);
    //gtk_table_attach(GTK_TABLE(main_table), frame_drift_flux_parameters, 2, 3, 0, 1, GTK_SHRINK, GTK_SHRINK, 5, 5); 
    //GtkWidget *table_drift_flux_parameters;
    //table_drift_flux_parameters = gtk_table_new(7, 4, FALSE);  
    //gtk_container_add(GTK_CONTAINER(frame_drift_flux_parameters), table_drift_flux_parameters);
    //GtkWidget *label_gas_liquid_profile_parameter ; label_gas_liquid_profile_parameter      = gtk_label_new("Gas-Liquid Profile Parameter");
    //GtkWidget *label_A ; label_A  = gtk_label_new("A");
    //GtkWidget *label_B ; label_B  = gtk_label_new("B");
    //GtkWidget *label_Fv ; label_Fv  = gtk_label_new("F_v");
    //GtkWidget *label_gas_liquid_drift_velocity ; label_gas_liquid_drift_velocity  = gtk_label_new("Gas-Liquid Drift Velocity");
    //GtkWidget *label_a1 ; label_a1  = gtk_label_new("a1");
    //GtkWidget *label_a2 ; label_a2  = gtk_label_new("a2");
    //GtkWidget *label_oil_water_profile_parameter ; label_oil_water_profile_parameter      = gtk_label_new("Oil-Water Profile Parameter");
    //GtkWidget *label_A_bar ; label_A_bar  = gtk_label_new("A'");
    //GtkWidget *label_B1 ; label_B1  = gtk_label_new("B'1");
    //GtkWidget *label_B2 ; label_B2  = gtk_label_new("B'2");

    //GtkObject *adjustment_A; adjustment_A   = gtk_adjustment_new(1.2, -1.0e50, 1.0e50, 1, 1, 0);
    //GtkObject *adjustment_B; adjustment_B   = gtk_adjustment_new(0.3, -1.0e50, 1.0e50, 1, 1, 0);
    //GtkObject *adjustment_Fv; adjustment_Fv = gtk_adjustment_new(1.0, -1.0e50, 1.0e50, 1, 1, 0);

    //GtkWidget *spin_button_A;  spin_button_A   = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_A), 0.1,1);
    //GtkWidget *spin_button_B;  spin_button_B   = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_B), 0.1,1);
    //GtkWidget *spin_button_Fv; spin_button_Fv  = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment_Fv), 0.1,1);  

    //gtk_table_attach(GTK_TABLE(table_drift_flux_parameters), label_gas_liquid_profile_parameter, 0, 2, 0, 1, GTK_SHRINK, GTK_SHRINK, 5, 5);
    //gtk_table_attach(GTK_TABLE(table_drift_flux_parameters), label_A, 0, 1, 1, 2, GTK_SHRINK, GTK_SHRINK, 5, 5);
    //gtk_table_attach(GTK_TABLE(table_drift_flux_parameters), label_B, 0, 1, 2, 3, GTK_SHRINK, GTK_SHRINK, 5, 5);
    //gtk_table_attach(GTK_TABLE(table_drift_flux_parameters), label_Fv, 0, 1, 3, 4, GTK_SHRINK, GTK_SHRINK, 5, 5);
    //gtk_table_attach(GTK_TABLE(table_drift_flux_parameters), spin_button_A, 1, 2, 1, 2, GTK_SHRINK, GTK_SHRINK, 5, 5);
    //gtk_table_attach(GTK_TABLE(table_drift_flux_parameters), spin_button_B, 1, 2, 2, 3, GTK_SHRINK, GTK_SHRINK, 5, 5);
    //gtk_table_attach(GTK_TABLE(table_drift_flux_parameters), spin_button_Fv, 1, 2, 3, 4, GTK_SHRINK, GTK_SHRINK, 5, 5);


    //gtk_table_attach(GTK_TABLE(table_drift_flux_parameters), label_gas_liquid_drift_velocity, 0, 4, 4, 5, GTK_SHRINK, GTK_SHRINK, 5, 5);
    //gtk_table_attach(GTK_TABLE(table_drift_flux_parameters), label_a1, 0, 1, 5, 6, GTK_SHRINK, GTK_SHRINK, 5, 5);
    //gtk_table_attach(GTK_TABLE(table_drift_flux_parameters), label_a2, 0, 1, 6, 7, GTK_SHRINK, GTK_SHRINK, 5, 5);


    //gtk_table_attach(GTK_TABLE(table_drift_flux_parameters), label_oil_water_profile_parameter, 2, 4, 0, 1, GTK_SHRINK, GTK_SHRINK, 5, 5);
    //gtk_table_attach(GTK_TABLE(table_drift_flux_parameters), label_A_bar, 2, 3, 1, 2, GTK_SHRINK, GTK_SHRINK, 5, 5);
    //gtk_table_attach(GTK_TABLE(table_drift_flux_parameters), label_B1, 2, 3, 2, 3, GTK_SHRINK, GTK_SHRINK, 5, 5);
    //gtk_table_attach(GTK_TABLE(table_drift_flux_parameters), label_B2, 2, 3, 3, 4, GTK_SHRINK, GTK_SHRINK, 5, 5);

    // ------------------------------------------------------------------------------

    GtkWidget *run_button;
    run_button = gtk_button_new_with_label("Run");
    gtk_box_pack_start(GTK_BOX(vbox), run_button, TRUE, TRUE, 0);

    GtkWidget *frame;
    frame = gtk_frame_new("Results");

    gtk_widget_set_size_request(GTK_WIDGET(frame), 300, 300);

    gtk_box_pack_start(GTK_BOX(vbox), frame, TRUE, TRUE, 0);

    GtkWidget *notebook; notebook = gtk_notebook_new();
    gtk_container_add(GTK_CONTAINER(frame), notebook);
    gtk_notebook_popup_enable(GTK_NOTEBOOK(notebook));
    gtk_notebook_set_scrollable(GTK_NOTEBOOK(notebook),true);


    GtkWidget *scheme; scheme = gtk_image_new_from_file(".\\images\\scheme.png");
    GtkWidget *frame_scheme; frame_scheme = gtk_frame_new ("Well Scheme"); 
    GtkWidget *frame_scheme_label; frame_scheme_label = gtk_label_new("Well Scheme");
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), frame_scheme, frame_scheme_label);
    gtk_container_add(GTK_CONTAINER(frame_scheme), scheme);

    GtkWidget *graph_pressure;
    graph_pressure = gtk_graph_new(XY);
    GtkWidget *frame_graph_pressure; frame_graph_pressure = gtk_frame_new ("Pressure vs Well Length [m]");  
    GtkWidget *frame_graph_label_pressure; frame_graph_label_pressure = gtk_label_new("Pressure");
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), frame_graph_pressure, frame_graph_label_pressure);
    gtk_container_add(GTK_CONTAINER(frame_graph_pressure), graph_pressure);

    GtkWidget *graph_gas_vol_frac;
    graph_gas_vol_frac = gtk_graph_new(XY);
    GtkWidget *frame_graph_gas_vol_frac; frame_graph_gas_vol_frac = gtk_frame_new ("Gas Volume Fraction vs Well Length [m]");  
    GtkWidget *frame_graph_label_gas_vol_frac; frame_graph_label_gas_vol_frac = gtk_label_new("Gas Volume Fraction");
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), frame_graph_gas_vol_frac, frame_graph_label_gas_vol_frac);
    gtk_container_add(GTK_CONTAINER(frame_graph_gas_vol_frac), graph_gas_vol_frac);

    GtkWidget *graph_oil_vol_frac;
    graph_oil_vol_frac = gtk_graph_new(XY);
    GtkWidget *frame_graph_oil_vol_frac; frame_graph_oil_vol_frac = gtk_frame_new ("Oil Volume Fraction vs Well Length [m]");  
    GtkWidget *frame_graph_label_oil_vol_frac; frame_graph_label_oil_vol_frac = gtk_label_new("Oil Volume Fraction");
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), frame_graph_oil_vol_frac, frame_graph_label_oil_vol_frac);
    gtk_container_add(GTK_CONTAINER(frame_graph_oil_vol_frac), graph_oil_vol_frac);

    GtkWidget *graph_water_vol_frac;
    graph_water_vol_frac = gtk_graph_new(XY);
    GtkWidget *frame_graph_water_vol_frac; frame_graph_water_vol_frac = gtk_frame_new ("Water Volume Fraction vs Well Length [m]");  
    GtkWidget *frame_graph_label_water_vol_frac; frame_graph_label_water_vol_frac = gtk_label_new("Water Volume Fraction");
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), frame_graph_water_vol_frac, frame_graph_label_water_vol_frac);
    gtk_container_add(GTK_CONTAINER(frame_graph_water_vol_frac), graph_water_vol_frac);

    GtkWidget *graph_volumetric_flux;
    graph_volumetric_flux = gtk_graph_new(XY);
    GtkWidget *frame_graph_volumetric_flux; frame_graph_volumetric_flux = gtk_frame_new ("Volumetric Flux vs Well Length [m]");  
    GtkWidget *frame_graph_label_volumetric_flux; frame_graph_label_volumetric_flux = gtk_label_new("Volumetric Flux");
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), frame_graph_volumetric_flux, frame_graph_label_volumetric_flux);
    gtk_container_add(GTK_CONTAINER(frame_graph_volumetric_flux), graph_volumetric_flux);

    GtkWidget *graph_mixture_velocity;
    graph_mixture_velocity = gtk_graph_new(XY);
    GtkWidget *frame_graph_mixture_velocity; frame_graph_mixture_velocity = gtk_frame_new ("Mixture Velocity vs Well Length [m]");  
    GtkWidget *frame_graph_label_mixture_velocity; frame_graph_label_mixture_velocity = gtk_label_new("Mixture Velocity");
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), frame_graph_mixture_velocity, frame_graph_label_mixture_velocity);
    gtk_container_add(GTK_CONTAINER(frame_graph_mixture_velocity), graph_mixture_velocity);

   

    
    
    //gtk_box_pack_start(GTK_BOX(vbox), aGraph, TRUE, TRUE, 0); 

    ////gtk_widget_show( aGraph );

    //gtk_notebook_set_page (GTK_NOTEBOOK(notebook), 0);
    gtk_widget_realize(graph_pressure);
    gtk_widget_realize(graph_gas_vol_frac);
    gtk_widget_realize(graph_oil_vol_frac); 
    gtk_widget_realize(graph_water_vol_frac);
    gtk_widget_realize(graph_mixture_velocity);
    gtk_widget_realize(graph_volumetric_flux);

     

    gtk_widget_set_no_show_all ( graph_pressure, TRUE );
    gtk_widget_set_no_show_all( graph_gas_vol_frac, TRUE );
    gtk_widget_set_no_show_all( graph_oil_vol_frac, TRUE );
    gtk_widget_set_no_show_all( graph_water_vol_frac, TRUE );
    gtk_widget_set_no_show_all( graph_mixture_velocity, TRUE );
    gtk_widget_set_no_show_all( graph_volumetric_flux, TRUE );

    //gtk_widget_set_no_show_all( notebook, TRUE );
    
    gtk_widget_show_all( window );            
    gtk_widget_show( GTK_WIDGET(frame_graph_pressure) );
    //show_about_info();
    g_signal_connect( G_OBJECT(window), "destroy" , G_CALLBACK(quit_program), NULL );
    g_signal_connect( G_OBJECT(quit  ), "activate", G_CALLBACK(quit_program), NULL );
    g_signal_connect( G_OBJECT(about  ), "activate", G_CALLBACK(show_about_info), NULL );

    /*Plot aPlot;
    aPlot.m_gtk_graph       = GTK_GRAPH(aGraph);
    aPlot.m_entry_x_min     = entry_x_min;
    aPlot.m_entry_x_max     = entry_x_max;
    aPlot.m_entry_y_min     = entry_y_min;
    aPlot.m_entry_y_max     = entry_y_max;
    aPlot.m_frame           = frame;
    aPlot.m_window          = window;*/
    SimulatorData simulator_data;
    simulator_data.m_notebook                       = notebook;
    simulator_data.m_graph_pressure                 = GTK_GRAPH(graph_pressure);
    simulator_data.m_graph_gas_vol_frac             = GTK_GRAPH(graph_gas_vol_frac);
    simulator_data.m_graph_oil_vol_frac             = GTK_GRAPH(graph_oil_vol_frac);
    simulator_data.m_graph_water_vol_frac           = GTK_GRAPH(graph_water_vol_frac);
    simulator_data.m_graph_mixture_velocity         = GTK_GRAPH(graph_mixture_velocity);
    simulator_data.m_graph_volumetric_flux          = GTK_GRAPH(graph_volumetric_flux);

    simulator_data.m_curve_id_pressure          = gtk_graph_trace_new(GTK_GRAPH(graph_pressure));
    simulator_data.m_curve_id_gas_vol_frac      = gtk_graph_trace_new(GTK_GRAPH(graph_gas_vol_frac));
    simulator_data.m_curve_id_oil_vol_frac      = gtk_graph_trace_new(GTK_GRAPH(graph_oil_vol_frac));
    simulator_data.m_curve_id_water_vol_frac    = gtk_graph_trace_new(GTK_GRAPH(graph_water_vol_frac));
    simulator_data.m_curve_id_mixture_velocity  = gtk_graph_trace_new(GTK_GRAPH(graph_mixture_velocity));
    simulator_data.m_curve_id_gas_volumetric_flux       = gtk_graph_trace_new(GTK_GRAPH(graph_volumetric_flux));
    simulator_data.m_curve_id_oil_volumetric_flux       = gtk_graph_trace_new(GTK_GRAPH(graph_volumetric_flux));
    simulator_data.m_curve_id_water_volumetric_flux     = gtk_graph_trace_new(GTK_GRAPH(graph_volumetric_flux));

    
    simulator_data.m_button_file_chooser_well_inflow = file_chooser_button_inflow;              
    simulator_data.m_button_number_of_nodes         = spin_button_number_of_nodes;
    simulator_data.m_button_well_length             = spin_button_well_length;
    simulator_data.m_button_well_diameter           = spin_button_well_diameter;
    simulator_data.m_button_well_inclination        = spin_button_well_inclination;
    simulator_data.m_button_initial_pressure        = spin_button_initial_pressure;
    simulator_data.m_button_initial_velocity        = spin_button_initial_velocity;
    simulator_data.m_button_initial_gas_vol_frac    = spin_button_initial_gas_vol_frac;
    simulator_data.m_button_initial_oil_vol_frac    = spin_button_initial_oil_vol_frac;
    simulator_data.m_button_toe_velocity            = spin_button_toe_velocity;
    simulator_data.m_button_heel_pressure           = spin_button_heel_pressure;
    simulator_data.m_button_timestep                = spin_button_timestep;
    simulator_data.m_button_max_delta_t             = spin_button_max_delta_t;
    simulator_data.m_button_final_time              = spin_button_final_time;
    simulator_data.m_button_tolerance               = spin_button_tolerance;
    simulator_data.m_well_initial_data              = well_initial_data;

    g_signal_connect(run_button, "clicked", G_CALLBACK(run), &simulator_data); 
    show_about_info();
    gtk_main ();    
   // std::cout.rdbuf(sb);   
    return 0;
}