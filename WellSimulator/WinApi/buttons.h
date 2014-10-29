#ifndef BUTTONS_H
#define BUTTONS_H

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
#include <simulator.h>
#include <algorithm>

GdkPixbuf *create_pixbuf(const gchar * filename);

//using namespace std;
class outbuf : public std::streambuf {
public:
    outbuf() {
        setp(0, 0);
    }

    virtual int_type overflow(int_type c = traits_type::eof()) {
        return fputc(c, stdout) == EOF ? traits_type::eof() : c;
    }
};

void quit_program(){
    gtk_main_quit();
    exit(0);
}

// This function needs to be implemented for each OS, e.g, MAC, WINDOWS, LINUX
static gboolean open_url_hook(GtkLabel* label, gchar* uri, gpointer user_data)
{      
    ShellExecute(NULL,"open",uri,NULL,NULL,SW_SHOWNORMAL);
    return true;    
}

void show_about_info(){ 
    GdkPixbuf *pixbuf = gdk_pixbuf_new_from_file(".\\images\\well_drift_logo.png", NULL);

    GtkWidget *dialog = gtk_about_dialog_new();
    gtk_window_set_icon(GTK_WINDOW(dialog), create_pixbuf(".\\images\\icon.png"));
    gtk_about_dialog_set_name(GTK_ABOUT_DIALOG(dialog), "WellDrift Simulator");
    gtk_about_dialog_set_version(GTK_ABOUT_DIALOG(dialog), "v0.1"); 
    gtk_about_dialog_set_copyright(GTK_ABOUT_DIALOG(dialog),"email: arthursoprano@gmail.com\nDepartament of Mechanical Engineering\nFederal University of Santa Catarina - UFSC\nFlorianopolis - SC - Brazil");
    
    const char *authors[]= {"Arthur B. Soprano\t-\tM.Sc. Mechanical Engineer",
                            "A. Fabio C. da Silva\t-\tProfessor",
                            "Clovis R. Maliska\t-\tProfessor", NULL};
   
    gtk_about_dialog_set_authors(GTK_ABOUT_DIALOG(dialog), (const gchar**)authors);
    gtk_about_dialog_set_comments(GTK_ABOUT_DIALOG(dialog), 
        "WellDrift is a three-phase, one-dimensional, drift-flux well simulator. The problem is solved using the Finite Volume Method and a Newton-Raphson algorithm to solve the coupled equations.");      
    
    g_signal_connect( G_OBJECT(dialog  ), "activate-link", G_CALLBACK(open_url_hook), NULL);                                             
    
    gtk_about_dialog_set_website(GTK_ABOUT_DIALOG(dialog),"http://www.sinmec.ufsc.br");
    gtk_about_dialog_set_logo(GTK_ABOUT_DIALOG(dialog), pixbuf);
    g_object_unref(pixbuf), pixbuf = NULL;
    gtk_dialog_run(GTK_DIALOG (dialog));
    gtk_widget_destroy(dialog);
}

GdkPixbuf *create_pixbuf(const gchar * filename)
{
    GdkPixbuf *pixbuf;
    GError *error = NULL;
    pixbuf = gdk_pixbuf_new_from_file(filename, &error);
    if(!pixbuf) {
        fprintf(stderr, "%s\n", error->message);
        g_error_free(error);
    }

    return pixbuf;
}

static gint counter = 0;

void increase(GtkWidget *widget, gpointer label)
{
    counter++;
    std::stringstream buf;
    buf << counter;
    gtk_label_set_text(GTK_LABEL(label), buf.str().c_str());
}

void decrease(GtkWidget *widget, gpointer label)
{
    counter--;     
    std::stringstream buf;
    buf << counter;
    gtk_label_set_text(GTK_LABEL(label), buf.str().c_str());
}

struct Plot{
    GtkGraph*   m_gtk_graph;
    GtkWidget*  m_entry_x_min;
    GtkWidget*  m_entry_x_max;
    GtkWidget*  m_entry_y_min;
    GtkWidget*  m_entry_y_max;
    GtkWidget*  m_frame;
    GtkWidget*  m_window;
    //GtkWidget*  m_entry_y_max;
};


struct SimulatorData{
    GtkWidget*  m_notebook;
    GtkGraph*   m_graph_pressure;
    GtkGraph*   m_graph_gas_vol_frac;
    GtkGraph*   m_graph_oil_vol_frac;
    GtkGraph*   m_graph_water_vol_frac;
    GtkGraph*   m_graph_mixture_velocity;
    GtkGraph*   m_graph_gas_velocity;
    GtkGraph*   m_graph_oil_velocity;
    GtkGraph*   m_graph_water_velocity;
    GtkGraph*   m_graph_volumetric_flux;

    int         m_curve_id_pressure;
    int         m_curve_id_gas_vol_frac;
    int         m_curve_id_oil_vol_frac;
    int         m_curve_id_water_vol_frac;
    int         m_curve_id_mixture_velocity;
    int         m_curve_id_gas_volumetric_flux;
    int         m_curve_id_oil_volumetric_flux;
    int         m_curve_id_water_volumetric_flux;
    GtkWidget*  m_button_file_chooser_well_inflow;
    GtkWidget*  m_button_number_of_nodes; 
    GtkWidget*  m_button_well_length;
    GtkWidget*  m_button_well_diameter;
    GtkWidget*  m_button_well_inclination; 
    GtkWidget*  m_button_initial_pressure;
    GtkWidget*  m_button_initial_velocity;
    GtkWidget*  m_button_initial_gas_vol_frac;
    GtkWidget*  m_button_initial_oil_vol_frac;
    GtkWidget*  m_button_toe_velocity;
    GtkWidget*  m_button_heel_pressure;
    GtkWidget*  m_button_timestep;
    GtkWidget*  m_button_max_delta_t;
    GtkWidget*  m_button_final_time;
    GtkWidget*  m_button_tolerance;
    SharedPointer<WellInitialData> m_well_initial_data;
};


void run(GtkWidget *widget, gpointer p_data)
{   
    SimulatorData* _data = reinterpret_cast< SimulatorData* >(p_data);   
    float64 degrees_to_radians = 0.0174532925; 
    //CurveID = 0;
    _data->m_well_initial_data->m_number_of_nodes           = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(_data->m_button_number_of_nodes     )); 
    _data->m_well_initial_data->m_well_length               = gtk_spin_button_get_value(GTK_SPIN_BUTTON(_data->m_button_well_length         ));
    _data->m_well_initial_data->m_well_diameter             = gtk_spin_button_get_value(GTK_SPIN_BUTTON(_data->m_button_well_diameter       ));
    _data->m_well_initial_data->m_well_inclination          = degrees_to_radians*gtk_spin_button_get_value(GTK_SPIN_BUTTON(_data->m_button_well_inclination    )); 
    _data->m_well_initial_data->m_initial_pressure          = gtk_spin_button_get_value(GTK_SPIN_BUTTON(_data->m_button_initial_pressure    ));
    _data->m_well_initial_data->m_initial_mixture_velocity  = gtk_spin_button_get_value(GTK_SPIN_BUTTON(_data->m_button_initial_velocity    ));
    _data->m_well_initial_data->m_initial_gas_vol_frac      = gtk_spin_button_get_value(GTK_SPIN_BUTTON(_data->m_button_initial_gas_vol_frac));
    _data->m_well_initial_data->m_initial_oil_vol_frac      = gtk_spin_button_get_value(GTK_SPIN_BUTTON(_data->m_button_initial_oil_vol_frac));
    _data->m_well_initial_data->m_toe_mixture_velocity      = gtk_spin_button_get_value(GTK_SPIN_BUTTON(_data->m_button_toe_velocity        ));
    _data->m_well_initial_data->m_heel_pressure             = gtk_spin_button_get_value(GTK_SPIN_BUTTON(_data->m_button_heel_pressure       ));
    _data->m_well_initial_data->m_delta_t                   = gtk_spin_button_get_value(GTK_SPIN_BUTTON(_data->m_button_timestep            ));
    _data->m_well_initial_data->m_max_delta_t               = gtk_spin_button_get_value(GTK_SPIN_BUTTON(_data->m_button_max_delta_t         ));
    _data->m_well_initial_data->m_final_time                = gtk_spin_button_get_value(GTK_SPIN_BUTTON(_data->m_button_final_time          ));
    _data->m_well_initial_data->m_tolerance                 = gtk_spin_button_get_value(GTK_SPIN_BUTTON(_data->m_button_tolerance           ));


    SharedPointer<vector_type> temp_oil_inflow     (new vector_type(_data->m_well_initial_data->m_number_of_nodes, 0.0));
    SharedPointer<vector_type> temp_water_inflow   (new vector_type(_data->m_well_initial_data->m_number_of_nodes, 0.0));
    SharedPointer<vector_type> temp_gas_inflow     (new vector_type(_data->m_well_initial_data->m_number_of_nodes, 0.0));

    _data->m_well_initial_data->m_oil_inflow    = temp_oil_inflow;
    _data->m_well_initial_data->m_water_inflow  = temp_water_inflow;
    _data->m_well_initial_data->m_gas_inflow    = temp_gas_inflow;

    char* temp = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER (_data->m_button_file_chooser_well_inflow)); 
    std::string well_name_path;       
    if( temp == NULL){
         well_name_path = "..\\WellData\\Inflow.txt";
    }
    else{
        well_name_path = temp;
    }
    std::ifstream inflow_file( well_name_path.c_str() );


    inflow_file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    int index = 0;
    while( !inflow_file.eof() && index < _data->m_well_initial_data->m_number_of_nodes){
        inflow_file >> (*(_data->m_well_initial_data->m_oil_inflow))[index];
        inflow_file >> (*(_data->m_well_initial_data->m_water_inflow))[index];
        inflow_file >> (*(_data->m_well_initial_data->m_gas_inflow))[index];
        ++index;
    }
    

    SharedPointer<DriftFluxWell> well = simulate(_data->m_well_initial_data);

    // For research only
    //SharedPointer<DriftFluxWell> well = simulate_provenzano(_data->m_well_initial_data);

    // Function definition

    int number_of_points = well->number_of_nodes();
    coord_type well_coordinates;
    std::vector<float> x_coord(number_of_points);
    std::vector<float> result_pressure(number_of_points);
    std::vector<float> result_gas_vol_frac(number_of_points);
    std::vector<float> result_oil_vol_frac(number_of_points);
    std::vector<float> result_water_vol_frac(number_of_points);
    std::vector<float> result_mixture_velocity(number_of_points);
    std::vector<float> result_gas_volumetric_flux(number_of_points);
    std::vector<float> result_oil_volumetric_flux(number_of_points);
    std::vector<float> result_water_volumetric_flux(number_of_points);

    std::vector<float> y_coord(number_of_points);
    for(int i = 0; i < number_of_points; ++i){
        x_coord[i] = well->coordinates(i)->getX();
        result_pressure[i]          = *(well->pressure(i));
        result_gas_vol_frac[i]      = well->get_gas_volume_fraction(i);
        result_oil_vol_frac[i]      = well->get_oil_volume_fraction(i); 
        result_water_vol_frac[i]    = well->get_water_volume_fraction(i); 
        result_mixture_velocity[i]  = well->get_mixture_velocity(i); 
        result_gas_volumetric_flux[i]   = well->get_gas_volumetric_flux(i);
        result_oil_volumetric_flux[i]   = well->get_oil_volumetric_flux(i);
        result_water_volumetric_flux[i] = well->get_water_volumetric_flux(i);
    }
    float x_coord_min = *( std::min_element(x_coord.begin(),x_coord.end()) ); 
    float x_coord_max = *( std::max_element(x_coord.begin(),x_coord.end()) );
    float result_pressure_min = *( std::min_element(result_pressure.begin(),result_pressure.end()) ); 
    float result_pressure_max = *( std::max_element(result_pressure.begin(),result_pressure.end()) );
    float result_gas_vol_frac_min   = *( std::min_element(result_gas_vol_frac.begin(),result_gas_vol_frac.end()) ); 
    float result_gas_vol_frac_max   = *( std::max_element(result_gas_vol_frac.begin(),result_gas_vol_frac.end()) );
    float result_oil_vol_frac_min   = *( std::min_element(result_oil_vol_frac.begin(),result_oil_vol_frac.end()) ); 
    float result_oil_vol_frac_max   = *( std::max_element(result_oil_vol_frac.begin(),result_oil_vol_frac.end()) );
    float result_water_vol_frac_min     = *( std::min_element(result_water_vol_frac.begin(),result_water_vol_frac.end()) ); 
    float result_water_vol_frac_max     = *( std::max_element(result_water_vol_frac.begin(),result_water_vol_frac.end()) );
    float result_mixture_velocity_min   = *( std::min_element(result_mixture_velocity.begin(),result_mixture_velocity.end()) ); 
    float result_mixture_velocity_max   = *( std::max_element(result_mixture_velocity.begin(),result_mixture_velocity.end()) );
    float result_gas_volumetric_flux_min        = *( std::min_element(result_gas_volumetric_flux.begin(),result_gas_volumetric_flux.end()) ); 
    float result_oil_volumetric_flux_min        = *( std::min_element(result_oil_volumetric_flux.begin(),result_oil_volumetric_flux.end()) ); 
    float result_water_volumetric_flux_min      = *( std::min_element(result_water_volumetric_flux.begin(),result_water_volumetric_flux.end()) ); 
    float result_gas_volumetric_flux_max        = *( std::max_element(result_gas_volumetric_flux.begin(),result_gas_volumetric_flux.end()) ); 
    float result_oil_volumetric_flux_max        = *( std::max_element(result_oil_volumetric_flux.begin(),result_oil_volumetric_flux.end()) ); 
    float result_water_volumetric_flux_max      = *( std::max_element(result_water_volumetric_flux.begin(),result_water_volumetric_flux.end()) ); 
    float result_volumetric_flux_min = std::min(result_gas_volumetric_flux_min, std::min(result_oil_volumetric_flux_min, result_water_volumetric_flux_min) );
    float result_volumetric_flux_max = std::max(result_gas_volumetric_flux_max, std::max(result_oil_volumetric_flux_max, result_water_volumetric_flux_max) );

    float x_min = x_coord_min;
    float x_max = x_coord_max; 

    
    gtk_graph_set_title (_data->m_graph_pressure, "Pressure", NULL);
    gtk_graph_legend_format(_data->m_graph_pressure, false, GTK_GRAPH_NORTH);
    gtk_graph_trace_set_data        (_data->m_graph_pressure, _data->m_curve_id_pressure, &x_coord[0], &result_pressure[0], number_of_points); 
    gtk_graph_trace_format_line     (_data->m_graph_pressure, _data->m_curve_id_pressure, SOLID, 2, &RED, 1);
    gtk_graph_trace_format_marker   (_data->m_graph_pressure, _data->m_curve_id_pressure, GTK_GRAPH_MARKER_NONE, 1, &RED, &RED, false);
    gtk_graph_axis_set_limits       (_data->m_graph_pressure, GTK_GRAPH_AXIS_INDEPENDANT, x_max, x_min);
    gtk_graph_axis_set_limits       (_data->m_graph_pressure, GTK_GRAPH_AXIS_DEPENDANT, result_pressure_max, result_pressure_min);
    gtk_graph_axis_format           (_data->m_graph_pressure, GTK_GRAPH_AXIS_INDEPENDANT, FLOATING_POINT, 1, "Well Distance [m]");
    gtk_graph_axis_format           (_data->m_graph_pressure, GTK_GRAPH_AXIS_DEPENDANT, FLOATING_POINT, 1, "[Pa]");
    gtk_graph_axis_set_crossing     (_data->m_graph_pressure, GTK_GRAPH_AXIS_INDEPENDANT, GTK_GRAPH_USERVALUE, 0.0);
    gtk_graph_axis_set_crossing     (_data->m_graph_pressure, GTK_GRAPH_AXIS_DEPENDANT, GTK_GRAPH_USERVALUE, 0.0);
    gtk_graph_axis_format_grid      (_data->m_graph_pressure, GTK_GRAPH_AXIS_INDEPENDANT, true);
    gtk_graph_axis_format_grid      (_data->m_graph_pressure, GTK_GRAPH_AXIS_DEPENDANT, true);


    gtk_graph_set_title (_data->m_graph_gas_vol_frac, "Gas Volume Fraction", NULL);
    gtk_graph_legend_format(_data->m_graph_gas_vol_frac, false, GTK_GRAPH_NORTH);
    gtk_graph_trace_set_data        (_data->m_graph_gas_vol_frac, _data->m_curve_id_gas_vol_frac, &x_coord[0], &result_gas_vol_frac[0], number_of_points); 
    gtk_graph_trace_format_line     (_data->m_graph_gas_vol_frac, _data->m_curve_id_gas_vol_frac, SOLID, 2, &ORANGE, 1);
    gtk_graph_trace_format_marker   (_data->m_graph_gas_vol_frac, _data->m_curve_id_gas_vol_frac, GTK_GRAPH_MARKER_NONE, 1, &RED, &RED, false);
    gtk_graph_axis_set_limits       (_data->m_graph_gas_vol_frac, GTK_GRAPH_AXIS_INDEPENDANT, x_max, x_min);
    gtk_graph_axis_set_limits       (_data->m_graph_gas_vol_frac, GTK_GRAPH_AXIS_DEPENDANT, result_gas_vol_frac_max, result_gas_vol_frac_min);
    gtk_graph_axis_format           (_data->m_graph_gas_vol_frac, GTK_GRAPH_AXIS_INDEPENDANT, FLOATING_POINT, 1, "Well Distance [m]");
    gtk_graph_axis_format           (_data->m_graph_gas_vol_frac, GTK_GRAPH_AXIS_DEPENDANT, SCIENTIFIC, 2, "[-]");
    gtk_graph_axis_set_crossing     (_data->m_graph_gas_vol_frac, GTK_GRAPH_AXIS_INDEPENDANT, GTK_GRAPH_USERVALUE, 0.0);
    gtk_graph_axis_set_crossing     (_data->m_graph_gas_vol_frac, GTK_GRAPH_AXIS_DEPENDANT, GTK_GRAPH_USERVALUE, 0.0);
    gtk_graph_axis_format_grid      (_data->m_graph_gas_vol_frac, GTK_GRAPH_AXIS_INDEPENDANT, true);
    gtk_graph_axis_format_grid      (_data->m_graph_gas_vol_frac, GTK_GRAPH_AXIS_DEPENDANT, true);

    gtk_graph_set_title (_data->m_graph_oil_vol_frac, "Oil Volume Fraction", NULL);
    gtk_graph_legend_format(_data->m_graph_oil_vol_frac, false, GTK_GRAPH_NORTH);
    gtk_graph_trace_set_data        (_data->m_graph_oil_vol_frac, _data->m_curve_id_oil_vol_frac, &x_coord[0], &result_oil_vol_frac[0], number_of_points); 
    gtk_graph_trace_format_line     (_data->m_graph_oil_vol_frac, _data->m_curve_id_oil_vol_frac, SOLID, 2, &BLACK, 1);
    gtk_graph_trace_format_marker   (_data->m_graph_oil_vol_frac, _data->m_curve_id_oil_vol_frac, GTK_GRAPH_MARKER_NONE, 1, &RED, &RED, false);
    gtk_graph_axis_set_limits       (_data->m_graph_oil_vol_frac, GTK_GRAPH_AXIS_INDEPENDANT, x_max, x_min);
    gtk_graph_axis_set_limits       (_data->m_graph_oil_vol_frac, GTK_GRAPH_AXIS_DEPENDANT, result_oil_vol_frac_max, result_oil_vol_frac_min);
    gtk_graph_axis_format           (_data->m_graph_oil_vol_frac, GTK_GRAPH_AXIS_INDEPENDANT, FLOATING_POINT, 1, "Well Distance [m]");
    gtk_graph_axis_format           (_data->m_graph_oil_vol_frac, GTK_GRAPH_AXIS_DEPENDANT, SCIENTIFIC, 2, "[-]");
    gtk_graph_axis_set_crossing     (_data->m_graph_oil_vol_frac, GTK_GRAPH_AXIS_INDEPENDANT, GTK_GRAPH_USERVALUE, 0.0);
    gtk_graph_axis_set_crossing     (_data->m_graph_oil_vol_frac, GTK_GRAPH_AXIS_DEPENDANT, GTK_GRAPH_USERVALUE, 0.0);
    gtk_graph_axis_format_grid      (_data->m_graph_oil_vol_frac, GTK_GRAPH_AXIS_INDEPENDANT, true);
    gtk_graph_axis_format_grid      (_data->m_graph_oil_vol_frac, GTK_GRAPH_AXIS_DEPENDANT, true);

    gtk_graph_set_title (_data->m_graph_water_vol_frac, "Water Volume Fraction", NULL);
    gtk_graph_legend_format(_data->m_graph_water_vol_frac, false, GTK_GRAPH_NORTH);
    gtk_graph_trace_set_data        (_data->m_graph_water_vol_frac, _data->m_curve_id_water_vol_frac, &x_coord[0], &result_water_vol_frac[0], number_of_points); 
    gtk_graph_trace_format_line     (_data->m_graph_water_vol_frac, _data->m_curve_id_water_vol_frac, SOLID, 2, &BLUE, 1);
    gtk_graph_trace_format_marker   (_data->m_graph_water_vol_frac, _data->m_curve_id_water_vol_frac, GTK_GRAPH_MARKER_NONE, 1, &RED, &RED, false);
    gtk_graph_axis_set_limits       (_data->m_graph_water_vol_frac, GTK_GRAPH_AXIS_INDEPENDANT, x_max, x_min);
    gtk_graph_axis_set_limits       (_data->m_graph_water_vol_frac, GTK_GRAPH_AXIS_DEPENDANT, result_water_vol_frac_max, result_water_vol_frac_min);
    gtk_graph_axis_format           (_data->m_graph_water_vol_frac, GTK_GRAPH_AXIS_INDEPENDANT, FLOATING_POINT, 1, "Well Distance [m]");
    gtk_graph_axis_format           (_data->m_graph_water_vol_frac, GTK_GRAPH_AXIS_DEPENDANT, SCIENTIFIC, 2, "[-]");
    gtk_graph_axis_set_crossing     (_data->m_graph_water_vol_frac, GTK_GRAPH_AXIS_INDEPENDANT, GTK_GRAPH_USERVALUE, 0.0);
    gtk_graph_axis_set_crossing     (_data->m_graph_water_vol_frac, GTK_GRAPH_AXIS_DEPENDANT, GTK_GRAPH_USERVALUE, 0.0);
    gtk_graph_axis_format_grid      (_data->m_graph_water_vol_frac, GTK_GRAPH_AXIS_INDEPENDANT, true);
    gtk_graph_axis_format_grid      (_data->m_graph_water_vol_frac, GTK_GRAPH_AXIS_DEPENDANT, true);


    gtk_graph_set_title (_data->m_graph_mixture_velocity, "Mixture Velocity", NULL);
    gtk_graph_legend_format(_data->m_graph_mixture_velocity, false, GTK_GRAPH_NORTH);
    gtk_graph_trace_set_data        (_data->m_graph_mixture_velocity, _data->m_curve_id_mixture_velocity, &x_coord[0], &result_mixture_velocity[0], number_of_points); 
    gtk_graph_trace_format_line     (_data->m_graph_mixture_velocity, _data->m_curve_id_mixture_velocity, SOLID, 2, &DARK_GREY, 1);
    gtk_graph_trace_format_marker   (_data->m_graph_mixture_velocity, _data->m_curve_id_mixture_velocity, GTK_GRAPH_MARKER_NONE, 1, &RED, &RED, false);
    gtk_graph_axis_set_limits       (_data->m_graph_mixture_velocity, GTK_GRAPH_AXIS_INDEPENDANT, x_max, x_min);
    gtk_graph_axis_set_limits       (_data->m_graph_mixture_velocity, GTK_GRAPH_AXIS_DEPENDANT, result_mixture_velocity_max, result_mixture_velocity_min);
    gtk_graph_axis_format           (_data->m_graph_mixture_velocity, GTK_GRAPH_AXIS_INDEPENDANT, FLOATING_POINT, 1, "Well Distance [m]");
    gtk_graph_axis_format           (_data->m_graph_mixture_velocity, GTK_GRAPH_AXIS_DEPENDANT, SCIENTIFIC, 3, "[m/s]");
    gtk_graph_axis_set_crossing     (_data->m_graph_mixture_velocity, GTK_GRAPH_AXIS_INDEPENDANT, GTK_GRAPH_USERVALUE, 0.0);
    gtk_graph_axis_set_crossing     (_data->m_graph_mixture_velocity, GTK_GRAPH_AXIS_DEPENDANT, GTK_GRAPH_USERVALUE, 0.0);
    gtk_graph_axis_format_grid      (_data->m_graph_mixture_velocity, GTK_GRAPH_AXIS_INDEPENDANT, true);
    gtk_graph_axis_format_grid      (_data->m_graph_mixture_velocity, GTK_GRAPH_AXIS_DEPENDANT, true);

    gtk_graph_set_title (_data->m_graph_volumetric_flux, "Volumetric Flux [m^3/s]", NULL);
    gtk_graph_legend_format(_data->m_graph_volumetric_flux, false, GTK_GRAPH_NORTH_EAST);
    
    // Gas flux curve
    gtk_graph_trace_format_title    (_data->m_graph_volumetric_flux, _data->m_curve_id_gas_volumetric_flux, "Gas");
    gtk_graph_trace_set_data        (_data->m_graph_volumetric_flux, _data->m_curve_id_gas_volumetric_flux, &x_coord[0], &result_gas_volumetric_flux[0], number_of_points); 
    gtk_graph_trace_format_line     (_data->m_graph_volumetric_flux, _data->m_curve_id_gas_volumetric_flux, SOLID, 2, &ORANGE, 1);
    gtk_graph_trace_format_marker   (_data->m_graph_volumetric_flux, _data->m_curve_id_gas_volumetric_flux, GTK_GRAPH_MARKER_NONE, 1, &RED, &RED, false);
    // Oil flux curve 
    gtk_graph_trace_format_title    (_data->m_graph_volumetric_flux, _data->m_curve_id_oil_volumetric_flux, "Oil"); 
    gtk_graph_trace_set_data        (_data->m_graph_volumetric_flux, _data->m_curve_id_oil_volumetric_flux, &x_coord[0], &result_oil_volumetric_flux[0], number_of_points); 
    gtk_graph_trace_format_line     (_data->m_graph_volumetric_flux, _data->m_curve_id_oil_volumetric_flux, SOLID, 2, &BLACK, 1);
    gtk_graph_trace_format_marker   (_data->m_graph_volumetric_flux, _data->m_curve_id_oil_volumetric_flux, GTK_GRAPH_MARKER_NONE, 1, &RED, &RED, false);
    // Water flux curve
    gtk_graph_trace_format_title    (_data->m_graph_volumetric_flux, _data->m_curve_id_water_volumetric_flux, "Water");
    gtk_graph_trace_set_data        (_data->m_graph_volumetric_flux, _data->m_curve_id_water_volumetric_flux, &x_coord[0], &result_water_volumetric_flux[0], number_of_points); 
    gtk_graph_trace_format_line     (_data->m_graph_volumetric_flux, _data->m_curve_id_water_volumetric_flux, SOLID, 2, &BLUE, 1);
    gtk_graph_trace_format_marker   (_data->m_graph_volumetric_flux, _data->m_curve_id_water_volumetric_flux, GTK_GRAPH_MARKER_NONE, 1, &RED, &RED, false);

    gtk_graph_axis_set_limits       (_data->m_graph_volumetric_flux, GTK_GRAPH_AXIS_INDEPENDANT, x_max, x_min);
    gtk_graph_axis_set_limits       (_data->m_graph_volumetric_flux, GTK_GRAPH_AXIS_DEPENDANT, result_volumetric_flux_max, result_volumetric_flux_min);
    gtk_graph_axis_format           (_data->m_graph_volumetric_flux, GTK_GRAPH_AXIS_INDEPENDANT, FLOATING_POINT, 1, "Well Distance [m]");
    gtk_graph_axis_format           (_data->m_graph_volumetric_flux, GTK_GRAPH_AXIS_DEPENDANT, SCIENTIFIC, 3, "[m^3/s]");
    gtk_graph_axis_set_crossing     (_data->m_graph_volumetric_flux, GTK_GRAPH_AXIS_INDEPENDANT, GTK_GRAPH_USERVALUE, 0.0);
    gtk_graph_axis_set_crossing     (_data->m_graph_volumetric_flux, GTK_GRAPH_AXIS_DEPENDANT, GTK_GRAPH_USERVALUE, 0.0);
    gtk_graph_axis_format_grid      (_data->m_graph_volumetric_flux, GTK_GRAPH_AXIS_INDEPENDANT, true);
    gtk_graph_axis_format_grid      (_data->m_graph_volumetric_flux, GTK_GRAPH_AXIS_DEPENDANT, true);
  
    gtk_widget_show( GTK_WIDGET(_data->m_graph_pressure) );
    gtk_widget_show( GTK_WIDGET(_data->m_graph_gas_vol_frac) );
    gtk_widget_show( GTK_WIDGET(_data->m_graph_oil_vol_frac) );
    gtk_widget_show( GTK_WIDGET(_data->m_graph_water_vol_frac) );
    gtk_widget_show( GTK_WIDGET(_data->m_graph_mixture_velocity) );
    gtk_widget_show( GTK_WIDGET(_data->m_graph_volumetric_flux) );

    gtk_widget_queue_draw_area(GTK_WIDGET(_data->m_graph_pressure), 0, 0, 800, 800);
    gtk_widget_queue_draw_area(GTK_WIDGET(_data->m_graph_gas_vol_frac), 0, 0, 800, 800);
    gtk_widget_queue_draw_area(GTK_WIDGET(_data->m_graph_oil_vol_frac), 0, 0, 800, 800);
    gtk_widget_queue_draw_area(GTK_WIDGET(_data->m_graph_water_vol_frac), 0, 0, 800, 800);
    gtk_widget_queue_draw_area(GTK_WIDGET(_data->m_graph_mixture_velocity), 0, 0, 800, 800);
    gtk_widget_queue_draw_area(GTK_WIDGET(_data->m_graph_volumetric_flux), 0, 0, 800, 800);

    
    std::cout << "done!\n";    
    gtk_notebook_set_current_page(GTK_NOTEBOOK(_data->m_notebook), 1);
}



#endif // BUTTONS_H