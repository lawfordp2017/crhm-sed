// 11/20/19
//---------------------------------------------------------------------------
#ifndef Sed_Overland
#define Sed_Overland
//---------------------------------------------------------------------------

#include "../../core/ClassModule.h"
#include "WQ_CRHM.h"

using namespace std;


struct sed_triple {
    double c;
    double z;
    double s;
};


class ClassSed_Overland : public ClassModule {
    public:

    long dayno{0};

    ClassSed_Overland(string Name, string Version = "undefined", 
                            LMODULE Lvl = LMODULE::PROTO) : ClassModule(Name, Version, Lvl) {};

    ClassSed_Overland* klone(string name) const;

    void decl(void);
    void init(void);
    void run(void);
    void finish(bool good);


/*****************************************************************
 * HYPE Routines
******************************************************************/

    void runoff_sed_by_erosion();


/*****************************************************************
 * Modified [or Daily] Morgan-Morgan-Finney Routines
******************************************************************/

    double calc_flowvel_manning(double depth, double n);
    double calc_flowvel_veg(double depth, double n);
    double calc_tillage_n();

    double calc_rainsplash_energy(double intensity);
    void calc_rainsplash_mobilization( sed_triple &rslt);
    void calc_flow_mobilization(double runoff, sed_triple &rslt);    // runoff: mm/int
    void calc_flow_deposition( bool mixed_sed, double flow_vel, sed_triple &rslt);
    void calc_delivered_to_transport(sed_triple &rslt);
    void calc_transport_capacity(sed_triple &rslt);
    void calc_mmf_sed_balance(sed_triple &rslt);
    double calc_erosion();
    void initialize_modMMF();


/*****************************************************************
 * VARIABLES
******************************************************************/

//    const double *runoff{ NULL };
    const double *scf{ NULL };

// Vars for the modified MMF formulation
    double *pct_sand{ NULL }; // ()

    double *DEP_immed_c{ NULL }; // (%)
    double *DEP_immed_z{ NULL }; // (%)
    double *DEP_immed_s{ NULL }; // (%)
    double *DEP_channel_c{ NULL }; // (%)
    double *DEP_channel_z{ NULL }; // (%)
    double *DEP_channel_s{ NULL }; // (%)

    double *v_std_bare{ NULL };  // m/s
    double *v_actual{ NULL };  // m/s
    double *v_veg{ NULL };  // m/s
    double *v_tillage{ NULL };  // m/s

    double *mob_rainsplash{ NULL };  // m/s
    double *mob_flow{ NULL };  // m/s
    double *sed_delivered_z{ NULL };  // m/s
    double *sed_transported_z{ NULL };  // m/s

// ******

    const double *soil_runoff{ NULL };   // mm (*km^2)/int
    double *soil_runoff_mWQ{ NULL };   // (g/int)
    double **soil_runoff_mWQ_lay{ NULL };   // (g/int)

    const double *net_rain{ NULL };

// Vars for the HYPE delay pool
    double *sedrelpool{ NULL }; // (Particulate P release from soil due to erosion (g/m^2)

// Dummy variables
    double *conc_soil_rechr{ NULL };   // concentration of contaminant (always zero, added to support interfacing to WQ_netroute)
    double **conc_soil_rechr_lay{ NULL };
    double *conc_soil_lower{ NULL };   // concentration of contaminant (always zero, added to support interfacing to WQ_netroute)
    double **conc_soil_lower_lay{ NULL };
   
/*******************
 * PARAMETERS
 *******************/

/* Parameters for the HYPE delay pool */
    const double *sedrelmax{ NULL }; // (mm)
    const double *sedrelexp{ NULL }; // ()
    const double *sreroexp{ NULL }; // ()

/* Parameters for modified MMF formulation */
    const double *mmf_mannings_n{ NULL }; // ()
    const double *mmf_repr_depth{ NULL }; // ()
    const double *sfcslope{ NULL }; // (m/m)
    const double *slopelen{ NULL }; // (m)

    const double *canopy_cvg{ NULL }; // ()
    const double *plant_height{ NULL }; // ()
    const double *ground_cover{ NULL }; // ()
    const double *stones_frac{ NULL }; // ()

    const double *stem_diam{ NULL }; // (m)
    const long   *stem_count{ NULL }; // (1/m^2)
    const double *RFR{ NULL }; // ()

    const double *mmf_d_bare{ NULL }; // ()
    const double *mmf_d_actual{ NULL }; // ()
    const double *mmf_d_tillage{ NULL }; // ()
    const double *mmf_n_bare{ NULL }; // ()
    const double *mmf_n_actual{ NULL }; // ()

    const double *erodibility_clay{ NULL }; // ()
    const double *erodibility_silt{ NULL }; // ()
    const double *erodibility_sand{ NULL }; // ()

    const double *detachibility_clay{ NULL }; // ()
    const double *detachibility_silt{ NULL }; // ()
    const double *detachibility_sand{ NULL }; // ()

    const double *pct_clay{ NULL }; // ()
    const double *pct_silt{ NULL }; // ()

    const double *channel_slope{ NULL }; // ()
    const double *hru_area{ NULL }; // (m)

};

#endif
