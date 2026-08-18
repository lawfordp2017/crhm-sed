// 11/20/19
//---------------------------------------------------------------------------

//#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <assert.h>     /* assert */

#include "ClassSed_overland.h"
#include "../newmodules/NewModules.h"
#include "../../core/GlobalDll.h"


using namespace std;
using namespace CRHM;


// ********* Extra constants:

constexpr double visc = 1.004e-6;  // kinematic viscosity coefficient of water (20C)
constexpr double g = 9.81;  // gravity (m/s2)
constexpr double rho_sed = 2650;  // kg/m3 for quartz
constexpr double rho_w = 1000;  // kg/m3 for water
constexpr double s = rho_sed / rho_w; 

// bed material sd could be calculated from 0.5 (D84/D50 + D16/D50) if the particle distribution profile was reliable enough
constexpr double bed_material_sd = 2.5;   // (sigma_s) From Van Rijn Suspended transport (84), p. 1630-1631
constexpr double von_Karman = 0.4;  // von Karman constant for clear flow
//double mannings_n = 0.104;  // check loch et al. 89 for values

// *******************

ClassSed_Overland* ClassSed_Overland::klone(string name) const{
  return new ClassSed_Overland(name);
}


void ClassSed_Overland::decl(void) {

  variation_set = VARIATION_ORG;

/*****************************************************************
 * VARIABLES
******************************************************************/

// NOTES: The MMF routines are originally written for SI units (I.e. m^2).
// These units are kept internally within this module and converted to km^2 (the default for CRHM) 
// as the final step in runoff_sed_by_erosion() for inter-module transfer.

//  declgetvar("*", "runoff", "(mm)", &runoff);
  declgetvar("*", "scf", "()", &scf);

  declgetvar("*", "soil_runoff", "(mm/int)", &soil_runoff);
  declputvar("*", "soil_runoff_mWQ", "(g/m^2/int)", &soil_runoff_mWQ,&soil_runoff_mWQ_lay);

  declstatvar("sedrelpool", TDim::NHRU, "sediment stored in conceptual storage pool", "(g/m^2)", &sedrelpool);

// This var is determined from silt and clay percentages, so not a parameter
  declstatvar("pct_sand", TDim::NHRU, "Calculated percentage of soil as sand", "(%)", &pct_sand);

  declstatvar("v_std_bare",TDim::NHRU, "surface velocity for the standard bare soil condition", "(m/s)", &v_std_bare);
  declstatvar("v_actual",  TDim::NHRU, "surface velocity for the 'actual' soil condition", "(m/s)", &v_actual);
  declstatvar("v_veg",     TDim::NHRU, "surface velocity for the vegetated soil condition", "(m/s)", &v_veg);
  declstatvar("v_tillage", TDim::NHRU, "surface velocity for the tilled soil condition", "(m/s)", &v_tillage);

  declstatvar("DEP_immed_c", TDim::NHRU, "Percentage of clay deposited immediately", "(%)", &DEP_immed_c);
  declstatvar("DEP_immed_z", TDim::NHRU, "Percentage of silt deposited immediately", "(%)", &DEP_immed_z);
  declstatvar("DEP_immed_s", TDim::NHRU, "Percentage of sand deposited immediately", "(%)", &DEP_immed_s);
  declstatvar("DEP_channel_c", TDim::NHRU, "Percentage of clay deposited in runoff channel", "(%)", &DEP_channel_c);
  declstatvar("DEP_channel_z", TDim::NHRU, "Percentage of silt deposited in runoff channel", "(%)", &DEP_channel_z);
  declstatvar("DEP_channel_s", TDim::NHRU, "Percentage of sand deposited in runoff channel", "(%)", &DEP_channel_s);

  declstatvar("mob_rainsplash",  TDim::NHRU, "Rainsplash mobilized sediment", "(g/m^2)", &mob_rainsplash);
  declstatvar("mob_flow",  TDim::NHRU, "Flow mobilized sediment", "(g/m^2)", &mob_flow);
  declstatvar("sed_delivered_z",  TDim::NHRU, "Flow mobilized sediment", "(g/m^2)", &sed_delivered_z);
  declstatvar("sed_transported_z",  TDim::NHRU, "Flow mobilized sediment", "(g/m^2)", &sed_transported_z);

  declstatvar("conc_soil_rechr", TDim::NDEFN, "Dummy variable to fulfill netroute requirements", "(mg/l)", 
                &conc_soil_rechr, &conc_soil_rechr_lay, numsubstances);
  declstatvar("conc_soil_lower", TDim::NDEFN, "Dummy variable to fulfill netroute requirements", "(mg/l)", 
                &conc_soil_lower, &conc_soil_lower_lay, numsubstances);

/*******************
 * PARAMETERS
 *******************/

// Parameters for the HYPE delay pool
// set this as maximum expected water runoff per timestep (mm)
  declparam("sedrelmax", TDim::NHRU, "[20]", "1","1000", "maximum release fraction from the sediment pool", "(mm)", &sedrelmax);
  declparam("sedrelexp", TDim::NHRU, "[1]", "0.1","10", "sedrelexp", "()", &sedrelexp);

//  declparam("cohesion", TDim::NHRU, "[100]", "0","100", "cohesion", "(kPa)", &cohesion);


// Parameters for the modified MMF formulation
  declparam("canopy_cvg", TDim::NHRU, "[0]", "0","100", "Percentage of canopy coverage", "(%)", &canopy_cvg);
  declparam("plant_height", TDim::NHRU, "[1.0]", "0","3", "Plant height", "(m)", &plant_height);
  declparam("ground_cover", TDim::NHRU, "[0.5]", "0","1", "ground cover fraction", "()", &ground_cover);
  declparam("stones_frac", TDim::NHRU, "[0.5]", "0","1", "stones fraction", "()", &stones_frac);

  declparam("stem_diam", TDim::NHRU, "[0.015]", "0","10", "Stem diameter", "(m)", &stem_diam);
  declparam("stem_count", TDim::NHRU, "[800]", "0","2000", "Stem count", "(1/m^2)", &stem_count);
  declparam("tillage_roughness", TDim::NHRU, "[15]", "5","50", "tillage roughness", "()", &RFR);

  declparam("sfcslope", TDim::NHRU, "[0.0001]", "0","1", "interrill slope angle as rise/run", "(m/m)", &sfcslope);
  declparam("slopelen", TDim::NHRU, "[50]", "0","1000", "slope length", "(m)", &slopelen);

  decldiagparam("mmf_d_bare",    TDim::BASIN, "[0.005]", "0","1", "mmf_d_bare",    "(m)", &mmf_d_bare);
// d_actual = 0.005 for unchanneled flow, 0.01 for shallow rills, and 0.25 for deeper rills
  decldiagparam("mmf_d_actual",  TDim::BASIN, "[0.010]", "0","1", "mmf_d_actual",  "(m)", &mmf_d_actual);
  decldiagparam("mmf_d_tillage", TDim::BASIN, "[0.005]", "0","1", "mmf_d_tillage", "(m)", &mmf_d_tillage);

  decldiagparam("mmf_n_bare",    TDim::BASIN, "[0.015]", "0","1", "mmf_n_bare",    "(s·m^{-1/3})", &mmf_n_bare);     // check loch et al. 89 for values of mannings n
  decldiagparam("mmf_n_actual",  TDim::BASIN, "[0.015]", "0","1", "mmf_n_actual",  "(s·m^{-1/3})", &mmf_n_actual);   // check loch et al. 89 for values of mannings n
// For tillage, n is calculated from tillage implement roughness table
//  decldiagparam("mmf_n_tillage", TDim::BASIN, "[0.3]", "0","1", "mmf_n_tillage", "(m)", &mmf_n_tillage);


// The default/suggested values for these parameters are based on Quansah 1982
// The value of Kc (clay erodibility) should be treated with care since, as shown by Poesen (1985) and Chisci et al. (1989), the detachability of clay particles has a very high variability
  decldiagparam("erodibility_clay", TDim::BASIN, "[0.1]", "0","1", "erodibility of clay", "(g/J)", &erodibility_clay);
  decldiagparam("erodibility_silt", TDim::BASIN, "[0.5]", "0","1", "erodibility of silt", "(g/J)", &erodibility_silt);
  decldiagparam("erodibility_sand", TDim::BASIN, "[0.3]", "0","1", "erodibility of sand", "(g/J)", &erodibility_sand);

  declparam("pct_clay", TDim::NHRU, "[33]", "0","100", "Percentage clay", "(%)", &pct_clay);
  declparam("pct_silt", TDim::NHRU, "[33]", "0","100", "Percentage silt", "(%)", &pct_silt);

// The default/suggested values for these parameters are based on Quansah 1982
  decldiagparam("detachibility_clay", TDim::BASIN, "[1.0]", "0","3", "detachibility of clay", "(g/J)", &detachibility_clay);
  decldiagparam("detachibility_silt", TDim::BASIN, "[1.6]", "0","3", "detachibility of silt", "(g/J)", &detachibility_silt);
  decldiagparam("detachibility_sand", TDim::BASIN, "[1.5]", "0","3", "detachibility of sand", "(g/J)", &detachibility_sand);

//  decldiagparam("mmf_mannings_n", TDim::NHRU, "[0.015]", "0","1", "MMF mannings n", "()", &mmf_mannings_n);
  decldiagparam("mmf_repr_depth", TDim::NHRU, "[0.005]", "0","0.3", "MMF representative depth", "(m)", &mmf_repr_depth);


  declparam("channel_slope", TDim::NHRU, "[0.0001]", "0","10", "slope of the rill channels", "(m/m)", &channel_slope);

// Parameters for the vanRijn formulation
//  declparam("channel_pct", TDim::NHRU, "[0]", "0","100", "fraction of area containing rills", "()", &channel_pct);
  declparam("hru_area", TDim::NHRU, "[1]", "1e-6", "1e+09", "hru area.", "(km^2)", &hru_area);

  declgetvar("*", "net_rain", "(mm/int)", &net_rain);

}


void ClassSed_Overland::init(void) {
  initialize_modMMF();

  for(hh=0; hh<nhru; ++hh) {
    for(int Sub=0; Sub < numsubstances; Sub++) {
      conc_soil_lower_lay[Sub][hh];
    }
  }
}

#define SED_CHANNEL 0

void ClassSed_Overland::run(void) {
    long step = getstep();
    long nstep = step% Global::Freq;

    if(step == 1) { // begining of run
    }

    dayno = julian("now");

    for(hh = 0; chkStruct(); ++hh){ // Using inhibit is dangerous
        runoff_sed_by_erosion();
    }

}


void ClassSed_Overland::finish(bool good) {
}


/*****************************************************************
 * 
 * HYPE Routines
 * 
******************************************************************/


void ClassSed_Overland::runoff_sed_by_erosion() {

    const double erodedSed = calc_erosion();       // total eroded sediment (g/m2)
    assert(erodedSed >= 0);

    if(soil_runoff[hh] > minFlow_WQ) { // eroded sed goes back to soil if no surface runoff
        sedrelpool[hh] += erodedSed;
    }

    // sediment released from delay pool
    double sedReleased;
    
    if ( sedrelmax[hh] > 0.0) {
      sedReleased = std::fmin(sedrelpool[hh],
                              sedrelpool[hh]* pow(soil_runoff[hh]/sedrelmax[hh], sedrelexp[hh])); // (g/m^2/int) export
    } else {
      sedReleased = sedrelpool[hh];
    }
    sedrelpool[hh] -= sedReleased;

    // sediment concentration for sediment released from delay pool
    if (soil_runoff[hh] > 0) {
//        double newSedConc = sedReleased / soil_runoff[hh];  // kg/(1000m^3) = g/m^3
//        assert(newSedConc >= 0);
//        soil_runoff_cWQ_lay[SED_CHANNEL][hh] = newSedConc;

//        double newSedMass = sedReleased*1e6; // g/m2 -> g/km2
        double newSedMass = sedReleased; 
        assert(newSedMass >= 0);
        soil_runoff_mWQ_lay[SED_CHANNEL][hh] = newSedMass;
    } else {
      soil_runoff_mWQ_lay[SED_CHANNEL][hh] = 0.0;
    }

}


/*****************************************************************
 * 
 * Modified [or Daily] Morgan-Morgan-Finney Routines
 * 
******************************************************************/


double ClassSed_Overland::calc_rainsplash_energy(double intensity_per_int) { 

// From Marshall and Palmer, suitable for North-western Europe
//    double Kintensity_A = 8.95;
//    double Kintensity_B = 8.44;

// From Laws and Parsons, suitable for North America east of the Rocky Mountains
// Also used in USLE (wischmeier and Smith, 1978)
  double Kintensity_A = 11.87;
  double Kintensity_B = 8.73;

  double intensity_per_hour = intensity_per_int * Global::Freq / 24;

  return max(0.0, Kintensity_A + Kintensity_B*log10(intensity_per_hour));
}


void ClassSed_Overland::calc_rainsplash_mobilization( sed_triple &rslt) {

  if(net_rain[hh] <= 5.0/Global::Freq) {  // 5.0 mm/day minimum for sediment detachment
    rslt.c = 0;
    rslt.z = 0;
    rslt.s = 0;

    mob_rainsplash[hh] = 0;  // DEBUGGING
    return;
  }

  const double KE_throughfall = calc_rainsplash_energy(net_rain[hh]);
  const double KE_leaf_drainage = max(0.0, (15.8*pow(plant_height[hh],0.5)-5.87) );
  const double cvg_frac = canopy_cvg[hh]/100;
  const double imed = (cvg_frac*KE_leaf_drainage + (1-cvg_frac)*KE_throughfall) * net_rain[hh] * (1.0-scf[hh]);
  assert(imed >= 0);
  rslt.c = erodibility_clay[0] * (pct_clay[hh]/100) * imed;  // (g/J) * (%/100) * (J/m2)
  rslt.z = erodibility_silt[0] * (pct_silt[hh]/100) * imed;  // (g/J) * (%/100) * (J/m2)
  rslt.s = erodibility_sand[0] * (pct_sand[hh]/100) * imed;  // (g/J) * (%/100) * (J/m2)

  mob_rainsplash[hh] = rslt.z;   // DEBUGGING

  return;
}



void ClassSed_Overland::calc_flow_mobilization(double runoff, sed_triple &rslt) {  // runoff: mm/int

// CRHM runoff is mm/int, but MMF expects mm/yr (annualized). 
// This must be accounted since the relationship is non-linear
  const double imed = pow(runoff*Global::Freq*365, 1.5) * (1-(ground_cover[hh]+stones_frac[hh])) * pow(sfcslope[hh], 0.3) * (1.0-scf[hh]) / (Global::Freq*365);
  rslt.c = detachibility_clay[0] * (pct_clay[hh]/100) * imed;
  rslt.z = detachibility_silt[0] * (pct_silt[hh]/100) * imed;
  rslt.s = detachibility_sand[0] * (pct_sand[hh]/100) * imed;

  mob_flow[hh] = rslt.z;   // DEBUGGING

}


// These functions rely on module parameters (not constexpr)
double ClassSed_Overland::calc_flowvel_manning(double depth, double n) {
  if (channel_slope[hh] <= 0.0) {
    CRHMException TExcept("Sed_overland: channel slope must be greater than zero", TExcept::TERMINATE);
    LogError(TExcept);
  }
  return pow(depth, 0.667) * pow(channel_slope[hh], 0.5) / n;
}


double ClassSed_Overland::calc_flowvel_veg(double depth, double n) {

  if (sfcslope[hh] <= 0.0) {
    CRHMException TExcept("Sed_overland: surface slope must be greater than zero", TExcept::TERMINATE);
    LogError(TExcept);
  }

// mMMF version of flow velocity with vegetation cover based on Jin et al. 2000
// stem_diam: (m), stem_count (1/m2)
//  return pow(2*g/(stem_diam[hh]*stem_count[hh]), 0.5) * pow(sfcslope[hh], 0.5);

// DMMF version of flow velocity with vegetation cover based on Petryk and Bosmajian 1975
  double modified_n = pow( (n*n + (stem_diam[hh]*stem_count[hh]*pow(depth, 1.33333))/(2*g)) , 0.5);
  return pow(depth, 0.666667) * pow(sfcslope[hh], 0.5) / modified_n;
}

double ClassSed_Overland::calc_tillage_n() {
  return exp(-2.1132 + 0.0349*RFR[hh]);
}



void ClassSed_Overland::calc_flow_deposition( bool mixed_sed, double flow_vel, sed_triple &rslt) {

// These values are from Modified MMF paper (Morgan and Duzant 2008)
  const double v_s_c = 2e-6;  // settling velocify for clay (m/s)
  const double v_s_z = 2e-3;  // settling velocity for silt (m/s)
  const double v_s_s = 2e-2;  // settling velocity for sand (m/s)

/* Explanation of the scale factor scl;
When applying these equations to deposition from runoff on hillslopes, the settling velocities used
earlier to describe deposition of particles immediately following their detachment by raindrop impact and runoff
need to be modified because the effective settling velocities for mixed particle sizes falling out of runoff in a
depositional environment are often much higher than those in sediment-free water (Zaneveld et al., 1982; Lovell
and Rose, 1991; Rose et al., 2003). They approach the maximum value of the distribution of settling velocities for
the particle size distribution of a given soil (Misra and Rose, 1991). A simple increase in value of the settling
velocities by an order of magnitude is used here to bring them in line with those obtained experimentally for
multi-particle sediments (Lovell and Rose, 1988; Hogarth et al., 2004). 
*/
  const double scl = (mixed_sed) ? 10.0 : 1.0;

  const double imed = slopelen[hh] / (flow_vel*mmf_repr_depth[hh]);

  rslt.c = min(100.0, 44.1 * pow( v_s_c*scl * imed , 0.29 ) );
  rslt.z = min(100.0, 44.1 * pow( v_s_z*scl * imed , 0.29 ) );
  rslt.s = min(100.0, 44.1 * pow( v_s_s*scl * imed , 0.29 ) );
}


void ClassSed_Overland::calc_delivered_to_transport(sed_triple &rslt) {   // g/m2
  static sed_triple detached_rain;
  static sed_triple detached_flow;
//  static sed_triple deposited_pct;

  calc_rainsplash_mobilization(detached_rain);
  calc_flow_mobilization(soil_runoff[hh], detached_flow);
//  calc_flow_deposition(false, deposited_pct);  // not mixed sed because it occurs at the detachment site
  rslt.c = (detached_rain.c + detached_flow.c) * ( 1.0 - (DEP_immed_c[hh]/100) ); // + UPSLOPE!
  rslt.z = (detached_rain.z + detached_flow.z) * ( 1.0 - (DEP_immed_z[hh]/100) ); // + UPSLOPE!
  rslt.s = (detached_rain.s + detached_flow.s) * ( 1.0 - (DEP_immed_s[hh]/100) ); // + UPSLOPE!

  sed_delivered_z[hh] = rslt.z;  // DEBUGGING
}


void ClassSed_Overland::calc_transport_capacity(sed_triple &rslt) {

// MMF expects annualized runoff, so scale as appropriate
  const double Q_annual = soil_runoff[hh]*Global::Freq*365;
// Use this if channel slope is in radians
//  const double imed = (v_actual[hh]*v_veg[hh]*v_tillage[hh]/ pow(v_std_bare[hh],3) ) * Q * Q * sin(channel_slope[hh]);
// Use this if channel slope is in m/m
  const double imed = (v_actual[hh]*v_veg[hh]*v_tillage[hh]/ pow(v_std_bare[hh],3) ) * channel_slope[hh] * Q_annual * Q_annual / (Global::Freq*365);

// Why is transport capacity dependent on percentage of sediment type?
  rslt.c = imed * (pct_clay[hh]/100);
  rslt.z = imed * (pct_silt[hh]/100);
  rslt.s = imed * (pct_sand[hh]/100);

}


void ClassSed_Overland::calc_mmf_sed_balance(sed_triple &rslt) {   // g/m2
  static sed_triple delivered;
  static sed_triple transport_cap;
//  static sed_triple deposited_pct;

  calc_delivered_to_transport(delivered);
  calc_transport_capacity(transport_cap);

  if (delivered.c <= transport_cap.c) {
    rslt.c = delivered.c;
  } else {
    rslt.c = max(transport_cap.c, delivered.c * (1 - DEP_channel_c[hh]/100) );
  }

  if (delivered.z <= transport_cap.z) {
    rslt.z = delivered.z;
  } else {
    rslt.z = max(transport_cap.z, delivered.z * (1 - DEP_channel_z[hh]/100) );
  }

  if (delivered.s <= transport_cap.s) {
    rslt.s = delivered.s;
  } else {
    rslt.s = max(transport_cap.s, delivered.s * (1 - DEP_channel_s[hh]/100) );
  }

  sed_transported_z[hh] = rslt.z; // DEBUGGING
}



double ClassSed_Overland::calc_erosion() {   // g/m2

  sed_triple rslt = {0,0,0};
  calc_mmf_sed_balance(rslt);
  return rslt.c + rslt.z + rslt.s;

}


void ClassSed_Overland::initialize_modMMF() {
  for (long hh = 0; hh < nhru; ++hh) {
    pct_sand[hh] = 100 - (pct_clay[hh] + pct_silt[hh]);
    assert(pct_clay[hh] >= 0);
    assert(pct_silt[hh] >= 0);
    assert(pct_sand[hh] >= 0);

    v_std_bare[hh] = calc_flowvel_manning(mmf_d_bare[0], mmf_n_bare[0]);
    v_actual[hh] = calc_flowvel_manning(mmf_d_actual[0], mmf_n_actual[0]);
    v_veg[hh] = calc_flowvel_veg(mmf_d_bare[0], mmf_n_bare[0]);
    double mmf_n_tillage = calc_tillage_n();
    v_tillage[hh] = calc_flowvel_manning(mmf_d_tillage[0], mmf_n_tillage);

    double v_active;
    if (stem_count[hh] > 0) {
      v_active = v_veg[hh];
    } else {
      v_active = v_std_bare[hh];
    }

    sed_triple dep_immed;
    calc_flow_deposition(false, v_active, dep_immed);
    DEP_immed_c[hh] = dep_immed.c;
    DEP_immed_z[hh] = dep_immed.z;
    DEP_immed_s[hh] = dep_immed.s;

    sed_triple dep_channel;
    calc_flow_deposition(true, v_actual[hh], dep_channel);
    DEP_channel_c[hh] = dep_channel.c;
    DEP_channel_z[hh] = dep_channel.z;
    DEP_channel_s[hh] = dep_channel.s;

  }
}
