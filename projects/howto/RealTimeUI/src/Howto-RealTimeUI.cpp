#include <biogears/cdm/compartment/SECompartmentManager.h>
#include <biogears/cdm/compartment/fluid/SELiquidCompartment.h>
#include <biogears/cdm/engine/PhysiologyEngineTrack.h>
#include <biogears/cdm/patient/SEPatient.h>
#include <biogears/cdm/patient/actions/SEAcuteStress.h>
#include <biogears/cdm/patient/actions/SEColdFaceTest.h>
#include <biogears/cdm/patient/actions/SEExercise.h>
#include <biogears/cdm/patient/actions/SESleep.h>
#include <biogears/cdm/patient/actions/SEConsumeNutrients.h>
#include <biogears/cdm/patient/actions/SESubstanceBolus.h>
#include <biogears/cdm/patient/SENutrition.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <biogears/cdm/substance/SESubstance.h>
#include <biogears/cdm/substance/SESubstanceFraction.h>
#include <biogears/cdm/substance/SESubstanceManager.h>
#include <biogears/cdm/system/environment/actions/SEEnvironmentChange.h>
#include <biogears/cdm/system/physiology/SEBloodChemistrySystem.h>
#include <biogears/cdm/system/physiology/SECardiovascularSystem.h>
#include <biogears/cdm/system/physiology/SEEnergySystem.h>
#include <biogears/cdm/system/physiology/SENervousSystem.h>
#include <biogears/cdm/system/physiology/SERenalSystem.h>
#include <biogears/cdm/system/physiology/SERespiratorySystem.h>
#include <biogears/cdm/system/equipment/ElectroCardioGram/SEElectroCardioGram.h>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/string/manipulation.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <queue>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

// POSIX sockets (macOS / Linux)
#include <arpa/inet.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <netdb.h>
#include <netinet/in.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>

using namespace biogears;
using clk = std::chrono::steady_clock;

namespace {
//  Constantes de configuracao
constexpr double kEngineHz             = 50.0;        // 1 dt = 0.02 s
constexpr int    kCsvDecimateTicks     = 50;          // 50 ticks @ 50Hz = 1 Hz CSV/MQTT
constexpr int    kSnapshotMaxSamples   = 6000;        // rolling buffer ~100 min @ 1 Hz
constexpr int    kHttpPort             = 8000;        // porta da UI

//  ECG (Lead III) — buffer dedicado, NAO decimado.
// O QRS tem duracao tipica de 80-100 ms; para revelar a sua morfologia
// precisamos de f_amostragem >> 1/0.08 Hz ~ 12.5 Hz (Nyquist). Amostrar ao
// mesmo ritmo dos restantes parametros (1 Hz, ver kCsvDecimateTicks) daria
// aliasing total: o "traçado" seria ruido, sem qualquer semelhanca com um
// ECG. Por isso a Lead III é gravada em CADA tick do motor (kEngineHz=50Hz,
// dt=20ms), muito acima do minimo de Nyquist para o complexo QRS.
constexpr int    kEcgBufferSeconds     = 10;
constexpr int    kEcgBufferMaxSamples  = static_cast<int>(kEngineHz) * kEcgBufferSeconds; // 50Hz*8s=400

//  CFT 
constexpr double kCFTGainReduction = 0.90;            // sev_AT_efetiva = sev * 0.90
constexpr double kCFTAfterglow_s   = 120.0;           // janela apos fim do CFT


constexpr double kStressEpiBolus_ug_per_severity = 1.5;   // ug por unidade de sev
constexpr double kStressEpiBolus_conc_ug_per_mL  = 0.1;   // = scenario padrao BG

static const bool        params_calculated = true;
static const char* const kExerciseParamsDir = "/Users/lucasteixeira/Desktop/PhysTwin/core_V5/exercise_params";
constexpr const char* kDefaultPatientXml = "patient/Male_25_Normal.xml";

constexpr const char* kPatientDirCandidates[] = {
  "patient",                                      
  "patients",                                      // ./patients
  "../../runtime/patients",                        // relativo a outputs/Release/bin/
  "../../../runtime/patients",
  "../../../../runtime/patients",
  "/Users/lucasteixeira/Desktop/PhysTwin/core_V5/build_V5/runtime/patients",

};

constexpr const char* kEnginePatientPrefix = "patient/";
constexpr const char* kLogPath      = "HowTo-RealTimeUI.log";
constexpr const char* kCsvPath      = "HowTo-RealTimeUI.csv";
constexpr const char* kDiaryJsonPath = "diary.json";

}  

class CommandQueue {
public:
  void push(std::string cmd) {
    std::lock_guard<std::mutex> lk(mu_);
    q_.push(std::move(cmd));
  }
  bool try_pop(std::string& out) {
    std::lock_guard<std::mutex> lk(mu_);
    if (q_.empty()) return false;
    out = std::move(q_.front());
    q_.pop();
    return true;
  }
private:
  std::mutex mu_;
  std::queue<std::string> q_;
};

struct ActiveActions {
  std::string exercise;
  std::string stress;
  std::string coldface;
  std::string sleep;

  std::string to_csv() const {
    std::vector<std::string> parts;
    if (!exercise.empty()) parts.push_back("exercise=" + exercise);
    if (!stress.empty())   parts.push_back("stress=" + stress);
    if (!coldface.empty()) parts.push_back("coldface=" + coldface);
    if (!sleep.empty())    parts.push_back("sleep=" + sleep);
    if (parts.empty()) return "baseline";
    std::string out;
    for (size_t i = 0; i < parts.size(); ++i) {
      if (i) out += ";";
      out += parts[i];
    }
    return out;
  }
};

//  Amostra timeseries (1 Hz)
struct Sample {
  double   t_s          = 0.0;
  long long wall_time_ms = 0;   // Unix epoch milliseconds 
  std::string actions;
  std::unordered_map<std::string, double> values; 
};

class Snapshot {
public:
  void push(const Sample& s) {
    std::lock_guard<std::mutex> lk(mu_);
    buf_.push_back(s);
    if (static_cast<int>(buf_.size()) > kSnapshotMaxSamples) {
      buf_.pop_front();
    }
  }
  std::vector<Sample> snapshot() const {
    std::lock_guard<std::mutex> lk(mu_);
    return std::vector<Sample>(buf_.begin(), buf_.end());
  }
  void clear() {
    std::lock_guard<std::mutex> lk(mu_);
    buf_.clear();
  }
private:
  mutable std::mutex mu_;
  std::deque<Sample> buf_;
};


//  Amostra de forma de onda do ECG (full engine rate, NAO decimada)
struct EcgSample {
  double t_s = 0.0;   // tempo de simulacao, segundos
  double mV  = 0.0;   // Lead III (unico lead com dado real no BioGears — ver ECG.cpp)
};

// Buffer circular de alta resolucao para a forma de onda do ECG. Ver o
// comentario junto a kEcgBufferMaxSamples para a justificacao (Nyquist) de
// por que este buffer nao pode partilhar a decimacao 1 Hz do Snapshot.
class EcgBuffer {
public:
  void push(const EcgSample& s) {
    std::lock_guard<std::mutex> lk(mu_);
    buf_.push_back(s);
    if (static_cast<int>(buf_.size()) > kEcgBufferMaxSamples) {
      buf_.pop_front();
    }
  }
  std::vector<EcgSample> snapshot() const {
    std::lock_guard<std::mutex> lk(mu_);
    return std::vector<EcgSample>(buf_.begin(), buf_.end());
  }
  void clear() {
    std::lock_guard<std::mutex> lk(mu_);
    buf_.clear();
  }
private:
  mutable std::mutex mu_;
  std::deque<EcgSample> buf_;
};


struct ScheduledStop {
  double      stop_at_s = 0.0;
  std::string what;
};


struct CFTTracker {
  double applied_at_s = -1.0;    
  double duration_s   =  0.0;
  bool   cleared      = false;
  double cleared_at_s = -1.0;

  void on_apply(double t, double dur) {
    applied_at_s = t;
    duration_s   = dur;
    cleared      = false;
    cleared_at_s = -1.0;
  }

  void on_clear(double t) {
    cleared      = true;
    cleared_at_s = t;
  }

  bool was_active_recently(double t_now) const {
    if (applied_at_s < 0.0) return false;
    const double t_window_end = cleared
        ? (cleared_at_s + kCFTAfterglow_s)
        : (applied_at_s + duration_s + kCFTAfterglow_s);
    return (t_now >= applied_at_s) && (t_now < t_window_end);
  }

  bool is_currently_applying(double t_now) const {
    if (applied_at_s < 0.0) return false;
    if (cleared)            return false;
    return (t_now >= applied_at_s) && (t_now < applied_at_s + duration_s);
  }
};

constexpr double kStressDecayTau_s        = 90.0;   
constexpr double kStressDecayUpdateDt_s   = 5.0;    
constexpr double kStressFloor             = 0.01;   

struct StressDecayTracker {
  bool   active           = false;
  double sev_start        = 0.0;     
  double start_t_s        = 0.0;     
  double tau_s            = kStressDecayTau_s;
  double last_update_t_s  = -1e9;    
  double last_applied_sev = 0.0;     

  void arm(double sev0, double t_now, double tau = kStressDecayTau_s) {
    active           = true;
    sev_start        = sev0;
    start_t_s        = t_now;
    tau_s            = tau;
    last_update_t_s  = t_now;
    last_applied_sev = sev0;
  }
  void disarm() {
    active           = false;
    last_applied_sev = 0.0;
  }
  double sev_at(double t_now) const {
    if (!active) return 0.0;
    const double dt = std::max(0.0, t_now - start_t_s);
    return sev_start * std::exp(-dt / tau_s);
  }
  
};

enum class ExerciseType { None, Generic, Running, Cycling };
static ExerciseType g_active_exercise_type = ExerciseType::None;

//  Globais (UI thread + engine thread)
static CommandQueue              g_queue;
static Snapshot                  g_snapshot;
static EcgBuffer                 g_ecg_buffer;   // forma de onda do ECG, alta resolucao
static std::atomic<bool>         g_running { true };
static std::atomic<double>       g_sim_time_s { 0.0 };
static CFTTracker                g_cft;            // engine thread escreve; UI thread le
static std::mutex                g_cft_mu;
static std::vector<ScheduledStop> g_scheduled;
static std::mutex                g_sched_mu;

// Aceleração do engine
static std::atomic<double>       g_speed_factor { 1.0 };


static StressDecayTracker        g_stress_decay;
static std::mutex                g_stress_decay_mu;


static std::atomic<double>       g_current_stress_sev { 0.0 };

//  Parametros selecionados pelo user para MQTT.
//  Defaut: os 6 parametros originais.
static std::vector<std::string> g_mqtt_params = {"hr","map","coreT","skinT","lac","met"};
static std::mutex                g_mqtt_params_mu;


static std::atomic<long long>    g_wall_epoch_origin_ms { 0LL };

//  Troca de paciente em runtime.

static std::atomic<bool>         g_reset_pending { false };
static std::string               g_pending_patient_xml;  
static std::mutex                g_pending_patient_mu;
static std::string               g_current_patient_xml = kDefaultPatientXml;
static std::mutex                g_current_patient_mu;

//  Helpers BioGears 

static void clear_exercise_by_type(PhysiologyEngine& bg, ExerciseType t) {
  switch (t) {
    case ExerciseType::Running: {
      SEExercise::SERunning r;
      r.SpeedRun.SetValue(0.0, LengthPerTimeUnit::m_Per_s);
      r.InclineRun.SetValue(0.0);
      SEExercise ex { r };
      bg.ProcessAction(ex);
      break;
    }
    case ExerciseType::Cycling: {
      SEExercise::SECycling c;
      c.PowerCycle.SetValue(0.0, PowerUnit::W);
      c.CadenceCycle.SetValue(0.0, FrequencyUnit::Per_min);
      SEExercise ex { c };
      bg.ProcessAction(ex);
      break;
    }
    default: {  // Generic ou None
      SEExercise::SEGeneric g;
      g.Intensity.SetValue(0.0);
      SEExercise ex { g };
      bg.ProcessAction(ex);
      break;
    }
  }
}

static void apply_exercise_generic(PhysiologyEngine& bg, double intensity) {

  if (g_active_exercise_type != ExerciseType::None &&
      g_active_exercise_type != ExerciseType::Generic) {
    clear_exercise_by_type(bg, g_active_exercise_type);
  }
  SEExercise::SEGeneric g;
  g.Intensity.SetValue(intensity);
  SEExercise ex {g};
  bg.ProcessAction(ex);
  g_active_exercise_type = (intensity > 0.0) ? ExerciseType::Generic : ExerciseType::None;
}

static void apply_exercise_run(PhysiologyEngine& bg, double speed_m_s, double incline) {

  if (g_active_exercise_type != ExerciseType::None &&
      g_active_exercise_type != ExerciseType::Running) {
    clear_exercise_by_type(bg, g_active_exercise_type);
  }
  SEExercise::SERunning r;
  r.SpeedRun.SetValue(speed_m_s, LengthPerTimeUnit::m_Per_s);
  r.InclineRun.SetValue(incline);
  SEExercise ex { r };
  bg.ProcessAction(ex);
  g_active_exercise_type = (speed_m_s > 0.0) ? ExerciseType::Running : ExerciseType::None;
}

static void apply_exercise_cycle(PhysiologyEngine& bg, double power_W, double cadence_per_min) {
  // Idem para Cycling
  if (g_active_exercise_type != ExerciseType::None &&
      g_active_exercise_type != ExerciseType::Cycling) {
    clear_exercise_by_type(bg, g_active_exercise_type);
  }
  SEExercise::SECycling c;
  c.PowerCycle.SetValue(power_W, PowerUnit::W);
  c.CadenceCycle.SetValue(cadence_per_min, FrequencyUnit::Per_min);
  SEExercise ex { c };
  bg.ProcessAction(ex);
  g_active_exercise_type = (power_W > 0.0) ? ExerciseType::Cycling : ExerciseType::None;
}

static void apply_meal(PhysiologyEngine& bg,
                       double carb_g    = 70.0,
                       double fat_g     = 20.0,
                       double protein_g = 20.0,
                       double sodium_g  =  0.5,
                       double water_L   =  0.35)
{
  SEConsumeNutrients meal;
  SENutrition& n = meal.GetNutrition();
  if (carb_g     > 0.0) n.GetCarbohydrate().SetValue(carb_g,     MassUnit::g);
  if (fat_g      > 0.0) n.GetFat()         .SetValue(fat_g,      MassUnit::g);
  if (protein_g  > 0.0) n.GetProtein()     .SetValue(protein_g,  MassUnit::g);
  if (sodium_g   > 0.0) n.GetSodium()      .SetValue(sodium_g,   MassUnit::g);
  if (water_L    > 0.0) n.GetWater()       .SetValue(water_L,    VolumeUnit::L);
  bg.ProcessAction(meal);
  std::ostringstream ss;
  ss << "[meal] carb=" << carb_g << "g fat=" << fat_g << "g prot=" << protein_g
     << "g Na=" << sodium_g << "g H2O=" << water_L << "L";
  bg.GetLogger()->Info(ss.str());
}

static void apply_sleep(PhysiologyEngine& bg) {
  SESleep sl;
  sl.SetSleepState(SEOnOff::On);
  bg.ProcessAction(sl);
}

static void clear_sleep(PhysiologyEngine& bg) {
  SESleep sl;
  sl.SetSleepState(SEOnOff::Off);
  bg.ProcessAction(sl);
}

static void clear_exercise(PhysiologyEngine& bg) {
  clear_exercise_by_type(bg, g_active_exercise_type);
  g_active_exercise_type = ExerciseType::None;
}

static void apply_stress(PhysiologyEngine& bg, double sev) {
  SEAcuteStress s;
  s.GetSeverity().SetValue(sev);
  bg.ProcessAction(s);
  g_current_stress_sev.store(sev);
}

static void clear_stress(PhysiologyEngine& bg) { apply_stress(bg, 0.0); }

static void apply_coldface(PhysiologyEngine& bg, double sev, double dur_s) {
  SEColdFaceTest c;
  c.GetSeverity().SetValue(sev);
  c.GetDuration().SetValue(dur_s, TimeUnit::s);
  bg.ProcessAction(c);
  std::ostringstream ss;
  ss << "[CFT] apply: severity=" << sev << ", duration=" << dur_s << " s";
  bg.GetLogger()->Info(ss.str());
}

static void clear_coldface(PhysiologyEngine& bg) {
  SEColdFaceTest c;
  c.GetSeverity().SetValue(0.0);
  bg.ProcessAction(c);
  bg.GetLogger()->Info("[CFT] clear: severity -> 0 (action removed)");
}

//  Lactato arterial (Aorta) em mmol/L
static double get_arterial_lactate_mmol_per_L(PhysiologyEngine& bg) {
  auto* aorta = bg.GetCompartments().GetLiquidCompartment(BGE::VascularCompartment::Aorta);
  if (!aorta) return std::nan("");
  auto* lactate = bg.GetSubstanceManager().GetSubstance(StandardSubstances::Lactate);
  if (!lactate) return std::nan("");
  auto* sq = aorta->GetSubstanceQuantity(*lactate);
  if (!sq) return std::nan("");
  return sq->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
}

static double calculate_MET_from_HR(double HR_bpm, int age_years, double weight_kg) {

    double EE = (-55.0969 + 0.6309 * HR_bpm + 0.1988 * weight_kg + 0.2017 * age_years) / 4.184;
    double MET = (EE * 60.0) / weight_kg;
    return std::max(1.0, std::min(20.0, MET));  // Limitando entre 1 e 20 METs
}

static double get_metabolic_equivalent(PhysiologyEngine& bg, double weight_kg) {

    double met = std::nan("");

    if (std::isnan(met) || met < 1.0) {
        const auto& cv = *bg.GetCardiovascularSystem();
        double hr = cv.GetHeartRate(FrequencyUnit::Per_min);
        int age = bg.GetPatient().GetAge(TimeUnit::yr);
        met = calculate_MET_from_HR(hr, age, weight_kg);
    }
    
    return met;
}


//  ECG, Lead III, potencial eletrico instantaneo (mV).
// Nota: esta funcao devolve UMA amostra por chamada — quando usada atraves
// do catalogo kParamCatalog (1 Hz, ver kCsvDecimateTicks) serve para CSV/MQTT
// como valor pontual, mas NAO reconstroi a morfologia do QRS (ver
// kEcgBufferMaxSamples). A forma de onda propriamente dita e servida por
// /ecg, alimentada por g_ecg_buffer a kEngineHz.
static double get_ecg_lead3_mV(PhysiologyEngine& bg) {
  const auto* ecg = bg.GetElectroCardioGram();
  if (ecg == nullptr) return std::nan("");
  if (!ecg->HasLead3ElectricPotential()) return std::nan("");
  return ecg->GetLead3ElectricPotential(ElectricPotentialUnit::mV);
}

//  Epinefrina plasmatica em ug/L
static double get_epinephrine_ug_per_L(PhysiologyEngine& bg) {
  auto* epi = bg.GetSubstanceManager().GetSubstance("Epinephrine");
  if (epi == nullptr) return std::nan("");
  if (!epi->HasPlasmaConcentration()) return std::nan("");
  return epi->GetPlasmaConcentration(MassPerVolumeUnit::ug_Per_L);
}

//  Catalogo de parametros fisiologicos disponíveis no BioGears:

struct ParamDef {
  std::string id;
  std::string label;
  std::string unit;
  std::string module;
  std::function<double(PhysiologyEngine&)> get;
};

// Helper: wraps a getter in try/catch, returning NaN on any exception
template<typename F>
static std::function<double(PhysiologyEngine&)> safe_get(F&& fn) {
  return [fn](PhysiologyEngine& e) -> double {
    try { return fn(e); }
    catch (...) { return std::numeric_limits<double>::quiet_NaN(); }
  };
}

static const std::vector<ParamDef> kParamCatalog = {

  // Cardiovascular 
  {"hr",    "Heart Rate",                   "bpm",        "Cardiovascular",
   safe_get([](PhysiologyEngine& e){ return e.GetCardiovascularSystem()->GetHeartRate(FrequencyUnit::Per_min); })},
  {"map",   "Mean Arterial Pressure",       "mmHg",       "Cardiovascular",
   safe_get([](PhysiologyEngine& e){ return e.GetCardiovascularSystem()->GetMeanArterialPressure(PressureUnit::mmHg); })},
  {"sbp",   "Systolic Pressure",            "mmHg",       "Cardiovascular",
   safe_get([](PhysiologyEngine& e){ return e.GetCardiovascularSystem()->GetSystolicArterialPressure(PressureUnit::mmHg); })},
  {"dbp",   "Diastolic Pressure",           "mmHg",       "Cardiovascular",
   safe_get([](PhysiologyEngine& e){ return e.GetCardiovascularSystem()->GetDiastolicArterialPressure(PressureUnit::mmHg); })},
  {"sv",    "Stroke Volume",                "mL",         "Cardiovascular",
   safe_get([](PhysiologyEngine& e){ return e.GetCardiovascularSystem()->GetHeartStrokeVolume(VolumeUnit::mL); })},
  {"co",    "Cardiac Output",               "L/min",      "Cardiovascular",
   safe_get([](PhysiologyEngine& e){ return e.GetCardiovascularSystem()->GetCardiacOutput(VolumePerTimeUnit::L_Per_min); })},
  {"svr",   "Systemic Vascular Resistance", "mmHg·s/mL",  "Cardiovascular",
   safe_get([](PhysiologyEngine& e){ return e.GetCardiovascularSystem()->GetSystemicVascularResistance(FlowResistanceUnit::mmHg_s_Per_mL); })},
  {"bv",    "Blood Volume",                 "mL",         "Cardiovascular",
   safe_get([](PhysiologyEngine& e){ return e.GetCardiovascularSystem()->GetBloodVolume(VolumeUnit::mL); })},
  {"ecg",   "ECG (Lead III)",               "mV",         "Cardiovascular",
   safe_get([](PhysiologyEngine& e){ return get_ecg_lead3_mV(e); })},


   //  Energy
  {"coreT", "Core Temperature",             "C",          "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetCoreTemperature(TemperatureUnit::C); })},
  {"skinT", "Skin Temperature",             "C",          "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetSkinTemperature(TemperatureUnit::C); })},
  {"tmr",   "Total Metabolic Rate",         "W",          "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetTotalMetabolicRate(PowerUnit::W); })},
  {"eed",   "Exercise Energy Demand",       "W",          "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetExerciseEnergyDemand(PowerUnit::W); })},
  {"ael",   "Achieved Exercise Level",      "",           "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetAchievedExerciseLevel(); })},
  {"fat",   "Fatigue Level",                "",           "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetFatigueLevel(); })},
  {"swr",   "Sweat Rate",                   "mg/s",       "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetSweatRate(MassPerTimeUnit::mg_Per_s); })},
  {"twrl",  "Total Work Rate Level",        "",           "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetTotalWorkRateLevel(); })},
  {"lpr",   "Lactate Production Rate",      "umol/min",   "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetLactateProductionRate(AmountPerTimeUnit::umol_Per_min); })},
  {"cpr",   "Creatinine Production Rate",    "umol/min",   "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetCreatinineProductionRate(AmountPerTimeUnit::umol_Per_min); })},
  {"emd",   "Exercise MAP Delta",           "mmHg",       "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetExerciseMeanArterialPressureDelta(PressureUnit::mmHg); })},
  {"sls",   "Sodium Lost to Sweat",         "mg",         "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetSodiumLostToSweat(MassUnit::mg); })},
  {"kls",   "Potassium Lost to Sweat",      "mg",         "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetPotassiumLostToSweat(MassUnit::mg); })},
  {"cls",   "Chloride Lost to Sweat",       "mg",         "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetChlorideLostToSweat(MassUnit::mg); })},
  {"edef",  "Energy Deficit",               "W",          "Energy",
   safe_get([](PhysiologyEngine& e){ return e.GetEnergySystem()->GetEnergyDeficit(PowerUnit::W); })},
  
  
   //  Respiratory 
  {"rr",    "Respiration Rate",             "bpm",        "Respiratory",
   safe_get([](PhysiologyEngine& e){ return e.GetRespiratorySystem()->GetRespirationRate(FrequencyUnit::Per_min); })},
  {"tv",    "Tidal Volume",                 "L",          "Respiratory",
   safe_get([](PhysiologyEngine& e){ return e.GetRespiratorySystem()->GetTidalVolume(VolumeUnit::L); })},
  {"tlv",   "Total Lung Volume",            "L",          "Respiratory",
   safe_get([](PhysiologyEngine& e){ return e.GetRespiratorySystem()->GetTotalLungVolume(VolumeUnit::L); })},
  {"etco2", "End-Tidal CO2 Fraction",       "",           "Respiratory",
   safe_get([](PhysiologyEngine& e){ return e.GetRespiratorySystem()->GetEndTidalCarbonDioxideFraction(); })},  
  

  // Renal 
  {"upr",   "Urine Production Rate",        "mL/min",     "Renal",
   safe_get([](PhysiologyEngine& e){ return e.GetRenalSystem()->GetUrineProductionRate(VolumePerTimeUnit::mL_Per_min); })},
  {"gfr",   "Glomerular Filtration Rate",   "mL/min",     "Renal",
   safe_get([](PhysiologyEngine& e){ return e.GetRenalSystem()->GetGlomerularFiltrationRate(VolumePerTimeUnit::mL_Per_min); })},
  
  
   //  Blood Chemistry 
  {"ph",    "Arterial Blood pH",            "",           "Blood Chemistry",
   safe_get([](PhysiologyEngine& e){ return e.GetBloodChemistrySystem()->GetArterialBloodPH(); })},
  {"paco2", "Arterial CO2 Pressure",        "mmHg",       "Blood Chemistry",
   safe_get([](PhysiologyEngine& e){ return e.GetBloodChemistrySystem()->GetArterialCarbonDioxidePressure(PressureUnit::mmHg); })},
  {"pao2",  "Arterial O2 Pressure",         "mmHg",       "Blood Chemistry",
   safe_get([](PhysiologyEngine& e){ return e.GetBloodChemistrySystem()->GetArterialOxygenPressure(PressureUnit::mmHg); })},
  {"hb",    "Hemoglobin Content",           "g",          "Blood Chemistry",
   safe_get([](PhysiologyEngine& e){ return e.GetBloodChemistrySystem()->GetHemoglobinContent(MassUnit::g); })},
 
 
   //  Derived 
  {"lac",   "Arterial Lactate",             "mmol/L",     "Derived",
   safe_get([](PhysiologyEngine& e){ return get_arterial_lactate_mmol_per_L(e); })},
  {"epi",   "Plasma Epinephrine",           "ug/L",       "Derived",
   safe_get([](PhysiologyEngine& e){ return get_epinephrine_ug_per_L(e); })},
  {"met",   "Metabolic Equivalent (MET)",   "MET",        "Derived",
   safe_get([](PhysiologyEngine& e){
     double w = e.GetPatient().GetWeight(MassUnit::kg);
     return get_metabolic_equivalent(e, w);
   })},
};

// DIARY — Rotina de 24h programavel que cicla indefinidamente

struct DiaryEntry {
  int         start_min = 0;          // [0, 1439] minutos desde meia-noite
  int         end_min   = 0;          // 0 = sem fim explicito (até ao próximo bloco)
  std::string action    = "baseline"; // "baseline","run","cycle","exercise","stress","cft","sleep"
  double      param1    = 0.0;        
  double      param2    = 0.0;        

  std::string start_str() const {
    char buf[6];
    std::snprintf(buf, sizeof(buf), "%02d:%02d", start_min / 60, start_min % 60);
    return buf;
  }
  std::string end_str() const {
    if (end_min <= 0) return "";
    char buf[6];
    std::snprintf(buf, sizeof(buf), "%02d:%02d", end_min / 60, end_min % 60);
    return buf;
  }
  bool has_end() const { return end_min > start_min; }
};

// Globals do diary
static std::vector<DiaryEntry>  g_diary_blocks;
static std::mutex               g_diary_mu;
static std::atomic<bool>        g_diary_active  { false };
static std::atomic<int>         g_diary_cur_idx { -2 };   // -2=uninit, -1=baseline

// parse "HH:MM" -> minutos [0, 1439]
static int parse_hhmm_to_min(const std::string& s) {
  int h = 0, m = 0;
  if (s.size() >= 5 && s[2] == ':') {
    h = std::atoi(s.substr(0, 2).c_str());
    m = std::atoi(s.substr(3, 2).c_str());
  }
  return std::max(0, std::min(1439, h * 60 + m));
}

//  Helpers JSON minimais (sem dependencias externas)

// Extrai string value de "key": "value"
static std::string mjson_str(const std::string& json, const std::string& key) {
  std::string pat = "\"" + key + "\"";
  auto pos = json.find(pat);
  if (pos == std::string::npos) return "";
  pos = json.find(':', pos + pat.size());
  if (pos == std::string::npos) return "";
  pos = json.find('"', pos + 1);
  if (pos == std::string::npos) return "";
  auto end = json.find('"', pos + 1);
  if (end == std::string::npos) return "";
  return json.substr(pos + 1, end - pos - 1);
}

// Extrai valor numerico de "key": number
static double mjson_dbl(const std::string& json, const std::string& key, double def = 0.0) {
  std::string pat = "\"" + key + "\"";
  auto pos = json.find(pat);
  if (pos == std::string::npos) return def;
  pos = json.find(':', pos + pat.size());
  if (pos == std::string::npos) return def;
  ++pos;
  while (pos < json.size() && (json[pos] == ' ' || json[pos] == '\t')) ++pos;
  if (pos >= json.size()) return def;
  return std::atof(json.c_str() + pos);
}

// Serializa DiaryEntry -> JSON object string
static std::string diary_entry_to_json(const DiaryEntry& e) {
  std::ostringstream o;
  o << std::fixed << std::setprecision(4);
  o << "{\"start\":\"" << e.start_str() << "\""
    << ",\"end\":\"" << e.end_str() << "\""
    << ",\"action\":\"" << e.action << "\""
    << ",\"param1\":" << e.param1
    << ",\"param2\":" << e.param2
    << "}";
  return o.str();
}

// Serializa vector<DiaryEntry> -> JSON array string
static std::string diary_to_json(const std::vector<DiaryEntry>& blocks) {
  std::string out = "[";
  for (std::size_t i = 0; i < blocks.size(); ++i) {
    if (i) out += ",";
    out += diary_entry_to_json(blocks[i]);
  }
  out += "]";
  return out;
}

// Desserializa JSON array/objeto com chave "blocks" -> vector<DiaryEntry>
static std::vector<DiaryEntry> diary_from_json_array(const std::string& json) {
  std::vector<DiaryEntry> out;
  std::string arr = json;
  // verifica se json e um objecto com chave "blocks"
  auto bpos = json.find("\"blocks\"");
  if (bpos != std::string::npos) {
    auto apos = json.find('[', bpos);
    if (apos != std::string::npos) {
      auto epos = json.rfind(']');
      if (epos != std::string::npos && epos > apos)
        arr = json.substr(apos, epos - apos + 1);
    }
  }
  // parseia cada objecto dentro do array
  int depth = 0;
  std::size_t obj_start = std::string::npos;
  for (std::size_t i = 0; i < arr.size(); ++i) {
    if (arr[i] == '{') {
      if (depth == 0) obj_start = i;
      ++depth;
    } else if (arr[i] == '}') {
      --depth;
      if (depth == 0 && obj_start != std::string::npos) {
        const std::string obj = arr.substr(obj_start, i - obj_start + 1);
        DiaryEntry e;
        e.start_min = parse_hhmm_to_min(mjson_str(obj, "start"));
        const std::string end_s = mjson_str(obj, "end");
        e.end_min   = end_s.empty() ? 0 : parse_hhmm_to_min(end_s);
        e.action    = mjson_str(obj, "action");
        if (e.action.empty()) e.action = "baseline";
        e.param1    = mjson_dbl(obj, "param1");
        e.param2    = mjson_dbl(obj, "param2");
        out.push_back(e);
        obj_start = std::string::npos;
      }
    }
  }
  return out;
}

// Ordena blocos por start_min (crescente)
static void diary_sort(std::vector<DiaryEntry>& blocks) {
  std::sort(blocks.begin(), blocks.end(),
    [](const DiaryEntry& a, const DiaryEntry& b){ return a.start_min < b.start_min; });
}

// Dado sim_day_s (segundos dentro do dia actual), retorna o indice do bloco activo (-1 se nenhum)
static int diary_find_idx(double sim_day_s, const std::vector<DiaryEntry>& blocks) {
  if (blocks.empty()) return -1;
  int cur_min = static_cast<int>(sim_day_s) / 60;
  int found = -1;
  for (int i = 0; i < static_cast<int>(blocks.size()); ++i) {
    if (blocks[i].start_min <= cur_min) found = i;
    else break;
  }
  // Se o bloco encontrado tem fim explicito e ja terminou -> baseline
  if (found >= 0 && blocks[found].has_end() && cur_min >= blocks[found].end_min)
    return -1;
  return found;
}

// Formata segundos dentro do dia -> "HH:MM:SS" (para logs do diary)
static std::string diary_hms(double sim_day_s) {
  int total_s = static_cast<int>(sim_day_s);
  int h = total_s / 3600;
  int m = (total_s % 3600) / 60;
  int s = total_s % 60;
  char buf[9];
  std::snprintf(buf, sizeof(buf), "%02d:%02d:%02d", h, m, s);
  return buf;
}

//  diary_save / load 
static void diary_save_json_file(const std::vector<DiaryEntry>& blocks) {
  std::ofstream f(kDiaryJsonPath);
  if (!f) {
    std::cerr << "[diary] nao foi possivel guardar " << kDiaryJsonPath << "\n";
    return;
  }
  f << "{\"blocks\":" << diary_to_json(blocks) << "}\n";
  std::cout << "[diary] guardado " << blocks.size()
            << " blocos em " << kDiaryJsonPath << "\n";
}

static std::vector<DiaryEntry> diary_load_json_file() {
  std::ifstream f(kDiaryJsonPath);
  if (!f) return {};
  std::string content((std::istreambuf_iterator<char>(f)),
                       std::istreambuf_iterator<char>());
  auto blocks = diary_from_json_array(content);
  diary_sort(blocks);
  std::cout << "[diary] carregado " << blocks.size()
            << " blocos de " << kDiaryJsonPath << "\n";
  return blocks;
}

//  Para todas as accoes activas (usado na transicao entre blocos do diary) 
static void diary_stop_all(PhysiologyEngine& bg, ActiveActions& act) {
  clear_exercise(bg);
  act.exercise.clear();
  clear_stress(bg);
  act.stress.clear();
  clear_coldface(bg);
  act.coldface.clear();
  clear_sleep(bg);
  act.sleep.clear();
}

//  Aplica o bloco new_idx do diary (ou baseline se new_idx == -1).
//
//  NOTA sobre refeicoes (meal) e sleep-off safety guard
//  ----------------------------------------------------
//  O BioGears consome substratos metabolicos (glucose, glicogenio, sodio) continuamente.
//  Em simulacoes longas (>2-3 h) sem ConsumeNutrients, os depositos depletam-se e o
//  paciente entra em hipoglicemia / choque antes de qualquer exercicio.
//
//  Por isso:
//    1. A accao "meal" chama apply_meal() com parametros configurados pelo utilizador
//       (param1 = carb_g, param2 = protein_g; fat/Na/H2O com defaults razoaveis).
//    2. Ao sair do sleep (i.e., o bloco ANTERIOR era sleep e o novo nao e),
//       injecta-se automaticamente um "pequeno-almoco de seguranca" (70 g carb,
//       20 g fat, 20 g protein, 0.5 g Na, 0.35 L H2O) para repor os depositos,
//       salvo se o novo bloco ja e "meal" (evita duplicacao).
//
//  Referencia nutritiva: FullDay.xml do BioGears (Scenarios/Miscellaneous/).

static bool g_diary_prev_was_sleep = false;  // track para o sleep-off guard

static void diary_apply_block(int new_idx, const std::vector<DiaryEntry>& blocks,
                               PhysiologyEngine& bg, ActiveActions& act) {

  const bool prev_was_sleep = g_diary_prev_was_sleep;

  diary_stop_all(bg, act);

  if (new_idx < 0 || new_idx >= static_cast<int>(blocks.size())) {
    std::cout << "[diary] -> baseline\n";
    // Sleep-off guard: se saimos do sleep para baseline, injecta pequeno-almoco
    if (prev_was_sleep) {
      std::cout << "[diary] sleep-off -> auto-meal (pequeno-almoco de seguranca)\n";
      apply_meal(bg);
    }
    g_diary_prev_was_sleep = false;
    return;
  }

  const DiaryEntry& e = blocks[new_idx];
  std::cout << "[diary] -> " << e.action
            << " (p1=" << e.param1 << " p2=" << e.param2 << ")\n";

  // Sleep-off guard: injecta refeicao ao acordar, salvo se o proprio bloco ja e "meal"
  if (prev_was_sleep && e.action != "sleep" && e.action != "meal") {
    std::cout << "[diary] sleep-off -> auto-meal (pequeno-almoco de seguranca)\n";
    apply_meal(bg);
  }

  if (e.action == "meal") {
    // param1 = carb (g), param2 = protein (g); fat/Na/H2O com defaults
    const double carb_g    = (e.param1 > 0.0) ? e.param1 : 70.0;
    const double protein_g = (e.param2 > 0.0) ? e.param2 : 20.0;
    apply_meal(bg, carb_g, /*fat_g=*/20.0, protein_g);
    // "meal" nao fica activo — e uma accao one-shot; ActiveActions fica em baseline
  } else if (e.action == "run") {
    const double speed_m_s = e.param1 / 3.6;  // km/h -> m/s
    apply_exercise_run(bg, speed_m_s, e.param2);
    std::ostringstream s;
    s << std::fixed << std::setprecision(2)
      << "run=" << speed_m_s << "m/s,i=" << e.param2;
    act.exercise = s.str();
  } else if (e.action == "cycle") {
    const double cadence = (e.param2 > 0) ? e.param2 : 60.0;
    apply_exercise_cycle(bg, e.param1, cadence);
    std::ostringstream s;
    s << std::fixed << std::setprecision(1)
      << "cyc=" << e.param1 << "W,c=" << cadence;
    act.exercise = s.str();
  } else if (e.action == "exercise") {
    const double intensity = std::clamp(e.param1, 0.0, 1.0);
    apply_exercise_generic(bg, intensity);
    std::ostringstream s;
    s << std::fixed << std::setprecision(2) << "ge=" << intensity;
    act.exercise = (intensity > 0.0) ? s.str() : "";
  } else if (e.action == "stress") {
    const double sev = std::clamp(e.param1, 0.0, 1.0);
    apply_stress(bg, sev);
    std::ostringstream s;
    s << std::fixed << std::setprecision(2) << sev;
    act.stress = s.str();
  } else if (e.action == "cft") {
    const double sev = std::clamp(e.param1, 0.0, 1.0);
    apply_coldface(bg, sev, 90.0);  // duracao sempre 90s conforme design
    std::ostringstream s;
    s << std::fixed << std::setprecision(2) << "sev=" << sev << "/90s";
    act.coldface = s.str();
  } else if (e.action == "sleep") {
    apply_sleep(bg);
    act.sleep = "on";
  }
  // "baseline" e qualquer outro valor: diary_stop_all ja tratou

  g_diary_prev_was_sleep = (e.action == "sleep");
}

// Wall-clock timestamp para string legivel, exmplo:  "2026-05-27 12:36:00.100"

static std::string format_wall_timestamp(long long epoch_ms) {
  long long ms  = epoch_ms % 1000LL;
  time_t    t   = static_cast<time_t>(epoch_ms / 1000LL);
  struct tm tm_info{};
#ifdef _WIN32
  localtime_s(&tm_info, &t);
#else
  localtime_r(&t, &tm_info);
#endif
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d %02d:%02d:%02d.%03lld",
    tm_info.tm_year + 1900, tm_info.tm_mon + 1, tm_info.tm_mday,
    tm_info.tm_hour, tm_info.tm_min, tm_info.tm_sec,
    static_cast<long long>(ms));
  return buf;
}


static std::vector<std::string> list_patient_xmls() {
  std::vector<std::string> out;
  std::string found_dir;
  for (const char* dir : kPatientDirCandidates) {
    DIR* d = ::opendir(dir);
    if (!d) continue;
    struct dirent* e;
    std::vector<std::string> in_this_dir;
    while ((e = ::readdir(d)) != nullptr) {
      std::string n = e->d_name;
      if (n.size() <= 4) continue;
      if (n.substr(n.size() - 4) != ".xml") continue;
      if (n[0] == '.') continue;
      in_this_dir.emplace_back(std::string(kEnginePatientPrefix) + n);
    }
    ::closedir(d);
    if (!in_this_dir.empty()) {
      out      = std::move(in_this_dir);
      found_dir = dir;
      break;
    }
  }
  if (!out.empty()) {

    static bool logged = false;
    if (!logged) {
      std::cout << "[patients] listing dir found: \"" << found_dir
                << "\" (" << out.size() << " XMLs)\n" << std::flush;
      logged = true;
    }
  } else {
    static bool warned = false;
    if (!warned) {
      std::cerr << "[patients] NENHUMA pasta candidata encontrada! Tentativas:\n";
      for (const char* dir : kPatientDirCandidates) {
        std::cerr << "           - " << dir << "\n";
      }
      std::cerr << "           (cwd actual + acrescente caminho em kPatientDirCandidates)\n";
      warned = true;
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

static std::string patient_display_name(const std::string& path) {
  auto slash = path.find_last_of('/');
  std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
  if (base.size() > 4 && base.substr(base.size() - 4) == ".xml") {
    base.resize(base.size() - 4);
  }
  return base;
}


static std::string build_exercise_params_cfg_path(const std::string& patient_xml) {
  const std::string base = patient_display_name(patient_xml);
  std::string out;
  out.reserve(std::strlen(kExerciseParamsDir) + 32 + base.size());
  out  = kExerciseParamsDir;
  out += "/exercise_params_active_";
  out += base;
  out += ".cfg";
  return out;
}


static void apply_exercise_params_cfg_for_patient(const std::string& patient_xml) {
  const std::string cfg = build_exercise_params_cfg_path(patient_xml);
  if (params_calculated) {
    std::ifstream check(cfg);
    if (!check.is_open()) {
      std::cerr << "[AVISO] params_calculated=true mas \"" << cfg
                << "\" nao encontrado.\n"
                << "        Continuando com parametros BioGears originais "
                << "(paciente=" << patient_display_name(patient_xml) << ").\n\n";
    } else {
      std::cout << "[params] using \"" << cfg << "\"\n" << std::flush;
    }
  } else {

    std::remove(cfg.c_str());
    std::cout << "[params] params_calculated=false -> removed \"" << cfg
              << "\" (if it existed)\n" << std::flush;
  }
}

//  MQTT minimal client (TCP)
//  A priori tem que ser localhost, rede INESC bloqueia port 1883
// O broker HiveMQ resolve a questão, é rede pública MAS funciona bem

namespace mqtt_impl {

// Configuracao (alterável via "mqtt start <broker> <port> <topic>") 
static std::string g_broker    = "broker.hivemq.com";
static int         g_port      = 1883;
static std::string g_topic     = "biogears/data/vitals";
static std::string g_client_id = "biogears-realtimeui";

static std::atomic<bool> g_active     { false };
static std::atomic<bool> g_connecting { false };  // connect em curso (background thread)
static std::mutex        g_sock_mu;
static int               g_sock   = -1;    // fd TCP; -1 = desligado
static bool              g_use_ws = false; // true = MQTT sobre WebSocket (porta 8080/8081)
static std::string       g_last_error;

// Helpers de codec

static std::vector<uint8_t> encode_remaining(size_t len) {
  std::vector<uint8_t> out;
  do {
    uint8_t b = static_cast<uint8_t>(len % 128);
    len /= 128;
    if (len > 0) b |= 0x80;
    out.push_back(b);
  } while (len > 0);
  return out;
}

static void append_string(std::vector<uint8_t>& buf, const std::string& s) {
  buf.push_back(static_cast<uint8_t>((s.size() >> 8) & 0xFF));
  buf.push_back(static_cast<uint8_t>( s.size()        & 0xFF));
  for (char c : s) buf.push_back(static_cast<uint8_t>(c));
}

static bool send_all(int fd, const uint8_t* data, size_t len) {
  while (len > 0) {
#ifdef MSG_NOSIGNAL
    ssize_t n = ::send(fd, data, len, MSG_NOSIGNAL);
#else
    ssize_t n = ::send(fd, data, len, 0);
#endif
    if (n <= 0) return false;
    data += n;
    len  -= static_cast<size_t>(n);
  }
  return true;
}

static bool recv_exactly(int fd, uint8_t* buf, size_t n, int timeout_ms = 4000) {
  size_t got = 0;
  while (got < n) {
    fd_set rset; FD_ZERO(&rset); FD_SET(fd, &rset);
    timeval tv; tv.tv_sec = timeout_ms / 1000; tv.tv_usec = (timeout_ms % 1000) * 1000;
    int sel = ::select(fd + 1, &rset, nullptr, nullptr, &tv);
    if (sel <= 0) return false;
    ssize_t r = ::recv(fd, buf + got, n - got, 0);
    if (r <= 0) return false;
    got += static_cast<size_t>(r);
  }
  return true;
}

// Conexao TCP + MQTT 

// Abre conexão TCP com timeout controlado.
// timeout_s: segundos maximos de espera pelo connect (10 s).
// Retorna fd válido em caso de sucesso, ou então -1 em caso de erro.

static int tcp_connect_fd(const std::string& host, int port, int timeout_s = 10) {
  struct addrinfo hints{}, *res = nullptr;
  hints.ai_family   = AF_UNSPEC;
  hints.ai_socktype = SOCK_STREAM;
  const std::string ps = std::to_string(port);
  if (::getaddrinfo(host.c_str(), ps.c_str(), &hints, &res) != 0 || !res)
    return -1;
  int fd = ::socket(res->ai_family, res->ai_socktype, res->ai_protocol);
  if (fd < 0) { ::freeaddrinfo(res); return -1; }
#ifdef SO_NOSIGPIPE
  { int one = 1; ::setsockopt(fd, SOL_SOCKET, SO_NOSIGPIPE, &one, sizeof(one)); }
#endif

  // Coloca o socket em modo nao-bloqueante (Isso pq da outra vez o código parava quando começavamos o MQTT)
  const int saved_flags = ::fcntl(fd, F_GETFL, 0);
  ::fcntl(fd, F_SETFL, saved_flags | O_NONBLOCK);

  int rc = ::connect(fd, res->ai_addr, res->ai_addrlen);
  ::freeaddrinfo(res);

  if (rc == 0) {
    ::fcntl(fd, F_SETFL, saved_flags);
    return fd;
  }
  if (errno != EINPROGRESS) {
    // Erro imediato (host invalido, porta bloqueada (INESC))
    ::close(fd); return -1;
  }

  fd_set wset, eset;
  FD_ZERO(&wset); FD_SET(fd, &wset);
  FD_ZERO(&eset); FD_SET(fd, &eset);
  timeval tv; tv.tv_sec = timeout_s; tv.tv_usec = 0;
  rc = ::select(fd + 1, nullptr, &wset, &eset, &tv);

  if (rc <= 0) {
    ::close(fd); return -1;
  }

  int sock_err = 0; socklen_t errlen = sizeof(sock_err);
  ::getsockopt(fd, SOL_SOCKET, SO_ERROR, &sock_err, &errlen);
  if (sock_err != 0) { ::close(fd); return -1; }

  ::fcntl(fd, F_SETFL, saved_flags);
  return fd;
}

// WebSocket helpers — MQTT over WebSocket (RFC 6455)
// Permite contornar firewalls que bloqueiam a porta 1883.
// Broker: test.mosquitto.org:8080 (ws) ou :8081 (wss, nao suportado aqui)

// Envia um frame WebSocket binario com masking obrigatorio (cliente → servidor)
static bool ws_send_frame(int fd, const uint8_t* data, size_t len) {
  std::vector<uint8_t> frame;
  // Byte 0: FIN=1, RSV=000, Opcode=0010 (binary)
  frame.push_back(0x82);
  // Byte 1..N: MASK=1 + payload length encoding
  if (len <= 125) {
    frame.push_back(0x80 | static_cast<uint8_t>(len));
  } else if (len <= 65535) {
    frame.push_back(0xFE);  // 0x80 | 126
    frame.push_back(static_cast<uint8_t>((len >> 8) & 0xFF));
    frame.push_back(static_cast<uint8_t>( len        & 0xFF));
  } else {
    frame.push_back(0xFF);  // 0x80 | 127
    for (int i = 7; i >= 0; --i)
      frame.push_back(static_cast<uint8_t>((len >> (8 * i)) & 0xFF));
  }
  // Masking key: 4 bytes derivados do timestamp (nao precisa de ser criptograficamente seguro)
  const auto now_ns = static_cast<uint32_t>(
      std::chrono::steady_clock::now().time_since_epoch().count() & 0xFFFFFFFFu);
  const uint8_t mask[4] = {
    static_cast<uint8_t>( now_ns        & 0xFF),
    static_cast<uint8_t>((now_ns >>  8) & 0xFF),
    static_cast<uint8_t>((now_ns >> 16) & 0xFF),
    static_cast<uint8_t>((now_ns >> 24) & 0xFF),
  };

  frame.push_back(mask[0]); frame.push_back(mask[1]);
  frame.push_back(mask[2]); frame.push_back(mask[3]);
  // Payload mascarado (XOR com chave rotativa)
  for (size_t i = 0; i < len; ++i)
    frame.push_back(data[i] ^ mask[i % 4]);
  return send_all(fd, frame.data(), frame.size());
}

// Le e desencapsula um frame WebSocket do servidor (servidor → cliente, sem masking)
// Retorna vector vazio em caso de erro ou timeout.
static std::vector<uint8_t> ws_recv_frame(int fd, int timeout_ms = 5000) {
  uint8_t hdr[2];
  if (!recv_exactly(fd, hdr, 2, timeout_ms)) return {};
  // hdr[0]: FIN + opcode  |  hdr[1]: MASK bit + payload length
  bool     masked = (hdr[1] & 0x80) != 0;
  uint64_t plen   = static_cast<uint64_t>(hdr[1] & 0x7F);
  if (plen == 126) {
    uint8_t ext[2];
    if (!recv_exactly(fd, ext, 2, timeout_ms)) return {};
    plen = (static_cast<uint64_t>(ext[0]) << 8) | ext[1];
  } else if (plen == 127) {
    uint8_t ext[8];
    if (!recv_exactly(fd, ext, 8, timeout_ms)) return {};
    plen = 0;
    for (int i = 0; i < 8; ++i) plen = (plen << 8) | ext[i];
  }
  uint8_t mkey[4] = {};
  if (masked && !recv_exactly(fd, mkey, 4, timeout_ms)) return {};
  std::vector<uint8_t> payload(static_cast<size_t>(plen));
  if (plen > 0) {
    if (!recv_exactly(fd, payload.data(), static_cast<size_t>(plen), timeout_ms)) return {};
    if (masked)
      for (size_t i = 0; i < static_cast<size_t>(plen); ++i) payload[i] ^= mkey[i % 4];
  }
  return payload;
}


static bool ws_handshake(int fd, const std::string& host, int port) {
  std::ostringstream req;
  req << "GET /mqtt HTTP/1.1\r\n"
      << "Host: " << host << ":" << port << "\r\n"
      << "Upgrade: websocket\r\n"
      << "Connection: Upgrade\r\n"
      << "Sec-WebSocket-Key: dGhlIHNhbXBsZSBub25jZQ==\r\n"
      << "Sec-WebSocket-Version: 13\r\n"
      << "Sec-WebSocket-Protocol: mqtt\r\n"
      << "\r\n";
  const std::string r = req.str();
  if (!send_all(fd, reinterpret_cast<const uint8_t*>(r.data()), r.size())) return false;
  std::string resp;
  resp.reserve(512);
  char c;
  while (resp.size() < 8192) {
    fd_set rset; FD_ZERO(&rset); FD_SET(fd, &rset);
    timeval tv; tv.tv_sec = 5; tv.tv_usec = 0;
    if (::select(fd + 1, &rset, nullptr, nullptr, &tv) <= 0) return false;
    if (::recv(fd, &c, 1, 0) != 1) return false;
    resp += c;
    const size_t n = resp.size();
    if (n >= 4 &&
        resp[n-4]=='\r' && resp[n-3]=='\n' &&
        resp[n-2]=='\r' && resp[n-1]=='\n') break;
  }
  return resp.find("101") != std::string::npos;
}

static bool do_mqtt_connect(int fd, const std::string& client_id, bool use_ws = false) {
  std::vector<uint8_t> var_hdr;
  append_string(var_hdr, "MQTT");
  var_hdr.push_back(4);
  var_hdr.push_back(0x02);
  var_hdr.push_back(0x00);
  var_hdr.push_back(0x3C);

  std::vector<uint8_t> payload;
  append_string(payload, client_id);

  std::vector<uint8_t> pkt;
  pkt.push_back(0x10);
  const size_t rem = var_hdr.size() + payload.size();
  auto rl = encode_remaining(rem);
  pkt.insert(pkt.end(), rl.begin(), rl.end());
  pkt.insert(pkt.end(), var_hdr.begin(), var_hdr.end());
  pkt.insert(pkt.end(), payload.begin(), payload.end());

  // Envio: raw TCP ou dentro de frame WebSocket
  if (use_ws) {
    if (!ws_send_frame(fd, pkt.data(), pkt.size())) return false;
    // CONNACK vem encapsulado num frame WebSocket
    const auto frame = ws_recv_frame(fd);
    if (frame.size() < 4) return false;
    return (frame[0] == 0x20 && frame[1] == 0x02 && frame[3] == 0x00);
  } else {
    if (!send_all(fd, pkt.data(), pkt.size())) return false;
    uint8_t connack[4];
    if (!recv_exactly(fd, connack, 4)) return false;
    return (connack[0] == 0x20 && connack[1] == 0x02 && connack[3] == 0x00);
  }
}

static bool do_mqtt_publish(int fd, const std::string& topic,
                             const std::string& payload, bool use_ws = false) {

  std::vector<uint8_t> var_hdr;
  append_string(var_hdr, topic);

  std::vector<uint8_t> pkt;
  pkt.push_back(0x30);
  const size_t rem = var_hdr.size() + payload.size();
  auto rl = encode_remaining(rem);
  pkt.insert(pkt.end(), rl.begin(), rl.end());
  pkt.insert(pkt.end(), var_hdr.begin(), var_hdr.end());
  for (char c : payload) pkt.push_back(static_cast<uint8_t>(c));

  if (use_ws) return ws_send_frame(fd, pkt.data(), pkt.size());
  return send_all(fd, pkt.data(), pkt.size());
}

static void do_mqtt_disconnect(int fd, bool use_ws = false) {
  uint8_t disc[2] = {0xE0, 0x00};
  if (use_ws) ws_send_frame(fd, disc, 2);
  else        send_all(fd, disc, 2);
}

// API pública

// Conecta ao broker e inicia publishing.
// Retorna string vazia em caso de sucesso ou mensagem de erro.

static std::string start(const std::string& broker, int port,
                          const std::string& topic, const std::string& client_id,
                          bool use_ws = false) {
  std::lock_guard<std::mutex> lk(g_sock_mu);
  if (g_active.load()) return "MQTT já ativo. Faça 'mqtt stop' primeiro.";

  int fd = tcp_connect_fd(broker, port);
  if (fd < 0) {
    return "Falha TCP: não conseguiu ligar a " + broker + ":" + std::to_string(port);
  }

  // Se WebSocket: HTTP Upgrade handshake antes do MQTT CONNECT
  if (use_ws) {
    if (!ws_handshake(fd, broker, port)) {
      ::close(fd);
      return "Falha WebSocket: handshake rejeitado por " + broker + ":" + std::to_string(port)
           + " (verifica se a porta suporta WS — tenta 8080)";
    }
    std::cout << "[mqtt] WebSocket upgrade OK\n";
  }

  if (!do_mqtt_connect(fd, client_id, use_ws)) {
    ::close(fd);
    return "Falha MQTT CONNECT/CONNACK (broker recusou ou timeout)";
  }

  g_broker    = broker;
  g_port      = port;
  g_topic     = topic;
  g_client_id = client_id;
  g_use_ws    = use_ws;
  g_sock      = fd;
  g_last_error.clear();
  g_active.store(true);
  return "";  // sucesso
}

// Para o publishing e fecha a conexão.
// Tambem cancela qualquer connect em curso 

static void stop() {
  g_connecting.store(false);
  std::lock_guard<std::mutex> lk(g_sock_mu);
  if (!g_active.load() && g_sock < 0) return;
  if (g_sock >= 0) {
    do_mqtt_disconnect(g_sock, g_use_ws);
    ::close(g_sock);
    g_sock = -1;
  }
  g_active.store(false);
  g_use_ws = false;
  g_last_error.clear();
}


static void publish(const std::string& payload) {
  std::lock_guard<std::mutex> lk(g_sock_mu);
  if (!g_active.load() || g_sock < 0) return;
  if (!do_mqtt_publish(g_sock, g_topic, payload, g_use_ws)) {
    ::close(g_sock);
    g_sock = -1;
    g_active.store(false);
    g_last_error = "publish falhou, broker desligou. Use 'mqtt start' para reconectar.";
    std::cerr << "[MQTT] " << g_last_error << "\n";
  }
}

}  

//  Dispatcher de comandos textuais

//      exercise <intensity> [for <dur_s>]
//      run      <speed_m/s> <incline>       [for <dur_s>]
//      cycle    <power_W>   <cadence_/min>  [for <dur_s>]
//      stress   <severity>                  [for <dur_s>]
//      coldface <severity>  <dur_s>           (dur_s é SEMPRE obrigatório)
//      stop     <exercise|stress|coldface>
//      status
//      help
//      quit

static double try_parse_for_duration(std::istringstream& is) {
  std::string for_kw;
  if (!(is >> for_kw)) return -1.0;
  if (for_kw != "for") return -1.0;
  double dur = -1.0;
  if (!(is >> dur)) return -1.0;
  return (dur > 0.0) ? dur : -1.0;
}

static void schedule_stop(const std::string& what, double dur_s) {
  if (dur_s <= 0.0) return;
  ScheduledStop ss { g_sim_time_s.load() + dur_s, what };
  std::lock_guard<std::mutex> lk(g_sched_mu);
  g_scheduled.push_back(ss);
}

static void cancel_scheduled(const std::string& kind) {
  std::lock_guard<std::mutex> lk(g_sched_mu);
  g_scheduled.erase(
      std::remove_if(g_scheduled.begin(), g_scheduled.end(), [&](const ScheduledStop& s) { return s.what == kind; }),
      g_scheduled.end());
}

//  Substitui o "clear_stress" abrupto por um decay exponencial:
//  sev(t) = sev_now * exp(-(t - t_now) / kStressDecayTau_s)

static void start_stress_decay(double sev_start) {
  if (sev_start <= kStressFloor) return;
  std::lock_guard<std::mutex> lk(g_stress_decay_mu);
  g_stress_decay.arm(sev_start, g_sim_time_s.load(), kStressDecayTau_s);
}

static bool tick_stress_decay(PhysiologyEngine& bg, std::string& act_stress_label) {
  std::lock_guard<std::mutex> lk(g_stress_decay_mu);
  if (!g_stress_decay.active) return false;
  const double t_now = g_sim_time_s.load();
  if (t_now - g_stress_decay.last_update_t_s < kStressDecayUpdateDt_s) return false;

  const double sev = g_stress_decay.sev_at(t_now);
  g_stress_decay.last_update_t_s  = t_now;
  if (sev < kStressFloor) {
    clear_stress(bg);             // sev=0 -> RemoveAcuteStress
    g_stress_decay.disarm();
    act_stress_label.clear();
    return true;
  }

  apply_stress(bg, sev);          // re-aplica com sev mais pequeno
  g_stress_decay.last_applied_sev = sev;
  std::ostringstream s;
  s << std::fixed << std::setprecision(3) << sev << " (decay tau="
    << kStressDecayTau_s << "s)";
  act_stress_label = s.str();
  return true;
}

static void print_help() {
  std::cout << "\n  Comandos:\n"
            << "    exercise <intensity 0-1> [for <dur_s>]\n"
            << "    run <speed_m/s> <incline 0-1> [for <dur_s>]\n"
            << "    cycle <power_W> <cadence_/min> [for <dur_s>]\n"
            << "    stress <severity 0-1> [for <dur_s>]\n"
            << "    coldface <severity 0-1> <dur_s>\n"
            << "    sleep                           (inicia sono)\n"
            << "    stop <exercise|stress|coldface|sleep>\n"
            << "    speed <factor>          (1=tempo real, 2=2x, 5=5x)\n"
            << "    patient <file.xml>      (reinicia o engine com novo paciente)\n"
            << "    mqtt start [broker port topic]  (publica dados via MQTT QoS 0)\n"
            << "    mqtt stop\n"
            << "    status | help | quit\n"
            << "  (UI: http://localhost:" << kHttpPort << ")\n\n";
}

//  Devolve true se o caller deve encerrar (quit).
static bool dispatch(const std::string& line,
                     PhysiologyEngine& bg,
                     ActiveActions& act)
{
  std::istringstream is(line);
  std::string cmd;
  is >> cmd;
  if (cmd.empty()) return false;

  // Quando o diary esta activo bloqueia tudo excepto: diary, mqtt, status, help, quit
  if (g_diary_active.load()) {
    static const std::array<const char*, 9> kManualCmds = {
      "exercise", "run", "cycle", "stress", "coldface", "sleep", "stop", "bolus", "patient"
    };
    for (const auto* mc : kManualCmds) {
      if (cmd == mc) {
        std::cout << "  [diary] bloqueado enquanto diary esta ativo."
                     " Use 'diary stop' para retomar controlo manual.\n";
        return false;
      }
    }
  }

  try {
    if (cmd == "exercise") {
      double i;
      if (!(is >> i)) { std::cout << "  uso: exercise <0..1> [for <s>]\n"; return false; }
      apply_exercise_generic(bg, i);
      std::ostringstream s; s << "ge=" << std::fixed << std::setprecision(2) << i;
      act.exercise = (i > 0.0) ? s.str() : "";
      cancel_scheduled("exercise");
      double dur = try_parse_for_duration(is);
      if (dur > 0.0) schedule_stop("exercise", dur);
      std::cout << "  [ok] exercise intensity=" << i
                << (dur > 0.0 ? " for " + std::to_string(dur) + " s" : "") << "\n";
    }

    else if (cmd == "run") {
      double v, inc;
      if (!(is >> v >> inc)) { std::cout << "  uso: run <m/s> <incline> [for <s>]\n"; return false; }
      apply_exercise_run(bg, v, inc);
      std::ostringstream s; s << "run=" << v << "m/s,i=" << inc;
      act.exercise = (v > 0.0) ? s.str() : "";
      cancel_scheduled("exercise");
      double dur = try_parse_for_duration(is);
      if (dur > 0.0) schedule_stop("exercise", dur);
      std::cout << "  [ok] run speed=" << v << " m/s, incline=" << inc
                << (dur > 0.0 ? " for " + std::to_string(dur) + " s" : "") << "\n";
    }
    else if (cmd == "cycle") {
      double p, c;
      if (!(is >> p >> c)) { std::cout << "  uso: cycle <W> <rpm> [for <s>]\n"; return false; }
      apply_exercise_cycle(bg, p, c);
      std::ostringstream s; s << "cyc=" << p << "W,c=" << c;
      act.exercise = (p > 0.0) ? s.str() : "";
      cancel_scheduled("exercise");
      double dur = try_parse_for_duration(is);
      if (dur > 0.0) schedule_stop("exercise", dur);
      std::cout << "  [ok] cycle power=" << p << " W, cadence=" << c << "/min"
                << (dur > 0.0 ? " for " + std::to_string(dur) + " s" : "") << "\n";
    }

    else if (cmd == "stress") {
      double sev;
      if (!(is >> sev)) { std::cout << "  uso: stress <0..1> [for <s>]\n"; return false; }

      const double t_now = g_sim_time_s.load();
      bool reduce;
      { std::lock_guard<std::mutex> lk(g_cft_mu);
        reduce = g_cft.was_active_recently(t_now); }
      const double sev_eff = reduce ? (sev * kCFTGainReduction) : sev;
      apply_stress(bg, sev_eff);
      { std::lock_guard<std::mutex> lk(g_stress_decay_mu);
        g_stress_decay.disarm(); }
      std::ostringstream s;
      s << std::fixed << std::setprecision(3) << sev_eff;
      if (reduce) {
        std::ostringstream s2;
        s2 << std::fixed << std::setprecision(2);
        s2 << "(<-CFT*" << kCFTGainReduction << " from " << sev << ")";
        s << s2.str();
      }
      act.stress = (sev_eff > 0.0) ? s.str() : "";
      cancel_scheduled("stress");
      double dur = try_parse_for_duration(is);
      if (dur > 0.0) schedule_stop("stress", dur);
      std::cout << "  [ok] stress sev_nominal=" << sev
                << " sev_eff=" << sev_eff
                << (reduce ? "  (CFT afterglow ON, x0.90)" : "")
                << (dur > 0.0 ? " for " + std::to_string(dur) + " s" : "") << "\n";
    }

    else if (cmd == "coldface") {
      double sev, dur;
      if (!(is >> sev >> dur)) { std::cout << "  uso: coldface <sev> <dur_s>\n"; return false; }
      apply_coldface(bg, sev, dur);
      std::ostringstream s; s << "sev=" << sev << "/" << dur << "s";
      act.coldface = (sev > 0.0) ? s.str() : "";
      cancel_scheduled("coldface");
      if (sev > 0.0) {
        { std::lock_guard<std::mutex> lk(g_cft_mu);
          g_cft.on_apply(g_sim_time_s.load(), dur); }
        // Auto-stop no fim da duracao nominal (engine ja para o efeito,
        // mas mantemos o tracker coerente).
        schedule_stop("coldface", dur);
      }
      std::cout << "  [ok] coldface sev=" << sev << " dur=" << dur << " s"
                << "  (afterglow window = " << kCFTAfterglow_s << " s)\n";
    }

    else if (cmd == "sleep") {
      // Sem parametros: liga o estado de sono. "stop sleep" para acordar.
      apply_sleep(bg);
      act.sleep = "on";
      cancel_scheduled("sleep");
      double dur = try_parse_for_duration(is);
      if (dur > 0.0) schedule_stop("sleep", dur);
      std::cout << "  [ok] sleep ON"
                << (dur > 0.0 ? " for " + std::to_string(dur) + " s" : "") << "\n";
    }

    else if (cmd == "stop") {
      std::string what; is >> what;
      if      (what == "exercise") { clear_exercise(bg); act.exercise.clear(); cancel_scheduled("exercise"); }
      else if (what == "sleep")    { clear_sleep(bg); act.sleep.clear(); cancel_scheduled("sleep");
                                     std::cout << "  [ok] sleep OFF (paciente acordado)\n"; return false; }
      else if (what == "stress")   {
        // Em vez de remover abruptamente, armamos o decay exponencial.
        // O ciclo de reapplicacao acontece no loop principal (tick_stress_decay).
        const double sev_now = g_current_stress_sev.load();
        cancel_scheduled("stress");
        if (sev_now > kStressFloor) {
          start_stress_decay(sev_now);
          std::ostringstream s;
          s << std::fixed << std::setprecision(3) << sev_now << " (decay tau="
            << kStressDecayTau_s << "s)";
          act.stress = s.str();
          std::cout << "  [ok] stress decay armed (sev_start=" << sev_now
                    << ", tau=" << kStressDecayTau_s << " s)\n";
        } else {
          // Sem stress activo, basta confirmar o estado limpo.
          clear_stress(bg);
          act.stress.clear();
          std::cout << "  [ok] stress parado (sem decay, ja estava em ~0)\n";
        }
      }
      else if (what == "coldface") { clear_coldface(bg); act.coldface.clear(); cancel_scheduled("coldface");
                                     { std::lock_guard<std::mutex> lk(g_cft_mu);
                                       g_cft.on_clear(g_sim_time_s.load()); } }
      else { std::cout << "  uso: stop <exercise|stress|coldface|sleep>\n"; return false; }
      if (what != "stress" && what != "sleep") std::cout << "  [ok] " << what << " parado\n";
    }

    else if (cmd == "speed") {
      // speed <factor>
      // Altera a velocidade do loop de wall-clock 
      double s;
      if (!(is >> s)) { std::cout << "  uso: speed <factor>  (ex: 1, 2, 5)\n"; return false; }
      if (s <= 0.0)   { std::cout << "  [erro] speed deve ser > 0\n"; return false; }
      if (s > 100.0)  { s = 100.0; std::cout << "  [aviso] speed limitado a 100x (engine atinge teto de wall-clock)\n"; }
      g_speed_factor.store(s);
      std::cout << "  [ok] speed factor = " << s << "x (dt_engine fixo em "
                << 1.0 / kEngineHz << " s)\n";
    }

    else if (cmd == "patient") {
      // Troca de paciente: reinicia o engine no loop principal.

      std::string xml;
      // Aceita até o resto da linha.
      std::getline(is, xml);
      
      while (!xml.empty() && (xml.front() == ' ' || xml.front() == '\t')) xml.erase(xml.begin());
      while (!xml.empty() && (xml.back()  == ' ' || xml.back()  == '\t' || xml.back() == '\r' || xml.back() == '\n')) xml.pop_back();
      if (xml.empty()) {
        std::cout << "  uso: patient <file.xml>  (ex: patient/Male_25_Normal.xml)\n";
        return false;
      }

      { std::lock_guard<std::mutex> lk(g_pending_patient_mu);
        g_pending_patient_xml = xml; }
      g_reset_pending.store(true);
      std::cout << "  [ok] patient switch requested -> " << xml
                << "  (engine sera reiniciado no proximo tick)\n";
    }

    else if (cmd == "mqtt") {
      std::string sub;
      if (!(is >> sub)) {
        std::cout << "  uso: mqtt start [broker port topic] | mqtt stop\n";
        return false;
      }
      
      if (sub == "start") {
        if (mqtt_impl::g_active.load()) {
          std::cout << "  [mqtt] ja ativo. Faz 'mqtt stop' primeiro.\n";
          return false;
        }
        if (mqtt_impl::g_connecting.load()) {
          std::cout << "  [mqtt] connect em curso, aguarda...\n";
          return false;
        }
        
        std::string broker = mqtt_impl::g_broker;
        int port           = mqtt_impl::g_port;
        std::string topic  = mqtt_impl::g_topic;
        bool use_ws        = false;
        std::string b;
        if (is >> b) {
          broker = b;
          int p = 0;
          if (is >> p) {
            port = p;
            std::string t;
            if (is >> t) {
              // t pode ser o topic ou a keyword "ws"
              if (t == "ws") {
                use_ws = true;
                if (port == 1883) port = 8080; // default WS port se nao alterado
              } else {
                topic = t;
                // ainda pode haver "ws" no fim
                std::string ws_kw;
                if (is >> ws_kw && ws_kw == "ws") {
                  use_ws = true;
                  if (port == 1883) port = 8080;
                }
              }
            }
          }
        }

        // IMPORTANTE: O connect TCP pode bloquear durante vários segundos, após alterações
        // O engine continua a correr normalmente enquanto isso acontece.

        mqtt_impl::g_connecting.store(true);
        { std::lock_guard<std::mutex> lk(mqtt_impl::g_sock_mu);
          mqtt_impl::g_last_error.clear(); }
        std::cout << "  [mqtt] connect em background -> " << broker << ":"
                  << port << "  topico=" << topic
                  << (use_ws ? "  [WebSocket]" : "") << "\n"
                  << "         (a simulacao continua a correr; estado disponivel na UI)\n"
                  << std::flush;
        std::thread([broker, port, topic, use_ws]() {
          // O engine thread nunca é bloqueado.
          if (!mqtt_impl::g_connecting.load()) return;
          const auto err = mqtt_impl::start(broker, port, topic, mqtt_impl::g_client_id, use_ws);
          mqtt_impl::g_connecting.store(false);
          if (err.empty()) {
            std::cout << "\n  [mqtt] ligado! Publishing ativo -> "
                      << broker << ":" << port << "  topico=" << topic
                      << (use_ws ? "  [WebSocket]" : "") << "\n> " << std::flush;
          } else {
            std::cout << "\n  [mqtt] falha: " << err << "\n> " << std::flush;
          }
        }).detach();
      } else if (sub == "stop") {
        mqtt_impl::stop();
        std::cout << "  [ok] MQTT parado.\n";
      } else {
        std::cout << "  uso: mqtt start [broker port topic] | mqtt stop\n";
      }
    }
    else if (cmd == "diary") {
      std::string sub;
      if (!(is >> sub)) {
        std::cout << "  uso: diary add HH:MM action param1 [param2]\n"
                  << "       diary clear\n"
                  << "       diary list\n"
                  << "       diary start\n"
                  << "       diary stop\n"
                  << "       diary save\n"
                  << "       diary load\n";
        return false;
      }
      if (sub == "add") {
        // Aceita "HH:MM" ou "HH:MM-HH:MM" como primeiro argumento
        std::string hhmm_range, action;
        double p1 = 0.0, p2 = 0.0;
        if (!(is >> hhmm_range >> action >> p1)) {
          std::cout << "  uso: diary add HH:MM[-HH:MM] action param1 [param2]\n"
                    << "       ex: diary add 08:00-08:30 run 10 0\n"
                    << "       acoes: run (km/h, incline)  cycle (W, rpm)\n"
                    << "              exercise (0..1)  stress (0..1)  cft (0..1)\n"
                    << "              baseline (sem params)\n";
          return false;
        }
        is >> p2; // param2 opcional
        // Parseia range de hora (start-end ou so start)
        std::string hhmm_start = hhmm_range, hhmm_end = "";
        auto dash = hhmm_range.find('-', 1);  // ignora possivel '-' inicial
        if (dash != std::string::npos) {
          hhmm_start = hhmm_range.substr(0, dash);
          hhmm_end   = hhmm_range.substr(dash + 1);
        }
        DiaryEntry e;
        e.start_min = parse_hhmm_to_min(hhmm_start);
        e.end_min   = hhmm_end.empty() ? 0 : parse_hhmm_to_min(hhmm_end);
        e.action    = action;
        e.param1    = p1;
        e.param2    = p2;
        { std::lock_guard<std::mutex> lk(g_diary_mu);
          g_diary_blocks.push_back(e);
          diary_sort(g_diary_blocks); }
        std::cout << "  [diary] bloco adicionado: " << e.start_str()
                  << (e.has_end() ? ("-" + e.end_str()) : "")
                  << "  " << e.action
                  << "  p1=" << e.param1 << " p2=" << e.param2 << "\n";
      } else if (sub == "clear") {
        { std::lock_guard<std::mutex> lk(g_diary_mu);
          g_diary_blocks.clear(); }
        g_diary_cur_idx.store(-2);
        std::cout << "  [diary] todos os blocos removidos.\n";
      } else if (sub == "list") {
        std::lock_guard<std::mutex> lk(g_diary_mu);
        if (g_diary_blocks.empty()) {
          std::cout << "  [diary] sem blocos definidos.\n";
        } else {
          std::cout << "  [diary] " << g_diary_blocks.size() << " blocos:\n";
          for (int i = 0; i < static_cast<int>(g_diary_blocks.size()); ++i) {
            const auto& e = g_diary_blocks[i];
            std::cout << "    [" << i << "]  " << e.start_str()
                      << (e.has_end() ? ("-" + e.end_str()) : "     ")
                      << "  " << e.action
                      << "  p1=" << e.param1
                      << "  p2=" << e.param2 << "\n";
          }
        }
      } else if (sub == "start") {
        if (g_diary_active.load()) {
          std::cout << "  [diary] ja ativo.\n"; return false;
        }
        { std::lock_guard<std::mutex> lk(g_diary_mu);
          if (g_diary_blocks.empty()) {
            std::cout << "  [diary] sem blocos. Usa 'diary add' primeiro.\n";
            return false;
          }
          diary_sort(g_diary_blocks); }
        // Activa diary ANTES do reset para que o 4.1b preserve g_diary_active=true
        g_diary_cur_idx.store(-2);
        g_diary_active.store(true);
        // Dispara reboot do engine com o paciente actual (como troca de paciente).
        // O diary arrancara a partir de sim_time=0 (i.e., 00:00 do horario do diary).
        {
          std::string cur_xml;
          { std::lock_guard<std::mutex> lk2(g_current_patient_mu);
            cur_xml = g_current_patient_xml; }
          { std::lock_guard<std::mutex> lk3(g_pending_patient_mu);
            g_pending_patient_xml = cur_xml; }
          g_reset_pending.store(true);
        }
        std::cout << "  [diary] INICIADO. Engine a reiniciar (sim_time -> 0)...\n"
                  << "          Controlo manual bloqueado (exceto MQTT).\n"
                  << "          Usa 'diary stop' para retomar controlo.\n";
      } else if (sub == "stop") {
        if (!g_diary_active.load()) {
          std::cout << "  [diary] nao esta ativo.\n"; return false;
        }
        g_diary_active.store(false);
        g_diary_cur_idx.store(-2);
        diary_stop_all(bg, act);
        std::cout << "  [diary] PARADO. Controlo manual disponivel.\n";
      } else if (sub == "save") {
        std::lock_guard<std::mutex> lk(g_diary_mu);
        diary_save_json_file(g_diary_blocks);
      } else if (sub == "load") {
        if (g_diary_active.load()) {
          std::cout << "  [diary] para o diary antes de carregar ('diary stop').\n";
          return false;
        }
        auto blocks = diary_load_json_file();
        { std::lock_guard<std::mutex> lk(g_diary_mu);
          g_diary_blocks = std::move(blocks); }
        g_diary_cur_idx.store(-2);
      } else {
        std::cout << "  subcomando desconhecido: '" << sub << "'\n";
      }
    }
    else if (cmd == "status") {
      const auto& cv = *bg.GetCardiovascularSystem();
      const auto& en = *bg.GetEnergySystem();
      std::cout << "  t=" << g_sim_time_s.load() << " s  "
                << "HR="   << cv.GetHeartRate(FrequencyUnit::Per_min)        << " bpm  "
                << "MAP="  << cv.GetMeanArterialPressure(PressureUnit::mmHg) << " mmHg  "
                << "CoreT="<< en.GetCoreTemperature(TemperatureUnit::C)      << " C  "
                << "SkinT="<< en.GetSkinTemperature(TemperatureUnit::C)      << " C  "
                << "Lac="  << get_arterial_lactate_mmol_per_L(bg)            << " mmol/L  "
                << "Epi="  << get_epinephrine_ug_per_L(bg)                   << " ug/L\n"
                << "  MQTT params: [";
      { std::lock_guard<std::mutex> lk(g_mqtt_params_mu);
        for (size_t i = 0; i < g_mqtt_params.size(); ++i) {
          if (i) std::cout << ",";
          std::cout << g_mqtt_params[i];
        }
      }
      std::cout << "]\n";
    }
    else if (cmd == "help") {
      print_help();
    }
    else if (cmd == "quit" || cmd == "exit") {
      return true;
    }
    else {
      std::cout << "  comando desconhecido: '" << cmd << "' (digite 'help')\n";
    }
  } catch (const std::exception& e) {
    std::cout << "  [erro] " << e.what() << "\n";
  }
  return false;
}

//  HTML embebido (single-page UI)
//  Vanilla JS + Chart.js 

constexpr const char* kHtmlPage = R"HTML(<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>BioGears RealTimeUI</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/luxon@3/build/global/luxon.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/chartjs-adapter-luxon@1.3.1/dist/chartjs-adapter-luxon.umd.min.js"></script>
<style>
  :root {
    --bg:#101418; --fg:#e6edf3; --muted:#8b949e; --accent:#58a6ff;
    --panel:#161b22; --border:#30363d; --btn:#21262d; --btn-h:#30363d;
    --ok:#3fb950; --warn:#d29922; --err:#f85149;
  }
  * { box-sizing: border-box; }
  html,body { margin:0; padding:0; height:100%; background:var(--bg); color:var(--fg);
              font-family: -apple-system, BlinkMacSystemFont, "SF Pro Text", Arial, sans-serif; font-size:13px; }

  /* ── Parameter Selection Overlay ──────────────────────────────────────── */
  .sel-overlay {
    display:none; position:fixed; inset:0; background:var(--bg);
    z-index:1000; overflow-y:auto; padding:32px 24px;
  }
  .sel-overlay.open { display:block; }
  .sel-header { text-align:center; margin-bottom:28px; }
  .sel-header h1 { font-size:22px; color:var(--accent); margin:0 0 6px; }
  .sel-header p  { color:var(--muted); font-size:13px; margin:0; }
  .sel-section { max-width:960px; margin:0 auto 28px; }
  .sel-section h2 { font-size:14px; font-weight:700; color:var(--fg); margin:0 0 6px; display:flex; align-items:center; gap:10px; }
  .sel-section h2 .badge {
    background:rgba(88,166,255,0.15); color:var(--accent);
    border:1px solid rgba(88,166,255,0.35); border-radius:20px;
    padding:2px 10px; font-size:11px; font-weight:600;
  }
  .sel-hint { font-size:11px; color:var(--muted); margin:0 0 12px; }
  .sel-module-group { margin-bottom:16px; }
  .sel-module-label { font-size:10px; text-transform:uppercase; letter-spacing:.08em;
                      color:var(--muted); margin-bottom:6px; }
  .sel-chips { display:flex; flex-wrap:wrap; gap:8px; }
  .sel-chip {
    display:flex; flex-direction:column; align-items:flex-start;
    background:var(--panel); border:1.5px solid var(--border); border-radius:8px;
    padding:8px 12px; cursor:pointer; min-width:140px; max-width:200px;
    transition: border-color .15s, background .15s;
  }
  .sel-chip:hover  { border-color:#58a6ff88; background:#1a2030; }
  .sel-chip.active { border-color:var(--accent); background:rgba(88,166,255,0.10); }
  .sel-chip.active .chip-name { color:var(--accent); }
  .sel-chip.disabled-chip { opacity:.35; cursor:not-allowed; pointer-events:none; }
  .chip-name { font-size:12px; font-weight:600; color:var(--fg); }
  .chip-unit { font-size:10px; color:var(--muted); margin-top:2px; }
  .sel-footer { max-width:960px; margin:0 auto; display:flex; justify-content:flex-end; align-items:center; gap:16px; padding-top:8px; }
  .sel-footer .hint { font-size:11px; color:var(--muted); }
  #sel-confirm-btn {
    background:var(--accent); color:#0d1117; border:none;
    padding:10px 28px; border-radius:6px; font-size:14px; font-weight:700;
    cursor:pointer;
  }
  #sel-confirm-btn:disabled { opacity:.35; cursor:not-allowed; }
  #sel-confirm-btn:not(:disabled):hover { filter:brightness(1.1); }

  /* ── Main App ─────────────────────────────────────────────────────────── */
  .app { display:grid; grid-template-columns: 1fr 360px; grid-template-rows: auto 1fr;
         height:100vh; gap:8px; padding:8px; }
  header { grid-column: 1 / 3; background:var(--panel); border:1px solid var(--border);
           border-radius:6px; padding:8px 12px; display:flex; gap:18px; align-items:center; flex-wrap:wrap; }
  header h1 { margin:0; font-size:14px; font-weight:600; color:var(--accent); }
  .vital { display:inline-flex; gap:4px; align-items:baseline; }
  .vital .k { color:var(--muted); }
  .vital .v { font-variant-numeric: tabular-nums; min-width:56px; display:inline-block; }
  .charts { background:var(--panel); border:1px solid var(--border); border-radius:6px;
            padding:8px; display:grid; grid-template-columns: 1fr 1fr; grid-template-rows: 1fr 1fr;
            gap:6px; min-height:0; }
  .chart-cell { background:#0d1117; border:1px solid var(--border); border-radius:4px; padding:6px;
                position:relative; min-height:0; min-width:0; }
  .chart-title { font-size:11px; color:var(--muted); position:absolute; top:6px; left:8px; z-index:2; }
  .chart-cell canvas { width:100% !important; height:100% !important; }
  .controls { background:var(--panel); border:1px solid var(--border); border-radius:6px;
              padding:10px; overflow-y:auto; display:flex; flex-direction:column; gap:10px; }
  .ctl-group { background:#0d1117; border:1px solid var(--border); border-radius:4px; padding:8px; }
  .ctl-group h3 { margin:0 0 6px 0; font-size:12px; font-weight:600; color:var(--accent); }
  .ctl-group .desc { font-size:11px; color:var(--muted); margin-bottom:6px; }
  .row { display:flex; gap:6px; align-items:center; margin-bottom:4px; flex-wrap:wrap; }
  .row label { min-width:64px; color:var(--muted); font-size:11px; }
  .row input[type=number] { width:60px; background:#010409; border:1px solid var(--border);
                             color:var(--fg); padding:3px 5px; border-radius:3px;
                             font-variant-numeric: tabular-nums; }
  button { background:var(--btn); color:var(--fg); border:1px solid var(--border);
           padding:5px 10px; border-radius:3px; font-size:11px; cursor:pointer; font-weight:600; }
  button:hover { background:var(--btn-h); }
  button.stop { color:var(--err); }
  button.go   { color:var(--ok); }
  button.warn { color:var(--warn); }
  #status { font-size:11px; color:var(--muted); padding:4px 0; word-break:break-all; }
  .pill { padding:4px 8px; border-radius:3px; font-size:11px; font-weight:600; }
  #cft-state { padding:4px 8px; border-radius:3px; font-size:11px; font-weight:600; }
  #cft-state.on  { background:rgba(63,185,80,0.15); color:var(--ok); }
  #cft-state.off { background:rgba(139,148,158,0.15); color:var(--muted); }
  #cft-applying.on  { background:rgba(248,81,73,0.18); color:var(--err); }
  #cft-applying.off { background:rgba(139,148,158,0.15); color:var(--muted); }
  #stress-decay.on  { background:rgba(210,153,34,0.18); color:var(--warn); }
  #stress-decay.off { background:rgba(139,148,158,0.15); color:var(--muted); }
  /* Botoes de speed: o ativo recebe .active (borda accent + fundo destacado) */
  .speed-btn { min-width:42px; }
  .speed-btn.active { background:var(--accent); color:#0d1117; border-color:var(--accent); }
  
  /* MET gauge style */
  .met-value { font-size:14px; font-weight:bold; color:var(--ok); }
  .met-category { font-size:11px; color:var(--muted); margin-left:8px; }

  /* Patient selector (header) */
  .patient-picker { position:relative; display:inline-block; }
  .patient-picker button.patient-btn {
    background:rgba(88,166,255,0.12); color:var(--accent);
    border:1px solid var(--border); padding:5px 10px; border-radius:4px;
    font-size:13px; font-weight:600; cursor:pointer; }
  .patient-picker button.patient-btn:hover { background:rgba(88,166,255,0.22); }
  .patient-picker button.patient-btn .chev { color:var(--muted); margin-left:6px; font-size:10px; }
  .patient-menu {
    display:none; position:absolute; top:calc(100% + 4px); left:0; z-index:50;
    background:var(--panel); border:1px solid var(--border); border-radius:4px;
    min-width:240px; max-height:340px; overflow-y:auto; padding:4px 0;
    box-shadow:0 6px 18px rgba(0,0,0,0.45);
  }
  .patient-menu.open { display:block; }
  .patient-menu .item {
    padding:6px 12px; color:var(--fg); font-size:12px; cursor:pointer;
    white-space:nowrap;
  }
  .patient-menu .item:hover { background:var(--btn-h); }
  .patient-menu .item.current { color:var(--accent); font-weight:600; }
  .patient-menu .item.current::before { content:"\2713  "; }   /* check mark */

  /* Diary modal */
  .diary-overlay {
    display:none; position:fixed; inset:0; background:rgba(0,0,0,0.72);
    z-index:200; align-items:flex-start; justify-content:center; padding-top:30px;
  }
  .diary-overlay.open { display:flex; }
  .diary-modal {
    background:var(--panel); border:1px solid var(--border); border-radius:8px;
    width:720px; max-width:96vw; max-height:88vh; overflow-y:auto;
    padding:18px; display:flex; flex-direction:column; gap:14px; box-shadow:0 8px 32px rgba(0,0,0,0.6);
  }
  .diary-modal h2 { margin:0; font-size:14px; color:var(--accent); }
  .diary-timeline {
    position:relative; height:38px; background:#010409;
    border:1px solid var(--border); border-radius:4px; overflow:hidden;
  }
  .diary-block-bar {
    position:absolute; top:4px; bottom:4px; border-radius:3px;
    font-size:9px; color:#fff; display:flex; align-items:center;
    justify-content:center; overflow:hidden; white-space:nowrap; min-width:18px;
    opacity:0.92; cursor:default;
  }
  .diary-table { width:100%; border-collapse:collapse; font-size:11px; }
  .diary-table th { background:#0d1117; color:var(--muted); padding:5px 6px;
                    text-align:left; border-bottom:1px solid var(--border); }
  .diary-table td { padding:5px 6px; border-bottom:1px solid #21262d; }
  .diary-table tr:hover td { background:rgba(255,255,255,0.03); }
  .diary-add-form { display:flex; gap:8px; align-items:center; flex-wrap:wrap;
                    background:#0d1117; border:1px solid var(--border);
                    border-radius:4px; padding:8px; }
  .diary-add-form label { font-size:11px; color:var(--muted); white-space:nowrap; }
  .diary-add-form input, .diary-add-form select {
    background:#010409; border:1px solid var(--border); color:var(--fg);
    padding:3px 6px; border-radius:3px; font-size:11px; }
  .diary-status-pill { padding:4px 10px; border-radius:4px; font-size:11px; font-weight:600; }
  .diary-status-pill.active   { background:rgba(63,185,80,0.15);  color:var(--ok);   }
  .diary-status-pill.inactive { background:rgba(139,148,158,0.12); color:var(--muted); }
  #diary-side-status.active   { color:var(--ok); }
  #diary-side-status.inactive { color:var(--muted); }

  /* Lock visual quando diary esta activo */
  .diary-locked { opacity:0.32; pointer-events:none; user-select:none; }
  #patient-btn.diary-locked  { opacity:0.32; pointer-events:none; }
  .speed-btn.diary-locked    { opacity:0.32; pointer-events:none; cursor:not-allowed; }
  .diary-active-banner {
    display:none; background:rgba(63,185,80,0.12); border:1px solid rgba(63,185,80,0.35);
    border-radius:4px; padding:5px 10px; font-size:11px; color:var(--ok);
    text-align:center; font-weight:600;
  }
  .diary-active-banner.show { display:block; }
</style>
</head>
<body>
<!-- ===== Parameter Selection Overlay ===== -->
<div class="sel-overlay" id="sel-overlay">
  <div class="sel-header">
    <h1>PhysTwin RealTimeUI &mdash; Configure Monitoring</h1>
    <p>Choose what to display and publish before starting the simulation</p>
  </div>

  <div class="sel-section">
    <h2>&#128202; Graph Parameters <span class="badge" id="chart-badge">0 / 4 selected</span></h2>
    <p class="sel-hint">Select <strong>exactly 4</strong> parameters to display in the real-time charts.</p>
    <div id="chart-chips"></div>
  </div>

  <div class="sel-section">
    <h2>&#128225; MQTT Parameters <span class="badge" id="mqtt-badge">0 selected</span></h2>
    <p class="sel-hint">Select any number of parameters to publish via MQTT (can be changed later).</p>
    <div id="mqtt-chips"></div>
  </div>

  <div class="sel-footer">
    <span class="hint" id="sel-hint-msg">Select 4 parameters for the graph to continue</span>
    <button id="sel-confirm-btn" disabled onclick="confirmSelection()">Confirm &amp; Start &#8594;</button>
  </div>
</div>

<!-- ===== Reconfigure button (shown in running UI) ===== -->

<div class="app">

  <header>
    <div class="patient-picker">
      <button class="patient-btn" id="patient-btn" onclick="togglePatientMenu(event)">
        PhysTwin / Patient: <span id="patient-name">--</span><span class="chev">&#9662;</span>
      </button>
      <div class="patient-menu" id="patient-menu"></div>
    </div>
    <div class="vital"><span class="k">t=</span><span class="v" id="v-t" style="min-width:148px">--</span></div>
    <div class="vital"><span class="k">HR=</span><span class="v" id="v-hr">--</span><span class="k">bpm</span></div>
    <div class="vital"><span class="k">MAP=</span><span class="v" id="v-map">--</span><span class="k">mmHg</span></div>
    <div class="vital"><span class="k">CoreT=</span><span class="v" id="v-ct">--</span><span class="k">C</span></div>
    <div class="vital"><span class="k">SkinT=</span><span class="v" id="v-st">--</span><span class="k">C</span></div>
    <div class="vital"><span class="k">Lac=</span><span class="v" id="v-lac">--</span><span class="k">mmol/L</span></div>
    <div class="vital"><span class="k">MET=</span><span class="v met-value" id="v-met">--</span><span class="k">MET</span>
    <span id="met-category" class="met-category"></span></div>
    <div id="cft-applying" class="pill off">CFT applying: OFF</div>
    <div id="cft-state" class="off">CFT afterglow: OFF</div>
    <div id="stress-decay" class="pill off">Stress decay: OFF</div>
    <div class="pill" style="background:rgba(88,166,255,0.15); color:var(--accent);">speed: <span id="v-speed">1</span>x</div>
    <div class="vital"><span class="k">actions:</span><span class="v" id="v-act" style="min-width:auto">baseline</span></div>
    <button onclick="openSelOverlay()" style="margin-left:auto;font-size:11px;color:var(--muted)">&#9881; Reconfigure</button>
  </header>

  <!-- Charts populated dynamically after param selection -->
  <div class="charts" id="charts-container"></div>

  <div class="controls">
    <div id="status">conectando...</div>
    <div class="diary-active-banner" id="diary-banner">
      &#9783; DIARY ATIVO &mdash; controlos bloqueados &mdash; apenas MQTT disponivel
    </div>

    <div class="ctl-group" id="grp-speed">
      <h3>Simulation Speed</h3>
      <div class="desc">
        Accelerate the wall-clock of loop.
      </div>
      <div class="row">
        <button id="spd-btn-1" class="speed-btn active" onclick="setSpeed(1)">1x</button>
        <button id="spd-btn-2" class="speed-btn"        onclick="setSpeed(2)">2x</button>
        <button id="spd-btn-5" class="speed-btn"        onclick="setSpeed(5)">5x</button>
      </div>
    </div>

    <div class="ctl-group" id="grp-run">
      <h3>Run (Treadmill)</h3>
      <div class="desc">SEExercise::SERunning -- speed + inclination.</div>
      <div class="row">
        <label>speed</label><input id="run-speed-kmh" type="number" step="0.1" min="0" value="10"><span>km/h</span>
        <label>incline</label><input id="run-inc" type="number" step="0.01" min="0" max="1" value="0">
      </div>
      <div class="row">
        <label>duration</label><input id="run-dur" type="number" step="1" min="0" value="20"><span>s</span>
        <button class="go" onclick="sendRun()">Run</button>
        <button class="stop" onclick="sendRaw('stop exercise')">Stop exercise</button>
      </div>
    </div>

    <div class="ctl-group" id="grp-cycle">
      <h3>Cycle (Bicycle)</h3>
      <div class="desc">SEExercise::SECycling .</div>
      <div class="row">
        <label>power</label><input id="cyc-pow" type="number" step="5" min="0" value="100"><span>W</span>
        <label>cadence</label><input id="cyc-cad" type="number" step="5" min="0" value="60"><span>/min</span>
      </div>
      <div class="row">
        <label>duration</label><input id="cyc-dur" type="number" step="1" min="0" value="30"><span>s</span>
        <button class="go" onclick="send('cycle',['cyc-pow','cyc-cad'],'cyc-dur')">Cycle</button>
      </div>
    </div>

    <div class="ctl-group" id="grp-exercise">
      <h3>Exercise (generic)</h3>
      <div class="desc">SEExercise::SEGeneric -- intensity.</div>
      <div class="row">
        <label>intensity</label><input id="ex-int" type="number" step="0.05" min="0" max="1" value="0.5">
        <label>duration</label><input id="ex-dur" type="number" step="1" min="0" value="60"><span>s</span>
        <button class="go" onclick="send('exercise',['ex-int'],'ex-dur')">Exercise</button>
      </div>
    </div>

    <div class="ctl-group" id="grp-cft">
      <h3>Cold Face Test</h3>
      <div class="desc">
        Trigeminocardiac Reflex.
      </div>
      <div class="row">
        <label>severity</label><input id="cf-sev" type="number" step="0.05" min="0" max="1" value="1.0">
        <label>duration</label><input id="cf-dur" type="number" step="5" min="0" value="90"><span>s</span>
        <button class="go" onclick="send('coldface',['cf-sev','cf-dur'])">Apply CFT</button>
        <button class="stop" onclick="sendRaw('stop coldface')">Stop CFT</button>
      </div>
    </div>

    <div class="ctl-group" id="grp-stress">
      <h3>Acute Stress</h3>
      <div class="desc">
      severity = sev x 0.90; If inside CFT Window.
      </div>
      <div class="row">
        <label>severity</label><input id="st-sev" type="number" step="0.05" min="0" max="1" value="1.0">
        <label>duration</label><input id="st-dur" type="number" step="1" min="0" value="60"><span>s</span>
        <button class="warn" onclick="send('stress',['st-sev'],'st-dur')">Stress</button>
        <button class="stop" onclick="sendRaw('stop stress')">Stop stress</button>
      </div>
    </div>

    <div class="ctl-group" id="grp-sleep">
      <h3>Sleep</h3>
      <div class="desc">
        SESleep — On/Off for Patient Sleep
      </div>
      <div class="row">
        <label>duration</label><input id="sl-dur" type="number" step="60" min="0" value="480"><span>s</span>
        <button id="sl-btn" class="go" onclick="sleepToggle()">&#128164; Sleep ON</button>
      </div>
      <div id="sl-status" style="font-size:11px;color:var(--muted);margin-top:4px">Estado: acordado</div>
    </div>

    <div class="ctl-group">
      <h3>MQTT Publish</h3>
      <div class="desc">
        Publish the data via MQTT.
        subscrive the topic:<br>
        <code style="font-size:10px;color:var(--accent)">mosquitto_sub -h &lt;broker&gt; -t &lt;topic&gt;</code>
      </div>
      <div class="row">
        <label>broker</label>
        <input id="mqtt-broker" type="text" style="width:150px" value="localhost">
        <label>port</label>
        <input id="mqtt-port" type="number" style="width:58px" value="1883">
      </div>
      <div class="row">
        <label>topic</label>
        <input id="mqtt-topic" type="text" style="width:210px" value="biogears/vitals">
      </div>
      <div class="row" style="align-items:center;gap:6px">
        <label style="cursor:pointer;display:flex;align-items:center;gap:4px;font-size:12px;font-weight:normal">
          <input id="mqtt-ws" type="checkbox" onchange="mqttWsChanged()" style="cursor:pointer">
          WebSocket
        </label>
        <span id="mqtt-ws-hint" style="font-size:10px;color:var(--muted)">(porta 8080, contorna firewall 1883)</span>
      </div>
      <div class="row">
        <button id="mqtt-btn" class="go" onclick="mqttToggle()">&#9654; Start MQTT</button>
        <span id="mqtt-status" style="font-size:11px;color:var(--muted);margin-left:6px">OFF</span>
      </div>
      <div id="mqtt-err" style="font-size:10px;color:var(--err);display:none;margin-top:4px"></div>
    </div>

    <div class="ctl-group">
      <h3>Patient Diary</h3>
      <div class="desc">24h Diary</div>
      <div class="row">
        <button class="go" style="font-size:12px" onclick="diaryOpenModal()">&#9783; Editar Diary</button>
        <span id="diary-side-status" class="inactive" style="font-size:11px;margin-left:8px">INATIVO</span>
      </div>
      <div class="row" style="margin-top:4px">
        <button id="diary-side-btn" class="go" onclick="diaryToggle()">&#9654; Start Diary</button>
        <span id="diary-side-block" style="font-size:10px;color:var(--muted);margin-left:8px">--</span>
      </div>
    </div>

    <div class="ctl-group">
      <h3>Sistem</h3>
      <div class="row">
        <button onclick="sendRaw('status')">Status (stdout)</button>
        <button class="stop" onclick="sendRaw('quit')">Quit engine</button>
      </div>
    </div>
  </div>
</div>

<!-- ===== Diary Modal ===== -->
<div class="diary-overlay" id="diary-overlay" onclick="diaryCloseOnBg(event)">
  <div class="diary-modal" onclick="event.stopPropagation()">

    <div style="display:flex;justify-content:space-between;align-items:center">
      <h2>&#9783; Patient Diary &mdash; Rotina 24h</h2>
      <button onclick="diaryCloseModal()" style="color:var(--muted);font-size:14px">&#10005;</button>
    </div>

    <!-- Timeline bar -->
    <div>
      <div style="font-size:10px;color:var(--muted);margin-bottom:4px">Linha do tempo (00:00 &rarr; 24:00)</div>
      <div class="diary-timeline" id="diary-tl"></div>
      <div style="display:flex;justify-content:space-between;font-size:9px;color:var(--muted);margin-top:3px">
        <span>00:00</span><span>06:00</span><span>12:00</span><span>18:00</span><span>24:00</span>
      </div>
    </div>

    <!-- Tabela de blocos -->
    <div>
      <div style="font-size:11px;color:var(--muted);margin-bottom:5px">Blocos definidos</div>
      <table class="diary-table">
        <thead><tr>
          <th>#</th><th>Inicio</th><th>Fim</th><th>Acao</th>
          <th>Param 1</th><th>Param 2</th><th></th>
        </tr></thead>
        <tbody id="diary-tbody"></tbody>
      </table>
    </div>

    <!-- Formulario de adicao -->
    <div class="diary-add-form">
      <label>Inicio</label>
      <input id="d-start" type="time" value="08:00" step="60" style="width:90px">
      <label>Fim (opt.)</label>
      <input id="d-end" type="time" value="" step="60" style="width:90px" placeholder="--:--">
      <label>Acao</label>
      <select id="d-action" onchange="diaryActionChanged()">
        <option value="baseline">baseline</option>
        <option value="run">run (km/h + incl)</option>
        <option value="cycle">cycle (W + rpm)</option>
        <option value="exercise">exercise (0-1)</option>
        <option value="stress">stress (0-1)</option>
        <option value="cft">CFT (sev 0-1)</option>
        <option value="sleep">sleep</option>
      </select>
      <label id="d-p1-lbl">Vel (km/h)</label>
      <input id="d-p1" type="number" step="0.5" min="0" value="8" style="width:65px">
      <label id="d-p2-lbl">Incline</label>
      <input id="d-p2" type="number" step="0.01" min="0" value="0" style="width:55px">
      <button class="go" onclick="diaryAddBlock()">+ Adicionar</button>
    </div>

    <!-- Botoes de controlo do modal -->
    <div style="display:flex;gap:8px;align-items:center;flex-wrap:wrap;padding-top:4px;border-top:1px solid var(--border)">
      <button id="diary-modal-btn" class="go" onclick="diaryToggle()">&#9654; Start Diary</button>
      <button onclick="diarySaveToServer()">&#128190; Guardar JSON</button>
      <button onclick="diaryLoadFromServer()">&#128196; Recarregar</button>
      <button class="stop" onclick="diaryClearAll()">&#128465; Limpar Tudo</button>
      <span id="diary-modal-status" class="diary-status-pill inactive">INATIVO</span>
    </div>

  </div>
</div>

<script>
'use strict';

// ══════════════════════════════════════════════════════════════════════════════
//  PARAMETER SELECTION OVERLAY
// ══════════════════════════════════════════════════════════════════════════════

// Module colour palette (same hues as diary, extended)
const MODULE_COLORS = {
  'Cardiovascular': '#e06c75',
  'Energy':         '#e5c07b',
  'Respiratory':    '#98c379',
  'Renal':          '#56b6c2',
  'Blood Chemistry':'#c678dd',
  'Derived':        '#61afef',
};

// Local state for selections
let selCatalog   = [];     // loaded from /params
let chartSel     = [];     // up to 4 ids
let mqttSel      = [];     // any ids

// Saved config (localStorage key)
const LS_KEY = 'bg_param_cfg_v1';

function selLoad() {
  try {
    const raw = localStorage.getItem(LS_KEY);
    if (raw) return JSON.parse(raw);
  } catch(_) {}
  return null;
}
function selSave(cfg) {
  localStorage.setItem(LS_KEY, JSON.stringify(cfg));
}

async function openSelOverlay() {
  document.getElementById('sel-overlay').classList.add('open');
  if (selCatalog.length === 0) {
    try {
      const r = await fetch('/params');
      const d = await r.json();
      selCatalog = d.params || [];
    } catch(e) {
      console.error('Failed to load params catalog', e);
    }
  }
  // Pre-fill from current config
  const saved = selLoad();
  chartSel = saved ? [...(saved.chart || [])] : [];
  mqttSel  = saved ? [...(saved.mqtt  || [])] : [];
  renderSelChips();
}

function closeSelOverlay() {
  document.getElementById('sel-overlay').classList.remove('open');
}

function renderSelChips() {
  renderChipGroup('chart-chips', chartSel, toggleChartSel, true);
  renderChipGroup('mqtt-chips',  mqttSel,  toggleMqttSel,  false);
  updateSelBadges();
}

function renderChipGroup(containerId, selected, toggleFn, limitTo4) {
  const container = document.getElementById(containerId);
  container.innerHTML = '';
  // Group by module
  const modules = {};
  selCatalog.forEach(p => {
    if (!modules[p.module]) modules[p.module] = [];
    modules[p.module].push(p);
  });
  Object.entries(modules).forEach(([mod, params]) => {
    const grp = document.createElement('div');
    grp.className = 'sel-module-group';
    const lbl = document.createElement('div');
    lbl.className = 'sel-module-label';
    lbl.style.color = MODULE_COLORS[mod] || '#888';
    lbl.textContent = mod;
    grp.appendChild(lbl);
    const chips = document.createElement('div');
    chips.className = 'sel-chips';
    params.forEach(p => {
      const chip = document.createElement('div');
      chip.className = 'sel-chip' + (selected.includes(p.id) ? ' active' : '');
      // Disable unchosen chips when 4 already chosen (chart mode only)
      if (limitTo4 && chartSel.length >= 4 && !selected.includes(p.id)) {
        chip.classList.add('disabled-chip');
      }
      chip.innerHTML =
        '<span class="chip-name">' + p.label + '</span>' +
        '<span class="chip-unit">' + (p.unit || 'unitless') + '</span>';
      if (p.id === 'ecg') {
        // "8s" deve acompanhar kEcgBufferSeconds no backend (C++).
        chip.title = 'Forma de onda de alta resolução (Lead III, ~50 Hz, janela de 8s) '
                   + '— gráfico dedicado tipo osciloscópio, independente do polling '
                   + 'a 1 Hz dos restantes parâmetros.';
      }
      chip.onclick = () => toggleFn(p.id);
      chips.appendChild(chip);
    });
    grp.appendChild(chips);
    container.appendChild(grp);
  });
}

function toggleChartSel(id) {
  const i = chartSel.indexOf(id);
  if (i >= 0) {
    chartSel.splice(i, 1);
  } else {
    if (chartSel.length >= 4) return; // already full
    chartSel.push(id);
  }
  renderSelChips();
}

function toggleMqttSel(id) {
  const i = mqttSel.indexOf(id);
  if (i >= 0) mqttSel.splice(i, 1);
  else         mqttSel.push(id);
  renderSelChips();
}

function updateSelBadges() {
  document.getElementById('chart-badge').textContent = chartSel.length + ' / 4 selected';
  document.getElementById('mqtt-badge').textContent  = mqttSel.length + ' selected';
  const btn  = document.getElementById('sel-confirm-btn');
  const hint = document.getElementById('sel-hint-msg');
  const ready = chartSel.length === 4;
  btn.disabled = !ready;
  hint.textContent = ready
    ? (mqttSel.length > 0 ? mqttSel.length + ' MQTT param(s) selected' : 'No MQTT params selected (optional)')
    : 'Select exactly 4 parameters for the graph to continue';
}

async function confirmSelection() {
  if (chartSel.length !== 4) return;
  // Persist locally
  selSave({ chart: chartSel, mqtt: mqttSel });
  // Push MQTT params to server
  try {
    await fetch('/config', {
      method: 'POST',
      headers: {'Content-Type':'application/json'},
      body: JSON.stringify({ mqtt_params: mqttSel })
    });
  } catch(_) {}
  closeSelOverlay();
  initChartsFromSelection(chartSel);
}

// ── Chart initialisation ────────────────────────────────────────────────────
const CHART_PALETTE = ['#a371f7','#58a6ff','#3fb950','#e5c07b','#e06c75','#56b6c2','#c678dd','#61afef'];
let charts = [];        // active Chart.js instances
let chartParamIds = []; // current chart param ids

// ECG: grafico + polling proprios (ver comentario junto a makeEcgChart).
// charts[idx] fica a null na posicao do slot 'ecg' para que o poll() geral
// (que corre a 1 Hz) o ignore automaticamente (if (charts[idx]) ...).
let ecgChart = null;
let ecgTimer = null;

function getCatalogParam(id) {
  return selCatalog.find(p => p.id === id) || {id, label:id, unit:'', module:''};
}

function initChartsFromSelection(ids) {
  chartParamIds = ids;
  const container = document.getElementById('charts-container');
  // Destroy old charts (charts[i] pode ser null no slot do ECG — ver abaixo)
  charts.forEach(c => c && c.destroy());
  charts = [];
  if (ecgChart) { ecgChart.destroy(); ecgChart = null; }
  if (ecgTimer) { clearInterval(ecgTimer); ecgTimer = null; }
  container.innerHTML = '';
  ids.forEach((id, idx) => {
    const p = getCatalogParam(id);
    const cell = document.createElement('div');
    cell.className = 'chart-cell';
    cell.innerHTML = '<div class="chart-title">' + p.label + (p.unit ? ' (' + p.unit + ')' : '') + '</div>' +
                     '<canvas id="c-' + id + '"></canvas>';
    container.appendChild(cell);
    if (id === 'ecg') {
      // Forma de onda de alta resolucao: grafico e polling dedicados
      // (o resto dos parametros usa poll() a 1 Hz, insuficiente para o QRS —
      // ver comentario em makeEcgChart / pollEcg).
      ecgChart = makeEcgChart('c-' + id);
      charts.push(null); // mantem charts[] alinhado por indice com chartParamIds
      ecgTimer = setInterval(pollEcg, 250);
      pollEcg();
    } else {
      const color = CHART_PALETTE[idx % CHART_PALETTE.length];
      charts.push(makeChart('c-' + id, color, p.label));
    }
  });
}

// ── ECG: grafico e poll dedicados ───────────────────────────────────────────
// Os outros parametros sao amostrados a 1 Hz (ver poll(), setInterval(poll,1000))
// porque sao grandezas fisiologicas lentas (HR, MAP, temperaturas, ...). O ECG
// e diferente: é a propria forma de onda eletrica do ciclo cardiaco, com um
// complexo QRS de ~80-100 ms — a 1 Hz nao ha sequer 1 amostra por batimento.
// O backend expõe /ecg com um buffer a kEngineHz=50Hz (nao decimado); aqui
// fazemos poll desse endpoint a 250 ms (4 Hz de refresco do traçado, mais do
// que suficiente para o olho, já que os dados em si vêm a 50 Hz por amostra).
function makeEcgChart(canvasId) {
  return new Chart(document.getElementById(canvasId), {
    type: 'line',
    data: { datasets: [{ label: 'ECG (Lead III)', data: [],
                          borderColor: '#e06c75', fill:false, tension:0.15 }] },
    options: {
      responsive: true, maintainAspectRatio: false, animation: false,
      layout: { padding: { top: 16 } },
      scales: {
        x: {
          type: 'linear',
          ticks: { color:'#8b949e', maxTicksLimit:6, callback: v => v.toFixed(1) + 's' },
          grid: { color:'#21262d' }
        },
        y: { ticks:{ color:'#8b949e' }, grid:{ color:'#21262d' } }
      },
      plugins: { legend: { display:false } },
      elements: { point: { radius:0 }, line: { borderWidth:1.3 } }
    }
  });
}

async function pollEcg() {
  if (!ecgChart) return;
  try {
    const r = await fetch('/ecg');
    if (!r.ok) return;
    const d = await r.json();
    const s = d.samples || [];
    if (s.length === 0) return;
    // Eixo X em segundos relativos a amostra mais recente (0 = "agora"),
    // como um osciloscopio: o traçado desliza da direita para a esquerda.
    const t0 = s[s.length - 1].t;
    ecgChart.data.datasets[0].data = s
      .filter(p => p.mv !== null && isFinite(p.mv))
      .map(p => ({ x: p.t - t0, y: p.mv }));
    ecgChart.update('none');
  } catch (_) { /* falha silenciosa: mantem o ultimo traçado valido no ecrã */ }
}

// ══════════════════════════════════════════════════════════════════════════════
//  STARTUP: check if selection exists, else show overlay
// ══════════════════════════════════════════════════════════════════════════════
async function startup() {
  const saved = selLoad();
  // Pre-load catalog (always needed for overlay)
  try {
    const r = await fetch('/params');
    const d = await r.json();
    selCatalog = d.params || [];
  } catch(_) {}

  if (saved && saved.chart && saved.chart.length === 4) {
    // Restore saved selection
    chartSel = saved.chart;
    mqttSel  = saved.mqtt || [];
    // Push MQTT params to server (in case of engine restart)
    try {
      await fetch('/config', {
        method:'POST', headers:{'Content-Type':'application/json'},
        body: JSON.stringify({ mqtt_params: mqttSel })
      });
    } catch(_) {}
    initChartsFromSelection(chartSel);
  } else {
    // Show selection overlay
    chartSel = [];
    mqttSel  = [];
    renderSelChips();
    document.getElementById('sel-overlay').classList.add('open');
  }
}

// ══════════════════════════════════════════════════════════════════════════════
//  CHARTS
// ══════════════════════════════════════════════════════════════════════════════

function getMETCategory(met) {
    if (met < 1.5) return "Sedentary";
    if (met < 3) return "Light activity";
    if (met < 6) return "Moderate activity";
    if (met < 9) return "Vigorous activity";
    return "Very vigorous activity";
}

// ---- Helpers de timestamp --------------------------------------------------
// Formata epoch ms para string legivel (sem dependencias externas)
function formatTs(ts_ms) {
  if (!ts_ms) return '--';
  const d = new Date(ts_ms);
  const p = n => String(n).padStart(2, '0');
  return d.getFullYear() + '-' + p(d.getMonth()+1) + '-' + p(d.getDate()) +
         ' ' + p(d.getHours()) + ':' + p(d.getMinutes()) + ':' + p(d.getSeconds());
}

// Formata segundos de simulacao como "HH:MM:SS" para o eixo de diary
function formatDiaryTime(t_s) {
  if (t_s == null || !isFinite(t_s)) return '--';
  const total = Math.floor(t_s);
  const h = Math.floor(total / 3600);
  const m = Math.floor((total % 3600) / 60);
  const s = total % 60;
  const p = n => String(n).padStart(2, '0');
  return p(h) + ':' + p(m) + ':' + p(s);
}

// ---- Charts ----------------------------------------------------------------
const commonOpts = (yLabel) => ({
  responsive: true, maintainAspectRatio: false, animation: false,
  layout: { padding: { top: 16 } },
  scales: {
    x: {
      type: 'time',
      time: {
        displayFormats: { second: 'HH:mm:ss', minute: 'HH:mm', hour: 'HH:mm' },
        tooltipFormat: 'yyyy-MM-dd HH:mm:ss'
      },
      adapters: { date: { locale: 'pt' } },
      ticks: { color:'#8b949e', maxTicksLimit:6, maxRotation:0 },
      grid: { color:'#21262d' }
    },
    y: { ticks:{color:'#8b949e'}, grid:{color:'#21262d'} }
  },
  plugins: { legend: { display:false } },
  elements: { point: { radius:0 }, line: { borderWidth:1.4 } }
});

function makeChart(canvasId, color, yLabel) {
  return new Chart(document.getElementById(canvasId), {
    type: 'line',
    data: { datasets: [{ label: yLabel, data: [], borderColor: color, fill:false, tension:0 }] },
    options: commonOpts(yLabel)
  });
}

// Helper: get a value from a sample's v map
function sv(sample, key) {
  if (!sample || !sample.v) return null;
  const v = sample.v[key];
  return (v !== undefined && v !== null && isFinite(v)) ? v : null;
}

// Polling
async function poll() {
  try {
    const r = await fetch('/data');
    if (!r.ok) throw new Error('HTTP ' + r.status);
    const d = await r.json();

    const n = d.samples.length;
    if (n > 0) {
      const last = d.samples[n-1];
      // Mostra diary time (HH:MM:SS) quando diary activo, senao timestamp real
      document.getElementById('v-t').textContent = diaryActive
        ? formatDiaryTime(last.t)
        : formatTs(last.ts_ms);
      const hr    = sv(last,'hr');
      const map_  = sv(last,'map');
      const coreT = sv(last,'coreT');
      const skinT = sv(last,'skinT');
      const lac   = sv(last,'lac');
      const met   = sv(last,'met');
      document.getElementById('v-hr').textContent  = hr    !== null ? hr.toFixed(1)    : '--';
      document.getElementById('v-map').textContent = map_  !== null ? map_.toFixed(1)  : '--';
      document.getElementById('v-ct').textContent  = coreT !== null ? coreT.toFixed(2) : '--';
      document.getElementById('v-st').textContent  = skinT !== null ? skinT.toFixed(2) : '--';
      document.getElementById('v-lac').textContent = lac   !== null ? lac.toFixed(3)   : '--';
      document.getElementById('v-met').textContent = met   !== null ? met.toFixed(2)   : '--';
      if (met !== null) document.getElementById('met-category').textContent = getMETCategory(met);

      document.getElementById('v-act').textContent = last.actions || 'baseline';
      syncSleepFromActions(last.actions || '');
    }

    const cftEl = document.getElementById('cft-state');
    if (d.cft_afterglow_on) {
      cftEl.className = 'on';
      cftEl.textContent = 'CFT afterglow: ON (stress x0.90)';
    } else {
      cftEl.className = 'off';
      cftEl.textContent = 'CFT afterglow: OFF';
    }

    // Indicador "CFT actualmente a actuar" (engine a modular vagal/simp.periferico).
    const cftApp = document.getElementById('cft-applying');
    if (d.cft_applying) {
      cftApp.className = 'pill on';
      const el = (d.cft_elapsed_s !== null && d.cft_dur_s) ? d.cft_elapsed_s.toFixed(0) : '?';
      const du = (d.cft_dur_s !== null) ? d.cft_dur_s.toFixed(0) : '?';
      cftApp.textContent = 'CFT applying: ON  (t=' + el + '/' + du + ' s)';
    } else {
      cftApp.className = 'pill off';
      cftApp.textContent = 'CFT applying: OFF';
    }

    // Indicador do decay exponencial do stress pos-acao.
    const stDec = document.getElementById('stress-decay');
    if (d.stress_decay_active) {
      stDec.className = 'pill on';
      const sNow = (typeof d.stress_decay_sev_now === 'number') ? d.stress_decay_sev_now.toFixed(3) : '?';
      const tau  = (typeof d.stress_decay_tau_s === 'number')   ? d.stress_decay_tau_s.toFixed(0)   : '?';
      const ela  = (typeof d.stress_decay_elapsed === 'number') ? d.stress_decay_elapsed.toFixed(0) : '?';
      stDec.textContent = 'Stress decay: ON (sev=' + sNow + ', τ=' + tau + 's, t=' + ela + 's)';
    } else {
      stDec.className = 'pill off';
      stDec.textContent = 'Stress decay: OFF';
    }

    // Sincroniza nome do paciente exibido no botao do header.
    if (typeof d.patient_name === 'string' && d.patient_name.length > 0) {
      const el = document.getElementById('patient-name');
      if (el && el.textContent !== d.patient_name) {
        el.textContent = d.patient_name;
      }
    }

    // Sincroniza display de velocidade (caso outro cliente a tenha mudado,
    // ou o factor inicial seja != 1).
    if (typeof d.speed_factor === 'number' && isFinite(d.speed_factor)) {
      const sStr = (Number.isInteger(d.speed_factor)) ? d.speed_factor.toFixed(0) : d.speed_factor.toFixed(2);
      document.getElementById('v-speed').textContent = sStr;
      // Marca o botao correspondente como activo, se for um valor preset.
      ['1','2','5'].forEach(v => {
        const b = document.getElementById('spd-btn-' + v);
        if (b) b.classList.toggle('active', parseFloat(v) === d.speed_factor);
      });
    }

    // Rolling window: ultimos 990 s de tempo de simulacao (t_s)
    const tMax = (n > 0) ? d.samples[n-1].t : 0;
    const tMin = Math.max(0, tMax - 990);
    const pts  = d.samples.filter(s => s.t >= tMin && s.ts_ms > 0);

    // Eixo X: diary activo → t*1000ms (UTC); senao → ts_ms (wall-clock)
    const xOf = (s) => diaryActive ? Math.round(s.t * 1000) : s.ts_ms;
    const xy = (key) => pts
      .filter(s => s.v && s.v[key] !== undefined && s.v[key] !== null && isFinite(s.v[key]))
      .map(s => ({x: xOf(s), y: s.v[key]}));

    // Update dynamic charts
    chartParamIds.forEach((id, idx) => {
      if (charts[idx]) {
        charts[idx].data.datasets[0].data = xy(id);
        charts[idx].update('none');
      }
    });

    // Estado MQTT na barra de status
    const lastTs = (n > 0 && d.samples[n-1].ts_ms > 0) ? formatTs(d.samples[n-1].ts_ms) : '--';
    document.getElementById('status').textContent =
      'conectado | ' + n + ' samples | ' + lastTs;

    // Sincroniza painel MQTT
    if (typeof d.mqtt_active === 'boolean') {
      const mqttOn   = d.mqtt_active;
      const mqttConn = d.mqtt_connecting;
      const btn  = document.getElementById('mqtt-btn');
      const st   = document.getElementById('mqtt-status');
      const errEl = document.getElementById('mqtt-err');
      if (btn) {
        if (mqttOn) {
          btn.innerHTML = '&#9632; Stop MQTT';
          btn.className = 'stop';
          btn.disabled = false;
          if (st) { st.style.color = 'var(--ok)';
                    st.textContent = 'ON → ' + d.mqtt_broker + ':' + d.mqtt_port + ' / ' + d.mqtt_topic; }
        } else if (mqttConn) {
          btn.innerHTML = '&#9711; Conectando...';
          btn.className = 'warn';
          btn.disabled = true;
          if (st) { st.style.color = 'var(--warn)';
                    st.textContent = 'Ligando a ' + d.mqtt_broker + ':' + d.mqtt_port + ' (timeout 10 s)'; }
        } else {
          btn.innerHTML = '&#9654; Start MQTT';
          btn.className = 'go';
          btn.disabled = false;
          if (st) { st.style.color = 'var(--muted)'; st.textContent = 'OFF'; }
        }
      }
      if (errEl) {
        if (d.mqtt_error && d.mqtt_error.length > 0) {
          errEl.textContent = '⚠ ' + d.mqtt_error; errEl.style.display = 'block';
        } else {
          errEl.style.display = 'none';
        }
      }
      // Preenche inputs com valores actuais do broker/topic (se nao modificados)
      const bi = document.getElementById('mqtt-broker');
      const pi = document.getElementById('mqtt-port');
      const ti = document.getElementById('mqtt-topic');
      if (!mqttOn && !mqttConn && bi && !bi._dirty) bi.value = d.mqtt_broker || bi.value;
      if (!mqttOn && !mqttConn && pi && !pi._dirty) pi.value = d.mqtt_port   || pi.value;
      if (!mqttOn && !mqttConn && ti && !ti._dirty) ti.value = d.mqtt_topic  || ti.value;
      // Reflecte modo WS no checkbox (enquanto ativo)
      const wsCb = document.getElementById('mqtt-ws');
      if (wsCb && typeof d.mqtt_ws === 'boolean' && (mqttOn || mqttConn)) {
        wsCb.checked = d.mqtt_ws;
      }
    }

    // Sincroniza estado do diary
    diaryUpdateFromPoll(d);

  } catch (e) {
    document.getElementById('status').textContent = 'erro: ' + e.message;
  }
}
// Startup: show selection overlay or restore saved config, then start polling
startup().then(() => {
  setInterval(poll, 1000);
  poll();
});

// Command dispatch 
async function sendRaw(cmd) {
  try {
    const r = await fetch('/cmd', {
      method: 'POST',
      headers: {'Content-Type':'application/json'},
      body: JSON.stringify({cmd})
    });
    const t = await r.text();
    document.getElementById('status').textContent = '> ' + cmd + ' [' + r.status + '] ' + t;
  } catch (e) {
    document.getElementById('status').textContent = 'erro envio: ' + e.message;
  }
}

function val(id) {
  const e = document.getElementById(id);
  return e ? e.value : '';
}

// send('run', ['run-speed','run-inc'], 'run-dur')
// constroi linha:  "run <speed> <inc> for <dur>"
// "coldface" usa <sev> <dur> diretamente (sem "for"), conforme o dispatcher.
function send(cmd, paramIds, durId) {
  const parts = [cmd, ...paramIds.map(val)];
  if (durId) {
    const d = val(durId);
    if (d && parseFloat(d) > 0) parts.push('for', d);
  }
  sendRaw(parts.join(' '));
}

// Run (treadmill): input em km/h, convertido para m/s antes de enviar para o engine.
// Conversao exata: v[m/s] = v[km/h] / 3.6  (1 km/h = 1000 m / 3600 s).
// O backend (SEExercise::SERunning) usa LengthPerTimeUnit::m_Per_s nativamente,
// entao mantemos o protocolo do dispatcher inalterado e fazemos a conversao aqui.
function sendRun() {
  const kmh = parseFloat(val('run-speed-kmh'));
  const inc = val('run-inc');
  const dur = val('run-dur');
  if (!isFinite(kmh) || kmh < 0) {
    document.getElementById('status').textContent = 'erro: speed (km/h) invalida';
    return;
  }
  const mps = (kmh / 3.6).toFixed(4);   // 4 casas: suficiente para 0.1 km/h de resolucao
  const parts = ['run', mps, inc];
  if (dur && parseFloat(dur) > 0) parts.push('for', dur);
  sendRaw(parts.join(' '));
}

/// Quando o checkbox WS muda: ajusta a porta default (1883 <-> 8080)
function mqttWsChanged() {
  const cb = document.getElementById('mqtt-ws');
  const pi = document.getElementById('mqtt-port');
  if (!cb || !pi) return;
  // Só altera se o utilizador não modificou manualmente (ou se ainda está no default)
  const curPort = parseInt(pi.value, 10);
  if (cb.checked) {
    if (curPort === 1883 || !pi._dirty) { pi.value = '8080'; pi._dirty = false; }
  } else {
    if (curPort === 8080 || !pi._dirty) { pi.value = '1883'; pi._dirty = false; }
  }
}

// MQTT toggle: envia "mqtt start broker port topic [ws]" ou "mqtt stop"
function mqttToggle() {
  const btn = document.getElementById('mqtt-btn');
  // Se o botão diz Stop => parar
  if (btn && btn.className === 'stop') {
    sendRaw('mqtt stop');
  } else {
    const broker = val('mqtt-broker') || 'localhost';
    const port   = val('mqtt-port')   || '1883';
    const topic  = val('mqtt-topic')  || 'biogears/vitals';
    const useWs  = document.getElementById('mqtt-ws') && document.getElementById('mqtt-ws').checked;
    const cmd    = 'mqtt start ' + broker + ' ' + port + ' ' + topic + (useWs ? ' ws' : '');
    sendRaw(cmd);
  }
}

// Marca os inputs MQTT como "sujos" (utilizador editou) para nao os sobrescrever no poll
(function() {
  ['mqtt-broker','mqtt-port','mqtt-topic'].forEach(id => {
    const el = document.getElementById(id);
    if (el) el.addEventListener('input', () => { el._dirty = true; });
  });
  // checkbox WS: marca porta como não-dirty quando o utilizador altera (mqttWsChanged cuida disso)
  const wsCb = document.getElementById('mqtt-ws');
  if (wsCb) wsCb.addEventListener('change', () => {
    const pi = document.getElementById('mqtt-port');
    if (pi) pi._dirty = false; // permite que mqttWsChanged ajuste a porta
  });
})();

// ---- Patient Diary UI -------------------------------------------------------

// Mapa de cores por tipo de accao
const DIARY_COLORS = {
  baseline: '#21262d', run: '#1f6feb', cycle: '#388bfd',
  exercise: '#58a6ff', stress: '#d29922', cft: '#0d419d',
  sleep: '#6e40c9',
};

// Estado local do diary (sincronizado com servidor via poll e fetch)
let diaryBlocks = [];   // array de {start, action, param1, param2, start_min}
let diaryActive = false;

// ---- helpers ----
function hhmm2min(s) {
  const [h,m] = s.split(':').map(Number);
  return (h||0)*60 + (m||0);
}
function min2hhmm(m) {
  return String(Math.floor(m/60)).padStart(2,'0')+':'+String(m%60).padStart(2,'0');
}

function diaryOpenModal() {
  document.getElementById('diary-overlay').classList.add('open');
  diaryRenderModal();
}
function diaryCloseModal() {
  document.getElementById('diary-overlay').classList.remove('open');
}
function diaryCloseOnBg(e) {
  if (e.target === document.getElementById('diary-overlay')) diaryCloseModal();
}

// Ajusta labels e visibilidade de param2 conforme a acao selecionada
function diaryActionChanged() {
  const a = document.getElementById('d-action').value;
  const p1l = document.getElementById('d-p1-lbl');
  const p1i = document.getElementById('d-p1');
  const p2l = document.getElementById('d-p2-lbl');
  const p2i = document.getElementById('d-p2');
  const cfg = {
    baseline: {l1:null,  l2:null,     d1:0,   d2:0},
    run:      {l1:'Vel (km/h)', l2:'Incline', d1:8,   d2:0},
    cycle:    {l1:'Power (W)', l2:'Cadence (rpm)', d1:100, d2:60},
    exercise: {l1:'Intensity (0-1)', l2:null, d1:0.5, d2:0},
    stress:   {l1:'Severity (0-1)', l2:null, d1:0.5, d2:0},
    cft:      {l1:'Severity (0-1)', l2:null, d1:0.5, d2:0},
    sleep:    {l1:null,             l2:null, d1:0,   d2:0},
  };
  const c = cfg[a] || {l1:'Param1',l2:'Param2',d1:0,d2:0};
  const show1 = c.l1 != null;
  const show2 = c.l2 != null;
  p1l.textContent = c.l1 || 'Param1';
  p1l.style.display = show1 ? '' : 'none';
  p1i.style.display = show1 ? '' : 'none';
  if (show1) p1i.value = c.d1;
  p2l.textContent = c.l2 || 'Param2';
  p2l.style.display = show2 ? '' : 'none';
  p2i.style.display = show2 ? '' : 'none';
  if (show2) p2i.value = c.d2;
}

function diaryAddBlock() {
  const startRaw = document.getElementById('d-start').value || '08:00';
  const endRaw   = document.getElementById('d-end').value || '';
  const action   = document.getElementById('d-action').value;
  const p1       = parseFloat(document.getElementById('d-p1').value) || 0;
  const p2       = parseFloat(document.getElementById('d-p2').value) || 0;
  const start_min = hhmm2min(startRaw);
  const end_min   = endRaw ? hhmm2min(endRaw) : 0;
  // Remove duplicado no mesmo minuto de inicio
  const i = diaryBlocks.findIndex(b => b.start_min === start_min);
  if (i >= 0) diaryBlocks.splice(i, 1);
  diaryBlocks.push({ start: min2hhmm(start_min), end: endRaw, action, param1: p1, param2: p2, start_min, end_min });
  diaryBlocks.sort((a,b) => a.start_min - b.start_min);
  diaryRenderModal();
}

function diaryRemoveBlock(idx) {
  diaryBlocks.splice(idx, 1);
  diaryRenderModal();
}

function diaryRenderModal() {
  // ---- Timeline bar ----
  const tl = document.getElementById('diary-tl');
  if (tl) {
    tl.innerHTML = '';
    diaryBlocks.forEach((b, i) => {
      const next = diaryBlocks[i+1];
      const pct  = (b.start_min / 1440) * 100;
      // Fim da barra: end_min (se definido) ou inicio do proximo bloco ou 24:00
      const barEnd = (b.end_min > b.start_min) ? b.end_min
                   : (next ? next.start_min : 1440);
      const wPct = (barEnd - b.start_min) / 1440 * 100;
      const div  = document.createElement('div');
      div.className = 'diary-block-bar';
      div.style.left = pct + '%';
      div.style.width = Math.max(wPct, 1) + '%';
      div.style.background = DIARY_COLORS[b.action] || '#888';
      div.title = b.start + (b.end_min > b.start_min ? '-'+b.end : '') + ' ' + b.action + ' p1=' + b.param1;
      div.textContent = b.action.substring(0,4);
      tl.appendChild(div);
    });
  }
  // ---- Table ----
  const tbody = document.getElementById('diary-tbody');
  if (tbody) {
    tbody.innerHTML = '';
    diaryBlocks.forEach((b, i) => {
      const tr = document.createElement('tr');
      const color = DIARY_COLORS[b.action] || '#fff';
      const endLabel = (b.end_min > b.start_min) ? b.end : '<span style="color:var(--muted)">→next</span>';
      tr.innerHTML =
        '<td>' + i + '</td>' +
        '<td>' + b.start + '</td>' +
        '<td>' + endLabel + '</td>' +
        '<td style="color:' + color + ';font-weight:600">' + b.action + '</td>' +
        '<td>' + b.param1 + '</td>' +
        '<td>' + b.param2 + '</td>' +
        '<td><button class="stop" style="padding:2px 7px;font-size:10px" ' +
          'onclick="diaryRemoveBlock(' + i + ')">&#10005;</button></td>';
      tbody.appendChild(tr);
    });
    if (!diaryBlocks.length) {
      const tr = document.createElement('tr');
      tr.innerHTML = '<td colspan="6" style="color:var(--muted);text-align:center;padding:10px">Sem blocos definidos</td>';
      tbody.appendChild(tr);
    }
  }
  // ---- Botoes Start/Stop (modal e side panel) ----
  const modalBtn  = document.getElementById('diary-modal-btn');
  const modalSt   = document.getElementById('diary-modal-status');
  const sideBtn   = document.getElementById('diary-side-btn');
  const sideSt    = document.getElementById('diary-side-status');
  if (diaryActive) {
    if (modalBtn) { modalBtn.className='stop'; modalBtn.textContent='⏹ Stop Diary'; }
    if (modalSt)  { modalSt.className='diary-status-pill active'; modalSt.textContent='ATIVO'; }
    if (sideBtn)  { sideBtn.className='stop'; sideBtn.textContent='⏹ Stop Diary'; }
    if (sideSt)   { sideSt.className='active'; sideSt.textContent='ATIVO'; }
  } else {
    if (modalBtn) { modalBtn.className='go'; modalBtn.innerHTML='&#9654; Start Diary'; }
    if (modalSt)  { modalSt.className='diary-status-pill inactive'; modalSt.textContent='INATIVO'; }
    if (sideBtn)  { sideBtn.className='go'; sideBtn.innerHTML='&#9654; Start Diary'; }
    if (sideSt)   { sideSt.className='inactive'; sideSt.textContent='INATIVO'; }
  }
}

// Envia os blocos actuais para o servidor (POST /diary)
async function diarySaveToServer() {
  const payload = { blocks: diaryBlocks.map(b => ({
    start: b.start, end: b.end || '', action: b.action, param1: b.param1, param2: b.param2
  }))};
  try {
    await fetch('/diary', { method:'POST',
      headers:{'Content-Type':'application/json'},
      body: JSON.stringify(payload) });
  } catch(e) { console.error('[diary] save error', e); }
}

// Carrega blocos do servidor (GET /diary)
async function diaryLoadFromServer() {
  try {
    const r = await fetch('/diary');
    if (!r.ok) return;
    const d = await r.json();
    diaryActive = !!d.active;
    diaryBlocks = (d.blocks||[]).map(b => ({
      ...b,
      start_min: hhmm2min(b.start || '00:00'),
      end_min:   (b.end && b.end.length >= 5) ? hhmm2min(b.end) : 0,
    }));
    diaryBlocks.sort((a,b)=>a.start_min-b.start_min);
    diaryRenderModal();
  } catch(e) { console.error('[diary] load error', e); }
}

// Toggle Start/Stop do diary
async function diaryToggle() {
  if (diaryActive) {
    sendRaw('diary stop');
  } else {
    if (!diaryBlocks.length) {
      alert('Adiciona pelo menos um bloco antes de iniciar o diary.');
      diaryOpenModal(); return;
    }
    // 1) Envia blocos ao servidor
    await diarySaveToServer();
    // 2) Carrega do ficheiro e arranca
    sendRaw('diary load');
    setTimeout(() => sendRaw('diary start'), 300);
  }
}

// Limpa todos os blocos (com confirmacao)
function diaryClearAll() {
  if (!confirm('Apagar todos os blocos do diary?')) return;
  diaryBlocks = [];
  diaryRenderModal();
}

// Actualiza o modo do eixo X dos graficos: diary-time (UTC, 00:00→) vs wall-clock (local)
function updateChartAxisMode(isDiary) {
  const zone  = isDiary ? 'UTC' : undefined;
  const ttFmt = isDiary ? 'HH:mm:ss' : 'yyyy-MM-dd HH:mm:ss';
  for (const ch of charts) {
    ch.options.scales.x.adapters.date.zone = zone;
    ch.options.scales.x.time.tooltipFormat = ttFmt;
    ch.update('none');
  }
}

// Bloqueia / desbloqueia os grupos de controlo fisiologico
// Unico controlo permitido durante diary: MQTT (e Diary stop + Quit)
const DIARY_LOCKED_GROUPS = ['grp-run','grp-cycle','grp-exercise','grp-cft','grp-stress','grp-sleep'];

function diaryLockControls(lock) {
  // Grupos de controlo fisiologico
  for (const id of DIARY_LOCKED_GROUPS) {
    const el = document.getElementById(id);
    if (el) {
      if (lock) el.classList.add('diary-locked');
      else      el.classList.remove('diary-locked');
    }
  }
  // Botao de paciente no header
  const pb = document.getElementById('patient-btn');
  if (pb) {
    if (lock) pb.classList.add('diary-locked');
    else      pb.classList.remove('diary-locked');
  }
  // Banner de aviso
  const banner = document.getElementById('diary-banner');
  if (banner) {
    if (lock) banner.classList.add('show');
    else      banner.classList.remove('show');
  }
}

// Actualiza estado do diary a partir dos dados do poll /data
function diaryUpdateFromPoll(d) {
  if (typeof d.diary_active !== 'boolean') return;
  const wasActive = diaryActive;
  diaryActive = d.diary_active;
  // Aplica lock e modo de eixo sempre que o estado muda
  if (wasActive !== diaryActive) {
    diaryLockControls(diaryActive);
    updateChartAxisMode(diaryActive);
  }
  diaryRenderModal();
}

// ---- Sleep toggle ----------------------------------------------------------
let slSleeping = false;

function sleepToggle() {
  const dur = parseFloat(val('sl-dur'));
  if (!slSleeping) {
    // Liga o sono
    const cmd = (dur > 0) ? ('sleep for ' + dur) : 'sleep';
    sendRaw(cmd);
    slSleeping = true;
  } else {
    // Acorda
    sendRaw('stop sleep');
    slSleeping = false;
  }
  _updateSleepUI();
}

function _updateSleepUI() {
  const btn = document.getElementById('sl-btn');
  const st  = document.getElementById('sl-status');
  if (!btn) return;
  if (slSleeping) {
    btn.innerHTML = '&#9728; Wake UP';
    btn.className = 'stop';
    if (st) { st.style.color = 'var(--accent)'; st.textContent = 'Estado: a dormir (SESleep ON)'; }
  } else {
    btn.innerHTML = '&#128164; Sleep ON';
    btn.className = 'go';
    if (st) { st.style.color = 'var(--muted)'; st.textContent = 'Estado: acordado'; }
  }
}

// Sincroniza o estado do sleep com o campo actions do poll
function syncSleepFromActions(actionsStr) {
  const wasSleeping = slSleeping;
  slSleeping = actionsStr && actionsStr.includes('sleep=on');
  if (wasSleeping !== slSleeping) _updateSleepUI();
}

// Speed control: ajusta o wall-clock do loop (dt do engine permanece 0.02 s).
// Optimista: actualiza UI imediatamente; o poll() seguinte reconcilia caso o
// backend tenha aplicado outro valor.
function setSpeed(s) {
  sendRaw('speed ' + s);
  document.getElementById('v-speed').textContent = s;
  ['1','2','5'].forEach(v => {
    const b = document.getElementById('spd-btn-' + v);
    if (b) b.classList.toggle('active', parseFloat(v) === s);
  });
}

// ---------- Patient picker ----------
// Mostra/oculta o menu de pacientes. Listagem e' obtida lazy via GET /patients.
async function togglePatientMenu(ev) {
  if (ev) ev.stopPropagation();
  const menu = document.getElementById('patient-menu');
  const willOpen = !menu.classList.contains('open');
  if (willOpen) {
    await loadPatients();
    menu.classList.add('open');
  } else {
    menu.classList.remove('open');
  }
}

// Fecha o menu se o utilizador clicar fora.
document.addEventListener('click', (ev) => {
  const picker = document.querySelector('.patient-picker');
  if (picker && !picker.contains(ev.target)) {
    document.getElementById('patient-menu').classList.remove('open');
  }
});

async function loadPatients() {
  const menu = document.getElementById('patient-menu');
  menu.innerHTML = '<div class="item" style="color:var(--muted)">a carregar...</div>';
  try {
    const r = await fetch('/patients');
    const d = await r.json();
    const current = d.current || '';
    menu.innerHTML = '';
    if (!d.patients || d.patients.length === 0) {
      menu.innerHTML = '<div class="item" style="color:var(--err)">nenhum XML encontrado</div>';
      return;
    }
    for (const p of d.patients) {
      const div = document.createElement('div');
      div.className = 'item' + ((p.path === current) ? ' current' : '');
      div.textContent = p.name;
      div.title = p.path;
      div.onclick = () => switchPatient(p.path, p.name, p.path === current);
      menu.appendChild(div);
    }
  } catch (e) {
    menu.innerHTML = '<div class="item" style="color:var(--err)">erro: ' + e.message + '</div>';
  }
}

// Solicita troca de paciente. BioGears nao permite trocar em runtime sem
// destruir + recriar o engine; portanto qualquer troca implica reset total
// da simulacao (t -> 0, snapshot limpo, scheduled stops zerados, trackers
// CFT e stress_decay zerados). O backend faz tudo isso atomicamente no
// inicio do tick seguinte ao receber "patient <file.xml>".
function switchPatient(path, name, isCurrent) {
  document.getElementById('patient-menu').classList.remove('open');
  if (isCurrent) {
    document.getElementById('status').textContent =
      'paciente "' + name + '" ja esta carregado.';
    return;
  }
  const ok = window.confirm(
    'Trocar para o paciente "' + name + '"?\n\n' +
    'A simulacao reiniciara: o tempo voltara para t=0 e todos os ' +
    'estados activos (exercicio, stress, CFT, decay) serao limpos.\n\n' +
    'O paciente actual sera substituido por:\n  ' + path 
  );
  if (!ok) return;
  sendRaw('patient ' + path);
}
</script>
</body>
</html>
)HTML";

//  HTTP server minimal

namespace {

static std::string http_response(int status, const std::string& reason,
                                  const std::string& ctype,
                                  const std::string& body)
{
  std::ostringstream o;
  o << "HTTP/1.1 " << status << " " << reason << "\r\n"
    << "Content-Type: " << ctype << "\r\n"
    << "Content-Length: " << body.size() << "\r\n"
    << "Access-Control-Allow-Origin: *\r\n"
    << "Connection: close\r\n\r\n"
    << body;
  return o.str();
}

static std::string json_escape(const std::string& in) {
  std::string out;
  out.reserve(in.size() + 2);
  for (char c : in) {
    switch (c) {
      case '"':  out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n";  break;
      case '\r': out += "\\r";  break;
      case '\t': out += "\\t";  break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) {
          char buf[8]; std::snprintf(buf, sizeof(buf), "\\u%04x", c);
          out += buf;
        } else out += c;
    }
  }
  return out;
}

static std::string num_or_null(double x) {
  if (std::isnan(x) || std::isinf(x)) return "null";
  std::ostringstream s; s << std::fixed << std::setprecision(4) << x;
  return s.str();
}

static std::string snapshot_to_json() {
  auto samples = g_snapshot.snapshot();
  const double t_now = g_sim_time_s.load();
  bool cft_on;          // afterglow window (stress*0.90)
  bool cft_applying;    // engine actualmente a modular vagal/simp.periferico
  double cft_elapsed_s = 0.0;
  double cft_dur_s     = 0.0;
  { std::lock_guard<std::mutex> lk(g_cft_mu);
    cft_on       = g_cft.was_active_recently(t_now);
    cft_applying = g_cft.is_currently_applying(t_now);
    if (g_cft.applied_at_s >= 0.0) {
      cft_elapsed_s = std::max(0.0, t_now - g_cft.applied_at_s);
      cft_dur_s     = g_cft.duration_s;
    }
  }
  bool   stress_decay_active = false;
  double stress_decay_sev_now = 0.0;
  double stress_decay_tau_s   = 0.0;
  double stress_decay_elapsed = 0.0;
  { std::lock_guard<std::mutex> lk(g_stress_decay_mu);
    stress_decay_active = g_stress_decay.active;
    if (stress_decay_active) {
      stress_decay_sev_now = g_stress_decay.sev_at(t_now);
      stress_decay_tau_s   = g_stress_decay.tau_s;
      stress_decay_elapsed = std::max(0.0, t_now - g_stress_decay.start_t_s);
    }
  }

  std::string current_patient;
  { std::lock_guard<std::mutex> lk(g_current_patient_mu);
    current_patient = g_current_patient_xml; }
  const std::string current_patient_display = patient_display_name(current_patient);

  // Estado MQTT para a UI
  const bool        mqtt_on         = mqtt_impl::g_active.load();
  const bool        mqtt_connecting = mqtt_impl::g_connecting.load();
  const std::string mqtt_broker     = mqtt_impl::g_broker;
  const int         mqtt_port       = mqtt_impl::g_port;
  const std::string mqtt_topic      = mqtt_impl::g_topic;
  bool              mqtt_ws         = false;
  std::string       mqtt_err;
  { std::lock_guard<std::mutex> lk(mqtt_impl::g_sock_mu);
    mqtt_ws  = mqtt_impl::g_use_ws;
    mqtt_err = mqtt_impl::g_last_error; }

  std::ostringstream o;
  o << "{\"sim_time\":"             << num_or_null(t_now)
    << ",\"speed_factor\":"         << num_or_null(g_speed_factor.load())
    << ",\"patient_xml\":\""        << json_escape(current_patient)         << "\""
    << ",\"patient_name\":\""       << json_escape(current_patient_display) << "\""
    << ",\"cft_afterglow_on\":"     << (cft_on ? "true" : "false")
    << ",\"cft_applying\":"         << (cft_applying ? "true" : "false")
    << ",\"cft_elapsed_s\":"        << num_or_null(cft_elapsed_s)
    << ",\"cft_dur_s\":"            << num_or_null(cft_dur_s)
    << ",\"stress_decay_active\":"  << (stress_decay_active ? "true" : "false")
    << ",\"stress_decay_sev_now\":" << num_or_null(stress_decay_sev_now)
    << ",\"stress_decay_tau_s\":"   << num_or_null(stress_decay_tau_s)
    << ",\"stress_decay_elapsed\":" << num_or_null(stress_decay_elapsed)
    << ",\"mqtt_active\":"           << (mqtt_on         ? "true" : "false")
    << ",\"mqtt_connecting\":"       << (mqtt_connecting ? "true" : "false")
    << ",\"mqtt_ws\":"               << (mqtt_ws         ? "true" : "false")
    << ",\"mqtt_broker\":\""         << json_escape(mqtt_broker) << "\""
    << ",\"mqtt_port\":"             << mqtt_port
    << ",\"mqtt_topic\":\""          << json_escape(mqtt_topic)  << "\""
    << ",\"mqtt_error\":\""          << json_escape(mqtt_err)    << "\""
    << ",\"diary_active\":"          << (g_diary_active.load() ? "true" : "false")
    << ",\"diary_cur_idx\":"         << g_diary_cur_idx.load()
    << ",\"samples\":[";

  for (size_t i = 0; i < samples.size(); ++i) {
    const auto& s = samples[i];
    if (i) o << ",";
    o << "{\"t\":"      << num_or_null(s.t_s)
      << ",\"ts_ms\":"  << s.wall_time_ms
      << ",\"actions\":\"" << json_escape(s.actions) << "\""
      << ",\"v\":{";
    bool first_v = true;
    for (const auto& kv : s.values) {
      if (!first_v) o << ",";
      o << "\"" << kv.first << "\":" << num_or_null(kv.second);
      first_v = false;
    }
    o << "}}";
  }
  o << "]}";
  return o.str();
}

// Serializa o buffer de alta resolucao do ECG (Lead III) para JSON.
// "t" = tempo de simulacao (s), "mv" = potencial eletrico (mV), amostrados
// a kEngineHz (NAO decimados) — ver EcgBuffer / kEcgBufferMaxSamples.
static std::string ecg_buffer_to_json() {
  auto samples = g_ecg_buffer.snapshot();
  std::ostringstream o;
  o << "{\"hz\":" << num_or_null(kEngineHz)
    << ",\"lead\":3"
    << ",\"samples\":[";
  for (size_t i = 0; i < samples.size(); ++i) {
    if (i) o << ",";
    o << "{\"t\":"  << num_or_null(samples[i].t_s)
      << ",\"mv\":" << num_or_null(samples[i].mV) << "}";
  }
  o << "]}";
  return o.str();
}


static bool parse_cmd_json(const std::string& body, std::string& out_cmd) {
  const std::string key = "\"cmd\"";
  auto p = body.find(key);
  if (p == std::string::npos) return false;
  p = body.find(':', p + key.size()); if (p == std::string::npos) return false;
  auto q1 = body.find('"', p + 1);    if (q1 == std::string::npos) return false;
  auto q2 = q1 + 1;
  std::string val;
  while (q2 < body.size() && body[q2] != '"') {
    if (body[q2] == '\\' && q2 + 1 < body.size()) {
      const char n = body[q2 + 1];
      switch (n) {
        case 'n':  val += '\n'; break;
        case 't':  val += '\t'; break;
        case 'r':  val += '\r'; break;
        case '"':  val += '"';  break;
        case '\\': val += '\\'; break;
        default:   val += n; break;
      }
      q2 += 2;
    } else {
      val += body[q2++];
    }
  }
  if (q2 >= body.size()) return false;
  out_cmd = val;
  return true;
}

static void handle_connection(int client_fd) {
  constexpr size_t kMaxReq = 16 * 1024;
  std::string req;
  req.reserve(2048);
  char buf[2048];
  ssize_t n;
  size_t headers_end = std::string::npos;
  size_t content_length = 0;

  while ((n = ::recv(client_fd, buf, sizeof(buf), 0)) > 0) {
    req.append(buf, buf + n);
    if (req.size() > kMaxReq) break;
    auto e = req.find("\r\n\r\n");
    if (e != std::string::npos) {
      headers_end = e + 4;
      const std::string hdrs = req.substr(0, e);
      auto cl = hdrs.find("Content-Length:");
      if (cl == std::string::npos) cl = hdrs.find("content-length:");
      if (cl != std::string::npos) {
        auto v = hdrs.find(':', cl) + 1;
        while (v < hdrs.size() && (hdrs[v] == ' ' || hdrs[v] == '\t')) ++v;
        content_length = static_cast<size_t>(std::atoll(hdrs.c_str() + v));
      }
      break;
    }
  }

  while (headers_end != std::string::npos
         && req.size() < headers_end + content_length
         && req.size() < kMaxReq)
  {
    n = ::recv(client_fd, buf, sizeof(buf), 0);
    if (n <= 0) break;
    req.append(buf, buf + n);
  }

  std::istringstream is(req);
  std::string method, path, version;
  is >> method >> path >> version;

  std::string resp;

  if (method == "GET" && (path == "/" || path == "/index.html")) {
    resp = http_response(200, "OK", "text/html; charset=utf-8", kHtmlPage);
  }
  else if (method == "GET" && path == "/data") {
    resp = http_response(200, "OK", "application/json; charset=utf-8",
                         snapshot_to_json());
  }
  else if (method == "GET" && path == "/ecg") {
    // Forma de onda do ECG (Lead III), alta resolucao — ver ecg_buffer_to_json().
    resp = http_response(200, "OK", "application/json; charset=utf-8",
                         ecg_buffer_to_json());
  }
  else if (method == "GET" && path == "/patients") {
    // Lista dinamica de XMLs disponiveis na pasta candidatas.

    auto files = list_patient_xmls();
    std::string current;
    { std::lock_guard<std::mutex> lk(g_current_patient_mu);
      current = g_current_patient_xml; }
    std::ostringstream o;
    o << "{\"current\":\"" << json_escape(current) << "\",\"patients\":[";
    for (size_t i = 0; i < files.size(); ++i) {
      if (i) o << ",";
      o << "{\"path\":\"" << json_escape(files[i]) << "\","
        <<  "\"name\":\"" << json_escape(patient_display_name(files[i])) << "\"}";
    }
    o << "]}";
    resp = http_response(200, "OK", "application/json; charset=utf-8", o.str());
  }
  else if (method == "POST" && path == "/cmd") {
    std::string body = (headers_end != std::string::npos)
                        ? req.substr(headers_end) : std::string{};
    std::string cmd;
    if (parse_cmd_json(body, cmd)) {
      g_queue.push(cmd);
      resp = http_response(200, "OK", "text/plain", "queued");
    } else {
      resp = http_response(400, "Bad Request", "text/plain",
                           "invalid JSON: expected {\"cmd\":\"...\"}");
    }
  }
  else if (method == "GET" && path == "/diary") {
    // Retorna estado actual do diary: blocos + active + cur_idx
    std::vector<DiaryEntry> snap;
    { std::lock_guard<std::mutex> lk(g_diary_mu); snap = g_diary_blocks; }
    std::ostringstream o;
    o << "{\"active\":"    << (g_diary_active.load() ? "true" : "false")
      << ",\"cur_idx\":"   << g_diary_cur_idx.load()
      << ",\"blocks\":"    << diary_to_json(snap)
      << "}";
    resp = http_response(200, "OK", "application/json; charset=utf-8", o.str());
  }
  else if (method == "POST" && path == "/diary") {
    // Corpo JSON: {"blocks":[{"start":"HH:MM","action":"...","param1":0,"param2":0}, ...]}
    // Substitui todos os blocos actuais. Nao inicia nem para o diary.
    std::string body = (headers_end != std::string::npos)
                        ? req.substr(headers_end) : std::string{};
    auto blocks = diary_from_json_array(body);
    diary_sort(blocks);
    { std::lock_guard<std::mutex> lk(g_diary_mu);
      g_diary_blocks = blocks; }
    g_diary_cur_idx.store(-2); // forca re-avaliacao no proximo tick se diary activo
    diary_save_json_file(blocks); // persiste automaticamente
    std::ostringstream o;
    o << "{\"ok\":true,\"count\":" << blocks.size() << "}";
    resp = http_response(200, "OK", "application/json; charset=utf-8", o.str());
  }
  else if (method == "GET" && path == "/params") {
    // Devolve catalogo completo de parametros fisiologicos
    std::ostringstream o;
    o << "{\"params\":[";
    for (size_t i = 0; i < kParamCatalog.size(); ++i) {
      const auto& p = kParamCatalog[i];
      if (i) o << ",";
      o << "{\"id\":\""     << json_escape(p.id)     << "\""
        << ",\"label\":\""  << json_escape(p.label)  << "\""
        << ",\"unit\":\""   << json_escape(p.unit)   << "\""
        << ",\"module\":\"" << json_escape(p.module) << "\"}";
    }
    o << "]}";
    resp = http_response(200, "OK", "application/json; charset=utf-8", o.str());
  }
  else if (method == "GET" && path == "/config") {
    // Devolve configuracao actual (mqtt_params)
    std::vector<std::string> mp;
    { std::lock_guard<std::mutex> lk(g_mqtt_params_mu); mp = g_mqtt_params; }
    std::ostringstream o;
    o << "{\"mqtt_params\":[";
    for (size_t i = 0; i < mp.size(); ++i) {
      if (i) o << ",";
      o << "\"" << json_escape(mp[i]) << "\"";
    }
    o << "]}";
    resp = http_response(200, "OK", "application/json; charset=utf-8", o.str());
  }
  else if (method == "POST" && path == "/config") {
    // Recebe {"mqtt_params":["hr","map",...]}
    std::string body = (headers_end != std::string::npos)
                        ? req.substr(headers_end) : std::string{};
    // Parseia array de strings para "mqtt_params"
    std::vector<std::string> new_params;
    auto bpos = body.find("\"mqtt_params\"");
    if (bpos != std::string::npos) {
      auto apos = body.find('[', bpos);
      auto epos = body.find(']', apos != std::string::npos ? apos : 0);
      if (apos != std::string::npos && epos != std::string::npos) {
        const std::string arr = body.substr(apos + 1, epos - apos - 1);
        std::string::size_type pos = 0;
        while (pos < arr.size()) {
          auto q1 = arr.find('"', pos);
          if (q1 == std::string::npos) break;
          auto q2 = arr.find('"', q1 + 1);
          if (q2 == std::string::npos) break;
          new_params.push_back(arr.substr(q1 + 1, q2 - q1 - 1));
          pos = q2 + 1;
        }
      }
    }
    { std::lock_guard<std::mutex> lk(g_mqtt_params_mu);
      g_mqtt_params = new_params; }
    std::ostringstream o;
    o << "{\"ok\":true,\"mqtt_params_count\":" << new_params.size() << "}";
    resp = http_response(200, "OK", "application/json; charset=utf-8", o.str());
  }
  else if (method == "OPTIONS") {
    std::ostringstream o;
    o << "HTTP/1.1 204 No Content\r\n"
      << "Access-Control-Allow-Origin: *\r\n"
      << "Access-Control-Allow-Methods: GET, POST, OPTIONS\r\n"
      << "Access-Control-Allow-Headers: Content-Type\r\n"
      << "Connection: close\r\n\r\n";
    resp = o.str();
  }
  else {
    resp = http_response(404, "Not Found", "text/plain", "not found");
  }

  ::send(client_fd, resp.data(), resp.size(), 0);
  ::close(client_fd);
}

static void http_server_thread() {
  int srv = ::socket(AF_INET, SOCK_STREAM, 0);
  if (srv < 0) {
    std::cerr << "[ERRO] socket() falhou\n";
    return;
  }
  int yes = 1;
  ::setsockopt(srv, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(yes));

  sockaddr_in addr{};
  addr.sin_family      = AF_INET;
  addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);   // 127.0.0.1 (so local)
  addr.sin_port        = htons(kHttpPort);

  if (::bind(srv, reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) < 0) {
    std::cerr << "[ERRO] bind() falhou na porta " << kHttpPort
              << " (ja em uso?). UI nao iniciada.\n";
    ::close(srv);
    return;
  }
  if (::listen(srv, 8) < 0) {
    std::cerr << "[ERRO] listen() falhou\n";
    ::close(srv);
    return;
  }

  std::cout << "[UI] HTTP server em http://localhost:" << kHttpPort << "\n"
            << "[UI] abra essa URL no browser.\n";

  while (g_running.load()) {
    timeval tv{}; tv.tv_sec = 1; tv.tv_usec = 0;
    fd_set rd; FD_ZERO(&rd); FD_SET(srv, &rd);
    int sel = ::select(srv + 1, &rd, nullptr, nullptr, &tv);
    if (sel < 0)  break;
    if (sel == 0) continue;
    sockaddr_in cli{}; socklen_t clen = sizeof(cli);
    int cfd = ::accept(srv, reinterpret_cast<sockaddr*>(&cli), &clen);
    if (cfd < 0) continue;
    handle_connection(cfd);
  }
  ::close(srv);
}

}  

//  main
// Cria + inicializa um novo engine BioGears com o paciente indicado.

static std::unique_ptr<PhysiologyEngine> make_engine_with_patient(const std::string& patient_xml) {
  auto eng = CreateBioGearsEngine(kLogPath);
  if (!eng) {
    std::cerr << "[ERRO] CreateBioGearsEngine devolveu nullptr.\n";
    return nullptr;
  }
  eng->GetLogger()->Info("HowTo-RealTimeUI engine create...");
  std::cout << "Carregando paciente: " << patient_xml
            << "  (estabilizacao 30-60 s)...\n" << std::flush;
  if (!eng->InitializeEngine(patient_xml)) {
    eng->GetLogger()->Error("Falha ao inicializar paciente.");
    std::cerr << "[ERRO] InitializeEngine falhou (\"" << patient_xml
              << "\"). Verifique " << kLogPath << "\n";
    return nullptr;
  }
  return eng;
}

int main()
{
  std::cout << "\n=========================================================\n"
            <<   " BioGears RealTimeUI Driver (HowTo-RealTimeUI, V5)\n"
            <<   "   - engine @ " << kEngineHz << " Hz (dt=" << 1.0 / kEngineHz << " s)\n"
            <<   "   - wall-clock speed factor: " << g_speed_factor.load() << "x (changeable: 1, 2, 5)\n"
            <<   "=========================================================\n";

  // 1. Engine init
  apply_exercise_params_cfg_for_patient(kDefaultPatientXml);
  std::unique_ptr<PhysiologyEngine> bg = make_engine_with_patient(kDefaultPatientXml);
  if (!bg) return 1;
  { std::lock_guard<std::mutex> lk(g_current_patient_mu);
    g_current_patient_xml = kDefaultPatientXml; }

  // 2. CSV (com Epinephrine)
  std::ofstream csv(kCsvPath);
  if (!csv) {
    std::cerr << "[ERRO] nao consegui abrir " << kCsvPath << "\n";
    return 1;
  }
  // CSV header dinamico baseado no catalogo
  csv << "timestamp";
  for (const auto& pd : kParamCatalog) {
    csv << "," << pd.id << "_" << pd.unit;
  }
  csv << ",active_actions\n";
  csv << std::fixed << std::setprecision(4);

  // 3. Threads auxiliares: stdin + HTTP UI
  std::thread input_thread([]() {
    std::string line;
    while (g_running.load(std::memory_order_relaxed)) {
      if (!std::getline(std::cin, line)) { g_queue.push("quit"); return; }
      g_queue.push(line);
    }
  });
  input_thread.detach();

  std::thread ui_thread(http_server_thread);

  // 4. Loop principal
  // Ancora o eixo temporal ao instante de arranque da simulacao.
  {
    const auto now = std::chrono::system_clock::now();
    g_wall_epoch_origin_ms.store(static_cast<long long>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()).count()));
  }

  ActiveActions actions;
  print_help();
  std::cout << "> " << std::flush;

  const double dt_s      = 1.0 / kEngineHz;
  // dt_chrono_base = duracao wall-clock de um tick em tempo real (1x).

  const auto   dt_chrono_base = std::chrono::microseconds(
                                   static_cast<long long>(1e6 / kEngineHz));
  auto  next_tick        = clk::now() + dt_chrono_base;
  long long tick         = 0;
  long long late_warn    = 0;
  bool quit              = false;

  // (cv, en, patient_weight_kg foram removidos: o catalogo kParamCatalog
  //  acede ao engine diretamente via PhysiologyEngine& em cada getter lambda.)


  while (!quit) {
    // 4.1 drena fila de comandos
    std::string cmd;
    while (g_queue.try_pop(cmd)) {
      if (dispatch(cmd, *bg, actions)) { quit = true; break; }
      std::cout << "> " << std::flush;
    }
    if (quit) break;

    // 4.1b  troca de paciente pedida pelo user

    if (g_reset_pending.load()) {
      std::string new_xml;
      { std::lock_guard<std::mutex> lk(g_pending_patient_mu);
        new_xml = g_pending_patient_xml;
        g_pending_patient_xml.clear(); }
      g_reset_pending.store(false);

      if (!new_xml.empty()) {
        std::cout << "\n[reset] reiniciando engine com paciente: " << new_xml
                  << "\n" << std::flush;

        // (a) Zera estados auxiliares ANTES de destruir o engine

        actions = ActiveActions{};
        { std::lock_guard<std::mutex> lk(g_sched_mu); g_scheduled.clear(); }
        { std::lock_guard<std::mutex> lk(g_cft_mu);    g_cft = CFTTracker{}; }
        { std::lock_guard<std::mutex> lk(g_stress_decay_mu);
          g_stress_decay.disarm(); }
          
        g_current_stress_sev.store(0.0);
        // Reinicia o diary para o inicio do ciclo (mantendo blocos e estado activo)
        g_diary_cur_idx.store(-2);
        g_snapshot.clear();
        g_ecg_buffer.clear();
        g_sim_time_s.store(0.0);
        tick = 0;
        late_warn = 0;
        // Reinicia a ancora temporal: novo t=0 corresponde ao instante actual.
        {
          const auto now = std::chrono::system_clock::now();
          g_wall_epoch_origin_ms.store(static_cast<long long>(
            std::chrono::duration_cast<std::chrono::milliseconds>(
              now.time_since_epoch()).count()));
        }

        // (b) Destroi o engine atual e cria um novo com o paciente pedido.

        const std::string previous_xml = [&] {
          std::lock_guard<std::mutex> lk(g_current_patient_mu);
          return g_current_patient_xml;
        }();
        bg.reset();
        apply_exercise_params_cfg_for_patient(new_xml);
        auto new_bg = make_engine_with_patient(new_xml);

        if (!new_bg) {
          std::cerr << "[ERRO] falha ao carregar \"" << new_xml
                    << "\". Tentando reverter para \"" << previous_xml << "\"...\n";
          apply_exercise_params_cfg_for_patient(previous_xml);
          new_bg = make_engine_with_patient(previous_xml);
          if (!new_bg) {
            std::cerr << "[FATAL] fallback tambem falhou. Encerrando simulacao.\n";
            quit = true;
            break;
          }
        } else {
          std::lock_guard<std::mutex> lk(g_current_patient_mu);
          g_current_patient_xml = new_xml;
        }
        bg = std::move(new_bg);

        // (c) Ponteiros cv/en foram removidos; catalogo acede engine diretamente.

        // (d) Re-baseline do wall-clock: comeca a contar a partir de "agora".
        next_tick = clk::now() + dt_chrono_base;

        // (e) Marca a quebra no CSV + reescreve o header (util para pos-processamento)
        {
          std::lock_guard<std::mutex> lk(g_current_patient_mu);
          const auto now_wall = std::chrono::system_clock::now();
          const long long ts_ms = static_cast<long long>(
            std::chrono::duration_cast<std::chrono::milliseconds>(
              now_wall.time_since_epoch()).count());
          csv << "# ---- patient reset: " << g_current_patient_xml
              << " @ " << format_wall_timestamp(ts_ms) << " (sim_t->0) ----\n";
          // Re-escreve header apos reset
          csv << "timestamp";
          for (const auto& pd : kParamCatalog) csv << "," << pd.id << "_" << pd.unit;
          csv << ",active_actions\n";
        }
        csv.flush();

        std::cout << "[reset] engine pronto. sim_time=0 s, paciente=";
        { std::lock_guard<std::mutex> lk(g_current_patient_mu);
          std::cout << g_current_patient_xml; }
        std::cout << "\n> " << std::flush;
        continue;   // pula o resto deste tick: comeca limpo no proximo
      }
    }

    // 4.2 varre scheduled-stops cujo prazo expirou
    {
      const double t_now = g_sim_time_s.load();
      std::vector<std::string> to_stop;
      {
        std::lock_guard<std::mutex> lk(g_sched_mu);
        for (auto it = g_scheduled.begin(); it != g_scheduled.end(); ) {
          if (it->stop_at_s <= t_now) {
            to_stop.push_back(it->what);
            it = g_scheduled.erase(it);
          } else ++it;
        }
      }
      for (const auto& w : to_stop) {
        if      (w == "exercise") { clear_exercise(*bg); actions.exercise.clear(); }
        else if (w == "stress")   {
          // Auto-stop do "stress sev for dur": em vez de remover abruptamente,
          // armamos o decay exponencial. Mesma logica do "stop stress" manual.
          const double sev_now = g_current_stress_sev.load();
          if (sev_now > kStressFloor) {
            start_stress_decay(sev_now);
            std::ostringstream s;
            s << std::fixed << std::setprecision(3) << sev_now << " (decay tau="
              << kStressDecayTau_s << "s)";
            actions.stress = s.str();
            std::cout << "\n  [auto-stop] stress decay armed (sev_start="
                      << sev_now << ", tau=" << kStressDecayTau_s
                      << " s, t=" << t_now << " s)\n> " << std::flush;
          } else {
            clear_stress(*bg);
            actions.stress.clear();
          }
        }
        else if (w == "coldface") { clear_coldface(*bg); actions.coldface.clear();
                                    { std::lock_guard<std::mutex> lk(g_cft_mu);
                                      g_cft.on_clear(t_now); } }
        else if (w == "sleep")    { clear_sleep(*bg); actions.sleep.clear(); }
        // "stress" ja imprimiu a sua propria mensagem dentro do branch (decay
        // armed); para os restantes mantemos a notificacao generica.
        if (w != "stress") {
          std::cout << "\n  [auto-stop] " << w << " (t=" << t_now << " s)\n> " << std::flush;
        }
      }
    }

    // 4.2b avanca o decay do stress 

    {
      std::string new_stress_label = actions.stress;
      if (tick_stress_decay(*bg, new_stress_label)) {
        actions.stress = new_stress_label;
      }
    }

    // 4.2c  diary: verifica transicao de bloco e ciclo de 24h
    if (g_diary_active.load()) {
      const double t_now_s   = g_sim_time_s.load();
      // sim_day_s: posicao dentro do dia actual (ciclo de 86400 s)
      const double sim_day_s = std::fmod(t_now_s, 86400.0);

      // Obtemos copia dos blocos (lock breve)
      std::vector<DiaryEntry> diary_snap;
      { std::lock_guard<std::mutex> lk(g_diary_mu);
        diary_snap = g_diary_blocks; }

      const int new_idx = diary_find_idx(sim_day_s, diary_snap);
      const int cur_idx = g_diary_cur_idx.load();

      if (new_idx != cur_idx) {
        // Transicao de bloco ou inicio de novo ciclo de 24h
        std::cout << "\n[diary] " << diary_hms(sim_day_s)
                  << "  bloco " << cur_idx << " -> " << new_idx << "\n> " << std::flush;
        diary_apply_block(new_idx, diary_snap, *bg, actions);
        g_diary_cur_idx.store(new_idx);
      }
    }

    // 4.3 advance engine
    //
    // Protecao contra NaN no circuito de energia do BioGears:
    //
    //   O BioGears usa um solver MNA (Modified Nodal Analysis) cujo vector de
    //   solucao persiste entre ticks. Se durante um Generic exercise em intensity=0
    //   o solver acumula NaN num no (e.g., divisao implícita no modelo de O2
    //   esqueletico), esse NaN mantem-se como "condicao inicial" nos ticks
    //   seguintes. Quando um exercise de tipo diferente (Running, Cycling) tenta
    //   usar esses nos, a excecao e lancada dentro de AdvanceModelTime — fora do
    //   try/catch do dispatcher, portanto terminava o processo.
    //
    //   Estrategia de recuperacao:
    //     1. catch(...) — apanha CommonDataModelException e qualquer outro tipo
    //        (inclui subclasses nao derivadas de std::exception).
    //     2. Envia zero para TODOS os tipos de exercicio para limpar os termos
    //        fonte do sistema MNA (b em Ax=b); se b e limpo, o solver converge.
    //     3. Tenta AdvanceModelTime novamente.
    //     4. Se ainda falhar, o estado do engine e irrecuperavel: forca reset
    //        do paciente (mesma logica do comando "patient <xml>").
    {
      bool advance_ok = false;
      try {
        bg->AdvanceModelTime(dt_s, TimeUnit::s);
        advance_ok = true;
      } catch (...) {
        bg->GetLogger()->Warning(
          "[engine] NaN/exception em AdvanceModelTime -- limpando exercicio e retentando");
        std::cout << "\n  [aviso] engine: NaN detectado -- exercicio cancelado, retentando\n> "
                  << std::flush;

        // Limpeza total: envia zero para todos os tipos de exercicio.
        // Isto limpa os termos fonte do sistema MNA independentemente de qual
        // tipo estava activo, pois o vector de solucao pode ter NaN em qualquer no.
        try {
          { SEExercise::SEGeneric g; g.Intensity.SetValue(0.0);
            SEExercise ex { g }; bg->ProcessAction(ex); }
          { SEExercise::SERunning r;
            r.SpeedRun.SetValue(0.0, LengthPerTimeUnit::m_Per_s);
            r.InclineRun.SetValue(0.0);
            SEExercise ex { r }; bg->ProcessAction(ex); }
          { SEExercise::SECycling c;
            c.PowerCycle.SetValue(0.0, PowerUnit::W);
            c.CadenceCycle.SetValue(0.0, FrequencyUnit::Per_min);
            SEExercise ex { c }; bg->ProcessAction(ex); }
        } catch (...) {}

        g_active_exercise_type = ExerciseType::None;
        actions.exercise.clear();
        cancel_scheduled("exercise");

        // Segunda tentativa com circuito limpo.
        try {
          bg->AdvanceModelTime(dt_s, TimeUnit::s);
          advance_ok = true;
        } catch (...) {
          // Estado do engine irrecuperavel — forca reset do paciente actual.
          std::cerr << "\n[FATAL] engine nao recuperou do NaN -- reiniciando paciente\n"
                    << std::flush;
          bg->GetLogger()->Error(
            "[engine] irrecuperavel -- a reiniciar com o mesmo paciente");
          std::string cur_xml;
          { std::lock_guard<std::mutex> lk(g_current_patient_mu);
            cur_xml = g_current_patient_xml; }
          { std::lock_guard<std::mutex> lk(g_pending_patient_mu);
            g_pending_patient_xml = cur_xml; }
          g_reset_pending.store(true);
          // Incrementa tick para nao perder o wall-clock, mas nao amostra.
          ++tick;
          g_sim_time_s.store(static_cast<double>(tick) * dt_s);
          next_tick += std::chrono::microseconds(
            static_cast<long long>(1e6 / (kEngineHz * g_speed_factor.load())));
          continue;
        }
      }
      (void)advance_ok;  // suprime warning de variavel nao usada
    }
    ++tick;
    g_sim_time_s.store(static_cast<double>(tick) * dt_s);

    // 4.3b ECG (Lead III): amostragem a full engine rate (kEngineHz=50Hz),
    // deliberadamente FORA do bloco decimado abaixo — ver comentario junto
    // a kEcgBufferMaxSamples (criterio de Nyquist para o complexo QRS).
    {
      EcgSample es;
      es.t_s = static_cast<double>(tick) * dt_s;
      es.mV  = get_ecg_lead3_mV(*bg);
      g_ecg_buffer.push(es);
    }

    // 4.4 amostragem (1 Hz @ kCsvDecimateTicks=50 ticks / 50 Hz)
    if ((tick % kCsvDecimateTicks) == 0) {
      Sample s;
      s.t_s   = static_cast<double>(tick) * dt_s;
      // Timestamp ancorado: T0_real + t_simulacao * 1000 ms.
      // A velocidade de simulação (1x/2x/5x) afeta t_s diretamente.
      s.wall_time_ms = g_wall_epoch_origin_ms.load()
                       + static_cast<long long>(s.t_s * 1000.0);
      s.actions = actions.to_csv();

      // Coleta todos os parametros do catalogo
      for (const auto& pd : kParamCatalog) {
        s.values[pd.id] = pd.get(*bg);
      }

      // CSV: header dinamico na primeira amostra (apos reset o CSV já tem header)
      // Linha: timestamp, todas as colunas do catalogo, actions
      csv << format_wall_timestamp(s.wall_time_ms);
      for (const auto& pd : kParamCatalog) {
        csv << "," << num_or_null(s.values.count(pd.id) ? s.values.at(pd.id)
                                                        : std::numeric_limits<double>::quiet_NaN());
      }
      csv << "," << s.actions << "\n";
      csv.flush();

      // MQTT: publica JSON com os parametros selecionados pelo user
      if (mqtt_impl::g_active.load()) {
        std::vector<std::string> mp;
        { std::lock_guard<std::mutex> lk(g_mqtt_params_mu); mp = g_mqtt_params; }
        if (!mp.empty()) {
          std::ostringstream mqtt_csv;
          mqtt_csv << std::fixed << std::setprecision(4);
          mqtt_csv << format_wall_timestamp(s.wall_time_ms);
          for (const auto& pid : mp) {
            auto it = s.values.find(pid);
            mqtt_csv << "," << (it != s.values.end() ? num_or_null(it->second) : "null");
          }
          mqtt_impl::publish(mqtt_csv.str());
        }
      }

      g_snapshot.push(s);
    }

    // 4.5 sleep_until, evita drift acumulado.
    const double speed = g_speed_factor.load();
    const auto   tick_us = std::chrono::microseconds(
                              static_cast<long long>(1e6 / (kEngineHz * speed)));
    const auto now = clk::now();
    if (now < next_tick) {
      std::this_thread::sleep_until(next_tick);
    } else {
      if ((++late_warn) % 500 == 0) {
        std::cerr << "[warn] engine atrasado em "
                  << std::chrono::duration<double>(now - next_tick).count()
                  << " s (cumulativo: " << late_warn << " ticks, speed="
                  << speed << "x)\n";
      }
      next_tick = now;
    }
    next_tick += tick_us;
  }

  // 5. Cleanup
  g_running.store(false);
  mqtt_impl::stop();   // DISCONNECT gracioso ao broker (se activo)
  if (ui_thread.joinable()) ui_thread.join();
  csv.close();

  std::cout << "\nSimulacao encerrada.\n"
            << "  CSV : " << kCsvPath << "\n"
            << "  Log : " << kLogPath << "\n";
  return 0;
}
