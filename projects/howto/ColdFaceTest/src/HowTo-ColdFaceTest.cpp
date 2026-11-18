/**************************************************************************************
 *  HowTo-ColdFaceTest.cpp
 *  ----------------------
 *  Protocolo D — Cold Face Test (CFT) com stress repetido, BioGears V5.
 *
 *  Sequência (StandardMale, single subject, fully agentic):
 *
 *      F1  Baseline                  0   - 600  s   (600 s)
 *      F2  CFT em baseline           600 - 780  s   (90 apply + 90 recov)
 *      F3  Stress ON #1              780 - 1080 s   (300 s, sev=0.80, step ON)
 *      F4  Stress ramp-down #1       1080- 1200 s   (120 s, sev 0.80 -> 0 LINEAR)
 *      F5  CFT após stress #1        1200- 1380 s   (90 apply + 90 recov)
 *      F6  Stress ON #2              1380- 1680 s   (300 s, sev=0.80, step ON)
 *      F7  Stress ramp-down #2       1680- 1800 s   (120 s, sev 0.80 -> 0 LINEAR)
 *      F8  CFT após stress #2        1800- 1980 s   (90 apply + 90 recov)
 *      F9  Post-baseline (washout)   1980- 2280 s   (300 s)
 *
 *  NOTA fisiológica: o desligamento do stress não é abrupto na vida real
 *  (clearance de NE/Epi tem τ ~ 90 s, Esler 1990). F4 e F7 fazem um ramp
 *  linear da severity entre 0.80 e 0 ao longo dos 120 s, atualizando o
 *  SEAcuteStress a 1 Hz. Para usar decaimento exponencial (mais fiel),
 *  ver run_stress_ramp_exp() comentado no fim do arquivo.
 *
 *  Total ≈ 38 min. O desenho permite estudar habituação/sensibilização do reflexo
 *  trigemino-cardíaco a stress repetido, comparando dHR/dSBP de CFT1 (baseline puro)
 *  contra CFT2 e CFT3 (em estados pós-stress parcialmente recuperados).
 *
 *  Outputs (na CWD do executável):
 *      HowTo-ColdFaceTest.csv          timeseries a 1 Hz
 *      HowTo-ColdFaceTest.log          log do engine
 *      HowTo-ColdFaceTest.summary.txt  sumário das métricas-chave
 *
 *  Referências:
 *      Heath ME, Downey JA. Clin Sci 1990; 78:139-147.
 *      Khurana RK, Wu R. Clin Auton Res 2006; 16:202-206.
 *      Brown CM, Sanya EO, Hilz MJ. Clin Auton Res 2003; 13:351-358.
 *      Schmidt H et al. Crit Care Med 2005; 33:1994-2002 (habituação simpática).
 *=====================================================================*/

#include <biogears/cdm/engine/PhysiologyEngineTrack.h>
#include <biogears/cdm/patient/SEPatient.h>
#include <biogears/cdm/patient/actions/SEAcuteStress.h>
#include <biogears/cdm/patient/actions/SEColdFaceTest.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <biogears/cdm/substance/SESubstance.h>
#include <biogears/cdm/substance/SESubstanceManager.h>
#include <biogears/cdm/system/physiology/SEBloodChemistrySystem.h>
#include <biogears/cdm/system/physiology/SECardiovascularSystem.h>
#include <biogears/cdm/system/physiology/SEEnergySystem.h>
#include <biogears/cdm/system/physiology/SERespiratorySystem.h>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/string/manipulation.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace biogears;

//  Configuração temporal — Protocolo D

namespace {
constexpr double kEngineHz = 50.0;             // 1 dt = 0.02 s
constexpr double kSampleHz = 1.0;              // CSV @ 1 Hz

// Janelas (s)
constexpr double kPhase1_Baseline_s     = 600.0;   // F1
constexpr double kCFT_apply_s           =  90.0;   // duração comum apply
constexpr double kCFT_recov_s           =  90.0;   // duração comum recovery
constexpr double kPhase3_StressOn1_s    = 300.0;   // F3
constexpr double kPhase4_StressOff1_s   = 60.0;   // F4
constexpr double kPhase6_StressOn2_s    = 300.0;   // F6
constexpr double kPhase7_StressOff2_s   = 120.0;   // F7
constexpr double kPhase9_PostBase_s     = 300.0;   // F9

// Severidades
constexpr double kAcuteStressSeverity   = 0.80;
constexpr double kCFTSeverity           = 0.5;

// Janela de tail-mean para definir baseline imediatamente antes de cada CFT
constexpr int    kTailWindow_s          = 60;     // últimos 60 s @ 1 Hz

// ---- Cinética de clearance do stress (propriedade do paciente, não do bloco)
// Linear:   sev(t) = max(sev_floor, sev_start - kStressDecayRate_per_s * t)
//   0.0067/s -> chega ao floor em ~118 s a partir de 0.80
// Exponencial: sev(t) = sev_start * exp(-t / kStressTau_s)
//   tau ~ 90 s consistente com clearance plasmática de NE (Esler 1990)
constexpr double kStressDecayRate_per_s = 0.0067; // unidades de severity por segundo
constexpr double kStressTau_s           = 90.0;
constexpr double kStressFloor           = 0.01;   // > IsZero() do BioGears


// Paths
static const char* const PATIENT_XML    = "patient/Male_25_Normal.xml";
constexpr const char* kCsvPath  = "HowTo-ColdFaceTest2.csv";
constexpr const char* kLogPath  = "HowTo-ColdFaceTest2.log";
constexpr const char* kSummary  = "HowTo-ColdFaceTest2.summary.txt";

} // namespace

//  Agregador de estatísticas com buffer p/ tail-mean
struct PhaseStats {
  double HR_min  = 1e9, HR_max = 0.0, HR_sum = 0.0;
  double SBP_min = 1e9, SBP_max = 0.0, SBP_sum = 0.0;
  double MAP_sum = 0.0;
  int    N = 0;
  double t_HRmin_rel_s = 0.0;
  std::vector<double> HR_samples;
  std::vector<double> SBP_samples;
  std::vector<double> MAP_samples;
};

static void update(PhaseStats& s, double hr, double sbp, double map, double t_rel) {
  if (hr  < s.HR_min)  { s.HR_min  = hr;  s.t_HRmin_rel_s = t_rel; }
  if (hr  > s.HR_max)    s.HR_max  = hr;
  if (sbp < s.SBP_min)   s.SBP_min = sbp;
  if (sbp > s.SBP_max)   s.SBP_max = sbp;
  s.HR_sum  += hr;
  s.SBP_sum += sbp;
  s.MAP_sum += map;
  ++s.N;
  s.HR_samples.push_back(hr);
  s.SBP_samples.push_back(sbp);
  s.MAP_samples.push_back(map);
}


static double tail_mean(const std::vector<double>& v, std::size_t k) {
  if (v.empty()) return 0.0;
  const std::size_t start = (v.size() > k) ? v.size() - k : 0;
  double sum = 0.0;
  for (std::size_t i = start; i < v.size(); ++i) sum += v[i];
  return sum / static_cast<double>(v.size() - start);
}


static void run_block(PhysiologyEngine& bg,
                      double duration_s,
                      const char* phase_label,
                      PhaseStats& st)
{
  const double dt = 1.0 / kEngineHz;
  const int    steps = static_cast<int>(duration_s / dt);
  const auto&  cv = *bg.GetCardiovascularSystem();
  const int    sample_every = static_cast<int>(kEngineHz / kSampleHz);

  bg.GetLogger()->Info(asprintf("=== %s : %.0f s ===", phase_label, duration_s));

  for (int i = 0; i < steps; ++i) {
    bg.AdvanceModelTime(dt, TimeUnit::s);
    if ((i % sample_every) != 0) continue;
    const double hr  = cv.GetHeartRate(FrequencyUnit::Per_min);
    const double sbp = cv.GetSystolicArterialPressure(PressureUnit::mmHg);
    const double map = cv.GetMeanArterialPressure(PressureUnit::mmHg);
    const double t_rel = static_cast<double>(i) * dt;
    update(st, hr, sbp, map, t_rel);
  }
}

//  Helpers compactos para acionar/desligar CFT e Stress
static void apply_cft(PhysiologyEngine& bg, double sev, double dur_s) {
  SEColdFaceTest cft;
  cft.GetSeverity().SetValue(sev);
  cft.GetDuration().SetValue(dur_s, TimeUnit::s);
  bg.ProcessAction(cft);
}
static void clear_cft(PhysiologyEngine& bg) {
  SEColdFaceTest off;
  off.GetSeverity().SetValue(0.0);
  bg.ProcessAction(off);
}
static void apply_stress(PhysiologyEngine& bg, double sev) {
  SEAcuteStress s;
  s.GetSeverity().SetValue(sev);
  bg.ProcessAction(s);
}
static void clear_stress(PhysiologyEngine& bg) {
  SEAcuteStress off;
  off.GetSeverity().SetValue(0.0);
  bg.ProcessAction(off);
}

// -----------------------------------------------------------------------
//  DECAIMENTO LINEAR com taxa fixa de clearance (propriedade do paciente).
//
//      sev(t) = max(kStressFloor, sev_start - rate_per_s * t)
//
//  A duração 'duration_s' apenas determina ATÉ ONDE no decaimento o bloco
//  vai amostrar — NÃO força o stress a chegar em zero. Isto é
//  fisiologicamente correto: a velocidade de clearance plasmática é
//  determinada pela cinética do paciente (Esler 1990), não pelo desenho
//  experimental.
//
//  Exemplos com sev_start=0.80, rate=0.0067/s:
//      duration =   1 s  ->  sev_final ≈ 0.793  (paciente ainda estressado)
//      duration =  60 s  ->  sev_final ≈ 0.398
//      duration = 118 s  ->  sev_final ≈ kStressFloor (clearance completa)
//      duration = 300 s  ->  sev_final ≈ kStressFloor (estabilizado)
//
//  Floor em kStressFloor=0.01 (> tolerância IsZero do SEScalar0To1)
//  evita que o engine chame RemoveAcuteStress() durante o ramp.
// -----------------------------------------------------------------------
static double run_stress_decay_linear(PhysiologyEngine& bg,
                                      double sev_start,
                                      double rate_per_s,
                                      double duration_s,
                                      const char* phase_label,
                                      PhaseStats& st)
{
  const double dt           = 1.0 / kEngineHz;
  const int    steps        = static_cast<int>(duration_s / dt);
  const auto&  cv           = *bg.GetCardiovascularSystem();
  const int    sample_every = static_cast<int>(kEngineHz / kSampleHz);
  const int    update_every = static_cast<int>(kEngineHz);   // 1 Hz update

  bg.GetLogger()->Info(asprintf(
      "=== %s : %.0f s (linear decay sev_start=%.2f, rate=%.4f /s) ===",
      phase_label, duration_s, sev_start, rate_per_s));

  double sev_last = sev_start;
  for (int i = 0; i < steps; ++i) {
    if ((i % update_every) == 0) {
      const double t = static_cast<double>(i) * dt;
      double sev = sev_start - rate_per_s * t;
      // Clamp defensivo (NaN-safe + bounds [floor, 1])
      if (!(sev > 0.0 && sev <= 1.0)) sev = kStressFloor;
      if (sev < kStressFloor)         sev = kStressFloor;
      apply_stress(bg, sev);
      sev_last = sev;
    }
    bg.AdvanceModelTime(dt, TimeUnit::s);
    if ((i % sample_every) != 0) continue;
    const double hr  = cv.GetHeartRate(FrequencyUnit::Per_min);
    const double sbp = cv.GetSystolicArterialPressure(PressureUnit::mmHg);
    const double map = cv.GetMeanArterialPressure(PressureUnit::mmHg);
    const double t_rel = static_cast<double>(i) * dt;
    update(st, hr, sbp, map, t_rel);
  }
  return sev_last;   // severity efetiva no fim do bloco
}

// -----------------------------------------------------------------------
//  DECAIMENTO EXPONENCIAL com tau fixo (propriedade do paciente).
//
//      sev(t) = max(kStressFloor, sev_start * exp(-t / tau_s))
//
//  Igual ao linear: a duração só determina até onde amostramos o
//  decaimento. tau ~ 90 s ≈ half-life × ln(2) ≈ 60 s × 1.5 da NE
//  plasmática (Esler 1990; Goldstein, Adv Pharmacol 1998).
//
//  Exemplos com sev_start=0.80, tau=90 s:
//      duration =   1 s  ->  sev ≈ 0.791
//      duration =  60 s  ->  sev ≈ 0.41   (half-life)
//      duration = 120 s  ->  sev ≈ 0.21
//      duration = 300 s  ->  sev ≈ 0.029  (~95% washout)
//
//  Update a 0.2 Hz (a cada 5 s) — Marshall/UnMarshall a 1 Hz é caro
//  e a dinâmica simpática tem constantes >> 1 s.
// -----------------------------------------------------------------------
static double run_stress_decay_exp(PhysiologyEngine& bg,
                                   double sev_start,
                                   double tau_s,
                                   double duration_s,
                                   const char* phase_label,
                                   PhaseStats& st)
{
  constexpr double kUpdateDt_s = 5.0;        // re-apply severity @ 0.2 Hz

  const double dt           = 1.0 / kEngineHz;
  const int    steps        = static_cast<int>(duration_s / dt);
  const int    sample_every = static_cast<int>(kEngineHz / kSampleHz);
  const int    update_every = static_cast<int>(kEngineHz * kUpdateDt_s);
  const auto&  cv           = *bg.GetCardiovascularSystem();

  bg.GetLogger()->Info(asprintf(
      "=== %s : %.0f s (exp decay sev_start=%.2f, tau=%.0f s) ===",
      phase_label, duration_s, sev_start, tau_s));

  double sev_last = sev_start;
  for (int i = 0; i < steps; ++i) {
    if ((i % update_every) == 0) {
      const double t = static_cast<double>(i) * dt;
      double sev = sev_start * std::exp(-t / tau_s);
      // Clamp defensivo (NaN-safe + bounds [floor, 1])
      if (!(sev > 0.0 && sev <= 1.0)) sev = kStressFloor;
      if (sev < kStressFloor)         sev = kStressFloor;
      apply_stress(bg, sev);
      sev_last = sev;
    }
    bg.AdvanceModelTime(dt, TimeUnit::s);
    if ((i % sample_every) != 0) continue;
    const double hr  = cv.GetHeartRate(FrequencyUnit::Per_min);
    const double sbp = cv.GetSystolicArterialPressure(PressureUnit::mmHg);
    const double map = cv.GetMeanArterialPressure(PressureUnit::mmHg);
    const double t_rel = static_cast<double>(i) * dt;
    update(st, hr, sbp, map, t_rel);
  }
  return sev_last;
}

//  Data requests
static void configure_data_requests(PhysiologyEngine& bg)
{
  auto& drm = bg.GetEngineTrack()->GetDataRequestManager();
  drm.SetSamplesPerSecond(kSampleHz);
  drm.CreatePhysiologyDataRequest().Set("HeartRate",                 FrequencyUnit::Per_min);
  drm.CreatePhysiologyDataRequest().Set("SystolicArterialPressure",  PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("DiastolicArterialPressure", PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("MeanArterialPressure",      PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("HeartStrokeVolume",         VolumeUnit::mL);
  drm.CreatePhysiologyDataRequest().Set("CardiacOutput",             VolumePerTimeUnit::L_Per_min);
  drm.CreatePhysiologyDataRequest().Set("SystemicVascularResistance", FlowResistanceUnit::mmHg_s_Per_mL);
  drm.CreatePhysiologyDataRequest().Set("OxygenSaturation");
  drm.CreatePhysiologyDataRequest().Set("RespirationRate",           FrequencyUnit::Per_min);
  drm.CreatePhysiologyDataRequest().Set("CoreTemperature",           TemperatureUnit::C);
  drm.CreatePhysiologyDataRequest().Set("OxygenConsumptionRate",       VolumePerTimeUnit::mL_Per_min);
  drm.CreatePhysiologyDataRequest().Set("CarbonDioxideProductionRate", VolumePerTimeUnit::mL_Per_min);
  drm.CreatePhysiologyDataRequest().Set("RespiratoryExchangeRatio");
  drm.CreatePhysiologyDataRequest().Set("TidalVolume",                 VolumeUnit::mL);
  drm.CreatePhysiologyDataRequest().Set("TotalPulmonaryVentilation",   VolumePerTimeUnit::L_Per_min);
  drm.CreatePhysiologyDataRequest().Set("TotalAlveolarVentilation",    VolumePerTimeUnit::L_Per_min);
  drm.CreatePhysiologyDataRequest().Set("TotalMetabolicRate",          PowerUnit::W);
  drm.CreatePhysiologyDataRequest().Set("SkinTemperature",       TemperatureUnit::C);

  SESubstance* epi = bg.GetSubstanceManager().GetSubstance("Epinephrine");
  if (epi != nullptr) {
    drm.CreateSubstanceDataRequest().Set(
        *epi, "PlasmaConcentration", MassPerVolumeUnit::ug_Per_L);
  } else {
    bg.GetLogger()->Warning("Epinephrine substance not found - skipping data request");
  }

  drm.SetResultsFilename(kCsvPath);
}

//  Programa principal

int HowToColdFaceTest()
{
  std::cout << "\n========================================================\n"
            << " Protocol D - CFT x Repeated Stress (TCR habituation)\n"
            << "========================================================\n";

  const auto t_wall_start = std::chrono::steady_clock::now();

  std::unique_ptr<PhysiologyEngine> bg = CreateBioGearsEngine(kLogPath);
  bg->GetLogger()->Info("HowTo-ColdFaceTest (Protocol D)");

  if (!bg->InitializeEngine(PATIENT_XML)) {
    bg->GetLogger()->Error("Could not load patient, check the error");
    return 1;
  }

  configure_data_requests(*bg);
  bg->SetAutoTrackFlag(true);
  bg->SetTrackStabilizationFlag(false);

  const double mass_kg = bg->GetPatient().GetWeight(MassUnit::kg);
  const double age_y   = bg->GetPatient().GetAge(TimeUnit::yr);
  bg->GetLogger()->Info(asprintf("Patient: StandardMale, age=%.1f y, mass=%.1f kg",
                                  age_y, mass_kg));

  //  F1
  PhaseStats f1;
  run_block(*bg, kPhase1_Baseline_s, "F1 Baseline", f1);
  const double HR_F1_mean  = f1.HR_sum  / std::max(1, f1.N);
  const double SBP_F1_mean = f1.SBP_sum / std::max(1, f1.N);
  const double MAP_F1_mean = f1.MAP_sum / std::max(1, f1.N);
  bg->GetLogger()->Info(asprintf("F1 baseline: HR=%.1f bpm, SBP=%.1f mmHg, MAP=%.1f mmHg",
                                  HR_F1_mean, SBP_F1_mean, MAP_F1_mean));

  // F2  CFT  baseline
  apply_cft(*bg, kCFTSeverity, kCFT_apply_s);
  bg->GetLogger()->Info(asprintf("CFT#1 ON  (sev=%.2f, dur=%.0f s) - baseline",
                                  kCFTSeverity, kCFT_apply_s));
  PhaseStats f2_apply;
  run_block(*bg, kCFT_apply_s, "F2 CFT#1 apply", f2_apply);
  clear_cft(*bg);
  PhaseStats f2_recov;
  run_block(*bg, kCFT_recov_s, "F2 CFT#1 recov", f2_recov);

  const double HR_F2_min   = f2_apply.HR_min;
  const double SBP_F2_max  = f2_apply.SBP_max;
  const double dHR_F2_pct  = (HR_F2_min  - HR_F1_mean ) / HR_F1_mean  * 100.0;
  const double dSBP_F2_mmHg=  SBP_F2_max - SBP_F1_mean;
  const double t_HRmin_F2  = f2_apply.t_HRmin_rel_s;

  //  F3  Stress ON #1
  apply_stress(*bg, kAcuteStressSeverity);
  bg->GetLogger()->Info(asprintf("Stress ON #1 (sev=%.2f)", kAcuteStressSeverity));
  PhaseStats f3;
  run_block(*bg, kPhase3_StressOn1_s, "F3 Stress ON #1", f3);
  const double HR_F3_mean  = f3.HR_sum  / std::max(1, f3.N);
  const double SBP_F3_mean = f3.SBP_sum / std::max(1, f3.N);
  bg->GetLogger()->Info(asprintf("F3 stress steady: HR=%.1f bpm (Δ=%+.1f), SBP=%.1f mmHg (Δ=%+.1f)",
                                  HR_F3_mean, HR_F3_mean - HR_F1_mean,
                                  SBP_F3_mean, SBP_F3_mean - SBP_F1_mean));

  // F4  Stress decay #1 (Usar exp (?))
  PhaseStats f4;
  const double sev_pos_F4 = run_stress_decay_linear(
      *bg, kAcuteStressSeverity, kStressDecayRate_per_s,
      kPhase4_StressOff1_s, "F4 Stress decay #1", f4);
  bg->GetLogger()->Info(asprintf(
      "F4 end: sev_residual=%.3f (paciente %s estressado)",
      sev_pos_F4,
      sev_pos_F4 > 0.30 ? "ainda" :
      sev_pos_F4 > 0.10 ? "parcialmente" : "praticamente nao"));
  const double HR_pre_CFT2  = tail_mean(f4.HR_samples,  kTailWindow_s);
  const double SBP_pre_CFT2 = tail_mean(f4.SBP_samples, kTailWindow_s);

  // F5  CFT após stress #1
  apply_cft(*bg, kCFTSeverity, kCFT_apply_s);
  bg->GetLogger()->Info("CFT#2 ON  (post-stress #1)");
  PhaseStats f5_apply;
  run_block(*bg, kCFT_apply_s, "F5 CFT#2 apply", f5_apply);
  clear_cft(*bg);
  PhaseStats f5_recov;
  run_block(*bg, kCFT_recov_s, "F5 CFT#2 recov", f5_recov);

  const double HR_F5_min    = f5_apply.HR_min;
  const double SBP_F5_max   = f5_apply.SBP_max;
  const double dHR_F5_pct   = (HR_F5_min  - HR_pre_CFT2 ) / HR_pre_CFT2  * 100.0;
  const double dSBP_F5_mmHg =  SBP_F5_max - SBP_pre_CFT2;
  const double t_HRmin_F5   = f5_apply.t_HRmin_rel_s;

  // F6  Stress ON #2
  apply_stress(*bg, kAcuteStressSeverity);
  bg->GetLogger()->Info(asprintf("Stress ON #2 (sev=%.2f)", kAcuteStressSeverity));
  PhaseStats f6;
  run_block(*bg, kPhase6_StressOn2_s, "F6 Stress ON #2", f6);
  const double HR_F6_mean  = f6.HR_sum  / std::max(1, f6.N);
  const double SBP_F6_mean = f6.SBP_sum / std::max(1, f6.N);

  //  F7  Stress decay #2 (taxa fixa)
  PhaseStats f7;
  const double sev_pos_F7 = run_stress_decay_linear(
      *bg, kAcuteStressSeverity, kStressDecayRate_per_s,
      kPhase7_StressOff2_s, "F7 Stress decay #2", f7);
  bg->GetLogger()->Info(asprintf(
      "F7 end: sev_residual=%.3f (paciente %s estressado)",
      sev_pos_F7,
      sev_pos_F7 > 0.30 ? "ainda" :
      sev_pos_F7 > 0.10 ? "parcialmente" : "praticamente nao"));
  const double HR_pre_CFT3  = tail_mean(f7.HR_samples,  kTailWindow_s);
  const double SBP_pre_CFT3 = tail_mean(f7.SBP_samples, kTailWindow_s);

  //  F8  CFT após stress #2
  apply_cft(*bg, kCFTSeverity, kCFT_apply_s);
  bg->GetLogger()->Info("CFT#3 ON  (post-stress #2)");
  PhaseStats f8_apply;
  run_block(*bg, kCFT_apply_s, "F8 CFT#3 apply", f8_apply);
  clear_cft(*bg);
  PhaseStats f8_recov;
  run_block(*bg, kCFT_recov_s, "F8 CFT#3 recov", f8_recov);

  const double HR_F8_min    = f8_apply.HR_min;
  const double SBP_F8_max   = f8_apply.SBP_max;
  const double dHR_F8_pct   = (HR_F8_min  - HR_pre_CFT3 ) / HR_pre_CFT3  * 100.0;
  const double dSBP_F8_mmHg =  SBP_F8_max - SBP_pre_CFT3;
  const double t_HRmin_F8   = f8_apply.t_HRmin_rel_s;

  //  F9  Post-baseline

  clear_stress(*bg);
  bg->GetLogger()->Info("F9 Post-baseline / washout (stress fully cleared)");
  const double dt  = 1.0 / kEngineHz;
  const int    sample_every = static_cast<int>(kEngineHz / kSampleHz);
  const int    steps = static_cast<int>(kPhase9_PostBase_s / dt);
  const double tol = 0.05 * HR_F1_mean;
  double tau_recovery_s = -1.0;
  for (int i = 0; i < steps; ++i) {
    bg->AdvanceModelTime(dt, TimeUnit::s);
    if ((i % sample_every) == 0 && tau_recovery_s < 0.0) {
      const double hr = bg->GetCardiovascularSystem()->GetHeartRate(FrequencyUnit::Per_min);
      if (std::abs(hr - HR_F1_mean) <= tol) {
        tau_recovery_s = static_cast<double>(i) * dt;
      }
    }
  }
  if (tau_recovery_s < 0.0) tau_recovery_s = kPhase9_PostBase_s; 

  const auto t_wall_end = std::chrono::steady_clock::now();
  const double wallclock_s = std::chrono::duration<double>(t_wall_end - t_wall_start).count();

  // Métricas comparativas

  std::cout << "Outputs:\n  " << kCsvPath << "\n  " << kLogPath << "\n  ";
  return 0;
}

int main()
{
  return HowToColdFaceTest();
}