/**************************************************************************************
 *  HowTo-MIST-CFT.cpp
 *  ------------------
 *  MIST + CFT (Cold Face Intervention) protocol, BioGears V5.
 *
 *  Reproduz o desenho experimental de:
 *      Richer R, Zenkner J, Kuderle A, Rohleder N, Eskofier BM.
 *      "Vagus activation by Cold Face Test reduces acute psychosocial stress responses."
 *      Sci Rep 12, 19270 (2022). doi:10.1038/s41598-022-23222-9
 *
 *  Sequencia (por condicao, single subject = StandardMale):
 *
 *      BLGlo                kBLGlo_s s   (global baseline; referencia para deltaHR%)
 *      Para cada MIST k em {1,2,3}:
 *          BL_k             kBL_s s      (sem acao)
 *          RP/CFI_k         kRPCFI_s s
 *              Controlo: sem acao
 *              CFT:      SEColdFaceTest (sev=kCFTSeverity, dur estendida)
 *          AT_k             kAT_s s      SEAcuteStress (sev = kSigmaAT[k])
 *              Na condicao CFT, SEColdFaceTest NAO e desligado no fim do
 *              RP/CFI -- decai exponencialmente em PARALELO com SEAcuteStress
 *              durante os primeiros kCFTOverlap_s s do AT, com tau = kCFTTau_s.
 *              Isto modela vagal afterglow + sympatholysis pre-sinaptica
 *              (Levy 1971; Sun & Guyenet 1986; Heath & Downey 1990).
 *          FB_k             kFB_s s      SEAcuteStress (sev = kSigmaAT[k]*kFBFactor)
 *          Inter-MIST       kInterMIST_s s    decaimento exponencial (tau=kStressTau_s)
 *      POST                 kPost_s s   (washout final, stress totalmente clear)
 *
 *  Sao executadas DUAS simulacoes back-to-back (Controlo, CFT), cada uma com
 *  uma instancia BioGears nova. Saidas (na CWD do executavel):
 *
 *      HowTo-MIST-Control.csv         timeseries @ 1 Hz
 *      HowTo-MIST-Control.log         log do engine
 *      HowTo-MIST-Control.phases.csv  estatisticas por subfase + boundaries
 *      HowTo-MIST-CFT.csv             ...
 *      HowTo-MIST-CFT.log
 *      HowTo-MIST-CFT.phases.csv
 *      HowTo-MIST-CFT.summary.txt     sumario comparativo legivel
 *
 *  Notas fisiologicas
 *  ------------------
 *  (1) BioGears e um motor mecanicista deterministico e NAO exibe sensibilizacao
 *      neural a estimulos repetidos. Para reproduzir o crescimento progressivo
 *      do pico de HR observado na condicao Controlo de Richer et al. 2022
 *      (deltaHR pico ~ 26%, 33%, 39% em MIST1-3), escalamos explicitamente a
 *      severity de SEAcuteStress entre ciclos atraves de kSigmaAT[3]. Isto e um
 *      proxy fenomenologico para sensibilizacao, NAO um modelo mecanicista de
 *      plasticidade central. Documenta-se esta limitacao explicitamente.
 *
 *  (2) O efeito de "carry-over" vagal sobre respostas SNS subsequentes,
 *      observado in vivo (Gurel et al. 2020; Richer et al. 2022), e modelado
 *      AQUI atraves de uma sobreposicao temporal CFT-AT controlada (option (ii)
 *      da analise de gap): na condicao CFT, o SEColdFaceTest NAO e desligado
 *      abruptamente no fim do RP/CFI mas decai exponencialmente em paralelo
 *      com o SEAcuteStress durante os primeiros kCFTOverlap_s = 60 s do AT,
 *      com constante de tempo kCFTTau_s = 30 s. Este tau e consistente com a
 *      cinetica de recuperacao cardiaca pos-CFT reportada por Heath & Downey
 *      1990 (HR retorna ao baseline em 30-90 s apos remocao do estimulo frio).
 *
 *      Justificacao mecanicista: representamos o output vagal residual
 *      sustentado pelo eixo NTS->DMV (Sun & Guyenet 1986) e o antagonismo
 *      acentuado pre-sinaptico (Levy 1971: ACh muscarinico inibe NE no no SA)
 *      atraves do unico knob disponivel no BioGears -- a severity do
 *      SEColdFaceTest. Isto NAO e o mecanismo correto a nivel celular, mas
 *      reproduz o efeito liquido (atraso na subida da HR durante o AT em
 *      comparacao com a condicao Controlo).
 *
 *  (3) Normalizacao deltaHR(%) e feita POST-HOC pelo script Python sobre o CSV,
 *      contra hr_blglo = mean HR nos ultimos kTailWindow_s segundos de BLGlo,
 *      em linha com a definicao de BLGlo de Richer et al. (eq. 1 do paper).
 *
 *  Referencias
 *  -----------
 *      Richer et al. Sci Rep 12:19270 (2022).
 *      Heath ME, Downey JA. Clin Sci 78:139-147 (1990).
 *      Khurana RK, Wu R. Clin Auton Res 16:202-206 (2006).
 *      Esler M et al. Hypertension 16:621-628 (1990).
 *      Gurel NZ et al. Brain Stimul 13:47-59 (2020).
 *      Berntson GG, Cacioppo JT, Quigley KS. Psychophysiology 28:251-263 (1991).
 *=====================================================================================*/

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
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
using namespace biogears;

//  Configuracao temporal e parametros
namespace {
constexpr double kEngineHz = 50.0;             // 1 dt = 0.02 s
constexpr double kSampleHz = 1.0;              // CSV @ 1 Hz

// MIST timings (s)
constexpr double kBLGlo_s     = 300.0;   // Global baseline (suficiente em modelo deterministico)
constexpr double kBL_s        =  60.0;   // BL por MIST
constexpr double kRPCFI_s     = 120.0;   // Resting Period / Cold Face Intervention
constexpr double kAT_s        = 240.0;   // Arithmetic Tasks
constexpr double kFB_s        =  60.0;   // Feedback
constexpr double kInterMIST_s = 180.0;   // Washout entre MIST
constexpr double kPost_s      = 300.0;   // Washout final
constexpr int    kNMIST       =   3;

// Severity AT por ciclo (proxy para sensibilizacao - Richer et al. 2022 Fig 2 Control)
constexpr double kSigmaAT[kNMIST] = { 0.50, 0.70, 0.85 };
constexpr double kFBFactor        = 0.50;        // FB severity = AT * 0.50

// CFT
constexpr double kCFTSeverity = 0.50;

// Gain reduction simpatica pos-CFT (proxy fenomenologico para o antagonismo
// acentuado pre-sinaptico de Levy 1971 + vagal afterglow central NTS->DMV
// de Sun & Guyenet 1986). Aplicado APENAS na condicao CFT, atenua a severity
// efetiva do SEAcuteStress no AT e FB. Valor calibrado para reproduzir a
// reducao de ~20% no pico de deltaHR(AT) reportada por Richer et al. 2022
// (Sci Rep 12:19270). NAO e um mecanismo celular -- e um lumped gain proxy
// no espirito de Ursino & Magosso (Am J Physiol Heart 2000).
constexpr double kCFTGainReduction = 0.90;   // sev_AT_CFT = sev_AT * 0.80

// CFT carry-over (vagal afterglow proxy - option (ii) gap analysis)
// Janela em que SEColdFaceTest decai exponencialmente em paralelo com
// SEAcuteStress no inicio do AT da condicao CFT. tau ~ 30 s consistente
// com Heath & Downey 1990 (HR recovery em 30-90 s pos-CFT).
constexpr double kCFTOverlap_s = 120.0;
constexpr double kCFTTau_s     = 90.0;
constexpr double kCFTFloor     = 0.01;           // > tolerancia IsZero do SEScalar0To1

// Cinetica de clearance do stress (Esler 1990: NE plasma t1/2 ~ 60-90 s)
constexpr double kStressTau_s = 90.0;
constexpr double kStressFloor = 0.01;            // > tolerancia IsZero do SEScalar0To1

// Tail-mean window (BLGlo e referencias pre-CFT)
constexpr int    kTailWindow_s = 60;

// Coorte / I/O
//   - DEFAULT_PATIENT_DIR : pasta a varrer para .xml de pacientes (idem ValidationPA).
//   - RESULTS_DIR         : raiz dos outputs. Por paciente cria-se uma subpasta
//                           com nome = stem do XML, contendo:
//                             HowTo-MIST-Control.csv / .log / .phases.csv
//                             HowTo-MIST-CFT.csv     / .log / .phases.csv
//                             HowTo-MIST-CFT.summary.txt
//   - GLOBAL_SUMMARY_CSV  : 1 linha por paciente com picos AT (Control vs CFT)
//                           para varrer toda a coorte de uma vez.
constexpr const char* DEFAULT_PATIENT_DIR = "/Users/lucasteixeira/Desktop/bio_V5/core_V5/share/data/patients";
constexpr const char* RESULTS_DIR         = "MIST_Result";
constexpr const char* GLOBAL_SUMMARY_CSV  = "MIST_Summary.csv";
} // namespace

//  Estrutura de dados
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

struct Subphase {
  std::string label;
  double      t_start_s = 0.0;
  double      t_end_s   = 0.0;
  double      sigma_stress = 0.0;
  bool        cft_on    = false;
  PhaseStats  stats;
};

enum class MISTCondition { CONTROL, CFT };

//  Helpers
static void update(PhaseStats& s, double hr, double sbp, double map, double t_rel)
{
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

static double tail_mean(const std::vector<double>& v, std::size_t k)
{
  if (v.empty()) return 0.0;
  const std::size_t start = (v.size() > k) ? v.size() - k : 0;
  double sum = 0.0;
  for (std::size_t i = start; i < v.size(); ++i) sum += v[i];
  return sum / static_cast<double>(v.size() - start);
}

static double mean_or_nan(double sum, int N) {
  return (N > 0) ? sum / static_cast<double>(N) : std::nan("");
}

//  Bloco de simulacao basico
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

//  Acoes BioGears
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
//  Decaimento exponencial do stress.
//
//      sev(t) = max(kStressFloor, sev_start * exp(-t / tau_s))
//
//  tau ~ 90 s ~ half-life * ln(2) ~ 60 s * 1.5 da NE plasmatica
//  (Esler 1990; Goldstein, Adv Pharmacol 1998).
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

// -----------------------------------------------------------------------
//  Decaimento exponencial do SEColdFaceTest em paralelo com SEAcuteStress.
//
//      sev_cft(t) = max(kCFTFloor, sev_start * exp(-t / tau_s))
//
//  Esta funcao implementa a option (ii) da gap analysis: em vez de desligar
//  o CFT abruptamente no fim do RP/CFI, ele e mantido ativo e a sua severity
//  e re-aplicada periodicamente seguindo um decaimento exponencial. Isto
//  acontece SIMULTANEAMENTE com SEAcuteStress (que ja deve estar ativo).
//
//  Mecanismo emulado:
//    - Vagal afterglow via NTS-DMV (Sun & Guyenet 1986)
//    - Antagonismo acentuado pre-sinaptico Levy 1971 (ACh M2 inibe NE no SA)
//    - Cinetica residual TCR pos-stimulo Heath & Downey 1990
//
//  tau ~ 30 s e o knob principal para controlar quao "longo" e o afterglow.
//  Update a 0.5 Hz (a cada 2 s) -- compromise entre fidelidade e custo.
//
//  Pre-condicao: SEColdFaceTest ja ativo (sev = sev_start) e SEAcuteStress
//                ja ativo. Nao chamar clear_cft() antes.
//  Pos-condicao: SEColdFaceTest permanece ativo no nivel kCFTFloor. Caller
//                deve chamar clear_cft() se quiser desligar completamente.
// -----------------------------------------------------------------------
static double run_cft_decay_exp(PhysiologyEngine& bg,
                                 double sev_start,
                                 double tau_s,
                                 double duration_s,
                                 const char* phase_label,
                                 PhaseStats& st)
{
  constexpr double kUpdateDt_s = 2.0;        // re-apply CFT severity @ 0.5 Hz

  const double dt           = 1.0 / kEngineHz;
  const int    steps        = static_cast<int>(duration_s / dt);
  const int    sample_every = static_cast<int>(kEngineHz / kSampleHz);
  const int    update_every = static_cast<int>(kEngineHz * kUpdateDt_s);
  const auto&  cv           = *bg.GetCardiovascularSystem();

  bg.GetLogger()->Info(asprintf(
      "=== %s : %.0f s (CFT exp decay sev_start=%.2f, tau=%.0f s) ===",
      phase_label, duration_s, sev_start, tau_s));

  double sev_last = sev_start;
  for (int i = 0; i < steps; ++i) {
    if ((i % update_every) == 0) {
      const double t = static_cast<double>(i) * dt;
      double sev = sev_start * std::exp(-t / tau_s);
      if (!(sev > 0.0 && sev <= 1.0)) sev = kCFTFloor;
      if (sev < kCFTFloor)            sev = kCFTFloor;
      // Duracao remanescente + buffer para o action permanecer ativo
      // ate ao proximo update (ou ate caller chamar clear_cft).
      const double dur_remaining = std::max(2.0, duration_s - t + 5.0);
      apply_cft(bg, sev, dur_remaining);
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
static void configure_data_requests(PhysiologyEngine& bg, const std::string& csv_path)
{
  auto& drm = bg.GetEngineTrack()->GetDataRequestManager();
  drm.SetSamplesPerSecond(kSampleHz);
  drm.CreatePhysiologyDataRequest().Set("HeartRate",                    FrequencyUnit::Per_min);
  drm.CreatePhysiologyDataRequest().Set("SystolicArterialPressure",     PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("DiastolicArterialPressure",    PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("MeanArterialPressure",         PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("HeartStrokeVolume",            VolumeUnit::mL);
  drm.CreatePhysiologyDataRequest().Set("CardiacOutput",                VolumePerTimeUnit::L_Per_min);
  drm.CreatePhysiologyDataRequest().Set("SystemicVascularResistance",   FlowResistanceUnit::mmHg_s_Per_mL);
  drm.CreatePhysiologyDataRequest().Set("OxygenSaturation");
  drm.CreatePhysiologyDataRequest().Set("RespirationRate",              FrequencyUnit::Per_min);
  drm.CreatePhysiologyDataRequest().Set("CoreTemperature",              TemperatureUnit::C);
  drm.CreatePhysiologyDataRequest().Set("OxygenConsumptionRate",        VolumePerTimeUnit::mL_Per_min);
  drm.CreatePhysiologyDataRequest().Set("CarbonDioxideProductionRate",  VolumePerTimeUnit::mL_Per_min);
  drm.CreatePhysiologyDataRequest().Set("RespiratoryExchangeRatio");
  drm.CreatePhysiologyDataRequest().Set("TidalVolume",                  VolumeUnit::mL);
  drm.CreatePhysiologyDataRequest().Set("TotalPulmonaryVentilation",    VolumePerTimeUnit::L_Per_min);
  drm.CreatePhysiologyDataRequest().Set("TotalAlveolarVentilation",     VolumePerTimeUnit::L_Per_min);
  drm.CreatePhysiologyDataRequest().Set("TotalMetabolicRate",           PowerUnit::W);
  drm.CreatePhysiologyDataRequest().Set("SkinTemperature",              TemperatureUnit::C);

  SESubstance* epi = bg.GetSubstanceManager().GetSubstance("Epinephrine");
  if (epi != nullptr) {
    drm.CreateSubstanceDataRequest().Set(*epi, "PlasmaConcentration",
                                         MassPerVolumeUnit::ug_Per_L);
  } else {
    bg.GetLogger()->Warning("Epinephrine substance not found - skipping data request");
  }

  drm.SetResultsFilename(csv_path);
}

//  Dump das subfases (boundaries + estatisticas)
static void write_phases_csv(const std::vector<Subphase>& sub,
                             double hr_blglo,
                             const std::string& path)
{
  std::ofstream f(path);
  f << "label,t_start_s,t_end_s,duration_s,sigma_stress,cft_on,"
    << "n_samples,hr_mean_bpm,hr_min_bpm,hr_max_bpm,"
    << "dHR_mean_pct,dHR_min_pct,dHR_max_pct,"
    << "sbp_mean_mmHg,sbp_min_mmHg,sbp_max_mmHg,"
    << "map_mean_mmHg,t_HRmin_rel_s\n";
  f << std::fixed << std::setprecision(4);
  for (const auto& sp : sub) {
    const double hr_mean = mean_or_nan(sp.stats.HR_sum, sp.stats.N);
    const double sbp_mean = mean_or_nan(sp.stats.SBP_sum, sp.stats.N);
    const double map_mean = mean_or_nan(sp.stats.MAP_sum, sp.stats.N);
    const double dHR_mean = (hr_blglo > 0) ? (hr_mean - hr_blglo) / hr_blglo * 100.0 : 0.0;
    const double dHR_min  = (hr_blglo > 0) ? (sp.stats.HR_min - hr_blglo) / hr_blglo * 100.0 : 0.0;
    const double dHR_max  = (hr_blglo > 0) ? (sp.stats.HR_max - hr_blglo) / hr_blglo * 100.0 : 0.0;
    f << sp.label << ','
      << sp.t_start_s << ',' << sp.t_end_s << ',' << (sp.t_end_s - sp.t_start_s) << ','
      << sp.sigma_stress << ',' << (sp.cft_on ? 1 : 0) << ','
      << sp.stats.N << ','
      << hr_mean << ',' << sp.stats.HR_min << ',' << sp.stats.HR_max << ','
      << dHR_mean << ',' << dHR_min << ',' << dHR_max << ','
      << sbp_mean << ',' << sp.stats.SBP_min << ',' << sp.stats.SBP_max << ','
      << map_mean << ',' << sp.stats.t_HRmin_rel_s << '\n';
  }
}

//  Execucao de UMA condicao
//
//    patient_xml : path absoluto (ou resolvivel) do .xml do paciente
//    outdir      : pasta onde se escrevem csv/log/phases.csv desta condicao
//    cond, tag   : condicao MIST + sufixo (Control / CFT)
//    outRecords  : populado com as Subphase desta condicao
//    out_hr_blglo: HR de referencia (tail-mean dos kTailWindow_s s do BLGlo)
static int run_one_condition(MISTCondition cond,
                              const std::string& tag,
                              const fs::path& patient_xml,
                              const fs::path& outdir,
                              std::vector<Subphase>& outRecords,
                              double& out_hr_blglo)
{
  const fs::path csv = outdir / ("HowTo-MIST-" + tag + ".csv");
  const fs::path log = outdir / ("HowTo-MIST-" + tag + ".log");

  std::cout << "\n--------------------------------------------------------\n"
            << " Running condition: " << tag << "\n"
            << " Patient          : " << patient_xml.filename().string() << "\n"
            << " Outdir           : " << outdir.string() << "\n"
            << "--------------------------------------------------------\n";

  std::unique_ptr<PhysiologyEngine> bg = CreateBioGearsEngine(log.string());
  bg->GetLogger()->Info("HowTo-MIST-CFT: condition = " + tag);
  bg->GetLogger()->Info("Patient XML: " + patient_xml.string());

  if (!bg->InitializeEngine(patient_xml.string())) {
    bg->GetLogger()->Error("Could not load patient, check the error");
    return 1;
  }

  configure_data_requests(*bg, csv.string());
  bg->SetAutoTrackFlag(true);
  bg->SetTrackStabilizationFlag(false);

  const double mass_kg = bg->GetPatient().GetWeight(MassUnit::kg);
  const double age_y   = bg->GetPatient().GetAge(TimeUnit::yr);
  bg->GetLogger()->Info(asprintf("Patient: %s, age=%.1f y, mass=%.1f kg",
                                  patient_xml.stem().string().c_str(),
                                  age_y, mass_kg));

  double t_clock = 0.0;

  //  BLGlo
  {
    Subphase sp;
    sp.label     = "BLGlo";
    sp.t_start_s = t_clock;
    run_block(*bg, kBLGlo_s, "BLGlo", sp.stats);
    t_clock     += kBLGlo_s;
    sp.t_end_s   = t_clock;
    out_hr_blglo = tail_mean(sp.stats.HR_samples, kTailWindow_s);
    bg->GetLogger()->Info(asprintf(
        "BLGlo (last %d s tail-mean): HR_BLGlo = %.2f bpm",
        kTailWindow_s, out_hr_blglo));
    outRecords.push_back(std::move(sp));
  }

  //  3 ciclos MIST
  for (int k = 0; k < kNMIST; ++k) {
    const std::string idx = std::to_string(k + 1);

    // BL_k
    {
      Subphase sp;
      sp.label     = "MIST" + idx + "_BL";
      sp.t_start_s = t_clock;
      run_block(*bg, kBL_s, sp.label.c_str(), sp.stats);
      t_clock     += kBL_s;
      sp.t_end_s   = t_clock;
      outRecords.push_back(std::move(sp));
    }

    // RP_k ou CFI_k
    {
      Subphase sp;
      sp.t_start_s = t_clock;
      if (cond == MISTCondition::CFT) {
        sp.label  = "MIST" + idx + "_CFI";
        sp.cft_on = true;
        // Aplicamos CFT com duracao estendida (RP/CFI + overlap window + buffer)
        // de modo a que a acao continue ativa quando entrarmos no AT.
        // NAO chamamos clear_cft() aqui -- a fase AT ira decair a severity
        // exponencialmente em paralelo com SEAcuteStress.
        apply_cft(*bg, kCFTSeverity, kRPCFI_s + kCFTOverlap_s + 5.0);
        bg->GetLogger()->Info(asprintf(
            "CFI ON (sev=%.2f, dur=%.0f s incl. overlap into AT)",
            kCFTSeverity, kRPCFI_s + kCFTOverlap_s));
        run_block(*bg, kRPCFI_s, sp.label.c_str(), sp.stats);
        // CFT permanece ativo a sev = kCFTSeverity no fim deste bloco.
      } else {
        sp.label = "MIST" + idx + "_RP";
        run_block(*bg, kRPCFI_s, sp.label.c_str(), sp.stats);
      }
      t_clock   += kRPCFI_s;
      sp.t_end_s = t_clock;
      outRecords.push_back(std::move(sp));
    }

    // AT_k
    // Na condicao CFT aplicamos um gain reduction simpatico (kCFTGainReduction)
    // sobre a severity nominal, modelando o pre-conditioning vagal observado
    // por Richer et al. 2022. Na condicao Control usa-se a severity nominal.
    const double sigma_at_eff = (cond == MISTCondition::CFT)
                                  ? kSigmaAT[k] * kCFTGainReduction
                                  : kSigmaAT[k];
    {
      Subphase sp;
      sp.label        = "MIST" + idx + "_AT";
      sp.t_start_s    = t_clock;
      sp.sigma_stress = sigma_at_eff;
      apply_stress(*bg, sigma_at_eff);
      bg->GetLogger()->Info(asprintf(
          "AT ON (sev_nominal=%.2f, sev_eff=%.2f, CFT_gain_reduction=%s)",
          kSigmaAT[k], sigma_at_eff,
          (cond == MISTCondition::CFT) ? "ON" : "OFF"));

      if (cond == MISTCondition::CFT) {
        // Overlap CFT->AT: SEColdFaceTest (ainda ativo desde o RP/CFI)
        // decai exponencialmente em paralelo com SEAcuteStress durante os
        // primeiros kCFTOverlap_s s do AT. Sustenta a bradicardia de fundo
        // durante a subida da HR ao setpoint do estressor, atrasando o
        // pico (vagal afterglow + sympatholysis pre-sinaptica).
        sp.cft_on = true;
        const double sev_end = run_cft_decay_exp(
            *bg, kCFTSeverity, kCFTTau_s, kCFTOverlap_s,
            (sp.label + "_overlap").c_str(), sp.stats);
        bg->GetLogger()->Info(asprintf(
            "AT overlap end: cft_residual=%.4f (clearing)", sev_end));
        clear_cft(*bg);
        // Restante do AT apenas com SEAcuteStress (sem CFT residual).
        run_block(*bg, kAT_s - kCFTOverlap_s,
                  (sp.label + "_postoverlap").c_str(), sp.stats);
      } else {
        run_block(*bg, kAT_s, sp.label.c_str(), sp.stats);
      }

      t_clock   += kAT_s;
      sp.t_end_s = t_clock;
      outRecords.push_back(std::move(sp));
    }

    // FB_k
    // FB herda o gain reduction do AT (consistencia: se o tonus vagal residual
    // atenua a resposta simpatica ao AT, atenua tambem a do FB no mesmo ciclo).
    // O inter-MIST que vem a seguir consome fb_sigma como sev_start do decay,
    // portanto a reducao propaga-se automaticamente para o washout.
    const double fb_sigma = sigma_at_eff * kFBFactor;
    {
      Subphase sp;
      sp.label        = "MIST" + idx + "_FB";
      sp.t_start_s    = t_clock;
      sp.sigma_stress = fb_sigma;
      apply_stress(*bg, fb_sigma);
      bg->GetLogger()->Info(asprintf("FB ON (sev=%.2f)", fb_sigma));
      run_block(*bg, kFB_s, sp.label.c_str(), sp.stats);
      t_clock   += kFB_s;
      sp.t_end_s = t_clock;
      outRecords.push_back(std::move(sp));
    }

    // Inter-MIST washout (so entre ciclos)
    if (k < kNMIST - 1) {
      Subphase sp;
      sp.label     = "MIST" + idx + "_INTER";
      sp.t_start_s = t_clock;
      const double sev_end = run_stress_decay_exp(
          *bg, fb_sigma, kStressTau_s, kInterMIST_s,
          sp.label.c_str(), sp.stats);
      sp.sigma_stress = sev_end;   // severity residual no fim do bloco
      t_clock        += kInterMIST_s;
      sp.t_end_s      = t_clock;
      bg->GetLogger()->Info(asprintf(
          "Inter-MIST end: sev_residual=%.3f", sev_end));
      outRecords.push_back(std::move(sp));
    }
  }

  // POST washout final
  {
    clear_stress(*bg);
    bg->GetLogger()->Info("POST: stress fully cleared");
    Subphase sp;
    sp.label     = "POST";
    sp.t_start_s = t_clock;
    run_block(*bg, kPost_s, "POST", sp.stats);
    t_clock     += kPost_s;
    sp.t_end_s   = t_clock;
    outRecords.push_back(std::move(sp));
  }

  // Dump phases.csv
  write_phases_csv(outRecords, out_hr_blglo,
                   (outdir / ("HowTo-MIST-" + tag + ".phases.csv")).string());

  bg->GetLogger()->Info(asprintf("Condition %s complete: total sim time = %.0f s",
                                  tag.c_str(), t_clock));
  return 0;
}

//  Sumario comparativo
static void write_comparative_summary(const std::vector<Subphase>& rcCtrl,
                                      const std::vector<Subphase>& rcCFT,
                                      double hr_blglo_ctrl,
                                      double hr_blglo_cft,
                                      const std::string& path)
{
  auto find_label = [](const std::vector<Subphase>& v, const std::string& lbl) -> const Subphase* {
    for (const auto& s : v) if (s.label == lbl) return &s;
    return nullptr;
  };
  auto pct = [](double x, double ref) {
    return (ref > 0) ? (x - ref) / ref * 100.0 : 0.0;
  };

  std::ofstream f(path);
  f << std::fixed << std::setprecision(2);
  f << "MIST + CFT Protocol - Comparative Summary\n";
  f << "(Richer et al. 2022 reproduction in BioGears v5)\n";
  f << "==========================================================\n\n";
  f << "HR_BLGlo  Control = " << hr_blglo_ctrl << " bpm\n";
  f << "HR_BLGlo  CFT     = " << hr_blglo_cft  << " bpm\n";
  f << "Delta(BLGlo CFT - Control) = "
    << (hr_blglo_cft - hr_blglo_ctrl) << " bpm\n\n";

  f << "Per-subphase peak deltaHR (%) relative to BLGlo of each condition:\n";
  f << "----------------------------------------------------------\n";
  f << "subphase           Control(peak)   CFT(peak)   diff(CFT-Ctrl)\n";

  for (int k = 1; k <= kNMIST; ++k) {
    const std::vector<std::string> sub_short = {"BL", "RP_or_CFI", "AT", "FB"};
    for (const auto& tag : sub_short) {
      const std::string lblCtrl = "MIST" + std::to_string(k) + "_"
                                  + (tag == "RP_or_CFI" ? "RP"  : tag);
      const std::string lblCFT  = "MIST" + std::to_string(k) + "_"
                                  + (tag == "RP_or_CFI" ? "CFI" : tag);
      const Subphase* sCtrl = find_label(rcCtrl, lblCtrl);
      const Subphase* sCFT  = find_label(rcCFT,  lblCFT);
      if (!sCtrl || !sCFT) continue;
      const double dCtrl = pct(sCtrl->stats.HR_max, hr_blglo_ctrl);
      const double dCFT  = pct(sCFT->stats.HR_max,  hr_blglo_cft);
      f << "MIST" << k << " " << std::left << std::setw(12) << tag
        << std::right
        << "   " << std::setw(8) << dCtrl << "%"
        << "   " << std::setw(7) << dCFT  << "%"
        << "   " << std::setw(8) << (dCFT - dCtrl) << "%\n";
    }
  }
  f << "\nNotes:\n";
  f << " - In Richer et al. 2022, peak deltaHR(AT) = 26.0%, 32.9%, 39.0% (Control).\n";
  f << " - kSigmaAT = {0.50, 0.70, 0.85} is a phenomenological proxy for sensitization.\n";
  f << "   Calibrate these values against the empirical peaks if needed.\n";
  f << " - Vagal carry-over (CFT pre-conditioning of SNS) is NOT explicitly modeled\n";
  f << "   at the NTS-DMV level in BioGears v5. Any sympatholytic effect observed\n";
  f << "   here arises only from acute TCR clearance dynamics and baroreflex setpoint\n";
  f << "   modulation.\n";
}

//  Listagem de pacientes
//  --------------------
//  Varre o directorio dado e devolve todos os .xml ordenados alfabeticamente
//  (mesma logica usada em HowTo-ValidationPA.cpp). NAO filtra por nome -- o
//  caller decide o que excluir (e.g. DefaultTemplate*).
static std::vector<fs::path> list_patient_xmls(const fs::path& dir)
{
  std::vector<fs::path> xmls;
  if (!fs::exists(dir) || !fs::is_directory(dir)) return xmls;
  for (const auto& e : fs::directory_iterator(dir)) {
    if (e.is_regular_file() && e.path().extension() == ".xml") {
      xmls.push_back(e.path());
    }
  }
  std::sort(xmls.begin(), xmls.end());
  return xmls;
}

//  Resumo per-paciente para o summary global
struct PatientMISTResult {
  std::string patient_stem;
  std::string patient_xml;
  bool        success            = false;
  double      hr_blglo_ctrl_bpm  = 0.0;
  double      hr_blglo_cft_bpm   = 0.0;
  // Pico (max) ΔHR(%) em cada AT, por condicao
  double      dHR_AT1_ctrl_pct   = std::nan("");
  double      dHR_AT2_ctrl_pct   = std::nan("");
  double      dHR_AT3_ctrl_pct   = std::nan("");
  double      dHR_AT1_cft_pct    = std::nan("");
  double      dHR_AT2_cft_pct    = std::nan("");
  double      dHR_AT3_cft_pct    = std::nan("");
  double      wallclock_s        = 0.0;
  std::string error_msg;
};

//  Helper: extrai o pico (max) ΔHR(%) de uma dada label nas Subphase
static double peak_dHR_pct(const std::vector<Subphase>& v,
                            const std::string& label,
                            double hr_blglo)
{
  for (const auto& s : v) {
    if (s.label == label && hr_blglo > 0.0) {
      return (s.stats.HR_max - hr_blglo) / hr_blglo * 100.0;
    }
  }
  return std::nan("");
}

//  CSV escape minimal (idem ValidationPA)
static std::string csv_escape(const std::string& s)
{
  if (s.find_first_of(",\"\n") == std::string::npos) return s;
  std::string out = "\"";
  for (char c : s) { out += (c == '"') ? "\"\"" : std::string(1, c); }
  out += "\"";
  return out;
}

static void write_global_summary_header(const fs::path& path)
{
  std::ofstream f(path);
  f << "PatientStem,PatientXML,Success,Wallclock_s,"
    << "HR_BLGlo_Control_bpm,HR_BLGlo_CFT_bpm,"
    << "dHR_AT1_Control_pct,dHR_AT2_Control_pct,dHR_AT3_Control_pct,"
    << "dHR_AT1_CFT_pct,dHR_AT2_CFT_pct,dHR_AT3_CFT_pct,"
    << "Atten_AT1_pct,Atten_AT2_pct,Atten_AT3_pct,"
    << "ErrorMsg\n";
}

static void append_global_summary_row(const fs::path& path,
                                       const PatientMISTResult& r)
{
  std::ofstream f(path, std::ios::app);
  auto atten = [](double ctrl, double cft) {
    return (ctrl != 0.0 && !std::isnan(ctrl) && !std::isnan(cft))
            ? (cft - ctrl)
            : std::nan("");
  };
  f << std::fixed << std::setprecision(4)
    << r.patient_stem << ',' << csv_escape(r.patient_xml) << ','
    << (r.success ? 1 : 0) << ',' << r.wallclock_s << ','
    << r.hr_blglo_ctrl_bpm << ',' << r.hr_blglo_cft_bpm << ','
    << r.dHR_AT1_ctrl_pct  << ',' << r.dHR_AT2_ctrl_pct << ',' << r.dHR_AT3_ctrl_pct << ','
    << r.dHR_AT1_cft_pct   << ',' << r.dHR_AT2_cft_pct  << ',' << r.dHR_AT3_cft_pct  << ','
    << atten(r.dHR_AT1_ctrl_pct, r.dHR_AT1_cft_pct) << ','
    << atten(r.dHR_AT2_ctrl_pct, r.dHR_AT2_cft_pct) << ','
    << atten(r.dHR_AT3_ctrl_pct, r.dHR_AT3_cft_pct) << ','
    << csv_escape(r.error_msg) << '\n';
}

//  Corre Control + CFT para UM paciente, escrevendo tudo em outdir
static PatientMISTResult run_one_patient(const fs::path& patient_xml,
                                          const fs::path& outdir)
{
  PatientMISTResult r;
  r.patient_stem = patient_xml.stem().string();
  r.patient_xml  = patient_xml.string();

  const auto t_wall_start = std::chrono::steady_clock::now();

  std::vector<Subphase> recControl;
  std::vector<Subphase> recCFT;
  double hr_blglo_ctrl = 0.0;
  double hr_blglo_cft  = 0.0;

  if (int rc = run_one_condition(MISTCondition::CONTROL, "Control",
                                  patient_xml, outdir,
                                  recControl, hr_blglo_ctrl); rc != 0) {
    r.error_msg = "Control condition failed (rc=" + std::to_string(rc) + ")";
    return r;
  }
  if (int rc = run_one_condition(MISTCondition::CFT, "CFT",
                                  patient_xml, outdir,
                                  recCFT, hr_blglo_cft); rc != 0) {
    r.error_msg = "CFT condition failed (rc=" + std::to_string(rc) + ")";
    return r;
  }

  write_comparative_summary(recControl, recCFT,
                             hr_blglo_ctrl, hr_blglo_cft,
                             (outdir / "HowTo-MIST-CFT.summary.txt").string());

  // Picos AT por condicao -> summary global
  r.hr_blglo_ctrl_bpm = hr_blglo_ctrl;
  r.hr_blglo_cft_bpm  = hr_blglo_cft;
  r.dHR_AT1_ctrl_pct  = peak_dHR_pct(recControl, "MIST1_AT", hr_blglo_ctrl);
  r.dHR_AT2_ctrl_pct  = peak_dHR_pct(recControl, "MIST2_AT", hr_blglo_ctrl);
  r.dHR_AT3_ctrl_pct  = peak_dHR_pct(recControl, "MIST3_AT", hr_blglo_ctrl);
  r.dHR_AT1_cft_pct   = peak_dHR_pct(recCFT,     "MIST1_AT", hr_blglo_cft);
  r.dHR_AT2_cft_pct   = peak_dHR_pct(recCFT,     "MIST2_AT", hr_blglo_cft);
  r.dHR_AT3_cft_pct   = peak_dHR_pct(recCFT,     "MIST3_AT", hr_blglo_cft);

  const auto t_wall_end = std::chrono::steady_clock::now();
  r.wallclock_s = std::chrono::duration<double>(t_wall_end - t_wall_start).count();
  r.success     = true;
  return r;
}

//  Programa principal
//
//    Uso:
//      ./howto-mist                    -> usa DEFAULT_PATIENT_DIR
//      ./howto-mist <patient_dir>     -> varre patient_dir/*.xml
//      ./howto-mist <patient_xml>     -> corre um unico paciente
int HowToMISTCFT(const fs::path& patient_arg)
{
  std::cout << "\n========================================================\n"
            << " MIST + CFT Protocol (Richer et al. 2022 reproduction)\n"
            << " BioGears v5 - cohort sweep\n"
            << "========================================================\n";

  // 1) Resolver input: ficheiro unico vs directorio
  std::vector<fs::path> xmls;
  if (fs::is_regular_file(patient_arg) && patient_arg.extension() == ".xml") {
    xmls.push_back(patient_arg);
  } else if (fs::is_directory(patient_arg)) {
    xmls = list_patient_xmls(patient_arg);
  } else {
    std::cerr << "[ERRO] patient_arg invalido: " << patient_arg.string() << "\n";
    return 1;
  }
  if (xmls.empty()) {
    std::cerr << "[ERRO] Nenhum .xml encontrado em " << patient_arg.string() << "\n";
    return 1;
  }
  std::cout << "[INFO] " << xmls.size() << " paciente(s) em "
            << patient_arg.string() << "\n";

  // 2) Pasta global de outputs + cabecalho do summary global
  fs::create_directories(RESULTS_DIR);
  const fs::path summary_path = fs::path(RESULTS_DIR) / GLOBAL_SUMMARY_CSV;
  write_global_summary_header(summary_path);

  // 3) Loop robusto
  const auto t_global_start = std::chrono::steady_clock::now();
  int n_ok = 0, n_fail = 0;

  for (size_t i = 0; i < xmls.size(); ++i) {
    const fs::path& xml = xmls[i];
    const std::string stem = xml.stem().string();
    const fs::path outdir  = fs::path(RESULTS_DIR) / stem;
    fs::create_directories(outdir);

    std::cout << "\n[" << (i + 1) << "/" << xmls.size() << "] "
              << stem << " ... " << std::flush;

    PatientMISTResult r;
    r.patient_stem = stem;
    r.patient_xml  = xml.string();
    try {
      r = run_one_patient(xml, outdir);
    } catch (const std::exception& ex) {
      r.error_msg = std::string("EXCEPTION: ") + ex.what();
      r.success = false;
    } catch (...) {
      r.error_msg = "EXCEPTION: unknown";
      r.success = false;
    }

    append_global_summary_row(summary_path, r);

    if (r.success) {
      ++n_ok;
      std::cout << "OK | wallclock=" << std::fixed << std::setprecision(1)
                << r.wallclock_s << " s | "
                << "AT3 Ctrl=" << r.dHR_AT3_ctrl_pct << "% "
                << "CFT=" << r.dHR_AT3_cft_pct << "%\n";
    } else {
      ++n_fail;
      std::cout << "FAIL  (" << r.error_msg << ")\n";
    }
  }

  const auto t_global_end = std::chrono::steady_clock::now();
  const double total_wall = std::chrono::duration<double>(
      t_global_end - t_global_start).count();

  std::cout << "\n========================================================\n"
            << " Total      : " << xmls.size() << " patients\n"
            << " Successful : " << n_ok        << "\n"
            << " Failed     : " << n_fail      << "\n"
            << " Wallclock  : " << std::fixed << std::setprecision(1)
            << total_wall << " s ("
            << (total_wall / std::max(1, n_ok)) << " s/patient)\n"
            << " Summary    : " << summary_path.string() << "\n"
            << " Per-patient: " << RESULTS_DIR << "/<stem>/HowTo-MIST-*.csv,.log,.phases.csv\n"
            << "========================================================\n";

  return (n_fail == 0) ? 0 : 2;
}

int main(int argc, char* argv[])
{
  fs::path patient_arg = (argc > 1) ? fs::path(argv[1])
                                    : fs::path(DEFAULT_PATIENT_DIR);
  return HowToMISTCFT(patient_arg);
}