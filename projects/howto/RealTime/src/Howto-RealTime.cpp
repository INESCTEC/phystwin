/**************************************************************************************
 *  Howto-RealTime.cpp
 *  ------------------
 *  Driver real-time interativo do BioGears V5.
 *
 *  Arquitetura:
 *    - Loop principal a 50 Hz (Δt = 0.02 s, passo nativo do engine), sincronizado
 *      com o relógio de parede via std::this_thread::sleep_until — evita drift
 *      acumulado e mantém a relação 1:1 entre tempo simulado e tempo real.
 *    - Thread separada lê comandos do stdin (std::getline bloqueante) e empurra
 *      para uma fila thread-safe; o loop principal drena a fila no início de
 *      cada tick e despacha as ações antes de chamar AdvanceModelTime.
 *
 *  Saída:
 *    - HowTo-RealTime.csv  (1 Hz):  time_s, HR_bpm, MAP_mmHg, CoreT_C, SkinT_C,
 *                                   Lactate_mmol_L, active_actions
 *    - HowTo-RealTime.log  (engine log padrão do BioGears).
 *
 *  Comandos disponíveis no terminal (durante a simulação):
 *      exercise <intensity 0-1>           ex: exercise 0.5
 *      run <speed_m/s> <incline 0-1>      ex: run 4 0.05
 *      cycle <power_W> <cadence_/min>     ex: cycle 100 60
 *      stress <severity 0-1>              ex: stress 0.7
 *      coldface <severity 0-1> <dur_s>    ex: coldface 0.5 120
 *      env <ambient_C> <RH 0-1>           ex: env 10 0.6
 *      stop <exercise|stress|coldface>    ex: stop exercise
 *      status                              imprime vitais atuais
 *      help                                lista comandos
 *      quit                                encerra
 *
 *  Notas físicas:
 *    - Lactato: medido como molaridade arterial (Aorta) em mmol/L, comparável
 *      com gasometria arterial em literatura de exercício (Brooks 2018).
 *    - "active_actions" no CSV permite reconstruir as janelas de cada estímulo
 *      sem cruzar logs de comando — essencial para análise off-line.
 *    - O engine é determinístico: rodar em "real-time" ou "fastest" dá
 *      exatamente o mesmo timeseries; o sleep só desacelera para UX.
 *=====================================================================*/

#include <biogears/cdm/compartment/SECompartmentManager.h>
#include <biogears/cdm/compartment/fluid/SELiquidCompartment.h>
#include <biogears/cdm/engine/PhysiologyEngineTrack.h>
#include <biogears/cdm/patient/SEPatient.h>
#include <biogears/cdm/patient/actions/SEAcuteStress.h>
#include <biogears/cdm/patient/actions/SEColdFaceTest.h>
#include <biogears/cdm/patient/actions/SEExercise.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <biogears/cdm/substance/SESubstance.h>
#include <biogears/cdm/substance/SESubstanceFraction.h>
#include <biogears/cdm/substance/SESubstanceManager.h>
#include <biogears/cdm/system/environment/actions/SEEnvironmentChange.h>
#include <biogears/cdm/system/physiology/SECardiovascularSystem.h>
#include <biogears/cdm/system/physiology/SEEnergySystem.h>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/string/manipulation.h>

#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <queue>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace biogears;
using clk = std::chrono::steady_clock;

namespace {

// Cadência do engine (fixa pelo BioGears: 50 Hz)
constexpr double kEngineHz             = 50.0;        // 1 dt = 0.02 s
constexpr int    kCsvDecimateTicks     = 50;          // 50 ticks @ 50Hz = 1 Hz CSV


// Caminhos relativos à CWD do executável (tipicamente runtime/)
static const bool        params_calculated = true;
static const char* const PARAMS_CFG_PATH   = "/Users/lucasteixeira/Desktop/bio_V5/core_V5/exercise_params/exercise_params_active_Male_25_Normal.cfg";
constexpr const char* kPatientXml = "patient/Male_25_Normal.xml";
constexpr const char* kLogPath    = "HowTo-RealTime.log";
constexpr const char* kCsvPath    = "HowTo-RealTime.csv";

}  // namespace

//  Fila de comandos thread-safe (produtor-consumidor clássico)
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

//  Tracker das ações ativas — strings curtas anotadas em cada linha do CSV
//  para que o pós-processamento saiba exatamente quando cada estímulo entrou
//  e saiu, sem precisar cruzar com logs de comando.

struct ActiveActions {
  std::string exercise;  // ex: "ge=0.50", "run=4m/s,i=0.05", "cyc=100W,c=60"
  std::string stress;    // ex: "0.70"
  std::string coldface;  // ex: "sev=0.5/120s"
  std::string env;       // ex: "10C/RH60"

  std::string to_csv() const {
    std::vector<std::string> parts;
    if (!exercise.empty()) parts.push_back("exercise=" + exercise);
    if (!stress.empty())   parts.push_back("stress=" + stress);
    if (!coldface.empty()) parts.push_back("coldface=" + coldface);
    if (!env.empty())      parts.push_back("env=" + env);
    if (parts.empty()) return "baseline";
    std::string out;
    for (size_t i = 0; i < parts.size(); ++i) {
      if (i) out += ";";
      out += parts[i];
    }
    return out;
  }
};

//  Helpers para aplicar / desligar ações

static void apply_exercise_generic(PhysiologyEngine& bg, double intensity) {
  SEExercise::SEGeneric g;
  g.Intensity.SetValue(intensity);
  SEExercise ex { g };
  bg.ProcessAction(ex);
}
static void apply_exercise_run(PhysiologyEngine& bg, double speed_m_s, double incline) {
  SEExercise::SERunning r;
  r.SpeedRun.SetValue(speed_m_s, LengthPerTimeUnit::m_Per_s);
  r.InclineRun.SetValue(incline);
  SEExercise ex { r };
  bg.ProcessAction(ex);
}
static void apply_exercise_cycle(PhysiologyEngine& bg, double power_W, double cadence_per_min) {
  SEExercise::SECycling c;
  c.PowerCycle.SetValue(power_W, PowerUnit::W);
  c.CadenceCycle.SetValue(cadence_per_min, FrequencyUnit::Per_min);
  SEExercise ex { c };
  bg.ProcessAction(ex);
}
static void clear_exercise(PhysiologyEngine& bg) {
  // Zerar o Exercise (setvalue(0.0))
  SEExercise::SEGeneric g;
  g.Intensity.SetValue(0.0);
  SEExercise ex { g };
  bg.ProcessAction(ex);
}

static void apply_stress(PhysiologyEngine& bg, double sev) {
  SEAcuteStress s;
  s.GetSeverity().SetValue(sev);
  bg.ProcessAction(s);
}
static void clear_stress(PhysiologyEngine& bg) { apply_stress(bg, 0.0); }

static void apply_coldface(PhysiologyEngine& bg, double sev, double dur_s) {
  SEColdFaceTest c;
  c.GetSeverity().SetValue(sev);
  c.GetDuration().SetValue(dur_s, TimeUnit::s);
  bg.ProcessAction(c);
}
static void clear_coldface(PhysiologyEngine& bg) {
  SEColdFaceTest c;
  c.GetSeverity().SetValue(0.0);
  bg.ProcessAction(c);
}

static void apply_env(PhysiologyEngine& bg, double ambient_C, double rh) {
  SEEnvironmentChange env(bg.GetSubstanceManager());
  SEEnvironmentalConditions& cd = env.GetConditions();
  cd.SetSurroundingType(SESurroundingType::Air);
  cd.GetAirVelocity().SetValue(0.1, LengthPerTimeUnit::m_Per_s);
  cd.GetAmbientTemperature().SetValue(ambient_C, TemperatureUnit::C);
  cd.GetAtmosphericPressure().SetValue(760.0, PressureUnit::mmHg);
  cd.GetClothingResistance().SetValue(0.5, HeatResistanceAreaUnit::clo);
  cd.GetEmissivity().SetValue(0.95);
  cd.GetMeanRadiantTemperature().SetValue(ambient_C, TemperatureUnit::C);
  cd.GetRelativeHumidity().SetValue(rh);
  cd.GetRespirationAmbientTemperature().SetValue(ambient_C, TemperatureUnit::C);


  auto& N2  = bg.GetSubstanceManager().GetSubstance(StandardSubstances::Nitrogen)->GetDefinition();
  auto& O2  = bg.GetSubstanceManager().GetSubstance(StandardSubstances::Oxygen)->GetDefinition();
  auto& CO2 = bg.GetSubstanceManager().GetSubstance(StandardSubstances::CarbonDioxide)->GetDefinition();
  cd.GetAmbientGas(N2 ).GetFractionAmount().SetValue(0.7896);
  cd.GetAmbientGas(O2 ).GetFractionAmount().SetValue(0.21);
  cd.GetAmbientGas(CO2).GetFractionAmount().SetValue(4.0e-4);

  bg.ProcessAction(env);
}

//  Lactato arterial em mmol/L 

static double get_arterial_lactate_mmol_per_L(PhysiologyEngine& bg) {
  auto* aorta = bg.GetCompartments().GetLiquidCompartment(BGE::VascularCompartment::Aorta);
  if (!aorta) return std::nan("");
  auto* lactate = bg.GetSubstanceManager().GetSubstance(StandardSubstances::Lactate);
  if (!lactate) return std::nan("");
  auto* sq = aorta->GetSubstanceQuantity(*lactate);
  if (!sq) return std::nan("");
  return sq->GetMolarity(AmountPerVolumeUnit::mmol_Per_L);
}

//  Help

static void print_help() {
  std::cout << "\n  Comandos disponiveis:\n"
            << "    exercise <intensity 0-1>           ex: exercise 0.5\n"
            << "    run <speed_m/s> <incline 0-1>      ex: run 4 0.05\n"
            << "    cycle <power_W> <cadence_/min>     ex: cycle 100 60\n"
            << "    stress <severity 0-1>              ex: stress 0.7\n"
            << "    coldface <severity 0-1> <dur_s>    ex: coldface 0.5 120\n"
            << "    env <ambient_C> <RH 0-1>           ex: env 10 0.6\n"
            << "    stop <exercise|stress|coldface>    ex: stop exercise\n"
            << "    status                              imprime vitais atuais\n"
            << "    help                                este menu\n"
            << "    quit                                encerra\n\n";
}

static bool dispatch(const std::string& line,
                     PhysiologyEngine& bg,
                     ActiveActions& act)
{
  std::istringstream is(line);
  std::string cmd;
  is >> cmd;
  if (cmd.empty()) return false;

  try {
    if (cmd == "exercise") {
      double i;
      if (!(is >> i)) { std::cout << "  uso: exercise <0..1>\n"; return false; }
      apply_exercise_generic(bg, i);
      std::ostringstream s; s << "ge=" << std::fixed << std::setprecision(2) << i;
      act.exercise = (i > 0.0) ? s.str() : "";
      std::cout << "  [ok] exercise generic intensity=" << i << "\n";
    }
    else if (cmd == "run") {
      double v, inc;
      if (!(is >> v >> inc)) { std::cout << "  uso: run <m/s> <incline 0-1>\n"; return false; }
      apply_exercise_run(bg, v, inc);
      std::ostringstream s; s << "run=" << v;
      act.exercise = (v > 0.0) ? s.str() : "";
      std::cout << "  [ok] run speed=" << v << " m/s, incline=" << inc << "\n";
    }
    else if (cmd == "cycle") {
      double p, c;
      if (!(is >> p >> c)) { std::cout << "  uso: cycle <W> <rpm>\n"; return false; }
      apply_exercise_cycle(bg, p, c);
      std::ostringstream s; s << "cyc=" << p << "W,c=" << c;
      act.exercise = (p > 0.0) ? s.str() : "";
      std::cout << "  [ok] cycle power=" << p << " W, cadence=" << c << "/min\n";
    }
    else if (cmd == "stress") {
      double sev;
      if (!(is >> sev)) { std::cout << "  uso: stress <0..1>\n"; return false; }
      apply_stress(bg, sev);
      std::ostringstream s; s << std::fixed << std::setprecision(2) << sev;
      act.stress = (sev > 0.0) ? s.str() : "";
      std::cout << "  [ok] stress severity=" << sev << "\n";
    }
    else if (cmd == "coldface") {
      double sev, dur;
      if (!(is >> sev >> dur)) { std::cout << "  uso: coldface <sev> <dur_s>\n"; return false; }
      apply_coldface(bg, sev, dur);
      std::ostringstream s; s << "sev=" << sev << "/" << dur << "s";
      act.coldface = (sev > 0.0) ? s.str() : "";
      std::cout << "  [ok] coldface severity=" << sev << ", duration=" << dur << " s\n";
    }
    else if (cmd == "env") {
      double T, rh;
      if (!(is >> T >> rh)) { std::cout << "  uso: env <C> <RH 0-1>\n"; return false; }
      apply_env(bg, T, rh);
      std::ostringstream s; s << T << "C/RH" << static_cast<int>(rh * 100);
      act.env = s.str();
      std::cout << "  [ok] env ambient=" << T << " C, RH=" << rh << "\n";
    }
    else if (cmd == "stop") {
      std::string what; is >> what;
      if      (what == "exercise") { clear_exercise(bg); act.exercise.clear(); }
      else if (what == "stress")   { clear_stress(bg);   act.stress.clear();   }
      else if (what == "coldface") { clear_coldface(bg); act.coldface.clear(); }
      else { std::cout << "  uso: stop <exercise|stress|coldface>\n"; return false; }
      std::cout << "  [ok] " << what << " parado\n";
    }
    else if (cmd == "status") {
      const auto& cv = *bg.GetCardiovascularSystem();
      const auto& en = *bg.GetEnergySystem();
      std::cout << "  HR="   << cv.GetHeartRate(FrequencyUnit::Per_min)        << " bpm  "
                << "MAP="    << cv.GetMeanArterialPressure(PressureUnit::mmHg) << " mmHg  "
                << "CoreT="  << en.GetCoreTemperature(TemperatureUnit::C)      << " C  "
                << "SkinT="  << en.GetSkinTemperature(TemperatureUnit::C)      << " C  "
                << "Lac="    << get_arterial_lactate_mmol_per_L(bg)            << " mmol/L\n";
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

//  main
int main()
{
  std::cout << "\n=========================================================\n"
            <<   " BioGears RealTime Driver  (Howto-RealTime, V5)\n"
            <<   "=========================================================\n";

  // 1. Inicialização do engine
  std::unique_ptr<PhysiologyEngine> bg = CreateBioGearsEngine(kLogPath);
  bg->GetLogger()->Info("HowTo-RealTime starting up...");

  std::cout << "Carregando paciente: " << kPatientXml
            << "  (a estabilizacao pode demorar 30-60 s)...\n" << std::flush;
  if (!bg->InitializeEngine(kPatientXml)) {
    bg->GetLogger()->Error("Falha ao inicializar paciente.");
    std::cerr << "[ERRO] InitializeEngine falhou. Verifique " << kLogPath << "\n";
    return 1;
  }

  if (params_calculated) {
    // Verifica se o ficheiro cfg existe
    std::ifstream check(PARAMS_CFG_PATH);
    if (!check.is_open()) {
      // Ficheiro não encontrado — avisa e continua com defaults
      // (Nervous.cpp fará o mesmo: sem cfg → usa valores originais BioGears)
      std::cerr << "[AVISO] params_calculated=true mas \""
                << PARAMS_CFG_PATH << "\" nao encontrado.\n"
                << "        Gera-o com: python generate_biogears_files.py <json>\n"
                << "        Continuando com parametros BioGears originais...\n\n";
    }
    // Se existe, Nervous.cpp vai carregá-lo automaticamente na 1ª iteração
  } else {
    // Garante que qualquer cfg anterior não interfere
    std::remove(PARAMS_CFG_PATH);
  }

  // 2. CSV header
  std::ofstream csv(kCsvPath);
  if (!csv) {
    std::cerr << "[ERRO] nao consegui abrir " << kCsvPath << "\n";
    return 1;
  }
  csv << "time_s,HR_bpm,MAP_mmHg,CoreT_C,SkinT_C,Lactate_mmol_L,active_actions\n";
  csv << std::fixed << std::setprecision(4);

  // 3. Thread de input — std::getline bloqueante; detach para não impedir saída
  CommandQueue queue;
  std::atomic<bool> running { true };
  std::thread input_thread([&]() {
    std::string line;
    while (running.load(std::memory_order_relaxed)) {
      if (!std::getline(std::cin, line)) {
        // EOF (Ctrl+D) ou stdin fechado: pede saída ordenada
        queue.push("quit");
        return;
      }
      queue.push(line);
    }
  });
  input_thread.detach();

  // 4. Loop principal
  ActiveActions actions;
  print_help();
  std::cout << "> " << std::flush;

  const double dt_s        = 1.0 / kEngineHz;
  const auto   dt_chrono   = std::chrono::microseconds(
                                static_cast<long long>(1e6 / kEngineHz));
  auto next_tick           = clk::now() + dt_chrono;
  long long tick           = 0;
  long long late_warn_cnt  = 0;
  bool quit                = false;

  // Refs persistentes — mais barato que chamar Get*System() a cada tick
  const auto& cv = *bg->GetCardiovascularSystem();
  const auto& en = *bg->GetEnergySystem();

  while (!quit) {
    // 4.1 drena fila de comandos
    std::string cmd;
    while (queue.try_pop(cmd)) {
      if (dispatch(cmd, *bg, actions)) { quit = true; break; }
      std::cout << "> " << std::flush;
    }
    if (quit) break;

    // 4.2 avança engine
    bg->AdvanceModelTime(dt_s, TimeUnit::s);
    ++tick;

    // 4.3 amostragem CSV (1 Hz)
    if ((tick % kCsvDecimateTicks) == 0) {
      const double t   = static_cast<double>(tick) * dt_s;
      const double hr  = cv.GetHeartRate(FrequencyUnit::Per_min);
      const double map = cv.GetMeanArterialPressure(PressureUnit::mmHg);
      const double Tc  = en.GetCoreTemperature(TemperatureUnit::C);
      const double Ts  = en.GetSkinTemperature(TemperatureUnit::C);
      const double lac = get_arterial_lactate_mmol_per_L(*bg);

      csv << t   << ","
          << hr  << ","
          << map << ","
          << Tc  << ","
          << Ts  << ","
          << lac << ","
          << actions.to_csv() << "\n";
      csv.flush();
    }

    // 4.4 sincronização wall-clock — sleep_until evita drift acumulado
    const auto now = clk::now();
    if (now < next_tick) {
      std::this_thread::sleep_until(next_tick);
    } else {
      // Engine atrasado: avisa raramente para não poluir stdout. Não
      // congelamos nem descartamos amostras — perde-se o real-time mas
      // a fidelidade física é preservada (engine é determinístico).
      if ((++late_warn_cnt) % 500 == 0) {
        std::cerr << "[warn] engine atrasado em "
                  << std::chrono::duration<double>(now - next_tick).count()
                  << " s (cumulativo: " << late_warn_cnt << " ticks)\n";
      }
    }
    next_tick += dt_chrono;
  }

  // 5. Cleanup
  running.store(false);
  csv.close();
  std::cout << "\nSimulacao encerrada.\n"
            << "  CSV : " << kCsvPath << "\n"
            << "  Log : " << kLogPath << "\n";
  return 0;
}
