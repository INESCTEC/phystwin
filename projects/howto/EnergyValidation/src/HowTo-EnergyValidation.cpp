/**************************************************************************************
Copyright 2015 Applied Research Associates, Inc.
Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the License
at:
http://www.apache.org/licenses/LICENSE-2.0
**************************************************************************************/

// -----------------------------------------------------------------------------
//  HowTo-EnergyValidation.cpp
//
//  Replica o protocolo de validação de Exercise (Tables 4a–d) da documentação
//  oficial do BioGears:
//      https://www.biogearsengine.com/documentation/_energy_methodology.html
//
//  Referência interna (método local): share/doc/methodology/EnergyMethodology.md
//
//  Protocolo (Bruce protocol variant, 30 s baseline + 4 estágios de 6 min com
//  2 min de rest entre eles; severidades 0.0375 / 0.075 / 0.1583 / 0.3583
//  ≈ 45 / 90 / 190 / 430 W).
//
//  Timeline scenario-time (s):
//      000 – 030  baseline
//      030 – 390  Severity 0.0375  (sample @ 389)
//      390 – 510  rest              (sample @ 509)
//      510 – 870  Severity 0.0750   (sample @ 869)
//      870 – 990  rest              (sample @ 989)
//      990 – 1350 Severity 0.1583   (sample @ 1349)
//     1350 – 1470 rest              (sample @ 1469)
//     1470 – 1830 Severity 0.3583   (sample @ 1829)
//     1830 – 1950 rest              (sample @ 1949)
//
//  Uso:
//      howto-energyvalidation [patient_xml] [csv_out]
//      howto-energyvalidation                              # defaults: StandardMale
//      howto-energyvalidation patients/StandardMale.xml    EnergyVal_StdMale.csv
//      howto-energyvalidation patients/Male_22_Fit_Soldier.xml EnergyVal_FitSoldier.csv
//
//  Comparar depois com:
//      /Users/lucasteixeira/Desktop/UI_App/EnergyValidation_V3/analysis.ipynb
// -----------------------------------------------------------------------------

#include <biogears/cdm/compartment/SECompartmentManager.h>
#include <biogears/cdm/engine/PhysiologyEngineTrack.h>
#include <biogears/cdm/patient/actions/SEExercise.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <biogears/cdm/system/physiology/SEBloodChemistrySystem.h>
#include <biogears/cdm/system/physiology/SECardiovascularSystem.h>
#include <biogears/cdm/system/physiology/SEEnergySystem.h>
#include <biogears/cdm/system/physiology/SERenalSystem.h>
#include <biogears/cdm/system/physiology/SERespiratorySystem.h>
#include <biogears/cdm/system/physiology/SETissueSystem.h>
#include <biogears/cdm/engine/PhysiologyEngineConfiguration.h>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/string/manipulation.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace biogears;

// -----------------------------------------------------------------------------
//  Estrutura que define um estágio do protocolo.
//  duration_s = duração da fase; intensity = 0 para rest; sample_time_s é o
//  scenario-time no qual a docs amostra o resultado (apenas para log).
// -----------------------------------------------------------------------------
struct Stage {
  const char* label;
  double intensity;
  double duration_s;
  double sample_time_s;   // informativo
};

static const std::vector<Stage> kProtocol = {
  // baseline
  {"baseline",          0.0,    30.0,  30.0 },
  // est. 1  (Table 4a)
  {"ex_0.0375_45W",     0.336, 360.0, 389.0},
  {"rest_after_4a",     0.0,    120.0, 509.0},
  // est. 2  (Table 4b)
  {"ex_0.0750_90W",     0.558, 360.0, 869.0},
  {"rest_after_4b",     0.0,    120.0, 989.0},
  // est. 3  (Table 4c)
  {"ex_0.1583_190W",    0.779, 360.0, 1349.0},
  {"rest_after_4c",     0.0,    120.0, 1469.0},
  // est. 4  (Table 4d)
  {"ex_0.3583_430W",    1.0, 360.0, 1829.0},
  {"rest_after_4d",     0.0,    120.0, 1950.0},  
};

static const bool        params_calculated = true;
static const char* const PARAMS_CFG_PATH   = "exercise_params_active.cfg";
static const char* const DEFAULT_PATIENT = "patients/StandardMale.xml";
static const char* const DEFAULT_CSV     = "EnergyValidation_StdMale.csv";
static const double      kSamplesPerSec  = 50.0;   // igual à docs

// -----------------------------------------------------------------------------
//  Regista todas as 15 métricas da Table 4 + algumas auxiliares.
//
//  Nomes seguem exactamente os usados em ExerciseStages.xml (runtime).
//  Para compatibilidade com a docs, unidades correspondem às da docs quando
//  a conversão é directa; `SystemicVascularResistance` em BioGears só tem
//  FlowResistanceUnit::mmHg_s_Per_mL — é convertido para mmHg·min/L no
//  notebook (factor 1000/60 ≈ 16.67).
// -----------------------------------------------------------------------------
static void AddValidationDataRequests(PhysiologyEngine& bg)
{
  auto& drm = bg.GetEngineTrack()->GetDataRequestManager();
  drm.SetSamplesPerSecond(kSamplesPerSec);

  // Table 4 — 15 medidas principais
  drm.CreatePhysiologyDataRequest().Set("CoreTemperature",            TemperatureUnit::C);
  drm.CreatePhysiologyDataRequest().Set("SkinTemperature",            TemperatureUnit::C);
  drm.CreatePhysiologyDataRequest().Set("HeartRate",                  FrequencyUnit::Per_min);
  drm.CreatePhysiologyDataRequest().Set("SystolicArterialPressure",   PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("DiastolicArterialPressure",  PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("MeanArterialPressure",       PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("HeartStrokeVolume",          VolumeUnit::mL);
  drm.CreatePhysiologyDataRequest().Set("CardiacOutput",              VolumePerTimeUnit::L_Per_min);
  drm.CreatePhysiologyDataRequest().Set("SystemicVascularResistance", FlowResistanceUnit::mmHg_s_Per_mL);
  drm.CreatePhysiologyDataRequest().Set("RespirationRate",            FrequencyUnit::Per_min);
  drm.CreatePhysiologyDataRequest().Set("TidalVolume",                VolumeUnit::mL);
  drm.CreatePhysiologyDataRequest().Set("RespiratoryExchangeRatio");
  drm.CreatePhysiologyDataRequest().Set("OxygenConsumptionRate",      VolumePerTimeUnit::mL_Per_min);
  drm.CreatePhysiologyDataRequest().Set("UrineProductionRate",        VolumePerTimeUnit::mL_Per_min);

  // Gut Blood Flow — medida 9 da Table 4 (inflow do compartimento gut)
  drm.CreateLiquidCompartmentDataRequest().Set("GutVasculature",      "InFlow", VolumePerTimeUnit::mL_Per_min);

  // Auxiliares úteis para a tese do Lucas (fluxo muscular, métrica do work rate
  // e trajectória metabólica — úteis para quantificar as correcções A/B/C/D).
  drm.CreateLiquidCompartmentDataRequest().Set("MuscleVasculature",   "InFlow", VolumePerTimeUnit::mL_Per_min);
  drm.CreatePhysiologyDataRequest().Set("AchievedExerciseLevel");
  drm.CreatePhysiologyDataRequest().Set("TotalWorkRateLevel");
  drm.CreatePhysiologyDataRequest().Set("FatigueLevel");
  drm.CreatePhysiologyDataRequest().Set("TotalMetabolicRate",         PowerUnit::W);
  drm.CreatePhysiologyDataRequest().Set("SweatRate",                  MassPerTimeUnit::g_Per_s);
  drm.CreatePhysiologyDataRequest().Set("MeanSkinFlow",               VolumePerTimeUnit::mL_Per_min);
}

// -----------------------------------------------------------------------------
//  Aplica uma `SEExercise::SEGeneric` com intensidade `i` (0 = rest).
// -----------------------------------------------------------------------------
static void ApplyIntensity(PhysiologyEngine& bg, double i)
{
  SEExercise::SEGeneric ge;
  ge.Intensity.SetValue(i);
  SEExercise act{ge};
  bg.ProcessAction(act);
}

// -----------------------------------------------------------------------------
//  Avança `duration_s` em incrementos de 1/kSamplesPerSec (igual à docs).
// -----------------------------------------------------------------------------
static void AdvanceFor(PhysiologyEngine& bg, double duration_s)
{
  const double dt = 1.0 / kSamplesPerSec;
  const int    N  = static_cast<int>(duration_s * kSamplesPerSec + 0.5);
  for (int k = 0; k < N; ++k) {
    bg.AdvanceModelTime(dt, TimeUnit::s);
  }
}

// -----------------------------------------------------------------------------
int RunEnergyValidation(const std::string& patient_xml, const std::string& csv_out)
{

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

  Logger logger("HowToEnergyValidation.log");
  std::unique_ptr<PhysiologyEngine> bg = CreateBioGearsEngine(&logger);

  AddValidationDataRequests(*bg);
  bg->GetEngineTrack()->GetDataRequestManager().SetResultsFilename(csv_out);

  bg->GetLogger()->Info(
    std::string("Energy Validation run — patient=") + patient_xml +
    "  csv=" + csv_out);

  if (!bg->InitializeEngine(patient_xml)) {
    bg->GetLogger()->Error("Could not load patient state; abort.");
    return 1;
  }

  // Cabeçalho no log com os tempos-alvo de amostragem (para referência humana)
  bg->GetLogger()->Info(
    "Protocol samples (scenario-time s): 30, 389, 509, 869, 989, 1349, 1469, 1829, 1949");

  double elapsed = 0.0;
  for (const Stage& s : kProtocol) {
    ApplyIntensity(*bg, s.intensity);
    bg->GetLogger()->Info(asprintf(
      "[t=%6.1f s] stage=%-20s intensity=%.4f  duration=%.0f s  (sample@%.0f s)",
      elapsed, s.label, s.intensity, s.duration_s, s.sample_time_s));
    AdvanceFor(*bg, s.duration_s);
    elapsed += s.duration_s;
  }

  bg->GetLogger()->Info(asprintf(
    "Finished. Total scenario duration = %.1f s (expected 1950).", elapsed));
  return 0;
}

// -----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  std::string patient = (argc >= 2) ? argv[1] : DEFAULT_PATIENT;
  std::string csv     = (argc >= 3) ? argv[2] : DEFAULT_CSV;

  std::cout << "[HowTo-EnergyValidation]\n"
            << "  patient : " << patient << "\n"
            << "  csv_out : " << csv     << "\n"
            << "  sps     : " << kSamplesPerSec << "\n"
            << "  stages  : " << kProtocol.size() << "  (expected total = 1950 s)\n";
  return RunEnergyValidation(patient, csv);
}
