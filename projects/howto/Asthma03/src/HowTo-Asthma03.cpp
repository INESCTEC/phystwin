/**************************************************************************************
  HowToAsthma3 - Simulação de crise asmática + registro ECG
  Autor: Lucas (ajustado para gerar dois CSVs separados)
**************************************************************************************/

#include <biogears/cdm/compartment/SECompartmentManager.h>
#include <biogears/cdm/engine/PhysiologyEngineTrack.h>
#include <biogears/cdm/patient/actions/SEAsthmaAttack.h>
#include <biogears/cdm/system/equipment/ElectroCardioGram/SEElectroCardioGram.h>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <fstream>

using namespace biogears;

int HowToAsthma3()
{
  Logger logger("HowToAsthma03.log");
  std::unique_ptr<PhysiologyEngine> bg = CreateBioGearsEngine(&logger);

  // --- Configuração dos dados fisiológicos padrão ---
  auto& drm = bg->GetEngineTrack()->GetDataRequestManager();
  drm.SetSamplesPerSecond(1);

  drm.CreatePhysiologyDataRequest().Set("HeartRate", FrequencyUnit::Per_min);
  drm.CreatePhysiologyDataRequest().Set("SkinTemperature", TemperatureUnit::C);
  drm.CreatePhysiologyDataRequest().Set("CoreTemperature", TemperatureUnit::C);
  drm.CreatePhysiologyDataRequest().Set("InspiratoryExpiratoryRatio");
  drm.CreatePhysiologyDataRequest().Set("CardiacOutput", VolumePerTimeUnit::mL_Per_min);
  drm.CreatePhysiologyDataRequest().Set("MeanArterialPressure", PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("SystolicArterialPressure", PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("DiastolicArterialPressure", PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("HemoglobinContent", MassUnit::g);
  drm.CreatePhysiologyDataRequest().Set("TotalVentilation", VolumePerTimeUnit::L_Per_min);
  drm.CreateGasCompartmentDataRequest().Set(BGE::PulmonaryCompartment::Trachea, "InFlow");

  drm.SetResultsFilename("HowToAsthma03.csv");   // ← Ficheiro fisiológico principal

  // --- Inicializa o paciente ---
  if (!bg->InitializeEngine("patients/Male_22_Fit_Soldier.xml")) {
    bg->GetLogger()->Error("Falha ao carregar paciente.");
    return 1;
  }

  // --- Cria ficheiro separado para ECG ---
  std::ofstream ecgCSV("AsthmaECG03.csv");
  ecgCSV << "Time_s,"
         << "Lead3_mV\n";

  // --- Simula baseline (2 minutos = 600s) ---
  bg->GetLogger()->Info("Simulando baseline (2min)...");
  for (int t = 0; t < 120; ++t) {
    bg->AdvanceModelTime(1, TimeUnit::s);

    const auto& ecg = *bg->GetElectroCardioGram();
    double time_s = bg->GetSimulationTime(TimeUnit::s);

    ecgCSV << time_s << ","
           << ecg.GetLead3ElectricPotential(ElectricPotentialUnit::mV) << "\n";
  }

  // --- Inicia crise asmática ---
  bg->GetLogger()->Info("Iniciando crise asmática...");
  SEAsthmaAttack asthma;
  asthma.GetSeverity().SetValue(0.5);
  bg->ProcessAction(asthma);

  // --- Simula durante o ataque (5 minutos = 300s) ---
  for (int t = 0; t < 300; ++t) {
    bg->AdvanceModelTime(1, TimeUnit::s);

    const auto& ecg = *bg->GetElectroCardioGram();
    double time_s = bg->GetSimulationTime(TimeUnit::s);

    ecgCSV << time_s << ","
           << ecg.GetLead3ElectricPotential(ElectricPotentialUnit::mV) << "\n";
  }

  // --- Recuperação ---
  bg->GetLogger()->Info("Recuperando após crise...");
  asthma.GetSeverity().SetValue(0.0);
  bg->ProcessAction(asthma);

  for (int t = 0; t < 300; ++t) {
    bg->AdvanceModelTime(1, TimeUnit::s);

    const auto& ecg = *bg->GetElectroCardioGram();
    double time_s = bg->GetSimulationTime(TimeUnit::s);

    ecgCSV << time_s << ","
           << ecg.GetLead3ElectricPotential(ElectricPotentialUnit::mV) << "\n";
  }
  //bg->AdvanceModelTime(300, TimeUnit::s); // 5 min recuperação

  // --- Finaliza ---
  ecgCSV.close();
  bg->GetLogger()->Info("Simulação concluída");


  return 0;
}

int main(int argc, char* argv[])
{
  return HowToAsthma3();
}
